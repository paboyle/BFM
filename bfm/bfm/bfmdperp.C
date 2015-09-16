#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include "bfm.h"
#include "bfm_cg_unprec.h"

// IroIro uses opposite sign gamma5
template <class Float>
void bfmbase<Float>::GeneralisedFiveDimEnd(void)
{
  if ( IsGeneralisedFiveDim() ) { 
    free(zdata -> a);
    free(zdata -> ap);
    free(zdata -> alpha);
    free(zdata -> beta);
    free(zdata -> gamma);
    free(zdata);
    freeFermion(g5d_tmp);
  }
}


template <class Float>
void bfmbase<Float>::GeneralisedFiveDimInit(void)
{
  if ( IsGeneralisedFiveDim() ) { 

    g5d_tmp = this->allocFermion();

    // Do we ever need Higham?
    // Zolotarev coeffs
    this->BossMessage("Generalised Five Dimensional Matrix for %s\n",this->SolverString(this->solver));
    this->BossMessage("\t Ls = %d\n",this->Ls);
    this->BossMessage("\t mq = %le\n",this->mass);
    this->BossMessage("\t M5 = %le\n",this->M5);

    double R=(1+this->mass)/(1-this->mass);

    ////////////////////////////
    //ContFrac:
    //  Ls always odd. Rational poly deg is either Ls or Ls-1
    //PartFrac 
    //  Ls always odd. Rational poly deg is either Ls or Ls-1
    //
    //Tanh: Ls always even, Rational poly deg is Ls
    // 
    // Just set nrational as Ls. Forget about Ls-1 cases.
    //
    // Require odd Ls for cont and part frac
    ////////////////////////////
    int  nrational = this->Ls;

    if ( this->Ls<=0 ) {
      this->Error("Ls must be positive %d\n",this->Ls);
      double d;
      double *pointer = (double *) NULL;
      d= *pointer;
      this->Error("I've just accessed a null pointer to get core dump and back trace to find out what moron called me like this %le\n",d);
    }


    if (   ( this->solver == HwPartFracTanh)
	 ||( this->solver == HwContFracTanh)
	 ||( this->solver == HtContFracTanh)
	   //	 ||( this->solver == HtPartFracTanh)
	   ){
      nrational=this->Ls-1;
      if ( (this->Ls&01) != 1 ) {
	this->BossMessage("Ls must be odd for Cont and Part Frac\n");
      }
      this->BossMessage("Cont/Part Frac Tanh approx:");
      this->BossMessage("Ls=%d but taking even RationalPoly=%d so matching Ls=%d DWF equivalent\n",Ls,nrational,nrational);
    }
    if (   ( this->solver == HwPartFracZolo)
	 ||( this->solver == HwContFracZolo) 
	   //	 ||( this->solver == HtPartFracZolo) 
	 ||( this->solver == HtContFracZolo) 
	   ){
      if ( (this->Ls&01) != 1 ) {
	this->BossMessage("Ls must be odd for Cont and Part Frac\n");
      }
      this->BossMessage("Cont/Part Frac Zolo approx:");
      this->BossMessage("Ls=%d and taking odd RationalPoly=%d for greatest accuracy\n",Ls,nrational);
    }

    if (   ( this->solver == HwPartFracZolo)
	 ||( this->solver == HwContFracZolo)
	   //||( this->solver == HtPartFracZolo)
	 ||( this->solver == HtContFracZolo)
	 ||( this->solver == HwCayleyZolo)
	 ||( this->solver == HtCayleyZolo)
	   ) {
      // Chroma uses a cast to "float" for the zolotarev range.
      // this creates a real difference in the operator which I do not like but we can replicate here
      // to demonstrate compatibility
      this->eps = (float)(this-> zolo_lo / this -> zolo_hi);
      // this->eps = (this-> zolo_lo / this -> zolo_hi);
      this->BossMessage("zolotarev range: (%f,%f)\n",this->zolo_lo,this->zolo_hi);
      this->BossMessage("zolotarev(eps=%f,Nrat=%d)\n",this->eps,nrational);
      zdata = bfm_zolotarev(this->eps,nrational,0);
    } else if (
 	         ( this->solver == HwPartFracTanh)
	       ||( this->solver == HwContFracTanh)
		 //	       ||( this->solver == HtPartFracTanh)
	       ||( this->solver == HtContFracTanh)
	       ||( this->solver == HwCayleyTanh)
	       ||( this->solver == HtCayleyTanh)
	       ||( this->solver == HmCayleyTanh) 
	       ) {
      this->zolo_hi=1;
      this->zolo_lo=1;
      this->eps = (float)(this-> zolo_lo / this -> zolo_hi);
      this->BossMessage("higham(%f,%d)\n",this->eps,nrational);
      zdata = bfm_higham(this->eps,nrational);
    } else if ( (this->solver == DWFTransfer) || (this->solver == DWFTransferInv)) {
      this->zolo_hi=1;
      this->zolo_lo=1;
      this->eps = (float)(this-> zolo_lo / this -> zolo_hi);
      this->BossMessage("higham(%f,%d)\n",this->eps,nrational);
      zdata = bfm_higham(this->eps,nrational);
      this->BossMessage("Setting up to do DWF transfer matrix\n");
      return;
    } else { 
      this->BossMessage("Don't know this action\n"); exit(-1);
    }

    if ( (solver == HwPartFracZolo) || (solver==HwPartFracTanh)
	 //|| (solver == HtPartFracZolo) 
	 //|| (solver==HtPartFracTanh)
	 ){ 
      if ( this->Ls != (2*zdata->da -1) ) { 
	this->Error("PartFrac Ls mismatch %d %d %d\n",zdata->da,this->Ls);
	exit(-1);
      }
    }

    if ( (solver == HwContFracZolo) 
	 ||(solver==HwContFracTanh) 
	 ||(solver==HtContFracZolo) 
	 ||(solver==HtContFracTanh) 
	 ){ 
      if ( this->Ls != zdata->db ) { 
	this->Error("ContFrac Ls mismatch %d %d\n",zdata->db,this->Ls);
	exit(-1);
      }
    }
    if (   ( this->solver == HwCayleyZolo)
	 ||( this->solver == HtCayleyZolo)
	 ||( this->solver == HwCayleyTanh)
	 ||( this->solver == HtCayleyTanh)
	 ||( this->solver == HmCayleyTanh) 
	   ) {
      if ( this->Ls != zdata->n ) { 
	this->Error("Cayley Ls mismatch %d %d\n",zdata->n,this->Ls);
	exit(-1);
      }
    }
    p.resize(zdata->da);
    q.resize(zdata->dd);


    if (   ( this->solver == HwCayleyZolo)
	 ||( this->solver == HtCayleyZolo)
	 ||( this->solver == HtContFracZolo)
	   //	 ||( this->solver == HtPartFracZolo)
	 ||( this->solver == HwContFracZolo)
	 ||( this->solver == HwPartFracZolo) ) {
      this->zolo_delta = zdata->Delta;
      this->BossMessage("zolo_delta = %20.12le\n",(double)this->zolo_delta);
    } else { 
      this->zolo_delta = 0;
    }

    //////////////////////////////////////////////////////////////
    // Part frac matrix coeffs
    // 
    // For Zolo case to scale the range from zolo_hi to 1.0 the 
    // present implementation is a function of "H" explicitly
    // (implicitly in Cayley case) and so it is easiest to simply scale 
    // H where it enters the 5d matrix. This is done in the Munprec_PartFrac routine.
    //
    // Here do nothing about approx max.
    //////////////////////////////////////////////////////////////

    for(int n=0;n<zdata->da;n++){
      p[n] = zdata -> alpha[n];
    }
    for(int n=0;n<zdata->dd;n++){
      q[n] = -zdata -> ap[n];
    }

    //////////////////////////////////////////////////////////////
    // Cont frac matrix coeffs
    // 
    // For Zolo case to scale the range from zolo_hi to 1.0 the 
    // present implementation is a function of "H" explicitly
    // (implicitly in Cayley case) and so it is easiest to simply scale 
    // H where it enters the 5d matrix. This is done in the Munprec_ContFrac routine.
    //
    // Here do nothing about approx max.
    //////////////////////////////////////////////////////////////

    Beta.resize(this->Ls);
    cc.resize(this->Ls);
    cc_d.resize(this->Ls);
    sqrt_cc.resize(this->Ls);
    for(int i=0; i < this->Ls ; i++){
      Beta[i] = zdata -> beta[i];
      cc[i] = 1.0/Beta[i];
      cc_d[i]=sqrt(cc[i]);
    }
    
    cc_d[this->Ls-1]=1.0;
    for(int i=0; i < this->Ls-1 ; i++){
      sqrt_cc[i]= sqrt(cc[i]*cc[i+1]);
    }    
    sqrt_cc[this->Ls-2]=sqrt(cc[this->Ls-2]);

    ///////////////////////////////////////////////////////////
    // The Cayley forms
    //
    // NB, Possible scale in Zolotarev case as spectrum does not end at "1".
    // Scale for Higham case is implicit in the Moebius scale of shamir (or whatever)
    // 
    // In the Zolo case, we have (see my notes) 
    // T_s^{-1} =  -(H_m+omega_s)/(H_m -omega_s  )
    //
    // To scale "H" in a product of "T_s" by a factor 1/zolo_hi, we can rewrite
    // as omega_s -> omega_s * scale. This mapps sgn(zolo_hi) to sgn(1.0) as required.
    //
    ///////////////////////////////////////////////////////////

    omega.resize(this->Ls);
    bs.resize(this->Ls);
    cs.resize(this->Ls);
    as.resize(this->Ls);

    double scale=1.0;
    if ( (solver==HwCayleyZolo) || (solver==HtCayleyZolo)  ){
      scale = this->zolo_hi;
    }

    if(solver==HmCayleyTanh) { 
      this->BossMessage("HmCayleyTanh scale %le\n",this->mobius_scale);
    } else if (solver== HwCayleyTanh) {
      this->BossMessage("HwCayleyTanh\n");
    } else if ( (solver == HtCayleyTanh) ||(solver== HtCayleyZolo)){  // Ht Tanh and Zolo
	this->BossMessage("HtCayleyTanh/HtCayleyZolo\n");
    }

    for(int i=0; i < zdata -> n; i++){
      as[i] = 1.0;
      omega[i] = ((double)zdata -> gamma[i]) * scale;
      if(solver==HmCayleyTanh) { 
	double bb=this->mobius_scale;
	bs[i] = 0.5*(bb/omega[i] +1.0); // Sanity check... zolo_hi>1 increases bb/omega term
	cs[i] = 0.5*(bb/omega[i] -1.0); // Same effect as  mobius_scale=1.7 for example
      } else if ( solver== HwCayleyZolo) {
	double bb=1.0;
	bs[i] = cs[i] = bb/(omega[i]);
      } else if (solver== HwCayleyTanh) {
	double bb=1.0;
	bs[i] = cs[i] = bb/(omega[i]);
      } else if ( (solver == HtCayleyTanh) ||(solver== HtCayleyZolo)){  // Ht Tanh and Zolo
	double bb=1.0;
	bs[i] = 0.5*(bb/(omega[i]) + 1.0);
	cs[i] = 0.5*(bb/(omega[i]) - 1.0);
      }
    }

  ////////////////////////////////////////////////////////
  // Constants for the preconditioned matrix Cayley form
  ////////////////////////////////////////////////////////
  if ( (solver == HmCayleyTanh)
    || (solver == HtCayleyTanh)
    || (solver == HwCayleyTanh)
    || (solver == HwCayleyZolo)
    || (solver == HtCayleyZolo)
       ) {

    bee.resize(this->Ls);
    cee.resize(this->Ls);
    beo.resize(this->Ls);
    ceo.resize(this->Ls);
    for(int i=0;i<this->Ls;i++){
      bee[i]=as[i]*(bs[i]*(4.0-this->M5) +1.0);
      cee[i]=as[i]*(1.0-cs[i]*(4.0-this->M5));
      beo[i]=as[i]*bs[i];
      ceo[i]=-as[i]*cs[i];
    }

    ///////////////////////////////////////////
    // Support for fifth dimension multigrid
    ///////////////////////////////////////////
    aee.resize(this->Ls);
    aeo.resize(this->Ls);
    for(int i=0;i<this->Ls;i++){
      aee[i]=cee[i];
      aeo[i]=ceo[i];
    }

    //////////////////////////////////////////
    // LDU decomposition of eeoo
    //////////////////////////////////////////
    dee.resize(this->Ls);
    lee.resize(this->Ls);
    leem.resize(this->Ls);
    uee.resize(this->Ls);
    ueem.resize(this->Ls);

    for(int i=0;i<this->Ls;i++){

      dee[i] = bee[i];

      if ( i < this->Ls-1 ) {

	lee[i] =-cee[i+1]/bee[i]; // sub-diag entry on the ith column

	leem[i]=this->mass*cee[this->Ls-1]/bee[0];
	for(int j=0;j<i;j++)  leem[i]*= aee[j]/bee[j+1];

	uee[i] =-aee[i]/bee[i];   // up-diag entry on the ith row

	ueem[i]=this->mass;
	for(int j=1;j<=i;j++) ueem[i]*= cee[j]/bee[j];
	ueem[i]*= aee[0]/bee[0];

      } else { 
	lee[i] =0.0;
	leem[i]=0.0;
	uee[i] =0.0;
	ueem[i]=0.0;
      }
    }

    { 
      double delta_d=this->mass*cee[this->Ls-1];
      for(int j=0;j<this->Ls-1;j++) delta_d *= cee[j]/bee[j];
      dee[this->Ls-1] += delta_d;
    }

  }


  if ( (solver == HwContFracZolo) 
       ||(solver==HwContFracTanh) 
       ||(solver==HtContFracZolo) 
       ||(solver==HtContFracTanh) 
       ){ 

    double scale =1.0;
    if ( this->solver==HwContFracZolo ) scale=1.0/this->zolo_hi;
    if ( this->solver==HtContFracZolo ) scale=1.0/this->zolo_hi;
    double dw_diag = (4.0-this->M5)*scale;
    
    See.resize(this->Ls);
    Aee.resize(this->Ls);
    int sign=1;
    for(int s=0;s<this->Ls;s++){
      Aee[s] = sign * Beta[s] * dw_diag;
      sign   = - sign;
    }
    Aee[this->Ls-1] += R;
    
    See[0] = Aee[0];
    for(int s=1;s<this->Ls;s++){
      See[s] = Aee[s] - 1.0/See[s-1];
    }
    for(int s=0;s<this->Ls;s++){
      this->BossDebug("s = %d Beta %le Aee %le See %le\n",s,Beta[s],Aee[s],See[s]);
    }
  }

  if(0){
    // Print out the matrix
    if ( (this->solver == HwPartFracTanh)||(this->solver == HwPartFracZolo) ){
      
    int nblock=(nrational-1)/2;
    for(int iblock=0;iblock<nblock;iblock++){
      for(int ii=0;ii<2;ii++){
	this->BossDebug("[");
	for(int jblock=0;jblock<nblock;jblock++){
	  if ( iblock==jblock ){
	    if ( ii==0) { 
	      this->BossDebug("        H %12.6e   ",-sqrt(q[nblock-iblock-1]));
	    } else {
	      this->BossDebug("%12.6e          -H ",-sqrt(q[nblock-iblock-1]));
	    }
	  } else {
  	      this->BossDebug("%12.6e   %12.6e   ",0.0,0.0);
	  }
	  }
	if ( ii==0){
	  this->BossDebug("%12.6e",sqrt(p[nblock-1-iblock]));
	} else { 
	  this->BossDebug("%12.6e",0);
	}
	this->BossDebug("]\n");

      }
    }
    this->BossDebug("[");
    for(int jblock=0;jblock<nblock;jblock++){
      this->BossDebug("%12.6e   %12.6e   ",sqrt(p[nblock-1-jblock]),0.0);
    }
    this->BossDebug("R g5 + %12.6e H",p[nblock]);
    this->BossDebug("]\n");

       
    } else if( (this->solver == HwContFracTanh)||(this->solver == HwContFracZolo) ) {
      
      this->BossDebug("Continued Fraction Representation:\n");
      int sign=1;
      for (int i = 0; i < this->Ls; i++){
	this->BossDebug("[");
	for (int j = 0; j < this->Ls; j++){

	  if ( (i==this->Ls-1) && (j==this->Ls-1) ) {
	    this->BossDebug("R g5 + %12.6e H",Beta[this->Ls-1]);
	  } else if (i==j){

	    if ( sign==1 ){
	      this->BossDebug("        H ");
	    } else {
	      this->BossDebug("     -  H ");
	    }
	    
	  }else if ( (j==(i-1)) ) {
	    this->BossDebug("%12.6e   ",sqrt_cc[j]);
	  }else if ( (i==(j-1)) ) {
	    this->BossDebug("%12.6e   ",sqrt_cc[i]);
	  }else{
	    this->BossDebug("%12.6e   ",0.0);
	  }
	}
	this->BossDebug("]\n");
	sign=-sign;
      }
    
    } else {
    
      this->BossDebug("Cayley representation\n");
      for (int i = 0; i < this->Ls; i++){
	this->BossDebug("[");
	for (int j = 0; j < this->Ls; j++){
	  if (i==j){
	    this->BossDebug("(%8.2e D + 1)   ",bs[i]);
	  }else if (j==(i+1)) {
	    this->BossDebug("-(1-%8.2e D) P- ",cs[i]);
	  }else if (j==(i-1)) {
	    this->BossDebug("-(1-%8.2e D) P+ ",cs[i]);
	  }else{


	    this->BossDebug("   %11.1f        ",0.0);
	  }
	}
	this->BossDebug("]\n");
      }
    }
  }
  }
}

//////////////////////////////////////////////
// Suport routines for mobius red black solver
//////////////////////////////////////////////
template <class Float>
double bfmbase<Float>::G5D_Mprec(Fermion_t psi, 
				 Fermion_t chi, 
				 Fermion_t tmp, 
				 int dag,int donrm) 
{
  int cb=this->rb_precondition_cb;
  int ocb=1-cb;

  if ( (solver == HmCayleyTanh)  || (solver == HtCayleyTanh) || (solver == HwCayleyTanh)  
    || (solver == HwCayleyZolo)  || (solver == HtCayleyZolo) || (solver == HwPartFracZolo)
    || (solver == HwPartFracTanh)||(solver == HwContFracZolo)||(solver == HwContFracTanh)
       ) {
    
    int me = this->thread_barrier();
    //Cayley form matrix
    double nrm;
    double m1=-1.0;

    // Moo + Moe Meeinv Meo
    G5D_Meo (psi,tmp,ocb,dag);   // odd->even 
    G5D_MooeeInv(tmp,chi,dag);    // even->even 
    G5D_Meo (chi,tmp,cb,dag);    // even->odd 
    G5D_Mooee   (psi,chi,dag);    // odd odd 
    nrm = axpy_norm (chi,tmp,chi,m1);     // Linear combination
    return nrm;
      
  } else if ( (solver == HtContFracZolo)||(solver == HtContFracTanh)){
    this->Error("Do not yet support Partial/Cont fraction as preconditioned system\n");
    exit(-1);
  } else if ((solver == DWFTransfer)
	   ||(solver == DWFTransferInv)){
    this->Error("Do not yet support Transfer matrix as preconditioned system\n");
    exit(-1);
  } else { 
    return Mprec(psi,chi,tmp,dag,donrm);
  }

}


template <class Float>
void bfmbase<Float>::G5D_Meo(Fermion_t psi, // Source
			     Fermion_t chi, // Result
			     int result_cb,
			     int dag)
{
  int Pminus=-1;
  int Pplus=1;
  double m=-0.5; // Factor of -0.5 in front of D_W=4+m-.5 dslash for hopping term
  if ( (solver == HmCayleyTanh)
    || (solver == HtCayleyTanh)
    || (solver == HwCayleyTanh)
    || (solver == HwCayleyZolo)
    || (solver == HtCayleyZolo)
       ) {

    int me = thread_barrier();


    Fermion_t tmp = g5d_tmp; // Hack.. add to bfm.h as meo_tmp
    
    if ( dag ) { 

      // Apply 4d dslash
    uint64_t t1=GetTimeBase();
    dslash(psi,tmp,result_cb,dag);
    uint64_t t2=GetTimeBase();
    this->G5D_CayleyMeoFiveD(tmp,chi,dag);
    uint64_t t3=GetTimeBase();
#if 0
      // Assemble the 5d matrix
      for(int s=0;s<Ls;s++){
      if ( s==0 ) {
	axpby_ssp_proj(chi,m*beo[s],tmp,   -m*ceo[s+1]  ,tmp,s,s+1,Pplus);
	axpby_ssp_proj(chi,   1.0,chi,m*mass*ceo[Ls-1],tmp,s,Ls-1,Pminus);
      } else if ( s==(Ls-1)) { 
	axpby_ssp_proj(chi,m*beo[s],tmp,m*mass*ceo[0],tmp,s,0,Pplus);
	axpby_ssp_proj(chi,1.0,chi,-m*ceo[s-1],tmp,s,s-1,Pminus);
      } else {
	axpby_ssp_proj(chi,m*beo[s],tmp,-m*ceo[s+1],tmp,s,s+1,Pplus);
	axpby_ssp_proj(chi,1.0   ,chi,-m*ceo[s-1],tmp,s,s-1,Pminus);
      }
      }
#endif

    if ( this->iter == this->time_report_iter+5 ) {
      this->ThreadBossPerformance("bfm:Cayley Meo: Dw dslash \t\t: %d cyc\n ",t2-t1);
      this->ThreadBossPerformance("bfm:Cayley Meo: s-direction    \t\t: %d cyc\n ",t3-t2);
    }

    } else { 

      this->G5D_CayleyMeoFiveD(psi,tmp,dag);
#if 0
      // Assemble the 5d matrix
      for(int s=0;s<Ls;s++){
	if ( s==0 ) {
	  //	tmp = bs psi[s] + cs[s] psi[s+1}
	  //      tmp+= -mass*cs[s] psi[s+1}
	  axpby_ssp_proj(tmp,m*beo[s],psi,-m*ceo[s],psi ,s, s+1,Pminus);
	  axpby_ssp_proj(tmp,1.0,tmp,m*mass*ceo[s],psi,s,Ls-1,Pplus);
	} else if ( s==(Ls-1)) { 
	  axpby_ssp_proj(tmp,m*beo[s],psi,m*mass*ceo[s],psi,s,0,Pminus);
	  axpby_ssp_proj(tmp,1.0,tmp,-m*ceo[s],psi,s,s-1,Pplus);
	} else {
	  axpby_ssp_proj(tmp,m*beo[s],psi,-m*ceo[s],psi,s,s+1,Pminus);
	  axpby_ssp_proj(tmp,1.0,tmp,-m*ceo[s],psi,s,s-1,Pplus);
	}
      }
#endif
      // Apply 4d dslash
      dslash(tmp,chi,result_cb,dag);

    }

  } else if ((solver == HwPartFracZolo)||(solver == HwPartFracTanh)){

    int sign = 1;
    double scale=1.0;
    double amax =1.0;
    
    if ( dag ) sign = -1;
    if (part_frac_chroma_convention) scale=2.0; // Chroma conventions annoy me
    if (this->solver==HwPartFracZolo) amax=this->zolo_hi;

    dslash(psi,chi,result_cb,0); // Dslash on diagonal. g5 Dslash is hermitian
    int nblock=(Ls-1)/2;
    for(int b=0;b<nblock;b++){
      int s = 2*b;
      ag5xpby_ssp(chi,-m*scale,chi,0.0,chi,s,s); 
      ag5xpby_ssp(chi, m*scale,chi,0.0,chi,s+1,s+1); 
    }
    ag5xpby_ssp(chi,p[nblock]*scale*m/amax,chi,0.0,chi,Ls-1,Ls-1);

  } else if ((solver == HwContFracZolo)||(solver == HwContFracTanh)){

    double scale =1.0;
    if ( this->solver==HwContFracZolo ) scale=1.0/this->zolo_hi;

    dslash(psi,chi,result_cb,0); // Dslash on diagonal. g5 Dslash is hermitian
    
    int sign=1;
    for(int s=0;s<Ls;s++){
      if ( s==(Ls-1) ){
	ag5xpby_ssp(chi,Beta[s]*m*scale,chi,0.0,chi,s,s);
      } else {
	ag5xpby_ssp(chi,cc[s]*Beta[s]*sign*m*scale,chi,0.0,chi,s,s);
      }
      sign=-sign; 
    }

  } else {
    this->Error("preconditioned unimplemented \n");
    exit(-1);
  }
}

template <class Float>
void bfmbase<Float>::G5D_Mooee(Fermion_t psi, 
			       Fermion_t chi, 
			       int dag) 
{

  int Pminus=-1;
  int Pplus=1;
  if ( (solver == HmCayleyTanh)
    || (solver == HtCayleyTanh)
    || (solver == HwCayleyTanh)
    || (solver == HwCayleyZolo)
    || (solver == HtCayleyZolo)
       ) {
    G5D_CayleyMooeeFiveD(psi,chi,dag);
#if 0
    Fermion_t tmp = this->g5d_tmp;
    if ( dag ) {
      for (int s=0;s<Ls;s++){
	// Assemble the 5d matrix
	if ( s==0 ) {
	  axpby_ssp_proj(chi,bee[s],psi,-cee[s+1]  ,psi,s,s+1,Pplus);
	  axpby_ssp_proj(chi,1.0,chi,mass*cee[Ls-1],psi,s,Ls-1,Pminus);
	} else if ( s==(Ls-1)) { 
	  axpby_ssp_proj(chi,bee[s],psi,mass*cee[0],psi,s,0,Pplus);
	  axpby_ssp_proj(chi,1.0,chi,-cee[s-1],psi,s,s-1,Pminus);
	} else {
	  axpby_ssp_proj(chi,bee[s],psi,-cee[s+1],psi,s,s+1,Pplus);
	  axpby_ssp_proj(chi,1.0   ,chi,-cee[s-1],psi,s,s-1,Pminus);
	}
      }
    } else { 
      for (int s=0;s<Ls;s++){
	if ( s==0 ) {
	  axpby_ssp_proj(chi,bee[s],psi ,-cee[s],psi,s,s+1,Pminus);
	  axpby_ssp_proj(chi,1.0,chi,mass*cee[s],psi,s,Ls-1,Pplus);
	} else if ( s==(Ls-1)) { 
	  axpby_ssp_proj(chi,bee[s],psi,mass*cee[s],psi,s,0,Pminus);
	  axpby_ssp_proj(chi,1.0,chi,-cee[s],psi,s,s-1,Pplus);
	} else {
	  axpby_ssp_proj(chi,bee[s],psi,-cee[s],psi,s,s+1,Pminus);
	  axpby_ssp_proj(chi,1.0,chi,-cee[s],psi,s,s-1,Pplus);
	}
      }
    }
#endif

  } else if ((solver == HwPartFracZolo)||(solver == HwPartFracTanh)){

    int sign = 1;
    double scale=1.0;
    double amax =1.0;
    double dw_diag = (4.0-this->M5);

    if ( dag ) sign = -1;
    if (part_frac_chroma_convention) scale=2.0; // Chroma conventions annoy me
    if (this->solver==HwPartFracZolo) amax=this->zolo_hi;

    int nblock=(Ls-1)/2;
    for(int b=0;b<nblock;b++){

      int s = 2*b;
      double pp = p[nblock-1-b];
      double qq = q[nblock-1-b];

      // Do each 2x2 block aligned at s and multiplies Dw site diagonal by G5 so Hw
      ag5xpby_ssp(chi,-dw_diag*scale,psi,amax*sqrt(qq)*scale,psi, s  ,s+1); 
      ag5xpby_ssp(chi, dw_diag*scale,psi,amax*sqrt(qq)*scale,psi, s+1,s);
      axpby_ssp  (chi, 1.0, chi,sqrt(amax*pp)*scale*sign,psi,s+1,Ls-1);
    }

    {
      double R=(1+this->mass)/(1-this->mass);
      //R g5 psi[Ls-1] + p[0] H
      ag5xpbg5y_ssp(chi,R*scale,psi,p[nblock]*scale*dw_diag/amax,psi,Ls-1,Ls-1);

      for(int b=0;b<nblock;b++){
	int s = 2*b+1;
	double pp = p[nblock-1-b];
	axpby_ssp(chi,1.0,chi,-sqrt(amax*pp)*scale*sign,psi,Ls-1,s);
      }
    }

  } else if ((solver == HwContFracZolo)||(solver == HwContFracTanh)){


    double scale =1.0;
    if ( this->solver==HwContFracZolo ) scale=1.0/this->zolo_hi;
    double dw_diag = (4.0-this->M5)*scale;
    
    int sign=1;
    for(int s=0;s<Ls;s++){
      if ( s==0 ) {
	ag5xpby_ssp(chi,cc[0]*Beta[0]*sign*dw_diag,psi,sqrt_cc[0],psi,s,s+1); // Multiplies Dw by G5 so Hw
      } else if ( s==(Ls-1) ){
	// Drop the CC here.
	double R=(1+this->mass)/(1-this->mass);
	ag5xpby_ssp(chi,Beta[s]*dw_diag,psi,sqrt_cc[s-1],psi,s,s-1);
	ag5xpby_ssp(chi,R,psi,1.0,chi,s,s);
      } else {
	ag5xpby_ssp(chi,cc[s]*Beta[s]*sign*dw_diag,psi,sqrt_cc[s],psi,s,s+1);
	axpby_ssp(chi,1.0,chi,sqrt_cc[s-1],psi,s,s-1);
      }
      sign=-sign; 
    }
    
  } else {
    this->Error("preconditioned unimplemented \n");
    exit(-1);
  }
}

template <class Float>
void bfmbase<Float>::G5D_MooeeInv(Fermion_t psi, 
				  Fermion_t chi, 
				  int dag) 
{
  int Pminus=-1;
  int Pplus=1;
  if ( (solver == HmCayleyTanh)
    || (solver == HtCayleyTanh)
    || (solver == HwCayleyTanh)
    || (solver == HwCayleyZolo)
    || (solver == HtCayleyZolo)
       ) {

    Fermion_t tmp = this->g5d_tmp;

    G5D_CayleyMooeeInvFiveD(psi,chi,dag);

#if 0
    if ( dag ) {
      // Apply (U^{\prime})^{-dagger}
      axpby_ssp (chi,1.0,psi,     0.0,psi,0,0);      // chi[0]=psi[0]
      for (int s=1;s<Ls;s++){
	axpby_ssp_proj(chi,1.0,psi,-uee[s-1],chi,s,s-1,Pminus);
      }
      // U_m^{-\dagger} 
      for (int s=0;s<Ls-1;s++){
	axpby_ssp_proj(chi,1.0,chi,-ueem[s],chi,Ls-1,s,Pplus);
      }
      // L_m^{-\dagger} D^{-dagger}
      for (int s=0;s<Ls-1;s++){
	axpby_ssp_proj(chi,1.0/dee[s],chi,-leem[s]/dee[Ls-1],chi,s,Ls-1,Pminus);
      }	
      axpby_ssp(chi,1.0/dee[Ls-1],chi,0.0,chi,Ls-1,Ls-1); // Modest avoidable 

      // Apply L^{-dagger}
      for (int s=Ls-2;s>=0;s--){
	axpby_ssp_proj (chi,1.0,chi,-lee[s],chi,s,s+1,Pplus);  // chi[Ls]
      }
    } else {  // NOT DAGGER
      // Apply (L^{\prime})^{-1}
      axpby_ssp (chi,1.0,psi,     0.0,psi,0,0);      // chi[0]=psi[0]
      for (int s=1;s<Ls;s++){
	axpby_ssp_proj(chi,1.0,psi,-lee[s-1],chi,s,s-1,Pplus);// recursion Psi[s] -lee P_+ chi[s-1]
      }
      // L_m^{-1} 
      for (int s=0;s<Ls-1;s++){ // Chi[ee] = 1 - sum[s<Ls-1] -leem[s]P_- chi
	axpby_ssp_proj(chi,1.0,chi,-leem[s],chi,Ls-1,s,Pminus);
      }
      // U_m^{-1} D^{-1}
      for (int s=0;s<Ls-1;s++){
	// Chi[s] + 1/d chi[s] 
	axpby_ssp_proj(chi,1.0/dee[s],chi,-ueem[s]/dee[Ls-1],chi,s,Ls-1,Pplus);
      }	
      axpby_ssp(chi,1.0/dee[Ls-1],chi,0.0,chi,Ls-1,Ls-1); // Modest avoidable 

      // Apply U^{-1}
      for (int s=Ls-2;s>=0;s--){
	axpby_ssp_proj (chi,1.0,chi,-uee[s],chi,s,s+1,Pminus);  // chi[Ls]
      }
    }
#endif
  } else if ((solver == HwPartFracZolo)||(solver == HwPartFracTanh)){

    int sign = 1;
    if ( dag) sign = -1;
    double scale=1.0;
    double amax =1.0;
    if (this->solver==HwPartFracZolo) amax=this->zolo_hi;

    double dw_diag = (4.0-this->M5);
    Fermion_t tmp = this->g5d_tmp;
    if (part_frac_chroma_convention) scale=2.0; // Chroma conventions annoy me

    ///////////////////////////////////////////////////////////////////////////////////////
    //Linv
    ///////////////////////////////////////////////////////////////////////////////////////
    int nblock=(Ls-1)/2;

    axpy(chi,psi,psi,0.0); // Identity piece

    for(int b=0;b<nblock;b++){
      int s = 2*b;
      double pp = p[nblock-1-b];
      double qq = q[nblock-1-b];
      double coeff1=sign*sqrt(amax*amax*amax*pp*qq) / ( dw_diag*dw_diag + amax*amax* qq);
      double coeff2=sign*sqrt(amax*pp)*dw_diag / ( dw_diag*dw_diag + amax*amax* qq); // Implicit g5 here
      axpby_ssp  (chi,1.0,chi,coeff1,psi,this->Ls-1,s);
      axpbg5y_ssp(chi,1.0,chi,coeff2,psi,this->Ls-1,s+1);
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////
    //Dinv (note D isn't really diagonal -- just diagonal enough that we can still invert)
    // Compute Seeinv (coeff of gamma5)
    ///////////////////////////////////////////////////////////////////////////////////////
    double R=(1+this->mass)/(1-this->mass);
    double Seeinv = R + p[nblock]*dw_diag/amax;
    for(int b=0;b<nblock;b++){
      Seeinv += p[nblock-1-b]*dw_diag/amax / ( dw_diag*dw_diag/amax/amax + q[nblock-1-b]);
    }    
    Seeinv = 1.0/Seeinv;

    for(int b=0;b<nblock;b++){
      int s = 2*b;
      double pp = p[nblock-1-b];
      double qq = q[nblock-1-b];
      double coeff1=dw_diag / ( dw_diag*dw_diag + amax*amax* qq); // Implicit g5 here
      double coeff2=amax*sqrt(qq) / ( dw_diag*dw_diag + amax*amax* qq);
      ag5xpby_ssp  (tmp,-coeff1,chi,coeff2,chi,s,s+1);
      ag5xpby_ssp  (tmp, coeff1,chi,coeff2,chi,s+1,s);
    }
    ag5xpby_ssp  (tmp, Seeinv,chi,0.0,chi,this->Ls-1,this->Ls-1);

    ///////////////////////////////////////////////////////////////////////////////////////
    // Uinv
    ///////////////////////////////////////////////////////////////////////////////////////
    for(int b=0;b<nblock;b++){
      int s = 2*b;
      double pp = p[nblock-1-b];
      double qq = q[nblock-1-b];
      double coeff1=-sign*sqrt(amax*amax*amax*pp*qq) / ( dw_diag*dw_diag + amax*amax* qq);
      double coeff2=-sign*sqrt(amax*pp)*dw_diag / ( dw_diag*dw_diag + amax*amax* qq); // Implicit g5 here
      axpby_ssp  (chi,1.0/scale,tmp,coeff1/scale,tmp,s,this->Ls-1);
      axpbg5y_ssp(chi,1.0/scale,tmp,coeff2/scale,tmp,s+1,this->Ls-1);
    }
    axpby_ssp  (chi, 1.0/scale,tmp,0.0,tmp,this->Ls-1,this->Ls-1);

  } else if ((solver == HwContFracZolo)||(solver == HwContFracTanh)){
#if 1
    // Apply Linv
    axpby_ssp(chi,1.0/cc_d[0],psi,0.0,psi,0,0);
    for(int s=1;s<Ls;s++){
      axpbg5y_ssp(chi,1.0/cc_d[s],psi,-1.0/See[s-1],chi,s,s-1);
    }
    // Apply Dinv
    for(int s=0;s<Ls;s++){
      ag5xpby_ssp(chi,1.0/See[s],chi,0.0,chi,s,s); //only appearance of See[0]
    }
    // Apply Uinv = (Linv)^T
    axpby_ssp(chi,1.0/cc_d[Ls-1],chi,0.0,chi,this->Ls-1,this->Ls-1);
    for(int s=Ls-2;s>=0;s--){
      axpbg5y_ssp(chi,1.0/cc_d[s],chi,-1.0*cc_d[s+1]/See[s]/cc_d[s],chi,s,s+1);
    }
#else
    // Debug of LDU -- apply LDU and directly compare to Mee
    // Apply U
    for(int s=0;s<Ls-1;s++){
      axpbg5y_ssp(chi,cc_d[s]*1.0,psi,cc_d[s+1]*1.0/See[s],psi,s,s+1);
    }
    axpby_ssp(chi,cc_d[Ls-1]*1.0,psi,0.0,psi,this->Ls-1,this->Ls-1);
    // Apply D
    for(int s=0;s<Ls;s++){
      ag5xpby_ssp(chi,See[s],chi,0.0,chi,s,s); //only appearance of See[0]
    }
    // Apply L
    for(int s=Ls-1;s>0;s--){
      axpbg5y_ssp(chi,1.0,chi,1.0/See[s-1],chi,s,s-1);
    }
    axpby_ssp(chi,1.0,chi,0.0,chi,0,0);
    //Dtilde
    for(int s=0;s<Ls;s++){
      axpby_ssp(chi,cc_d[s],chi,0.0,chi,s,s);
    }

#endif


  } else {
    this->Error("preconditioned unimplemented \n");
    exit(-1);
  }
}


template <class Float>
void bfmbase<Float>::G5D_Munprec(Fermion_t psi[2], 
				 Fermion_t chi[2], 
				 int dag)
{

  if ( solveMobiusDminus ) {
    G5D_Dminus(psi,chi,dag);
  } else if ( (solver == HmCayleyTanh)
    || (solver == HtCayleyTanh)
    || (solver == HwCayleyTanh)
    || (solver == HwCayleyZolo)
    || (solver == HtCayleyZolo)
       ) {
    G5D_Munprec_Cayley(psi,chi,dag);
  } else if ((solver == HwPartFracZolo)||(solver == HwPartFracTanh)){
    G5D_Munprec_PartFrac(psi,chi,dag);
  } else if ((solver == HwContFracZolo)||(solver == HwContFracTanh)
	   ||(solver == HtContFracZolo)||(solver == HtContFracTanh)){
    G5D_Munprec_ContFrac(psi,chi,dag);
  } else if (solver == DWFTransfer) {
 // T      =  [ 2g5 +(1+g5) D_W ]^{-1} [ 2g5-(1-g5) D_W ]
    int sgn=1; // Denominator has plus sign
    G5D_TransferNumDen   (psi,chi,sgn,dag);
  } else if (solver == DWFTransferInv) { 
 // T^{-1} =  [ 2g5-(1-g5) D_W ]^{-1} [ 2g5+(1+g5) D_W ]
    int sgn=-1;// Denominator has minus sign
    G5D_TransferNumDen   (psi,chi,sgn,dag);

  } else { 
    Fermion_t tmp = threadedAllocFermion();
    Munprec(psi,chi,tmp,dag);
    threadedFreeFermion(tmp);
  }
}

// Needed tools
// Z[s] = A X[s]+ B P+/- Y[s']  // axpby_ssp_proj
// Z[s] = A X[s]+ B Y[s']       // axpby_ssp
// Z[s] = A G5 X[s] + B Y[s']   // ag5xpby_ssp
// Z[s] = A G5 X[s] + B G5 Y[s']// ag5xpbg5y_ssp
void  bfm_ptr_check (void *ptr,int size=0);
template <class Float>
void bfmbase<Float>::G5D_Munprec_Cayley(Fermion_t psi[2], 
					Fermion_t chi[2], 
					int dag)
{
  Fermion_t Din[2];
  Din[0]  = threadedAllocFermion(); 
  Din[1]  = threadedAllocFermion(); 
  int Pminus=-1;
  int Pplus=1;

  int words = 24 * nsimd * simd_cbvol *cbLs;
  int bytes = words*sizeof(Float);

  bfm_ptr_check(Din[0],bytes);
  bfm_ptr_check(Din[1],bytes);
  printf("Din initial check done\n");fflush(stdout);
  if ( dag ) { 
    // Under adjoint
    //D1+        D1- P-    ->   D1+^dag   P+ D2-^dag
    //D2- P+     D2+            P-D1-^dag D2+dag
    for(int cb=0;cb<2;cb++){
      // Apply Dw
      G5D_DW(psi,Din[cb],cb,dag); 
    }
    printf("dperp dag checking after Dw\n");fflush(stdout);
    bfm_ptr_check(Din[0],bytes);
    bfm_ptr_check(Din[1],bytes);
    printf("checking after Dw done\n");fflush(stdout);

    for(int cb=0;cb<2;cb++){
      for(int s=0;s<Ls;s++){
        // Collect the terms in DW
	//	Chi = bs Din[s] + cs[s] Din[s+1}
	//    Chi+= -mass*cs[s] psi[s+1}
	if ( s==0 ) {
	  axpby_ssp_proj(chi[cb],bs[s],Din[cb],cs[s+1],Din[cb],s,s+1,Pplus);
	  axpby_ssp_proj(chi[cb],1.0,chi[cb],-mass*cs[Ls-1],Din[cb],s,Ls-1,Pminus);
	} else if ( s==(Ls-1)) { 
	  axpby_ssp_proj(chi[cb],bs[s],Din[cb],-mass*cs[0],Din[cb],s,0,Pplus);
	  axpby_ssp_proj(chi[cb],1.0,chi[cb],cs[s-1],Din[cb],s,s-1,Pminus);
	} else {
	  axpby_ssp_proj(chi[cb],bs[s],Din[cb],cs[s+1],Din[cb],s,s+1,Pplus);
	  axpby_ssp_proj(chi[cb],1.0,chi[cb],cs[s-1],Din[cb],s,s-1,Pminus);
	}
	// Collect the terms indept of DW
	if ( s==0 ){
	  axpby_ssp_proj(chi[cb],1.0,chi[cb],-1.0,psi[cb],s,s+1,Pplus);
	  axpby_ssp_proj(chi[cb],1.0,chi[cb],mass,psi[cb],s,Ls-1,Pminus);
	} else if ( s==(Ls-1)) {
	  axpby_ssp_proj(chi[cb],1.0,chi[cb],mass,psi[cb],s,0,Pplus);
	  axpby_ssp_proj(chi[cb],1.0,chi[cb],-1.0,psi[cb],s,s-1,Pminus);
	} else {
	  axpby_ssp_proj(chi[cb],1.0,chi[cb],-1.0,psi[cb],s,s+1,Pplus);
	  axpby_ssp_proj(chi[cb],1.0,chi[cb],-1.0,psi[cb],s,s-1,Pminus);
	}
      }
      // ((b D_W + D_w hop terms +1) on s-diag
      axpby (chi[cb],chi[cb],psi[cb],1.0,1.0); 
      
    }
    printf("dperp dag checking after ssp\n");fflush(stdout);
    bfm_ptr_check(Din[0],bytes);
    bfm_ptr_check(Din[1],bytes);
    printf("checking after ssp done\n");fflush(stdout);


  } else { 
    for(int cb=0;cb<2;cb++){
      // Assemble Din
      for(int s=0;s<Ls;s++){
	if ( s==0 ) {
	  //	Din = bs psi[s] + cs[s] psi[s+1}
	  axpby_ssp_proj(Din[cb],bs[s],psi[cb],cs[s],psi[cb],s,s+1,Pminus);
	  //      Din+= -mass*cs[s] psi[s+1}
	  axpby_ssp_proj(Din[cb],1.0,Din[cb],-mass*cs[s],psi[cb],s,Ls-1,Pplus);
	} else if ( s==(Ls-1)) { 
	  axpby_ssp_proj(Din[cb],bs[s],psi[cb],-mass*cs[s],psi[cb],s,0,Pminus);
	  axpby_ssp_proj(Din[cb],1.0,Din[cb],cs[s],psi[cb],s,s-1,Pplus);
	} else {
	  axpby_ssp_proj(Din[cb],bs[s],psi[cb],cs[s],psi[cb],s,s+1,Pminus);
	  axpby_ssp_proj(Din[cb],1.0,Din[cb],cs[s],psi[cb],s,s-1,Pplus);
	}
      }
    }     

    printf("dperp checking after ssp\n");fflush(stdout);
    bfm_ptr_check(Din[0],bytes);
    bfm_ptr_check(Din[1],bytes);
    printf("checking after ssp done\n");fflush(stdout);


    for(int cb=0;cb<2;cb++){
      // Apply Dw
      G5D_DW(Din,chi[cb],cb,dag); 
      // ((b D_W + D_w hop terms +1) on s-diag
      axpby (chi[cb],chi[cb],psi[cb],1.0,1.0); 
    }    

    printf("dperp checking after Dw\n");fflush(stdout);
    bfm_ptr_check(Din[0],bytes);
    bfm_ptr_check(Din[1],bytes);
    printf("checking after Dw done\n");fflush(stdout);


    for(int cb=0;cb<2;cb++){
      // Dmobius
      for(int s=0;s<Ls;s++){
	if ( s==0 ){
	  axpby_ssp_proj(chi[cb],1.0,chi[cb],-1.0,psi[cb],s,s+1,Pminus);
	  axpby_ssp_proj(chi[cb],1.0,chi[cb],mass,psi[cb],s,Ls-1,Pplus);
	} else if ( s==(Ls-1)) {
	  axpby_ssp_proj(chi[cb],1.0,chi[cb],mass,psi[cb],s,0,Pminus);
	  axpby_ssp_proj(chi[cb],1.0,chi[cb],-1.0,psi[cb],s,s-1,Pplus);
	} else {
	  axpby_ssp_proj(chi[cb],1.0,chi[cb],-1.0,psi[cb],s,s+1,Pminus);
	  axpby_ssp_proj(chi[cb],1.0,chi[cb],-1.0,psi[cb],s,s-1,Pplus);
	}
      }
    }

  }
    printf("checking after alldone \n");fflush(stdout);
    bfm_ptr_check(Din[0],bytes);
    bfm_ptr_check(Din[1],bytes);
    printf("checking after Dw done\n");fflush(stdout);


  printf("Freeing Din\n");fflush(stdout);
  threadedFreeFermion(Din[0]);
  threadedFreeFermion(Din[1]);
}

template <class Float>
void bfmbase<Float>::G5D_DW(Fermion_t psi[2],Fermion_t chi,int cb,int dag,double scale)
{
  this->dslash(psi[1-cb],chi,cb,dag);
  double coeff_psi = ( 4.0 - this->M5 )*scale;
  double coeff_chi = -0.5*scale;
  this->axpby(chi,psi[cb],chi,coeff_psi,coeff_chi);
}

template <class Float>
void bfmbase<Float>::G5D_Dminus(Fermion_t psi[2],Fermion_t chi[2],int dag)
{
  for (int cb=0;cb<2;cb++){
      G5D_DW(psi,chi[cb],cb,dag);       // chi = D_W psi
      for(int s=0;s<Ls;s++){
	axpby_ssp(chi[cb],1.0,psi[cb],-cs[s],chi[cb],s,s);// chi = (1-c[s] D_W) psi
      }
  }
}
template <class Float>
void bfmbase<Float>::G5D_Dplus(Fermion_t psi[2],Fermion_t chi[2],int dag)
{
  for (int cb=0;cb<2;cb++){
      G5D_DW(psi,chi[cb],cb,dag);       // chi = D_W psi
      for(int s=0;s<Ls;s++){
	axpby_ssp(chi[cb],1.0,psi[cb],bs[s],chi[cb],s,s);// chi = (1+b[s] D_W) psi
      }
  }
}




template <class Float>
void bfmbase<Float>::G5D_HW(Fermion_t psi[2],Fermion_t chi,int cb,int dag)
{
  G5D_DW(psi,chi,cb,dag);
  this->Error("g5 mult not implemented\n");
  exit(-1);
}

template <class Float>
  void bfmbase<Float>::G5D_Munprec_PartFrac(Fermion_t psi[2], 
				  Fermion_t chi[2], 
				  int dag)
{
  Fermion_t D= g5d_tmp;


  // For partial frac Hw case (b5=c5=1) chroma quirkily computes
  
  int sign = 1;
  if ( dag ) sign = -1;
  double scale=1.0;
  if (part_frac_chroma_convention) scale=2.0;

  //
  // Conventions for partfrac appear to be a mess.
  // Tony's Nara lectures have
  //
  // BlockDiag(  H/p_i  1             | 1       )    
  //          (  1      p_i H / q_i^2 | 0       )  
  //           ---------------------------------
  //           ( -1      0                | R  +p0 H  )
  //
  //Chroma     ( -2H    2sqrt(q_i)    |   0         )
  //           (2 sqrt(q_i)   2H      |  2 sqrt(p_i) )
  //           ---------------------------------
  //           ( 0     -2 sqrt(p_i)   |  2 R gamma_5 + p0 2H
  //
  // Edwards/Joo/Kennedy/Wenger
  //
  //
  // Here, the "beta's" selected by chroma to scale the unphysical bulk constraint fields
  // incorporate the approx scale factor. This is obtained by propagating the
  // scale on "H" out to the off diagonal elements as follows:
  //
  // BlockDiag(  H/p_i  1             | 1       ) 
  //          (  1      p_i H / q_i^2 | 0       )  
  //           ---------------------------------
  //          ( -1      0                | R  + p_0 H  )
  //
  // becomes:
  // BlockDiag(  H/ sp_i  1               | 1             ) 
  //          (  1      sp_i H / s^2q_i^2 | 0             )  
  //           ---------------------------------
  //           ( -1      0                | R + p_0/s H   )
  //
  //
  // This is implemented in Chroma by
  //           p0' = p0/approxMax
  //           p_i' = p_i*approxMax
  //           q_i' = q_i*approxMax*approxMax
  //
  // After the equivalence transform is applied the matrix becomes
  // 
  //Chroma     ( -2H    sqrt(q'_i)    |   0         )
  //           (sqrt(q'_i)   2H       |   sqrt(p'_i) )
  //           ---------------------------------
  //           ( 0     -sqrt(p'_i)    |  2 R gamma_5 + p'0 2H
  //
  //     =     ( -2H    sqrt(q_i)amax    |   0              )
  //           (sqrt(q_i)amax   2H       |   sqrt(p_i*amax) )
  //           ---------------------------------
  //           ( 0     -sqrt(p_i)*amax   |  2 R gamma_5 + p0/amax 2H
  //
  for(int cb=0;cb<2;cb++){

    double amax =1.0;
    if ( this->solver==HwPartFracZolo ) amax=this->zolo_hi;

    G5D_DW(psi,D,cb,0,1.0); 

    int nblock=(Ls-1)/2;
    for(int b=0;b<nblock;b++){

      int s = 2*b;
      double pp = p[nblock-1-b];
      double qq = q[nblock-1-b];

      // Do each 2x2 block aligned at s and
      ag5xpby_ssp(chi[cb],-1.0*scale,D,amax*sqrt(qq)*scale,psi[cb], s  ,s+1); // Multiplies Dw by G5 so Hw
      ag5xpby_ssp(chi[cb], 1.0*scale,D,amax*sqrt(qq)*scale,psi[cb], s+1,s);

      // Pick up last column
      axpby_ssp  (chi[cb], 1.0, chi[cb],sqrt(amax*pp)*scale*sign,psi[cb],s+1,Ls-1);

    }

    {
      double R=(1+this->mass)/(1-this->mass);
      //R g5 psi[Ls] + p[0] H
      ag5xpbg5y_ssp(chi[cb],R*scale,psi[cb],p[nblock]*scale/amax,D,Ls-1,Ls-1);

      for(int b=0;b<nblock;b++){
	int s = 2*b+1;
	double pp = p[nblock-1-b];
	axpby_ssp(chi[cb],1.0,chi[cb],-sqrt(amax*pp)*scale*sign,psi[cb],Ls-1,s);
      }
    }
  }

}

template <class Float>
  void bfmbase<Float>::G5D_Munprec_ContFrac(Fermion_t psi[2], 
					    Fermion_t chi[2], 
					    int dag)
{
  // This matrix is already hermitian. (g5 Dw) = Dw dag g5 = (g5 Dw)dag
  // The rest of matrix is symmetric.
  // Can ignore "dag".

  if ( (this->solver==HwContFracZolo) ||
       (this->solver==HwContFracTanh) ){

    Fermion_t D = g5d_tmp;

    for(int cb=0;cb<2;cb++){

      double scale =1.0;
      if ( this->solver==HwContFracZolo ) scale=1.0/this->zolo_hi;
      G5D_DW(psi,D,cb,0,scale); 

      int sign=1;
      for(int s=0;s<Ls;s++){
	if ( s==0 ) {
	  ag5xpby_ssp(chi[cb],cc[0]*Beta[0]*sign,D,sqrt_cc[0],psi[cb],s,s+1); // Multiplies Dw by G5 so Hw
	} else if ( s==(Ls-1) ){
	  // Drop the CC here.
	  double R=(1+this->mass)/(1-this->mass);
	  ag5xpby_ssp(chi[cb],Beta[s],D,sqrt_cc[s-1],psi[cb],s,s-1);
	  ag5xpby_ssp(chi[cb],R,psi[cb],1.0,chi[cb],s,s);
	  
	} else {
	  ag5xpby_ssp(chi[cb],cc[s]*Beta[s]*sign,D,sqrt_cc[s],psi[cb],s,s+1);
  	  axpby_ssp(chi[cb],1.0,chi[cb],sqrt_cc[s-1],psi[cb],s,s-1);
	}
	sign=-sign; 
      }
    }

  } else if (   (this->solver == HtContFracTanh) 
		||(this->solver == HtContFracZolo)  ) {
    // Hw/2+Dw kernel
    // Multipy through by (2+D2)
    //

    Fermion_t D = g5d_tmp;

    for(int cb=0;cb<2;cb++){

      double scale =1.0;
      if ( this->solver==HwContFracZolo ) scale=1.0/this->zolo_hi;
      G5D_DW(psi,D,cb,0,scale); 

      int sign=1;
      for(int s=0;s<Ls;s++){
	if ( s==0 ) {
	  ag5xpby_ssp(chi[cb],cc[0]*Beta[0]*sign,D,sqrt_cc[0]*2.0,psi[cb],s,s+1); // Multiplies Dw by G5 so Hw
	  axpby_ssp(chi[cb],1.0,chi[cb],sqrt_cc[0],D,s,s+1); 
	} else if ( s==(Ls-1) ){
	  // Drop the CC here.
	  double R=(1+this->mass)/(1-this->mass);
	  ag5xpby_ssp(chi[cb],Beta[s],D,sqrt_cc[s-1]*2.0,psi[cb],s,s-1);
	  ag5xpby_ssp(chi[cb],R,psi[cb],1.0,chi[cb],s,s);
	  axpby_ssp(chi[cb],1.0,chi[cb],sqrt_cc[s-1],D,s,s-1); 
	} else {
	  ag5xpby_ssp(chi[cb],cc[s]*Beta[s]*sign,D,sqrt_cc[s]*2.0,psi[cb],s,s+1);
  	  axpby_ssp(chi[cb],1.0,chi[cb],sqrt_cc[s-1]*2.0,psi[cb],s,s-1);
	  axpby_ssp(chi[cb],1.0,chi[cb],sqrt_cc[s-1],D,s,s-1); 
	  axpby_ssp(chi[cb],1.0,chi[cb],sqrt_cc[s  ],D,s,s+1); 
	}
	sign=-sign; 
      }
    }
    
  } else { 
    this->Error("ContFrac: bad solver type\n");
    exit(-1);
  }

}


template <class Float>
void bfmbase<Float>::DperpDWFcompat_cpp (Fermion_t psi_t,Fermion_t chi_t,int cb,int dag) 
{
#if 0
  Float *psi = (Float *) psi_t;
  Float *chi = (Float *) chi_t;
  
  //Interesting dependence on cb because of the BC's
  int me = this->thread_barrier();
  if ( me == 0 ) {

  int nspinco = 12;
  int x[4];

  for ( int s=0;s<this->Ls;s++ ) { 
  for ( x[3]=0;x[3]<this->simd_latt[3];x[3]++ ) { 
  for ( x[2]=0;x[2]<this->simd_latt[2];x[2]++ ) { 
  for ( x[1]=0;x[1]<this->simd_latt[1];x[1]++ ) { 
  for ( x[0]=0;x[0]<this->simd_latt[0];x[0]++ ) { 

    int p = (s+x[3]+x[2]+x[1]+x[0])& 0x1;

    if ( p == cb ) { 

      Float c_p,  c_m ;
      Float *p_m,*p_p,*ch;
      int s_p,s_m;

      int site = this->psite(x,this->simd_latt);

      c_m = c_p = 2.0;

      if ( s== 0    ) c_m = -2*this->mass; // mf boundary wrap
      if ( s== this->Ls-1 ) c_p = -2*this->mass;

      s_p = (s+1)%this->Ls;
      s_m = (s+this->Ls-1)%this->Ls;

      p_m = &psi[this->bagel_idx5d(site,s_m,nspinco)];
      p_p = &psi[this->bagel_idx5d(site,s_p,nspinco)];
      ch  = &chi[this->bagel_idx5d(site,s,nspinco)];

      for(int spinco=0;spinco<nspinco/2;spinco++){
      for(int simd_el=0;simd_el<this->nsimd;simd_el++){
      for(int reim=0;reim<2;reim++){

        int idx = reim + simd_el*2 + spinco*this->nsimd*2;
        int hdx = reim + simd_el*2 + (spinco+6)*this->nsimd*2;

        if ( dag ) {
          ch[idx]+= c_p*p_p[idx]; // 1+gamma_5
          ch[hdx]+= c_m*p_m[hdx]; // 1-gamma_5
        }else{
          ch[idx]+= c_m*p_m[idx]; // 1+gamma_5
          ch[hdx]+= c_p*p_p[hdx]; // 1-gamma_5
        }

      }}}


    }
  }}}}}
  }
  this->thread_barrier();
#endif
  exit(-1);
}


template <class Float>
void bfmbase<Float>::G5D_TransferNumDen   (Fermion_t psi[2],Fermion_t chi[2],int sgn,int dag)
{
  
  if ( (this->Ls != 1) ) { 
    this->Error("Oops bad Ls for G5D Transfer\n");
    exit(-1);
  }
  if ( (this->solver!= DWFTransfer) && (this->solver!=DWFTransferInv ) ) {
    this->Error("Oops bad solver for G5D Transfer\n");
    exit(-1);
  }
  if ( abs(sgn) != 1 ) { 
    this->Error("Oops bad sign for G5D Transfer\n");
    exit(-1);
  }

  Fermion_t D[2];
  D[0]  = this->threadedAllocFermion(); 
  D[1]  = this->threadedAllocFermion(); 

  int s=0;

  // b-c =1 ; b+c = 2b-1 = m
  double b,c;
  double m = mobius_scale;
  b=(m+1.0)/2.0;
  c=b-1.0;
 
  double cc;
  if ( sgn==1 )    cc=1.0;
  else             cc=-1.0;

  this->BossLog("G5D_Transfer: mob=%f b=%f c=%f cc=%f dag=%d\n",m,b,c,cc,dag);

#if 1
  if ( dag==1 ) {
  // Generalised pieces needed for mobius:  
  // [ 2g5+ [cc*(b+c)+(b-c)g5] D_W^{\dagger}  ] 
    for(int cb=0;cb<2;cb++)   G5D_DW(psi,D[cb],cb,dag); 
    for(int cb=0;cb<2;cb++)   ag5xpbg5y_ssp(chi[cb],2.0,psi[cb],b-c,D[cb],s,s);
    for(int cb=0;cb<2;cb++)   axpby_ssp(chi[cb],1.0,chi[cb],cc*(b+c),D[cb],s,s);
  } else { 
  // [ 2g5+ D_W [ cc*(b+c)+(b-c)g5 ] ]
    for(int cb=0;cb<2;cb++)   axpbg5y_ssp(D[cb],cc*(b+c),psi[cb],(b-c),psi[cb],s,s);
    for(int cb=0;cb<2;cb++)   G5D_DW(D,chi[cb],cb,dag); 
    for(int cb=0;cb<2;cb++)   ag5xpby_ssp(chi[cb],2.0,psi[cb], 1.0,chi[cb],s,s);
  }
#else
  int Pminus=-1;
  int Pplus=1;
  if ( dag==1 ) {

    for(int cb=0;cb<2;cb++){
      G5D_DW(psi,D[cb],cb,dag); 
      if ( sgn==1 ) 
	ag5xpby_ssp_proj(chi[cb],2.0,psi[cb], 2.0,D[cb],s,s,Pplus);
      else 
	ag5xpby_ssp_proj(chi[cb],2.0,psi[cb],-2.0,D[cb],s,s,Pminus);
    }
  } else { 
    for(int cb=0;cb<2;cb++){
      if ( sgn==1 ) 
	axpby_ssp_proj(D[cb],0.0,psi[cb],2.0,psi[cb],s,s,Pplus);
      else 
	axpby_ssp_proj(D[cb],0.0,psi[cb],-2.0,psi[cb],s,s,Pminus);
    }
    for(int cb=0;cb<2;cb++){
      G5D_DW(D,chi[cb],cb,dag); 
      ag5xpby_ssp(chi[cb],2.0,psi[cb], 1.0,chi[cb],s,s);
    }
  }

#endif
  this->threadedFreeFermion(D[0]);
  this->threadedFreeFermion(D[1]);
}


 // Transfer matrix is 
 // 
 // T^{-1} =  [ 2g5- D_W(1-g5) ]^{-1} [ 2g5+ D_W(1+g5) ]
 // T      =  [ 2g5+ D_W(1+g5) ]^{-1} [ 2g5- D_W(1-g5) ]
 //
 // Applying this requires solving
 //
 // chi = T^{-1}  psi => [ 2g5+ D_W (1+g5) ] chi =[ 2g5 - D_W (1-g5)]psi = eta
 // 

template <class Float>
void bfmbase<Float>::G5D_Transfer   (Fermion_t psi[2],Fermion_t chi[2])
{

  int sgn;
  int dag=0;


  Fermion_t eta[2];
  eta[0] = this->threadedAllocFermion();
  eta[1] = this->threadedAllocFermion();

  if ( this->solver==DWFTransfer ) {
    sgn=-1; // Numerator has minus sign
  } else if ( this->solver==DWFTransferInv ) {
    sgn=1;  // Numerator has plus sign
  } else { 
    this->Error("G5D_Transfer: bad solver type %d\n",this->solver);
    exit(-1);
  }
  G5D_TransferNumDen(psi,eta,sgn,dag);

  // Invert the denominator
  this->CGNE_unprec(chi,eta);

  this->threadedFreeFermion(eta[0]);    
  this->threadedFreeFermion(eta[1]);    
}
#define ALIGNIT(A) (double *)( (((uint64_t)A) + 31)& (~0x1FUL) );

template <class Float>
void bfmbase<Float>::G5D_CayleyMooeeFiveD(Fermion_t psi,
					  Fermion_t chi,
					  int dag)
{
  int me,thrlen,throff;
  this->thread_work(simd_cbvol,me,thrlen,throff);

  Float *pp = ((Float * )psi) + throff*this->Ls*24*this->nsimd;
  Float *cc = ((Float * )chi) + throff*this->Ls*24*this->nsimd;

  double coeffs[3*this->Ls+1];
  double * ptr = coeffs;

  integer neighbours[this->Ls*2];

  integer args[8];
  args[0] = (integer)pp; // input
  args[1] = (integer)cc; //output
  args[2] = thrlen;
  args[3] = this->Ls;
  args[4] = (integer)coeffs;
  args[5] = (integer)neighbours;


  if ( this -> cbLs != this->Ls ) {
    this->Error("CayleyMeFiveD : Oops ... supposed to support 4d even odd\n");
    exit(-1);
  }

  for(int s=0;s<this->Ls;s++){

    int s_m = s-1;
    int s_p = s+1;

    if ( s == 0 ) s_m = this->Ls-1;
    if ( s == this->Ls-1 ) s_p = 0;

    if ( dag ) {

      neighbours[2*s]   = s_m;
      if ( s == 0 )  *(ptr++) = mass*aee[Ls-1];
      else           *(ptr++) = -aee[s-1];

      *(ptr++) = bee[s];

      neighbours[2*s+1] = s_p;
      if ( s == (Ls-1) )  *(ptr++) = mass*cee[0];
      else                *(ptr++) = -cee[s+1];

    } else {

      neighbours[2*s]   = s_p; // P_- 
      if ( s == (Ls-1) )  *(ptr++) = mass*aee[s];
      else                *(ptr++) = -aee[s];

      *(ptr++) = bee[s];
      
      neighbours[2*s+1] = s_m; // P_+
      if ( s == 0 )  *(ptr++) = mass*cee[s];
      else           *(ptr++) = -cee[s];

    }

  }
  
  if(sizeof(Float) == sizeof(double))
    vmx_cayley((integer)args);
  else
    vmx_cayley_s((integer)args);

  this->thread_barrier();
}

template <class Float>
void bfmbase<Float>::G5D_CayleyMeoFiveD(Fermion_t psi,
					Fermion_t chi,
					int dag)
{
  int me,thrlen,throff;
  this->thread_work(simd_cbvol,me,thrlen,throff);

  Float *pp = ((Float * )psi) + throff*this->Ls*24*this->nsimd;
  Float *cc = ((Float * )chi) + throff*this->Ls*24*this->nsimd;

  double coeffs[3*this->Ls+1];
  double * ptr = coeffs;

  integer neighbours[this->Ls*2];

  integer args[8];
  args[0] = (integer)pp; // input
  args[1] = (integer)cc; //output
  args[2] = thrlen;
  args[3] = this->Ls;
  args[4] = (integer)coeffs;
  args[5] = (integer)neighbours;


  if ( this -> cbLs != this->Ls ) {
    this->Error("CayleyMeoFiveD : Oops ... supposed to support 4d even odd\n");
    exit(-1);
  }

  for(int s=0;s<this->Ls;s++){

    int s_m = s-1;
    int s_p = s+1;

    double m = -0.5;

    if ( s == 0 ) s_m = this->Ls-1;
    if ( s == this->Ls-1 ) s_p = 0;

    if ( dag ) {

      neighbours[2*s]   = s_m;
      if ( s == 0 )  *(ptr++) = m*mass*aeo[Ls-1];
      else           *(ptr++) = -m*aeo[s-1];

       *(ptr++) = m*beo[s];

      neighbours[2*s+1] = s_p;
      if ( s == (Ls-1) )  *(ptr++) = m*mass*ceo[0];
      else                *(ptr++) = -m*ceo[s+1];

    } else {

      neighbours[2*s]   = s_p; // P_-                                              
      if ( s == (Ls-1) )  *(ptr++) = m*mass*aeo[s];

      else                *(ptr++) = -m*aeo[s];

       *(ptr++) = m*beo[s];
      
      neighbours[2*s+1] = s_m; // P_+                                        
      if ( s == 0 )  *(ptr++) = m*mass*ceo[s];
      else           *(ptr++) = -m*ceo[s];

    }

  }
  
  if(sizeof(Float) == sizeof(double))
    vmx_cayley((integer)args);
  else
    vmx_cayley_s((integer)args);

  this->thread_barrier();
}


template <class Float>
void bfmbase<Float>::G5D_CayleyMooeeInvFiveD(Fermion_t psi,
					     Fermion_t chi,
					     int dag)
{
  int me,thrlen,throff;
  this->thread_work(simd_cbvol,me,thrlen,throff);

  Float *pp = ((Float * )psi) + throff*this->Ls*24*this->nsimd;
  Float *cc = ((Float * )chi) + throff*this->Ls*24*this->nsimd;

  double coeffs[5*this->Ls]; // uee ueem dee leem lee
  double *ptr = coeffs;

  integer args[8];
  args[0] = (integer)pp; // input
  args[1] = (integer)cc; //output
  args[2] = thrlen;
  args[3] = this->Ls;
  args[4] = (integer)coeffs;

  if ( this -> cbLs != this->Ls ) {
    this->Error("CayleyMeeInvFiveD : Oops ... supposed to support 4d even odd\n");
    exit(-1);
  }


  if ( dag ) {

    for(int s=0;s<this->Ls;s++) {
      if(s==0)  *(ptr++) = 1.0;
      else      *(ptr++) = -uee[s-1];

      if(s==this->Ls-1)  *(ptr++) = 1.0;
      else               *(ptr++) = -ueem[s];

    }

     *(ptr++) = 1.0/dee[this->Ls-1];

    for(int s=1;s<this->Ls;s++) {
       *(ptr++) = 1.0/dee[this->Ls-1-s];
       *(ptr++) = -lee[this->Ls-1-s];
       *(ptr++) = -leem[this->Ls-1-s];
    }

  } else {

    for(int s=0;s<this->Ls;s++) {
      if(s==0)  *(ptr++) = 1.0;
      else      *(ptr++) = -lee[s-1];

      if(s==this->Ls-1)  *(ptr++) = 1.0;
      else               *(ptr++) = -leem[s];

    }

     *(ptr++) = 1.0/dee[this->Ls-1];

    for(int s=1;s<this->Ls;s++) {
       *(ptr++) = 1.0/dee[this->Ls-1-s];
       *(ptr++) = -uee[this->Ls-1-s];
       *(ptr++) = -ueem[this->Ls-1-s];
    }
  }

  if(dag){ 
    if(sizeof(Float) == sizeof(double))
      vmx_cayley_inv_dag((integer)args);
    else
      vmx_cayley_inv_dag_s((integer)args);
  } else {
    if(sizeof(Float) == sizeof(double))
      vmx_cayley_inv((integer)args);
    else
      vmx_cayley_inv_s((integer)args);
  }
  this->thread_barrier();
}

//template <class Float>
//void bfmbase<Float>::G5D_MprecDeriv(Fermion_t chi, 
//			Fermion_t psi, 
//			Matrix_t force[2], 
//			int dag)
//{
//}
template <class Float>
void bfmbase<Float>::G5D_MeoDeriv(Fermion_t X, 
				  Fermion_t Y,
				  Matrix_t force,
				  int cb,
				  int dag)
{
  int Pminus=-1;
  int Pplus=1;
  double m=-0.5; // Factor of -0.5 in front of D_W=4+m-.5 dslash for hopping term
  if ( (solver == HmCayleyTanh)
    || (solver == HtCayleyTanh)
    || (solver == HwCayleyTanh)
    || (solver == HwCayleyZolo)
    || (solver == HtCayleyZolo)
       ) {


    int me = thread_barrier();

    Fermion_t tmp = g5d_tmp; // Hack.. add to bfm.h as meo_tmp
    int mdag = 1-dag;

    if ( dag ) { 

      //      X d/du [D_w D5]^dag Y = X D5^dag d/du DW^dag Y
      G5D_CayleyMeoFiveD(X,tmp,mdag);
      DslashDeriv(tmp,Y,force,cb,dag);
      
    } else { 

      //      X d/du [D_w D5] Y = X d/du DW D5 Y
      this->G5D_CayleyMeoFiveD(Y,tmp,dag);
      DslashDeriv(X,tmp,force,cb,dag);

    }

  } else { 
    this->Error("Oops not implemented G5D_MprecDeriv\n");
    exit(-1);
  }
}
template class bfmbase<double>;
template class bfmbase<float>;
