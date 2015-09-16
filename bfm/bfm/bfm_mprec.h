#ifndef BFM_MPREC_H
#define BFM_MPREC_H
#include "bfm.h"
#include <stdio.h>
#include <stdlib.h>

/****************************************************************************
 * 2007
 *  Calls PAB's assembler routines to give an implementation
 *  of the dwf dslash. Aim is to give very high performance
 *  in a reasonably portable/retargettable library.
 ****************************************************************************
 * 2011
 *  Revamp for multiple actions.
 *  Probably want to move to a class heirarchy
 *
 *  bfmarg  <--- multiple initialiser functions?
 *     |
 *  bfmbase <--- merge these into single "bfm" class.
 *     |
 *  bfmdop
 *     |
 *  _____________________________________________......______
 *  |          |          |             |                    |
 *bfmDWF   bfmWilson   bfmWilsonTM    bfmDWFrb4d           bfmClover etc...
 * 
 * The routine in this file should be virtual
 ******************************************************************
 * Conventions:
 *
 *-----------------------------------------------------------
 * Munprec  =  (Nd+M)               - 1/2 D
 *          =   (4+m)               - 1/2 Dw   [ Wilson ]
 *          =   (4+m)  + i tm g5    - 1/2 Dw   [ WilsonTM ]
 *          =   (5-M5)              - 1/2 Ddwf [ DWF ]
 *-----------------------------------------------------------
 * Red black Schur decomposition
 *
 *  M = (Mee Meo) =  (1             0 )   (Mee   0               )  (1 Mee^{-1} Meo)
 *      (Moe Moo)    (Moe Mee^-1    1 )   (0   Moo-Moe Mee^-1 Meo)  (0   1         )
 *                =         L                     D                     U
 *
 * Mprec = Doo = Moo-MoeMee^{-1}Meo
 *
 * L^-1 = (1              0 )
 *        (-MoeMee^{-1}   1 )   
 * L^{dag} = ( 1       Mee^{-dag} Moe^{dag} )
 *           ( 0       1                    )
 * L^{-d}  = ( 1      -Mee^{-dag} Moe^{dag} )
 *           ( 0       1                    )
 *
 * U^-1 = (1   -Mee^{-1} Meo)
 *        (0    1           )
 * U^{dag} = ( 1                 0)
 *           (Meo^dag Mee^{-dag} 1)
 * U^{-dag} = (  1                 0)
 *            (-Meo^dag Mee^{-dag} 1)
 * 
 * To solve M psi = eta
 *
 *Odd
 * i)   (D_oo)^{\dag} D_oo psi_o = (D_oo)^\dag L^{-1}  eta_o
 *                        eta_o' = D_oo (eta_o - Moe Mee^{-1} eta_e)
 *Even
 * ii)  Mee psi_e + Meo psi_o = src_e
 *
 *   => sol_e = M_ee^-1 * ( src_e - Meo sol_o )...
 *
 *
 *******************************************************************
 * Case by case:
 *
 * DWF : Mee = Moo = (5-M5)
 *       Meo = -1/2 Ddwf_eo (5d hopping term)
 *       Moe = -1/2 Ddwf_oe (5d hopping term)
 *
 * DWFrb4d : Mee = Moo = (5-M5) - 1/2 Dperp(m_f)
 *           Meo = -1/2 Dw_eo (4d hopping term)
 *           Moe = -1/2 Dw_oe (4d hopping term)
 *
 * Wilson : Mee = Moo = (4+m)
 *          Meo = -1/2 Ddwf_eo (5d hopping term)
 *          Moe = -1/2 Ddwf_oe (5d hopping term)
 *
 * WilsonTM : Mee = Moo = 4+m + i tm g5
 *            Meeinv    = (4+m - i tm g5)/ ( (4+m)^2 + tm^2 ) 
 *                      =    (4+m)
 *                        -----------------   .  ( 1 -i tm/(4+m) g5 )
 *                        ( (4+m)^2 + tm^2 )
 *
 *
 */

template <class Float>
void bfmbase<Float>::MooeeInv(Fermion_t psi, 
			       Fermion_t chi, 
			      int dag, int cb)
{
  if ( IsGeneralisedFiveDim() ) {
    this->G5D_MooeeInv(psi,chi,dag);
    return;
  }


  if ( this->solver == DWF) {
    if ( ! this->precon_5d ) {
      printf("logic bomb -- dwf but not precon5d\n");
      exit(-1);
    }
    MooeeInv5dprec(psi,chi,dag); // In practice unused
  } else if ( this->solver == DWFrb4d ) {
    if ( this->precon_5d ) {
      printf("logic bomb -- dwfrb4d but precon5d\n");
      exit(-1);
    }
    MooeeInv4dprec(psi,chi,dag);
  } else if ( this->solver == WilsonFermion ) {
    double coeff = 1.0/ ( 4.0 + this->mass );
    double zero  = 0.0;
    axpby(chi,psi,psi,coeff,zero);
  } else if (this->solver == CloverFermion ) {
    applyClover (psi,chi,cb,1);
  } else if ( this->solver == WilsonTM ) {
    double m    = this->mass;
    double tm   = this->twistedmass;
    double mtil = 4.0+this->mass;

    double sq = mtil*mtil + tm*tm;

    double a = mtil/sq;
    double b = -tm /sq;
    if(dag) b=-b;
    axpibg5x(chi,psi,a,b);
  }
}
template <class Float>
void bfmbase<Float>::Mooee(Fermion_t psi, 
			   Fermion_t chi, 
			   int dag,int cb)
{

  if ( IsGeneralisedFiveDim() ) {
    this->G5D_Mooee(psi,chi,dag);
    return;
  }

  if ( this->solver == DWF) {
    Mooee5dprec(psi,chi,dag);
  } else if ( this->solver == DWFrb4d ) {
    Mooee4dprec(psi,chi,dag);
  } else if ( this->solver == WilsonFermion ) {
    double coeff = ( 4.0 + this->mass );
    double zero  = 0.0;
    axpby(chi,psi,psi,coeff,zero);
  } else if (this->solver == CloverFermion ) {
    applyClover (psi,chi,cb,0);
  } else if ( this->solver == WilsonTM ) {
    double a = 4.0+this->mass;
    double b = this->twistedmass;
    if(dag) b=-b;
    axpibg5x(chi,psi,a,b);
  }
}
template <class Float>
void bfmbase<Float>::Meo(Fermion_t psi,
			 Fermion_t chi, 
			 int result_cb,
			 int dag)
{
  if ( IsGeneralisedFiveDim() ) {
    this->G5D_Meo(psi,chi,result_cb,dag);
    return;
  }
  double m=-0.5;
  dslash(psi,chi,result_cb,dag);
  scale (chi,m); // Low performance alert. Never use this in CG.
}


template <class Float>
void bfmbase<Float>::MprecDeriv(Fermion_t X, 
				Fermion_t Y, 
				Matrix_t  Force[2],
				int dag,
				double coeff)
{
  int me = this->thread_barrier();

  if ( CGdiagonalMee ) { 
    printf("MprecDeriv CGdiagonalMee not implemented\n");
    exit(-1);
  }

  Fermion_t tmp1=threadedAllocFermion();
  Fermion_t tmp2=threadedAllocFermion();

  //  X^dag Der_oe MeeInv Meo Y
  // Use Mooee as nontrivial but gauge field indept
  Meo      (Y,tmp1,Even,dag);      // odd->even -- implicit -0.5 factor to be applied
  MooeeInv (tmp1,tmp2,dag,Even);   // even->even 
  axpby(tmp2,tmp2,tmp2,-coeff,0.0);  // Flip the sign
  MeoDeriv (X,tmp2,Force[Odd],Odd,dag);
  
  //  Accumulate X^dag M_oe MeeInv Der_eo Y
  int mdag = 1-dag;
  Meo     (X,tmp1,Even,mdag);    // even->odd -- implicit -0.5 factor to be applied
  MooeeInv(tmp1,tmp2,mdag,Even); // even->even 
  axpby(tmp2,tmp2,tmp2,-coeff,0.0);// Flip the sign
  MeoDeriv(tmp2,Y,Force[Even],Even,dag);

  threadedFreeFermion(tmp1);
  threadedFreeFermion(tmp2);
}


template <class Float>
void bfmbase<Float>::MeoDeriv(Fermion_t X, 
			      Fermion_t Y, 
			      Matrix_t  Force,
			      int cb,
			      int dag)
{
  if ( IsGeneralisedFiveDim() ) {
    this->G5D_MeoDeriv(X,Y,Force,cb,dag);
    return;
  }
  
  Fermion_t tmp=threadedAllocFermion();
  double m=-0.5;
  axpby(tmp,Y,Y,m,0.0);
  DslashDeriv(X,tmp,Force,cb,dag);
  threadedFreeFermion(tmp);
  return;
}


template <class Float>
double bfmbase<Float>::MprecTilde(Fermion_t psi, 
				  Fermion_t chi, 
				  Fermion_t tmp, 
				  int dag,
				  int donrm)
{
  int cb=this->rb_precondition_cb;
  //
  // Left multiple system of equations by diag(MeeInv MooInv)
  // Question: how does this work for HMC?
  //
  // Mee Meo              1 MeeInv Meo
  // Moe Moo         ->   MooInvMoe 1
  //
  // eta' = MeeInv eta_e
  //        MooInv eta_o
  //
  //
 
  /*
   *  Mt= (1 MeeInv Meo) =  (1             0 )   (1   0               )         (1 Mee^{-1} Meo)
   *      (MooInv Moe 1)    (MooInv Moe    1 )   (0   1-MooInv Moe Mee^-1 Meo)  (0   1         )
   *                =         L                     D                     U
   *
   * MprecTilde = Doo = 1-Moo^{-1}MoeMee^{-1}Meo
   *
   *
   *
   * L^-1 = (1              0 )
   *        (-MooInv Moe      1 )   
   * L^{dag} = ( 1       Moe^{dag} MooInv^{dag} )
   *           ( 0       1                      )
   * L^{-d}  = ( 1      -Moe^{dag} MooInv^{dag} )
   *           ( 0       1                      )
   *
   * U^-1 = (1   -Mee^{-1} Meo)
   *        (0    1           )
   * U^{dag} = ( 1                 0)
   *           (Meo^dag Mee^{-dag} 1)
   * U^{-dag} = (  1                 0)
   *            (-Meo^dag Mee^{-dag} 1)
   * 
   * To solve Mt psi = eta'
   * <=>
   *    Ddag L^-1 L D U psi = Ddag L^-1 eta'
   * <=>
   *    (1                   )(psi_e+Mee^{-1} Meo psi_o)       = [ eta_e'                                   ]
   *    (     D_oo^{dag} D_oo)(psi_o                   )         [D_oo^{dag} ( eta_o' - Moo^{-1}Moe eta_e') ]
   *
   *                                                           = [ Mee^{-1} eta_e                                             ]
   *                                                             [ D_oo^{dag} ( Moo^{-1} eta_o - Moo^{-1} Moe Mee^{-1} eta_e) ]
   *Odd
   * i)   (D_oo)^{\dag} D_oo psi_o = (D_oo)^\dag ( Moo^{-1} eta_o - Moo^{-1} Moe Mee^{-1} eta_e) 
   * (WAS)(D_oo)^{\dag} D_oo psi_o = (D_oo)^\dag (eta_o - Moe Mee^{-1} eta_e)
   *                       
   *Even
   * ii)  psi_e + Mee^{-1} Meo psi_o = src_e
   *
   *   => sol_e =  Mee^{-1} eta_e - M_ee^-1 Meo sol_o ...
   */
  double nrm;
  double m1=-1.0;

  int ocb=1-cb; // other checker board
  if ( dag ) {
    //   * MprecTilde = Doo = 1-(Moo^{-1}MoeMee^{-1}Meo)^dag
    MooeeInv (psi,tmp,dag,cb);    // odd odd 
    Meo (tmp,chi,ocb,dag);        // odd->even 
    MooeeInv(chi,tmp,dag,ocb);    // even->even 
    Meo (tmp,chi,cb,dag);        // even->odd 
    nrm = axpy_norm (chi,chi,psi,m1);     // Linear combination
  } else { 
    //   * MprecTilde = Doo = 1-Moo^{-1}MoeMee^{-1}Meo
    Meo (psi,tmp,ocb,dag);   // odd->even 
    MooeeInv(tmp,chi,dag,ocb);    // even->even 
    Meo (chi,tmp,cb,dag);    // even->odd 
    MooeeInv (tmp,chi,dag,cb);   // odd odd 
    nrm = axpy_norm (chi,chi,psi,m1);     // Linear combination
  }
  return nrm;
}

template <class Float>
double bfmbase<Float>::CompactMprec(Fermion_t compact_psi, 
				    Fermion_t compact_chi, 
				    Fermion_t psi,
				    Fermion_t chi,
				    Fermion_t tmp, 
				    int dag,
				    int donrm)
{
	copy(psi, compact_psi);
  
	double result = Mprec(psi, chi, tmp, dag, donrm);
  
	copy(compact_chi, chi);	

	return result;
}
template <class Float>
double bfmbase<Float>::Mprec(Fermion_t psi, 
			     Fermion_t chi, 
			     Fermion_t tmp, 
			     int dag,
			     int donrm)
{
  int cb=this->rb_precondition_cb;

  if ( CGdiagonalMee ) { 
     return MprecTilde(psi,chi,tmp,dag,donrm);
  }

  if ( IsGeneralisedFiveDim() ) {
    return this->G5D_Mprec(psi,chi,tmp,dag,donrm);
  }

  double nrm = 0.0;

  // Use Mooee as nontrivial
  // Wilson, DWF are faster to apply as above as Mooee is a scalar mult
  int ocb=1-cb;// other checkerboard

  dslash  (psi,tmp,ocb,dag);   // odd->even -- implicit -0.5 factor to be applied
  MooeeInv(tmp,chi,dag,ocb);   // even->even 
  dslash  (chi,tmp,cb,dag);    // even->odd -- implicit -0.5 factor to be applied
  Mooee   (psi,chi,dag,cb);    // odd odd 
  double m1=-0.25;
  nrm = axpy_norm (chi,tmp,chi,m1);     // Linear combination

  return nrm;

}

template <class Float>
void bfmbase<Float>::CompactMunprec(Fermion_t compact_psi[2], 
			     Fermion_t compact_chi[2],
					 Fermion_t psi[2],
					 Fermion_t chi[2],
			     Fermion_t tmp,
   		             int dag)
{
	copy(psi[0], compact_psi[0]);
	copy(psi[1], compact_psi[1]);

	Munprec(psi, chi, tmp, dag);

	copy(compact_chi[0], chi[0]);
	copy(compact_chi[1], chi[1]);
}
template <class Float>
void bfmbase<Float>::Munprec(Fermion_t psi[2], 
			     Fermion_t chi[2], 
			     Fermion_t tmp,
   		             int dag) 
{

  if ( IsGeneralisedFiveDim() ) {
    this->G5D_Munprec(psi,chi,dag);
    return;
  }

  double m=-0.5;

  Mooee(psi[0],chi[0],dag,0);
  Mooee(psi[1],chi[1],dag,1);

  dslash(psi[0],tmp,1,dag);
  axpy  (chi[1],tmp,chi[1],m);     // Linear combination

  dslash(psi[1],tmp,0,dag);
  axpy  (chi[0],tmp,chi[0],m);     // Linear combination

}



// CPS style .. these are site diagonal and trivial
template <class Float>
void bfmbase<Float>::Mooee5dprec(Fermion_t psi_t, 
				  Fermion_t chi_t, 
				  int dag)
{
  // This is just a reference implementation
  Float *psi = (Float *) psi_t;
  Float *chi = (Float *) chi_t;
  int nspinco = 12;
  Float invtwokappa = 5.-M5;

  Float *xx;
  Float *rr;

  int me,thrlen,throff;
  thread_work(simd_cbvol,me,thrlen,throff);
  for(int site=throff;site<throff+thrlen;site++){
   for(int s=0;s<cbLs;s++) {

    xx = &psi[s*nspinco*nsimd*2   + site*nspinco*nsimd*2*cbLs];
    rr = &chi[s*nspinco*nsimd*2   + site*nspinco*nsimd*2*cbLs];

    for(int spinco=0;spinco<nspinco;spinco++){
    for(int simd_el=0;simd_el<nsimd;simd_el++){
    for(int reim=0;reim<2;reim++){
      int idx = reim + simd_el*2 + spinco*nsimd*2;
      rr[idx] = invtwokappa*xx[idx];
    }}}
   }
  }  
  thread_barrier();
  return;
}
template <class Float>
void bfmbase<Float>::MooeeInv5dprec(Fermion_t psi_t, 
				     Fermion_t chi_t, 
				     int dag)
{

  // This is just a reference implementation
  Float *psi = (Float *) psi_t;
  Float *chi = (Float *) chi_t;
  int nspinco = 12;
  Float invtwokappa = 5.-M5;
  Float twokappa = 1.0/invtwokappa;

  Float *xx;
  Float *rr;

  int me,thrlen,throff;
  thread_work(simd_cbvol,me,thrlen,throff);
  for(int site=throff;site<throff+thrlen;site++){

   for(int s=0;s<cbLs;s++) {

    xx = &psi[s*nspinco*nsimd*2   + site*nspinco*nsimd*2*cbLs];
    rr = &chi[s*nspinco*nsimd*2   + site*nspinco*nsimd*2*cbLs];

    for(int spinco=0;spinco<nspinco;spinco++){
    for(int simd_el=0;simd_el<nsimd;simd_el++){
    for(int reim=0;reim<2;reim++){
      int idx = reim + simd_el*2 + spinco*nsimd*2;
      rr[idx] = twokappa*xx[idx];
    }}}
   }
  }  
  thread_barrier();
  return;
}



// CHROMA style oo/ee part of preconditioned matrix
template <class Float>
void bfmbase<Float>::Mooee4dprec(Fermion_t psi_t, 
				  Fermion_t chi_t, 
				  int dag)
{
  // There is not yet any optimised version of this code
  // Clearly this code is DWFrb4d specific.
  Float *psi = (Float *) psi_t;
  Float *chi = (Float *) chi_t;
  int nspinco = 12;

  int me,thrlen,throff;
  thread_work(simd_cbvol,me,thrlen,throff);
  for(int site=throff;site<throff+thrlen;site++){
  for(int s=0;s<Ls;s++) {

    Float invtwokappa = 5.-M5;

    Float  c,  c_p,  c_m ;
    Float *p, *p_p, *p_m, *ch ;
    integer    s_p,  s_m ;

    c_m = c_p = -1.0;
    c   = invtwokappa;

    if ( s== 0   ) c_m = mass; // mf boundary wrap
    if ( s== Ls-1 ) c_p = mass;

    s_p = (s+1)%Ls;
    s_m = (s+Ls-1)%Ls;

    // Real/Imaginary innermost
    // Can I use the virtual bagel_idx function????

    p   = &psi[bagel_idx5d(site,s,nspinco)];
    p_m = &psi[bagel_idx5d(site,s_m,nspinco)];
    p_p = &psi[bagel_idx5d(site,s_p,nspinco)];
    ch  = &chi[bagel_idx5d(site,s,nspinco)];
    /*
     * For bluegene this can be made very efficient
     * by each time around the s-loop loading 
     * L 6-regs psi   (+6 regs reusing last psi_p)
     * L 6-regs psi_p, 
     *   psi_m -> old psi regs 
     *  6/12-regs chi are assigned and stored.
     * Footprint is 24-30 registers.
     *
     * Need to unroll by two to avoid register copies.
     * Modulo scheduling would have been nice...
     * But unrolling & icaches are the reason modulo scheduling 
     * is always avoidable I guess....
     */

    for(int spinco=0;spinco<nspinco/2;spinco++){
    for(int simd_el=0;simd_el<nsimd;simd_el++){
    for(int reim=0;reim<2;reim++){

      int idx = reim + simd_el*2 + spinco*nsimd*2;
      int hdx = reim + simd_el*2 + (spinco+6)*nsimd*2;

      if ( dag ) {
        ch[idx]= c*p[idx]+c_p*p_p[idx];
        ch[hdx]= c*p[hdx]+c_m*p_m[hdx];
      }else{
        ch[idx]= c*p[idx]+c_m*p_m[idx];
        ch[hdx]= c*p[hdx]+c_p*p_p[hdx];
      }
    }}}

    
  }
  }

  thread_barrier();

}

template <class Float>
void bfmbase<Float>::MooeeInv4dprec(Fermion_t psi_t, 
				    Fermion_t chi_t, 
				    int dag)
{

  Float *psi = (Float *) psi_t;
  Float *chi = (Float *) chi_t;

  int nspinco = 12;
  int idx;
  int hdx;
  int adx;

  Float invtwokappa = 5.-M5;
  Float twokappa = 1.0/(5.-M5);
  Float kappa    = 0.5*twokappa;
  Float invDfactor = 1.0/(1.0+mass/pow(invtwokappa,Ls));

  Float *p,*ch,*chm,*chp;

  int me,thrlen,throff;
  thread_work(simd_cbvol,me,thrlen,throff);

  for(int site=throff;site<throff+thrlen;site++){

   Float ChiAdjust[12*simd()];
   bzero((void *)ChiAdjust,12*simd()*sizeof(Float));

   Float fact=mass*twokappa;

   for(int s=0;s<Ls;s++) {

    p   = &psi[bagel_idx5d(site,s,nspinco)];
    ch  = &chi[bagel_idx5d(site,s,nspinco)];

    for(int spinco=0;spinco<nspinco/2;spinco++){
    for(int simd_el=0;simd_el<simd();simd_el++){
    for(int reim=0;reim<2;reim++){

      if ( dag ) {
        hdx = reim + simd_el*2 + spinco*simd()*2;
        idx = reim + simd_el*2 + (spinco+6)*simd()*2;
	adx = hdx;
      } else {  
        idx = reim + simd_el*2 + spinco*simd()*2;
        hdx = reim + simd_el*2 + (spinco+6)*simd()*2;
        adx = idx;
      }

      // Copy and scale by two kappa
      ch[idx] = p[idx]*twokappa;
      ch[hdx] = p[hdx]*twokappa;

      // inverse of Lm
      if ( s<(Ls-1) ) {
	ChiAdjust[adx] -= fact*ch[hdx]; // P- takes lower components
      } else { 
	ch[hdx] = ch[hdx]+ChiAdjust[adx];
      }

      //inverse of L -- forward elimination
      if ( s>0 ) { 
        chm = &chi[bagel_idx5d(site,s-1,nspinco)];
	ch[idx]+=twokappa*chm[idx];
      }

      //Inverse of D
      if ( s==(Ls-1) ) { 
	ch[idx] *= invDfactor;
	ch[hdx] *= invDfactor;
      }

    }}}    

    fact*=twokappa;

   }

   fact=mass*pow(twokappa,Ls-1);
   // mass*pow(twokappa,s+1)
   for(int s=Ls-2;s>=0;s--) {
     Float *chl;
     ch  = &chi[bagel_idx5d(site,s,nspinco)];
     chp = &chi[bagel_idx5d(site,s+1,nspinco)];
     chl = &chi[bagel_idx5d(site,Ls-1,nspinco)];

     for(int spinco=0;spinco<nspinco/2;spinco++){
     for(int simd_el=0;simd_el<simd();simd_el++){
     for(int reim=0;reim<2;reim++){

       if ( dag ) {
	 hdx = reim + simd_el*2 + spinco*simd()*2;
	 idx = reim + simd_el*2 + (spinco+6)*simd()*2;
       } else {  
	 idx = reim + simd_el*2 + spinco*simd()*2;
	 hdx = reim + simd_el*2 + (spinco+6)*simd()*2;
       }

       ch[hdx] = ch[hdx]+twokappa*chp[hdx];
       ch[idx] -= fact*chl[idx];
     }}}
     fact /= twokappa;     
   }
  }
  thread_barrier();

}

/***********************************************
		Taku PV 
************************************************/


template <class Float>
double bfmbase<Float>::MprecTW(Fermion_t psi, 
			     Fermion_t chi, 
			     Fermion_t tmp, 
			     int dag, int n, int L5,
			     int donrm)
{

  double nrm = 0.0;
  
  if ( Ls == 1 ) {
    // This is the Chroma style preconditioning
    // Use Mooee as nontrivial
    // Wilson, DWF are faster to apply as below    
    dslash  (psi,tmp,Even,dag);   // odd->even -- implicit -0.5 factor to be applied    
    MooeeInv5dprec_TW(tmp,chi,n,L5,dag);     // even->even * 1/(5-M5)
    dslash  (chi,tmp,Odd,dag);   // even->odd -- implicit -0.5 factor to be applied
    Mooee5dprec_TW   (psi,chi,n,L5,dag);     // odd odd * (5-M5)

    double m1=-0.25;
    nrm = axpy_norm (chi,tmp,chi,m1);     // Linear combination

  } else {


exit(1);

  }
  return nrm;
}

template <class Float>
void bfmbase<Float>::MunprecTW(Fermion_t psi[2], 
			     Fermion_t chi[2], 
			     Fermion_t tmp,
			     int n, int L5,
   		             int dag) 
{			     	       
  double m=-0.5;
  Mooee5dprec_TW(psi[0],chi[0],n, L5,dag);
  Mooee5dprec_TW(psi[1],chi[1],n, L5, dag);

  dslash(psi[0],tmp,1,dag);
  axpy  (chi[1],tmp,chi[1],m);     // Linear combination
  dslash(psi[1],tmp,0,dag);
  axpy  (chi[0],tmp,chi[0],m);     // Linear combination  
  
}


template <class Float>
void bfmbase<Float>::Mooee5dprec_TW(Fermion_t psi_t, 
				    Fermion_t chi_t, int n, int L5, 
				    int dag)
{

  Float *psi = (Float *) psi_t;
  Float *chi = (Float *) chi_t;
  int nspinco = 12;
  
  int dg = 1;
  if(dag) dg = -1;
  
  Float invtwokappa = 5.-M5 - cos( 2.0 * M_PI * (n + 0.5) / L5 );
  Float sn = dg*sin( 2.0 * M_PI * (n + 0.5) / L5 );
  Float *xx;
  Float *rr;

  int me,thrlen,throff, sign;
  thread_work(simd_cbvol,me,thrlen,throff);
  for(int site=throff;site<throff+thrlen;site++){
   for(int s=0;s<cbLs;s++) {

    xx = &psi[s*nspinco*nsimd*2   + site*nspinco*nsimd*2*cbLs];
    rr = &chi[s*nspinco*nsimd*2   + site*nspinco*nsimd*2*cbLs];

    for(int spinco=0;spinco<nspinco;spinco++){
    sign = 1;
    if(spinco > 5) sign = -1;
    for(int simd_el=0;simd_el<nsimd;simd_el++){
    
    //real part reim = 0, imaginary part reim = 1;
    //Lower chirality co > 5
    
      int idx_re = 0 + simd_el*2 + spinco*nsimd*2;
      int idx_im = 1 + idx_re;// + simd_el*2 + spinco*nsimd*2;
      Float tmpr, tmpi;
      tmpr = xx[idx_re]; tmpi = xx[idx_im];
      rr[idx_re] = invtwokappa*tmpr - sign*sn*tmpi;
      rr[idx_im] = invtwokappa*tmpi + sign*sn*tmpr;
    
    }}
   }
  }  
  thread_barrier();
  return;
  
}


template <class Float>
void bfmbase<Float>::MooeeInv5dprec_TW(Fermion_t psi_t, 
				     Fermion_t chi_t, int n, int L5, 
				     int dag)
{

  Float *psi = (Float *) psi_t;
  Float *chi = (Float *) chi_t;
  int nspinco = 12;
  //Float invtwokappa = 5.-M5;
  //Float twokappa = 1.0/invtwokappa;

  int dg = -1;
  if(dag) dg = 1;
  
  Float invtwokappa = 5.-M5 - cos( 2.0 * M_PI * (n + 0.5) / L5 );
  Float sn = sin( 2.0 * M_PI * (n + 0.5) / L5 );
  
  Float invx = invtwokappa / (invtwokappa*invtwokappa + sn*sn);
  Float invy = dg*sn*invx / invtwokappa;
  
  Float *xx;
  Float *rr;

  int me,thrlen,throff,sign;
  thread_work(simd_cbvol,me,thrlen,throff);
  for(int site=throff;site<throff+thrlen;site++){

   for(int s=0;s<cbLs;s++) {

    xx = &psi[s*nspinco*nsimd*2   + site*nspinco*nsimd*2*cbLs];
    rr = &chi[s*nspinco*nsimd*2   + site*nspinco*nsimd*2*cbLs];

    for(int spinco=0;spinco<nspinco;spinco++){
    sign = 1;
    if(spinco > 5) sign = -1;
    for(int simd_el=0;simd_el<nsimd;simd_el++){
    
    //real part reim = 0, imaginary part reim = 1;
    //Lower chirality co > 5
    
      int idx_re = 0 + simd_el*2 + spinco*nsimd*2;
      int idx_im = 1 + simd_el*2 + spinco*nsimd*2;
      Float tmpr, tmpi;
      tmpr = xx[idx_re]; tmpi = xx[idx_im];
      rr[idx_re] = invx*tmpr - sign*invy*tmpi;
      rr[idx_im] = invx*tmpi + sign*invy*tmpr;
    
    
    }}
   }
  }  
  thread_barrier();
  


  return;
}

/***********************************************/



#endif 
