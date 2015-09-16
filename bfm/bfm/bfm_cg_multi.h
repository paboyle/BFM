#ifndef BFM_MULTI_H
#define BFM_MULTI_H

template <class Float>
int bfmbase<Float>::CGNE_prec_MdagM_multi_shift(Fermion_t psi[], 
						Fermion_t src,
						double    mass[],
						double    alpha[],
						int       nshift,
						double mresidual[],
						int single,
						int forcecontinue)
{
  int me = thread_barrier();
 
  if ( this->isBoss() && (!me) ) { 
    this->InverterEnter();
  }
  
  // Per shift fields
  Fermion_t ps [nshift]; // search directions 
  double    bs [nshift];
  double    rsq[nshift];
  double    z[nshift][2];
  int       converged[nshift];

  const int       primary =0;

  //Primary shift fields CG iteration
  double a,b,c,d;
  double cp,bp; //prev
  Fermion_t p ;
  Fermion_t r ;

  struct timeval start,stop;

  // Matrix mult fields
  Fermion_t tmp = threadedAllocFermion(mem_fast); 
  Fermion_t mp  = threadedAllocFermion(mem_fast); 
  Fermion_t mmp = threadedAllocFermion(mem_fast); 

  // Check lightest mass
  for(int s=0;s<nshift;s++){
    if ( mass[s] < mass[primary] ) {
      Error("First shift not lightest - oops\n");
      exit(-1);
    }
  }

  //Allocate per shift fields
  r = threadedAllocFermion(mem_slow);
  p = threadedAllocFermion(mem_slow);
  for(int i=0;i<nshift;i++){
    ps[i] = threadedAllocFermion(mem_slow);
    converged[i]=0;
  }

  // Wire guess to zero
  // Residuals "r" are src
  // First search direction "p" is also src
  cp = norm(src);
  for(int s=0;s<nshift;s++){
    rsq[s] = cp * mresidual[s] * mresidual[s];
    this->ThreadBossDebug("CGNE_prec_multi: shift %d target resid %le %le %le\n",s,rsq[s],cp,mresidual[s]);
    copy(ps[s],src);
  }
  // r and p for primary
  copy(r,src);
  copy(p,src);
  
  //MdagM+m[0]
#if 0
   d= Mprec(p,mp,tmp,DaggerNo,1); 
      Mprec(mp,mmp,tmp,DaggerYes);
#else
//not correct yet
   d= Mprec(p,mp,tmp,DaggerNo,1); 
      Mprec(mp,mmp,tmp,DaggerYes);
   axpy(mmp,p,mmp,mass[0]);
   double rn = norm(p);
   d += rn*mass[0];
#endif

  b = -cp /d;


  // Set up the various shift variables
  int       iz=0;
  z[0][1-iz] = 1.0;
  z[0][iz]   = 1.0;
  bs[0]      = b;
  for(int s=1;s<nshift;s++){
    z[s][1-iz] = 1.0;
    z[s][iz]   = 1.0/( 1.0 - b*(mass[s]-mass[0]));
    bs[s]      = b*z[s][iz]; // Sign relative to Mike - FIXME
  }
  
  // r += b[0] A.p[0]
  // c= norm(r)
  c=axpy_norm(r,mmp,r,b);


  // Linear algebra overhead can be minimised if summing results
  // psi-= b[0]p[0] = -b[0]src

  //  if ( single ) {
  //    double coeff=0;
  //    for(int s=0;s<nshift;s++) coeff+=-bs[s]*alpha[s];
  //    axpby(psi[0],src,src,0.0,coeff);
  //  }  else {
  //  }
  for(int s=0;s<nshift;s++) {
    axpby(psi[s],src,src,0.,-bs[s]*alpha[s]);
  }
  

  // Iteration loop
  int k;

  gettimeofday(&start,NULL);
  for (k=1;k<=max_iter;k++){

    a = c /cp;
    axpy(p,p,r,a);

    // Note to self - direction ps is iterated seperately
    // for each shift. Does not appear to have any scope
    // for avoiding linear algebra in "single" case.
    // 
    // However SAME r is used. Could load "r" and update
    // ALL ps[s]. 2/3 Bandwidth saving
    // New Kernel: Load r, vector of coeffs, vector of pointers ps
    for(int s=0;s<nshift;s++){
      if ( (! converged[s]) || forcecontinue ) { 
        if (s==0){
          axpy(ps[s],ps[s],r,a);
        } else{
	  double as =a *z[s][iz]*bs[s] /(z[s][1-iz]*b);
	  axpby(ps[s],r,ps[s],z[s][iz],as);
        }
      }
    }

    cp=c;
    
   d= Mprec(p,mp,tmp,DaggerNo,1); 
      Mprec(mp,mmp,tmp,DaggerYes);
   axpy(mmp,p,mmp,mass[0]);
   double rn = norm(p);
   d += rn*mass[0];

    bp=b;
    b=-cp/d;
    
    c=axpy_norm(r,mmp,r,b);

    // Toggle the recurrence history
    bs[0] = b;
    iz = 1-iz;
    for(int s=1;s<nshift;s++){
      if((!converged[s])||forcecontinue){
	double z0 = z[s][1-iz];
	double z1 = z[s][iz];
	z[s][iz] = z0*z1*bp
	  / (b*a*(z1-z0) + z1*bp*(1- (mass[s]-mass[0])*b)); 
	bs[s] = b*z[s][iz]/z0; // NB sign  rel to Mike
      }
    }
    
    // Fixme - Mike has a variable "at" that is missing here.
    // Understand this. NB sign of bs rel to Mike
    for(int s=0;s<nshift;s++){
      int ss = s;
      // Scope for optimisation here in case of "single".
      // Could load psi[0] and pull all ps[s] in.
      //      if ( single ) ss=primary;
      // Bandwith saving in single case is Ls * 3 -> 2+Ls, so ~ 3x saving
    // Pipelined CG gain:
    //
    // New Kernel: Load r, vector of coeffs, vector of pointers ps
    // New Kernel: Load psi[0], vector of coeffs, vector of pointers ps
    // If can predict the coefficient bs then we can fuse these and avoid write reread cyce
    //  on ps[s].
    // Before:  3 x npole  + 3 x npole
    // After :  2 x npole (ps[s])        => 3x speed up of multishift CG.

      if( (!converged[s])|| forcecontinue ) { 
	axpy(psi[ss],ps[s],psi[ss],-bs[s]*alpha[s]);
      }
    }

    // Convergence checks
    int all_converged = 1;
    if ( (k%500)==0) {
      this->ThreadBossMessage("CGNE_prec_multi: k=%d c=%g\n",k,c);
    }
    for(int s=0;s<nshift;s++){

      if ( (!converged[s]) || forcecontinue ){

	double css  = c * z[s][iz]* z[s][iz];
	
	if(css<rsq[s]){
	  if ( ! converged[s] )
	    this->ThreadBossLog("CGNE_prec_multi: k=%d Shift %d has converged\n",k,s);
	  converged[s]=1;
	} else {
	  all_converged=0;
	}

      }
    }

    if ( all_converged ){

      gettimeofday(&stop,NULL);
      struct timeval diff;
      timersub(&stop,&start,&diff);
      double t = diff.tv_sec*1.0E6 + diff.tv_usec;

      
      this->ThreadBossMessage("CGNE_prec_multi: k=%d All shifts have converged after %le s\n",k,t*1.0e-6);
      this->ThreadBossDebug("CGNE_prec_multi: k=%d Checking solutions\n",k);
      // Check answers 
      for(int s=0; s < nshift; s++) { 
	Mprec(psi[s],mp,tmp,DaggerNo);
	Mprec(mp,mmp,tmp,DaggerYes);
        axpy(tmp,psi[s],mmp,mass[s]);
	axpy(r,tmp,src,-1);
	double rn = norm(r);
	double cn = norm(src);
	ThreadBossLog("CGNE_prec_multi: shift[%d] true residual %le \n",s,sqrt(rn/cn));
      }

      if ( single ) {
	for(int s=1; s < nshift; s++) { 
	  axpy(psi[0],psi[s],psi[0],1.0);
	}      
      }

      threadedFreeFermion(tmp);
      threadedFreeFermion(p);
      threadedFreeFermion(mp);
      threadedFreeFermion(mmp);
      threadedFreeFermion(r);
      for(int i=0;i<nshift;i++){
	threadedFreeFermion(ps[i]);
      }  

      if ( this->isBoss() && (!me) ) { 
	this->InverterExit();
      }
      return k;
    }
  }

  this->ThreadBossDebug("CGNE_prec_multi: CG not converged after %d iterations\n",k);
  for(int s=0; s < nshift; s++) { 
    Mprec(psi[s],mp,tmp,DaggerNo);
    Mprec(mp,mmp,tmp,DaggerYes);
    axpy(tmp,psi[s],mmp,mass[s]);
    axpy(r,tmp,src,-1);
    double rn = norm(r);
    double cn = norm(src);
    this->ThreadBossDebug("CGNE_prec_multi: shift[%d] true residual %le \n",s,sqrt(rn/cn));
  }

  threadedFreeFermion(tmp);
  threadedFreeFermion(p);
  threadedFreeFermion(mp);
  threadedFreeFermion(mmp);
  threadedFreeFermion(r);
  for(int i=0;i<nshift;i++){
    threadedFreeFermion(ps[i]);
  }  
  if ( this->isBoss() && (!me) ) { 
    this->InverterExit();
  }
}

#endif
