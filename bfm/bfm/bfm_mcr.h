#ifndef BFM_MCR_H
#define BFM_MCR_H

/*  2007
 *  Calls PAB's assembler routines to give an implementation
 *  of the dwf dslash. Aim is to give very high performance
 *  in a reasonably portable/retargettable library.
 */
#include "bfm.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <complex>
template <class Float>
int bfmbase<Float>::MCR_prec(Fermion_t psi, Fermion_t src)
{
  return MCR_prec(psi,src,(Fermion_t*)NULL,0);
}

template <class Float>
int bfmbase<Float>::MCR_prec(Fermion_t psi, Fermion_t src,
			     Fermion_t *residuals,int nresid) // Internal
{

  double f;
  double cp,c,a,d,b,ci;
  std::complex<double> num;
  int me = thread_barrier();

  int stored_resid=0;

  if ( this->isBoss() && (!me) ) { 
    this->InverterEnter();
  }

  Fermion_t p   = threadedAllocFermion(mem_fast); 
  Fermion_t tmp = threadedAllocFermion(mem_fast); 
  Fermion_t mp  = threadedAllocFermion(mem_fast); 
  Fermion_t r   = threadedAllocFermion(mem_fast); 
  Fermion_t Ap  = threadedAllocFermion(mem_fast); 
  Fermion_t Ar  = threadedAllocFermion(mem_fast); 

  //Initial residual computation & set up
  double guess = norm(psi);
  Mprec(psi,mp,tmp,DaggerNo);
  Mprec(mp,Ar,tmp,DaggerYes);

  axpy (r, Ar, src,-1.0);
  axpy (p, Ar, src,-1.0);

  a =norm(p);
  cp=norm(r);

  Float ssq =  norm(src);
  Float rsq =  residual* residual*ssq;

  //Check if guess is really REALLY good :)
  if ( cp <= rsq ) {
    ThreadBossMessage("MCR_prec k=0 converged - suspiciously nice guess %le %le\n",cp,rsq);
    threadedFreeFermion(tmp);
    threadedFreeFermion(p);
    threadedFreeFermion(mp);
    threadedFreeFermion(r);
    threadedFreeFermion(Ap);
    threadedFreeFermion(Ar);
    if ( this->isBoss() && (!me) ) { 
      this->InverterExit();
    }
    return 0;
  }
  ThreadBossMessage("MCR_prec k=0 residual %le rsq %le\n",cp,rsq);

  if ( this->watchfile ) {
    ThreadBossMessage("MCR_prec watching file \"%s\"\n",this->watchfile);
  }
  struct timeval start,stop;
  gettimeofday(&start,NULL);


  c = Mprec(p,mp,tmp,DaggerNo);
  d = Mprec(mp,Ap,tmp,DaggerYes);  //Ap
  axpy(Ar,Ap,Ap,0);	       //Ar

  //c = real( dot(r,Ar) );
  int k;
  for (k=1;k<=max_iter;k++){

    this->iter=k;

  //c = real( dot(r,Ar) );		//c = rAr
  a = c/d;

  axpy(psi,p,psi,a);		//x = x + ap
  cp = axpy_norm(r,Ap,r,-a);		//r = r - aAp

    if ( ((k%10 == 0) && (verbose!=0)) || (verbose > 10)){
      ThreadBossMessage("MCR_prec: k= %d r^2= %le %le\n",k,cp,sqrt(cp/ssq));
    }

    // Stopping condition
    if ( cp <= rsq ) { 

      if ( stored_resid >= nresid ) { 

	gettimeofday(&stop,NULL);
	struct timeval diff;
	timersub(&stop,&start,&diff);
	
	ThreadBossMessage("MCR_prec converged in %d iterations\n",k);
	ThreadBossMessage("MCR_prec converged in %d.%6.6d seconds\n",diff.tv_sec,diff.tv_usec);


	double flops = mprecFlops()*2.0 + 2.0*axpyNormFlops() + axpyFlops()*2.0;
	flops = flops * k;

	double t = diff.tv_sec*1.0E6 + diff.tv_usec;
	ThreadBossMessage("MCR_prec: %d mprec flops/site\n",mprecFlopsPerSite());
	ThreadBossMessage("MCR_prec: %le flops\n",flops);
	ThreadBossMessage("MCR_prec: %le mflops per node\n",flops/t);

	Mprec(psi,mp,tmp,0);
	Mprec(mp,Ap,tmp,1); 
	axpy(tmp,src,Ap,-1.0);
	double true_residual = sqrt(norm(tmp)/norm(src));
	ThreadBossMessage("MCR_prec: true residual is %le \n",true_residual);

	threadedFreeFermion(tmp);
	threadedFreeFermion(p);
	threadedFreeFermion(mp);
	threadedFreeFermion(r);
	threadedFreeFermion(Ap);
	threadedFreeFermion(Ar);
	if ( this->isBoss() && (!me) ) { 
	  this->InverterExit();
	}
	return k;
	
      } else { 
	
	//	Mprec(psi,mp,tmp,0);
	//	Mprec(mp,Ap,tmp,1); 
	//	axpy(residuals[stored_resid],src,Ap,-1.0);
	//
	axpy(residuals[stored_resid],r,r,0.0);
	ThreadBossMessage("MCR_prec: storing residual %d/%d\n", stored_resid,nresid);
	//	if ( (k % 20) == 0 ) { 
	  stored_resid++;
	  //	}
      }
    }
    ci = d; //(Ap,Ap)
   //d = real( dot(r,Ar) );
    d = Mprec(r,mp,tmp,DaggerNo); //(r,Ar)
        Mprec(mp,Ar,tmp,DaggerYes);
    b = d/c;
    c = d;
    axpy(p,p,r,b);
    d = axpy_norm(Ap,Ap,Ar,b);

  }
  
  {
    Mprec(psi,mp,tmp,0);
    Mprec(mp,Ap,tmp,1); 
    axpy(tmp,src,Ap,-1.0);
    double true_residual = sqrt(norm(tmp)/norm(src));
    ThreadBossMessage("MCR_prec: true residual is %le \n",true_residual);

    ThreadBossMessage("MCR_prec: CG not converged after %d iterations\n",k);
  }

  threadedFreeFermion(tmp);
  threadedFreeFermion(p);
  threadedFreeFermion(mp);
  threadedFreeFermion(r);
  threadedFreeFermion(Ap);
  threadedFreeFermion(Ar);
  if ( this->isBoss() && (!me) ) { 
    this->InverterExit();
  }
  return k;

}

template <class Float>
int bfmbase<Float>::MCR_prec_G5R(Fermion_t psi, Fermion_t src)
{

  double f;
  double cp,c,a,d,b,ci,dold, aold,g,del;
  std::complex<double> num;
  int me = thread_barrier();

  if ( this->isBoss() && (!me) ) { 
    this->InverterEnter();
  }

  Fermion_t p   = threadedAllocFermion(mem_fast); 
  Fermion_t pold   = threadedAllocFermion(mem_fast); 
  Fermion_t Apold   = threadedAllocFermion(mem_fast); 
  Fermion_t tmp = threadedAllocFermion(mem_fast); 
  Fermion_t tmp1 = threadedAllocFermion(mem_fast); 
  Fermion_t mp  = threadedAllocFermion(mem_fast); 
  Fermion_t r   = threadedAllocFermion(mem_fast); 
  Fermion_t Ap  = threadedAllocFermion(mem_fast); 
  Fermion_t Ar  = threadedAllocFermion(mem_fast); 

  //Initial residual computation & set up
  double guess = norm(psi);
  Mprec(psi,mp,tmp,DaggerNo);
  //Mprec(mp,Ar,tmp,DaggerYes);
  G5R(mp,Ar);

  axpy (r, Ar, src,-1.0);
  axpy (p, Ar, src,-1.0);
  axpby (pold, pold,pold,0,0);
  axpby (Apold,Apold,Apold,0,0);
  aold = 0;
  dold = 0;
  d = 0;
  a = 0;
  cp= norm(r);

  Float ssq =  norm(src);
  Float rsq =  residual* residual*ssq;


  //Check if guess is really REALLY good :)
  if ( cp <= rsq ) {
    ThreadBossMessage("MCR_prec k=0 converged - suspiciously nice guess %le %le\n",cp,rsq);
    threadedFreeFermion(tmp);
    threadedFreeFermion(tmp1);
    threadedFreeFermion(p);
    threadedFreeFermion(pold);
    threadedFreeFermion(mp);
    threadedFreeFermion(r);
    threadedFreeFermion(Ap);
    threadedFreeFermion(Apold);
    threadedFreeFermion(Ar);
    if ( this->isBoss() && (!me) ) { 
      this->InverterExit();
    }
    return 0;
  }
  ThreadBossMessage("MCR_prec k=0 residual %le rsq %le\n",cp,rsq);

  if ( this->watchfile ) {
    ThreadBossMessage("MCR_prec watching file \"%s\"\n",this->watchfile);
  }
  struct timeval start,stop;
  gettimeofday(&start,NULL);


  Mprec(p,mp,tmp,DaggerNo);
  G5R(mp,Ap);

  d = axpy_norm(Ar,Ap,Ap,0);	       //Ar
  a = 0;

  double small = 1e-2;

  for (int k=1;k<=max_iter;k++){

  this->iter=k;

  c = real( dot(r,Ap) );
  aold = a;
  a = c/d;

  axpy(psi,p,psi,a);			//x = x + ap
  cp = axpy_norm(r,Ap,r,-a);		//r = r - aAp


    if ( ((k%100 == 0) && (verbose!=0)) || (verbose > 10)){
      ThreadBossMessage("MCR_prec: k= %d r^2= %le %le\n",k,cp,sqrt(cp/ssq));
    }

    // Stopping condition
    if ( cp <= rsq ) { 

      gettimeofday(&stop,NULL);
      struct timeval diff;
      timersub(&stop,&start,&diff);

      ThreadBossMessage("MCR_prec converged in %d iterations\n",k);
      ThreadBossMessage("MCR_prec converged in %d.%6.6d seconds\n",diff.tv_sec,diff.tv_usec);


      double flops = mprecFlops()*2.0 + 2.0*axpyNormFlops() + axpyFlops()*2.0;
             flops = flops * k;

      double t = diff.tv_sec*1.0E6 + diff.tv_usec;
      ThreadBossMessage("MCR_prec: %d mprec flops/site\n",mprecFlopsPerSite());
      ThreadBossMessage("MCR_prec: %le flops\n",flops);
      ThreadBossMessage("MCR_prec: %le mflops per node\n",flops/t);

      Mprec(psi,mp,tmp,0);
      Mprec(mp,Ap,tmp,1); 
      axpy(tmp,src,Ap,-1.0);
      double true_residual = sqrt(norm(tmp)/norm(src));
      ThreadBossMessage("MCR_prec: true residual is %le \n",true_residual);

      threadedFreeFermion(tmp);
    threadedFreeFermion(tmp1);
      threadedFreeFermion(p);
    threadedFreeFermion(Apold);
    threadedFreeFermion(pold);
      threadedFreeFermion(mp);
      threadedFreeFermion(r);
      threadedFreeFermion(Ap);
    threadedFreeFermion(Ar);
      if ( this->isBoss() && (!me) ) { 
	this->InverterExit();
      }
      return k;
    }

if( abs(a)>small || k==1){

  Mprec(r,mp,tmp,DaggerNo);
  G5R(mp,Ar);
  c = real( dot(Ar,Ap) );
  b = -c/d; //d = (Ap,Ap)

  axpy(pold,p,p,0);
  axpy(p,p,r,b);
  axpy(Apold,Ap,Ap,0);

  dold = d;
  d = axpy_norm(Ap,Ap,Ar,b);
}
else{

  Mprec(Ap,mp,tmp,DaggerNo); //(r,Ar)
  G5R(mp,tmp); // tmp = AAp
  g = real(dot(tmp,Ap))/d;
  del =  d / dold;
  if( abs(aold) > small ){
  	del *= (-1./aold);
  }

  axpy(tmp1,p,p,0);
  axpy(mp,p,Ap,-g);
  axpy(p,pold,mp,-del); //p = Ap - gp - delpold
  axpy(pold,tmp1,tmp1,0);

  axpy(tmp1,Ap,Ap,0);
  axpy(mp,Ap,tmp,-g);
  dold = d;
  d = axpy_norm(Ap,Apold,mp,-del); //  Ap = AAp - gAp - del Apold
  axpy(Apold,tmp1,tmp1,0);

  axpy(Ar,tmp,Ar,-a);

}

  }
  ThreadBossMessage("MCR_prec: CG not converged \n");
  threadedFreeFermion(tmp);
    threadedFreeFermion(tmp1);
  threadedFreeFermion(p);
    threadedFreeFermion(Apold);
    threadedFreeFermion(pold);
  threadedFreeFermion(mp);
  threadedFreeFermion(r);
  threadedFreeFermion(Ap);
    threadedFreeFermion(Ar);
  if ( this->isBoss() && (!me) ) { 
    this->InverterExit();
  }

  return -1;

}
template <class Float>
int bfmbase<Float>::MCR(Fermion_t solution[2], Fermion_t source[2])
{
  return MCR(solution,source,0,(bfm_fermion *)NULL,(double *)NULL);
}

template <class Float>
int bfmbase<Float>::MCR(Fermion_t solution[2], Fermion_t source[2], int Nev,bfm_fermion *evecs, double *evals)
{
  int me = thread_barrier();
  Fermion_t src = threadedAllocFermion(); 
  Fermion_t hsrc = threadedAllocFermion(); 
  Fermion_t tmp = threadedAllocFermion(); 
  Fermion_t Mtmp= threadedAllocFermion(); 
  Fermion_t lsol= threadedAllocFermion(); 

  ThreadBossMessage("MCR: Solver type is %d\n",solver);

  double f;
  f = norm(source[0]);
  f+= norm(source[1]);

  ThreadBossMessage("MCR: Source norm is %le\n",f);
  f = norm(solution[0]);
  f+= norm(solution[1]);
  ThreadBossMessage("MCR: Guess norm is %le\n",f);


  // src_o = Mdag * (source_o - Moe MeeInv source_e)
  //
  // When using the CGdiagonalMee, we left multiply system of equations by 
  // diag(MeeInv,MooInv), and must multiply by an extra MooInv
  //


    MooeeInv(source[Even],tmp,DaggerNo);
    Meo     (tmp,src,Odd,DaggerNo);
    axpy    (tmp,src,source[Odd],-1.0);
    Mprec(tmp,src,Mtmp,DaggerYes);  
  
    //G5R(src,hsrc);
    //double zo = axpy_norm(hsrc,hsrc,tmp,-1.0);
    //ThreadBossMessage("MCR: G5R twice %le\n", sqrt( zo ) );

  if(Nev > 0){

    ThreadBossMessage("MCR: deflating with %d evecs\n",Nev);

    axpy    (hsrc,src,src,0.0);
    axpby   (lsol,src,src,0.0,0.0);
  
    for(int n = 0; n<Nev; n++){
      std::complex<double> cn = dot (evecs[n][1], src) ;
      caxpy(lsol, evecs[n][1], lsol, real(cn)/ evals[n], imag(cn)/ evals[n]);
    }

    axpy(solution[Odd], lsol,lsol, 0.0);  

    f = norm(lsol);
    ThreadBossMessage("MCR: deflated guess %le\n",f);

  }
  
  int iter = MCR_prec(solution[Odd],src);
 
  threadedFreeFermion(Mtmp);
  threadedFreeFermion(lsol);
  threadedFreeFermion(hsrc);

  // sol_e = M_ee^-1 * ( src_e - Meo sol_o )...
  Meo(solution[Odd],tmp,Even,DaggerNo);
  axpy(src,tmp,source[Even],-1.0);
  MooeeInv(src,solution[Even],DaggerNo);
  
  f =norm(solution[0]);
  f+=norm(solution[1]);
  ThreadBossMessage("MCR: unprec sol norm is %le\n",f);
  
  threadedFreeFermion(tmp);
  threadedFreeFermion(src);

  return iter;

}

template <class Float>
int bfmbase<Float>::MCR_G5R(Fermion_t solution[2], Fermion_t source[2])
{
  int me = thread_barrier();
  Fermion_t src = threadedAllocFermion(); 
  Fermion_t hsrc = threadedAllocFermion(); 
  Fermion_t tmp = threadedAllocFermion(); 
  Fermion_t Mtmp= threadedAllocFermion(); 
  Fermion_t lsol= threadedAllocFermion(); 

  ThreadBossMessage("MCR: Solver type is %d\n",solver);

  double f;
  f = norm(source[0]);
  f+= norm(source[1]);

  ThreadBossMessage("MCR: Source norm is %le\n",f);
  f = norm(solution[0]);
  f+= norm(solution[1]);
  ThreadBossMessage("MCR: Guess norm is %le\n",f);


  // src_o = Mdag * (source_o - Moe MeeInv source_e)
  //
  // When using the CGdiagonalMee, we left multiply system of equations by 
  // diag(MeeInv,MooInv), and must multiply by an extra MooInv
  //
{
  ThreadBossMessage("MCR: check \n");
  G5R(src,hsrc);
  Mprec(hsrc,tmp,Mtmp,DaggerYes);

  Mprec(src,hsrc,Mtmp,DaggerNo);
  G5R(hsrc,lsol);

  
  f = axpy_norm(lsol,hsrc,lsol,-1.);
  ThreadBossMessage("MCR: diff is %le\n",f);

}


    MooeeInv(source[Even],tmp,DaggerNo);
    Meo     (tmp,src,Odd,DaggerNo);
    axpy    (tmp,src,source[Odd],-1.0);
    //Mprec(tmp,src,Mtmp,DaggerYes);  
    G5R(tmp,src);

  int iter = MCR_prec_G5R(solution[Odd],src);
 
  threadedFreeFermion(Mtmp);
  threadedFreeFermion(lsol);
  threadedFreeFermion(hsrc);

  // sol_e = M_ee^-1 * ( src_e - Meo sol_o )...
  Meo(solution[Odd],tmp,Even,DaggerNo);
  axpy(src,tmp,source[Even],-1.0);
  MooeeInv(src,solution[Even],DaggerNo);
  
  f =norm(solution[0]);
  f+=norm(solution[1]);
  ThreadBossMessage("MCR_G5R: unprec sol norm is %le\n",f);
  
  threadedFreeFermion(tmp);
  threadedFreeFermion(src);

  return iter;

}



template <class Float>
int bfmbase<Float>::MCR_single_shift(Fermion_t psi, Fermion_t src,double shift)
{

  double f;
  double cp,c,a,d,b,ci;
  std::complex<double> num;
  int me = thread_barrier();

  int stored_resid=0;

  if ( this->isBoss() && (!me) ) { 
    this->InverterEnter();
  }

  Fermion_t p   = threadedAllocFermion(mem_fast); 
  Fermion_t tmp = threadedAllocFermion(mem_fast); 
  Fermion_t mp  = threadedAllocFermion(mem_fast); 
  Fermion_t r   = threadedAllocFermion(mem_fast); 
  Fermion_t Ap  = threadedAllocFermion(mem_fast); 
  Fermion_t Ar  = threadedAllocFermion(mem_fast); 

  //Initial residual computation & set up
  fill(psi,0.0);
  double guess = norm(psi);
  Mprec(psi,mp,tmp,DaggerNo);
  Mprec(mp,Ar,tmp,DaggerYes);
  axpy(Ar,psi,Ar,shift);

  axpy (r, Ar, src,-1.0);
  axpy (p, Ar, src,-1.0);

  a =norm(p);
  cp=norm(r);

  Float ssq =  norm(src);
  Float rsq =  residual* residual*ssq;

  //Check if guess is really REALLY good :)
  if ( cp <= rsq ) {
    ThreadBossMessage("MCR_single_shift k=0 converged - suspiciously nice guess %le %le\n",cp,rsq);
    threadedFreeFermion(tmp);
    threadedFreeFermion(p);
    threadedFreeFermion(mp);
    threadedFreeFermion(r);
    threadedFreeFermion(Ap);
    threadedFreeFermion(Ar);
    if ( this->isBoss() && (!me) ) { 
      this->InverterExit();
    }
    return 0;
  }
  ThreadBossMessage("MCR_single_shift k=0 residual %le rsq %le\n",cp,rsq);

  if ( this->watchfile ) {
    ThreadBossMessage("MCR_single_shift watching file \"%s\"\n",this->watchfile);
  }
  struct timeval start,stop;
  gettimeofday(&start,NULL);


  c = Mprec(p,mp,tmp,DaggerNo);    // c = pdag[Mdag M] p
  d = Mprec(mp,Ap,tmp,DaggerYes);  // Ap   ; d=pdag Adag A p
  d=axpy_norm(Ap,p,Ap,shift);
  c+= shift*norm(p);
  
  axpy(Ar,Ap,Ap,0);	       //Ar

  //c = real( dot(r,Ar) );
  int k;
  for (k=1;k<=max_iter;k++){

    this->iter=k;

  //c = real( dot(r,Ar) );		//c = rAr
  a = c/d;

  axpy(psi,p,psi,a);		//x = x + ap
  cp = axpy_norm(r,Ap,r,-a);		//r = r - aAp

    if ( ((k%10 == 0) && (verbose!=0)) || (verbose > 10)){
      ThreadBossMessage("MCR_single_shift: k= %d r^2= %le %le\n",k,cp,sqrt(cp/ssq));
    }

    // Stopping condition
    if ( cp <= rsq ) { 

      gettimeofday(&stop,NULL);
      struct timeval diff;
      timersub(&stop,&start,&diff);
	
      ThreadBossMessage("MCR_single_shift converged in %d iterations %d.6.6d seconds\n",k,diff.tv_sec,diff.tv_usec);
      

      double flops = mprecFlops()*2.0 + 2.0*axpyNormFlops() + axpyFlops()*2.0;
      flops = flops * k;
      
      double t = diff.tv_sec*1.0E6 + diff.tv_usec;
      ThreadBossPerformance("MCR_single_shift: %le mflops per node\n",flops/t);
      
      Mprec(psi,mp,tmp,0);
      Mprec(mp,Ap,tmp,1); 
      axpy(tmp,src,Ap,-1.0);
      double true_residual = sqrt(norm(tmp)/norm(src));
      ThreadBossMessage("MCR_single_shift: true residual is %le \n",true_residual);
      
      threadedFreeFermion(tmp);
      threadedFreeFermion(p);
      threadedFreeFermion(mp);
      threadedFreeFermion(r);
      threadedFreeFermion(Ap);
      threadedFreeFermion(Ar);
      if ( this->isBoss() && (!me) ) { 
	this->InverterExit();
      }
      return k;
    }
    ci = d; //(Ap,Ap)
   //d = real( dot(r,Ar) );
    d = Mprec(r,mp,tmp,DaggerNo); //(r,Ar)
        Mprec(mp,Ar,tmp,DaggerYes);
    axpy(Ar,r,Ar,shift);
    d+=norm(r);
    b = d/c;
    c = d;
    axpy(p,p,r,b);
    d = axpy_norm(Ap,Ap,Ar,b);

  }
  
  {
    Mprec(psi,mp,tmp,0);
    Mprec(mp,Ap,tmp,1); 
    axpy(Ap,psi,Ap,shift);
    axpy(tmp,src,Ap,-1.0);
    double true_residual = sqrt(norm(tmp)/norm(src));
    ThreadBossMessage("MCR_single_shift: true residual is %le \n",true_residual);
    ThreadBossMessage("MCR_single_shift: CG not converged after %d iterations\n",k);
  }

  threadedFreeFermion(tmp);
  threadedFreeFermion(p);
  threadedFreeFermion(mp);
  threadedFreeFermion(r);
  threadedFreeFermion(Ap);
  threadedFreeFermion(Ar);
  if ( this->isBoss() && (!me) ) { 
    this->InverterExit();
  }
  return k;

}


#endif
