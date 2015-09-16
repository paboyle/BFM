
#ifndef BFM_CG_H
#define BFM_CG_H

/*  2007
 *  Calls PAB's assembler routines to give an implementation
 *  of the dwf dslash. Aim is to give very high performance
 *  in a reasonably portable/retargettable library.
 */
#include "bfm.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>


template <class Float>
int bfmbase<Float>::mprecFlopsPerSite ( void ) 
{
  int fourdwilson = 1320;
  int fivedhop    = 24;

  int flops;

  if ( precon_5d && (Ls > 1) ) {
    // A + B d d
    flops =  2*(fourdwilson + fivedhop) + 72;
  } else if ( Ls > 1 ) {
    flops = 2*(fourdwilson + 24*5);
  } if ( this->solver == CloverFermion ) {
    flops =  2*(fourdwilson + 576);
  }else {
    // A+ B d d 
    flops =  2*fourdwilson + 72;
  }
  return flops;
}

template <class Float>
int bfmbase<Float>::cbSites( void ) 
{
      int lx = node_latt[0];
      int ly = node_latt[1];
      int lz = node_latt[2];
      int lt = node_latt[3];

      int sites = (lx*ly*lz*lt*Ls)/2;

      return sites;
}

template <class Float>
int bfmbase<Float>::mprecFlops( void ) 
{
  return mprecFlopsPerSite() * cbSites();
}
template <class Float>
int bfmbase<Float>::axpyFlops ( void ) 
{
  return 2*24*cbSites();
}
template <class Float>
int bfmbase<Float>::axpyNormFlops ( void ) 
{
  return 4*24*cbSites();
}
template <class Float>
int bfmbase<Float>::axpyBytes ( void ) 
{
  return 8*24*3*cbSites();
}


  /*
   * Red black Schur decomposition
   *
   *  M = (Mee Meo) =  (1             0 )   (Mee   0               )  (1 Mee^{-1} Meo)
   *      (Moe Moo)    (Moe Mee^-1    1 )   (0   Moo-Moe Mee^-1 Meo)  (0   1         )
   *                =         L                     D                     U
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
   ***********************************************************
   */

  /*
   *
   * Routines:
   *
   *  bfm::CGNE_M 
   *  bfm::CGNE (historical name for CGNE_M)
   ***********************
   *     M psi = eta
   ***********************
   *Odd
   * i)   (D_oo)^{\dag} D_oo psi_o = (D_oo)^\dag L^{-1}  eta_o
   *                        eta_o' = D_oo (eta_o - Moe Mee^{-1} eta_e)
   *Even
   * ii)  Mee psi_e + Meo psi_o = src_e
   *
   *   => sol_e = M_ee^-1 * ( src_e - Meo sol_o )...
   *
   */

template <class Float>
int bfmbase<Float>::CGNE(Fermion_t solution[2], Fermion_t source[2])
{
  return CGNE_M(solution,source);
}
template <class Float>
int bfmbase<Float>::CGNE_M(Fermion_t solution[2], Fermion_t source[2])
{
  return CGNE_M(solution,source,0,(bfm_fermion *)NULL,(double *)NULL);
}
#ifdef BFM_QDP_BINDING
template <class Float>
int bfmbase<Float>::CGNE(Fermion_t solution[2], Fermion_t source[2], multi1d<bfm_fermion> &evecs, multi1d<double> &evals)
{
  return CGNE_M(solution,source, evecs, evals);
}

template <class Float>
int bfmbase<Float>::CGNE_M(Fermion_t solution[2], Fermion_t source[2], multi1d<bfm_fermion> &evecs, multi1d<double> &evals)
{
  return CGNE_M(solution,source,evecs.size(),&evecs[0],&evals[0]);
}
#endif 
template <class Float>
int bfmbase<Float>::CGNE_M(Fermion_t solution[2], Fermion_t source[2], int Nev, bfm_fermion *evecs, double *evals)
{
  int me = thread_barrier();
  Fermion_t src = threadedAllocFermion(); 
  Fermion_t hsrc = threadedAllocFermion(); 
  Fermion_t tmp = threadedAllocFermion(); 
  Fermion_t Mtmp= threadedAllocFermion(); 
  Fermion_t lsol= threadedAllocFermion(); 

  ThreadBossMessage("CGNE: Solver type is %s\n",SolverString(solver));
  ThreadBossMessage("CGNE: Float  size is %d\n",sizeof(Float));

  double f;
  f = norm(source[0]);
  f+= norm(source[1]);

  ThreadBossMessage("CGNE: Source norm is %le\n",f);
  f = norm(solution[0]);
  f+= norm(solution[1]);
  ThreadBossMessage("CGNE: Guess norm is %le\n",f);

  // src_o = Mdag * (source_o - Moe MeeInv source_e)
  //
  // When using the CGdiagonalMee, we left multiply system of equations by 
  // diag(MeeInv,MooInv), and must multiply by an extra MooInv
  //

  if ( CGdiagonalMee ) { 
    MooeeInv(source[Even],tmp,DaggerNo,Even);
    Meo     (tmp,src,Odd,DaggerNo);
    axpy    (src,src,source[Odd],-1.0);
    MooeeInv(src,tmp,DaggerNo);
    Mprec(tmp,src,Mtmp,DaggerYes,Even);  
  } else { 
    MooeeInv(source[Even],tmp,DaggerNo,Even);
    Meo     (tmp,src,Odd,DaggerNo);
    axpy    (tmp,src,source[Odd],-1.0);
    Mprec(tmp,src,Mtmp,DaggerYes);  
  }
  /*
  MooeeInv(source[Even],tmp,DaggerNo);
  f = norm(tmp);
  ThreadBossMessage("CGNE: prec tmp norm (MooeeInv) is %le\n",f);
  Meo     (tmp,src,Odd,DaggerNo);
  f = norm(src);
  ThreadBossMessage("CGNE: prec src norm (Meo) is %le\n",f);
  axpy    (tmp,src,source[Odd],-1.0);
  f = norm(tmp);
  ThreadBossMessage("CGNE: prec tmp norm (axpy) is %le\n",f);
  Mprec(tmp,src,Mtmp,DaggerYes);  
  f = norm(src);
  ThreadBossMessage("CGNE: prec src norm (MprecDag) is %le\n",f);
  */  
  if(Nev > 0){

    ThreadBossMessage("CGNE: deflating with %d evecs\n",Nev);

    axpy    (hsrc,src,src,0.0);
    axpby   (lsol,src,src,0.0,0.0);
  
    for(int n = 0; n<Nev; n++){
      std::complex<double> cn = dot (evecs[n][1], src) ;
      ///deflate src
      //axpy    (Mtmp,evecs[n][1],evecs[n][1],0.0);
      //scale   (Mtmp,real(cn),imag(cn) );
      //axpy    (hsrc,Mtmp,hsrc,-1.0);
      ///low solution
      //axpy    (lsol, Mtmp, lsol, 1.0/evals[n]);
	caxpy(lsol, evecs[n][1], lsol, real(cn)/ evals[n], imag(cn)/ evals[n]);
    }
  
    // Do we want to be able to control whether
    // i)   guess becomes hsrc.
    // ii)  guess becomes zero
    // iii) guess becomes orthogonal subspace projection of guess (costs a bit more)
    //      but means a perfect guess finishes and is useful to restarted CG.
    //      Importance is reduced once I implement defect correction single/double CG.

    axpy(solution[Odd], lsol,lsol, 0.0);  

  }

  int iter = CGNE_prec_MdagM(solution[Odd],src);
  
  if ( Nev > 0 ) { 
    //axpy(solution[Odd], lsol, solution[Odd], 1.0);  
  }
 
  threadedFreeFermion(lsol);
  threadedFreeFermion(hsrc);

  // sol_e = M_ee^-1 * ( src_e - Meo sol_o )...
  Meo(solution[Odd],tmp,Even,DaggerNo);
  axpy(src,tmp,source[Even],-1.0);
  MooeeInv(src,solution[Even],DaggerNo,Even);
  
  f =norm(solution[0]);
  f+=norm(solution[1]);
  ThreadBossMessage("CGNE: unprec sol norm is %le\n",f);

  // Even defect
  Meo(solution[Odd],tmp,Even,DaggerNo);
  Mooee(solution[Even],Mtmp,DaggerNo,Even);
  axpy(Mtmp,Mtmp,tmp,1.0);
  axpy(Mtmp,Mtmp,source[Even],-1.0);
  double rf=norm(Mtmp);
  // Odd defect
  Meo(solution[Even],tmp,Odd,DaggerNo);
  Mooee(solution[Odd],Mtmp,DaggerNo,Odd);
  axpy(Mtmp,Mtmp,tmp,1.0);
  axpy(Mtmp,Mtmp,source[Odd],-1.0);
  rf+=norm(Mtmp);
  ThreadBossMessage("CGNE: unprec sol true residual is %le\n",sqrt(rf/f));
  
  threadedFreeFermion(Mtmp);
  threadedFreeFermion(tmp);
  threadedFreeFermion(src);


  return iter;
}

  /*
   ***********************
   * Routine:
   *
   *  bfm::CGNE_Mdag
   *
   ***********************
   *     Mdag psi = U^dag, D^dag L^{dag} = eta
   ***********************
   *
   *Odd
   * i)  D_oo^dagD_oo phi_o = U^{-dag} eta_o
   *
   *     eta_o' = eta_o - Moe^{dag}  Mee^{-dag} eta_e
   *
   *     phi_o = D_oo^{-1} Doo^{-dag} eta_o'
   *
   *     psi_o = D_oo phi_o 
   *
   *Even
   *ii) 
   *     Mee^{d} psi_e + Meo^{d} psi_o = src_e
   * =>  psi_e = Mee^{-d} ( src_e - Meo^{d} psi_o)
   *
   */

template <class Float>
int bfmbase<Float>::CGNE_Mdag(Fermion_t solution[2], Fermion_t source[2])
{
  Fermion_t src = threadedAllocFermion(); 
  Fermion_t tmp = threadedAllocFermion(); 

  // src_o = source_o - Moe^d MeeInv^d source_e
  MooeeInv(source[Even],src,DaggerYes,Even);
  Meo     (src,tmp,Odd,DaggerYes);
  axpy    (src,tmp,source[Odd],-1.0);

  //  phi_o = D_oo^{-1} Doo^{-dag} eta_o'
  int iter = CGNE_prec_MdagM(tmp,src);

  //     psi_o = D_oo phi_o 
  Mprec(tmp,solution[Odd],src,DaggerNo);  
  
  // psi_e = Mee^{-d} ( src_e - Meo^{d} psi_o)
  Meo(solution[Odd],tmp,Even,DaggerYes);
  axpy(src,tmp,source[Even],-1.0);
  MooeeInv(src,solution[Even],DaggerYes,Even);

  threadedFreeFermion(tmp);
  threadedFreeFermion(src);

  return iter;
}

  /*
   *
   ***********************
   * Routine:
   *  bfm::CGNE_prec_MdagM
   ***********************
   *     MdagM psi = eta
   ***********************
   *   Simply solves
   *   D_oo^dagD_oo phi_o = L^{-1} eta_o
   *
   */


template <class Float>
int bfmbase<Float>::CGNE_prec_MdagM(Fermion_t psi, Fermion_t src)
{

  int iter;


  int me= thread_barrier();


  if ( reproduce ) { 

    uint64_t c1,c2;

    Fermion_t tmp_r = threadedAllocFermion(); 
    axpy    (tmp_r,psi,psi,0.0);

    c1=checksum(tmp_r);
    c2=checksum(psi);

    repro_control(ReproRecord);
    iter = CGNE_prec(tmp_r,src);
    ThreadBossMessage(" Checking reproducibility\n");
    repro_control(ReproCheck);
    iter = CGNE_prec(psi,src);
    repro_control(ReproIdle);

    c1=checksum(tmp_r);
    c2=checksum(psi);
    if ( c1!=c2  ) {
       Error(" result vector checksum miscompare in reproducing CG %16.16lx %16.16lx\n",c1,c2);
       exit(-1);
    }
    ThreadBossMessage(" reproducing cg passed -- checksums are %16.16lx %16.16lx\n",c1,c2);

    threadedFreeFermion(tmp_r);

  } else {

    double guess = norm (psi);

    int me = this->thread_barrier();

    ThreadBossDebug("Not reproducing -- guess is %le\n",guess);
    iter = CGNE_prec(psi,src);
    repro_control(ReproIdle);

  }
  return iter;
}
template <class Float>
int bfmbase<Float>::CGNE_prec(Fermion_t psi, Fermion_t src)
{
  return CGNE_prec(psi,src,(Fermion_t*)NULL,0);
}
template <class Float>
int bfmbase<Float>::CGNE_prec(Fermion_t psi, Fermion_t src,
			     Fermion_t *residuals,int nresid) // Internal
{

  double f;
  double cp,c,a,d,b;

  int me = thread_barrier();

  int stored_resid=0;

  if ( this->isBoss() && (!me) ) { 
    this->InverterEnter();
  }

  Fermion_t p   = threadedAllocFermion(mem_fast); 
  Fermion_t tmp = threadedAllocFermion(mem_fast); 
  Fermion_t mp  = threadedAllocFermion(mem_fast); 
  Fermion_t mmp = threadedAllocFermion(mem_fast); 
  Fermion_t r   = threadedAllocFermion(mem_fast); 

  //Initial residual computation & set up
  double guess = norm(psi);

  d=Mprec(psi,mp,tmp,DaggerNo);
  b=Mprec(mp,mmp,tmp,DaggerYes);

  axpy (r, mmp, src,-1.0);
  axpy (p, mmp, src,-1.0);

  a =norm(p);
  cp= norm(r);

  Float ssq =  norm(src);
  ThreadBossLog("CGNE_prec gues %le \n",guess);
  ThreadBossLog("CGNE_prec src  %le \n",ssq);
  ThreadBossLog("CGNE_prec  mp  %le \n",d);
  ThreadBossLog("CGNE_prec  mmp %le \n",b);
  ThreadBossLog("CGNE_prec   r  %le \n",cp);
  ThreadBossLog("CGNE_prec   p  %le \n",a);

  Float rsq =  residual* residual*ssq;

  //Check if guess is really REALLY good :)
  if ( cp <= rsq ) {
    ThreadBossMessage("CGNE_prec k=0 converged - suspiciously nice guess %le %le\n",cp,rsq);
    threadedFreeFermion(tmp);
    threadedFreeFermion(p);
    threadedFreeFermion(mp);
    threadedFreeFermion(mmp);
    threadedFreeFermion(r);
    if ( this->isBoss() && (!me) ) { 
      this->InverterExit();
    }
    return 0;
  }

  ThreadBossDebug("CGNE_prec k=0 residual %le rsq %le\n",cp,rsq);

  if ( this->watchfile ) {
    ThreadBossDebug("CGNE_prec watching file \"%s\"\n",this->watchfile);
  }
  struct timeval start,stop;
  gettimeofday(&start,NULL);

  int k;
  for (k=1;k<=max_iter;k++){

    this->iter=k;
    uint64_t t_iter_1=GetTimeBase();

    c=cp;

    uint64_t t_mprec_1=GetTimeBase();
    d = Mprec(p,mp,tmp,0,1);// Dag no
    uint64_t t_mprec_2=GetTimeBase();
    a = c/d;

    //NB: No longer need mp, overwrite
    //    With optimised version need to take
    //            mp-> Dmp
    //        Dmp,mp-> amp+DDmp -> mp inplace.
    // Implement "Mprec_inplace(mp,mmp,1); 
    uint64_t t_mprec_3=GetTimeBase();
    double qq=Mprec(mp,mmp,tmp,1); // Dag yes
    uint64_t t_mprec_4=GetTimeBase();
    double b_pred = a*(a*qq-d)/c;

#if 0
    // Can fuse these together. 
    // Further, producing "p" in a cache-blocked fashion could allow no reload
    // of p from main memory (~8 timeslices at a time?)

    //Fused axpy-norm required.    
    //As know "a" in advance could save bandwidth by
    //not writing out mmp, but doing linear combination.
    //in Mprec.
    uint64_t tr1=GetTimeBase();
    cp = axpy_norm(r,mmp,r,-a);
    b = cp/c;
    uint64_t tr2=GetTimeBase();

    uint64_t tpsi1=GetTimeBase();
    axpy(psi,p,psi,a);
    uint64_t tpsi2=GetTimeBase();
    // New (conjugate/M-orthogonal) search direction
    uint64_t tp1=GetTimeBase();
    axpy(p,p,r,b);
    uint64_t tp2=GetTimeBase();
#else
    uint64_t tpsi1=GetTimeBase();
    uint64_t tp1=GetTimeBase();
    uint64_t tp2=tp1+1;
    uint64_t tr1=tp1;
    uint64_t tr2=tp2;
    b  = b_pred;
    cp = cg_psi_p(psi,p,r,mmp,a,b);
    uint64_t tpsi2=GetTimeBase();
#endif
    uint64_t t_iter_2=GetTimeBase();

    
    // verbose nonsense
    if ( (this->iter==this->time_report_iter) ) {
      
      int lx = node_latt[0];
      int ly = node_latt[1];
      int lz = node_latt[2];
      int lt = node_latt[3];

      int cb4dsites = (lx*ly*lz*lt)/2;
      ThreadBossDebug("fermionCacheFootprint: %ld \n",7*axpyBytes()/3);
      ThreadBossDebug("gauge  CacheFootprint: %ld \n",2*18*8*cb4dsites*2);
      ThreadBossDebug("fermionVecBytes      : %ld \n",axpyBytes()/3);
      //      ThreadBossDebug("axpyBytes            : %ld \n",axpyBytes());
      //ThreadBossDebug("axpy      (soln)     : %ld cyc %le MB/s\n",(tpsi2-tpsi1),(double)axpyBytes()*1600./(tpsi2-tpsi1));
      //ThreadBossDebug("axpy_norm (residual) : %ld cyc %le MB/s\n",(tr2-tr1),(double)axpyBytes()*1600./(tr2-tr1));
      // ThreadBossDebug("axpy      (search)   : %ld cyc %le MB/s\n",(tp2-tp1),(double)axpyBytes()*1600./(tp2-tp1));
      ThreadBossDebug("Iter time            : %ld cyc\n",t_iter_2-t_iter_1);
      ThreadBossDebug("linalg time          : %ld cyc\n",t_iter_2-t_iter_1-(t_mprec_2-t_mprec_1)-(t_mprec_4-t_mprec_3));
      ThreadBossDebug("Mprec time           : %ld cyc\n",t_mprec_2-t_mprec_1);
      ThreadBossDebug("Mprec time           : %ld cyc\n",t_mprec_4-t_mprec_3);
      fflush(stdout);
    }


    if ( ((k%1 == 0) && (verbose!=0)) || (verbose > 10) ){
      double tr=0;
      if(0){
	Mprec(psi,mp,tmp,0);
	Mprec(mp,mmp,tmp,1); 
	axpy(tmp,src,mmp,-1.0);
        tr = norm(tmp);
      }
      ThreadBossMessage("CGNE_prec: k=%d r^2=%le rel %le true %le \n",k,cp,sqrt(cp/ssq),sqrt(tr/ssq));
    }

    {
      double rr=sqrt(cp/ssq);
      this->InverterLogIteration(k, rr,src,psi,mp,mmp,tmp);
    }


    // Stopping condition
    if ( cp <= rsq ) { 

      if ( stored_resid >= nresid ) { 

	gettimeofday(&stop,NULL);
	struct timeval diff;
	timersub(&stop,&start,&diff);
	


	double flops = mprecFlops()*2.0 + 2.0*axpyNormFlops() + axpyFlops()*2.0;
	flops = flops * k;

	double t = diff.tv_sec*1.0E6 + diff.tv_usec;
	ThreadBossMessage("CGNE_prec: %d mprec flops/site\n",mprecFlopsPerSite());
	ThreadBossMessage("CGNE_prec: %le flops\n",flops);
	ThreadBossMessage("CGNE_prec: %le mflops per node\n",flops/t);
	
	Mprec(psi,mp,tmp,0);
	Mprec(mp,mmp,tmp,1); 
	axpy(tmp,src,mmp,-1.0);
	
	double  mpnorm = sqrt(norm(mp));
	double mmpnorm = sqrt(norm(mmp));
	double psinorm = sqrt(norm(psi));
	double srcnorm = sqrt(norm(src));
	double tmpnorm = sqrt(norm(tmp));
	double true_residual = tmpnorm/srcnorm;
	ThreadBossMessage("CGNE_prec converged in %d iterations %d.%6.6d s, resid=%le, mass=%le\n",
			  k,diff.tv_sec,diff.tv_usec,true_residual,this->mass);
	
	threadedFreeFermion(tmp);
	threadedFreeFermion(p);
	threadedFreeFermion(mp);
	threadedFreeFermion(mmp);
	threadedFreeFermion(r);
	if ( this->isBoss() && (!me) ) { 
	  this->InverterExit();
	}
	return k;
      } else { 
	axpy(residuals[stored_resid],r,r,0.0);
	ThreadBossLog("CG_prec: storing residual %d/%d\n",stored_resid,nresid);
	//	if ( (k % 20) == 0 ) { 
	  stored_resid++;
	//	}
      }
    }

  }
  ThreadBossMessage("CGNE_prec: CG not converged after %d iter resid %le\n",k,sqrt(cp/ssq));
  threadedFreeFermion(tmp);
  threadedFreeFermion(p);
  threadedFreeFermion(mp);
  threadedFreeFermion(mmp);
  threadedFreeFermion(r);
  if ( this->isBoss() && (!me) ) { 
    this->InverterExit();
  }

  return -1;

}


/**********************************************************
*                                                         *
*           PV taku trick. Can probably do this better    *
*                                                         *
***********************************************************/

template <class Float>
int bfmbase<Float>::CGNE_PV_M(Fermion_t solution[2], Fermion_t source[2],int n, int L5)
{
  Fermion_t src = threadedAllocFermion(); 
  Fermion_t tmp = threadedAllocFermion(); 
  Fermion_t Mtmp= threadedAllocFermion(); 

  double f;
  f = norm(source[0]);
  f+= norm(source[1]);

  int me = thread_barrier();
  //ThreadBossMessage("CGNE: Source norm is %le\n",f);

  f = norm(solution[0]);
  f+= norm(solution[1]);

  //ThreadBossMessage("CGNE: Guess norm is %le\n",f);

  // src_o = Mdag * (source_o - Moe MeeInv source_e)
  MooeeInv5dprec_TW(source[Even],tmp,n,L5,DaggerNo);
  f = norm(tmp);

  //ThreadBossMessage("CGNE: tmp norm is %le\n",f);
  Meo     (tmp,src,Odd,DaggerNo);
  f = norm(src);

  //ThreadBossMessage("CGNE: dsl norm is %le\n",f*4);
  axpy    (tmp,src,source[Odd],-1.0);
  
  MprecTW(tmp,src,Mtmp,DaggerYes,n,L5,0);  
  threadedFreeFermion(Mtmp);
  
  int iter = CGNE_PV_prec_M(solution[Odd],src,n,L5);
  axpby(solution[Even],tmp,tmp,0.0,0.0);
  
  // sol_e = M_ee^-1 * ( src_e - Meo sol_o )...
  Meo(solution[Odd],tmp,Even,DaggerNo);
  axpy(src,tmp,source[Even],-1.0);
  MooeeInv5dprec_TW(src,solution[Even],n,L5,DaggerNo);


  threadedFreeFermion(tmp);
  threadedFreeFermion(src);

  return iter;
}

template <class Float>
int bfmbase<Float>::CGNE_PV_Mdag(Fermion_t solution[2], Fermion_t source[2],int n, int L5)
{
  Fermion_t src = threadedAllocFermion(); 
  Fermion_t tmp = threadedAllocFermion(); 

  // src_o = source_o - Moe^d MeeInv^d source_e
  MooeeInv5dprec_TW(source[Even],src,n,L5,DaggerYes);
  Meo     (src,tmp,Odd,DaggerYes);
  
  axpy    (src,tmp,source[Odd],-1.0);

  //  phi_o = D_oo^{-1} Doo^{-dag} eta_o'
  int iter = CGNE_PV_prec_M(tmp,src, n, L5);

  //     psi_o = D_oo phi_o 
  MprecTW(tmp,solution[Odd],src,DaggerNo,n,L5,0);  
  
  // psi_e = Mee^{-d} ( src_e - Meo^{d} psi_o)
  Meo(solution[Odd],tmp,Even,DaggerYes);
  
  axpy(src,tmp,source[Even],-1.0);
  MooeeInv5dprec_TW(src,solution[Even],n,L5,DaggerYes);

  threadedFreeFermion(tmp);
  threadedFreeFermion(src);

  return iter;
}

template <class Float>
int bfmbase<Float>::CGNE_PV_prec_M(Fermion_t psi, Fermion_t src,int n, int L5)
{
  int iter;
  int me= thread_barrier();

  if ( reproduce ) { 

    uint64_t c1,c2;

    Fermion_t tmp_r = threadedAllocFermion(); 
    axpy    (tmp_r,psi,psi,0.0);


    c1=checksum(tmp_r);
    c2=checksum(psi);

    repro_control(ReproRecord);
    iter = CGNE_PV_prec(tmp_r,src,n,L5);
    ThreadBossMessage(" Checking reproducibility\n");
    repro_control(ReproCheck);
    iter = CGNE_PV_prec(psi,src,n,L5);
    repro_control(ReproIdle);

    c1=checksum(tmp_r);
    c2=checksum(psi);
    if ( c1!=c2  ) {
       Error(" result vector checksum miscompare in reproducing CG %16.16lx %16.16lx\n",c1,c2);
       exit(-1);
    }
    ThreadBossMessage(" reproducing cg passed -- checksums are %16.16lx %16.16lx\n",c1,c2);

    threadedFreeFermion(tmp_r);

  } else {
    iter = CGNE_PV_prec(psi,src,n,L5);
    repro_control(ReproIdle);

  }
  return iter;
}


template <class Float>
int bfmbase<Float>::CGNE_PV_prec(Fermion_t psi, Fermion_t src,int n, int L5)
{

  double f;
  double cp,c,a,d,b;
  int me = thread_barrier();

  Fermion_t p   = threadedAllocFermion(mem_fast); 
  Fermion_t tmp = threadedAllocFermion(mem_fast); 
  Fermion_t mp  = threadedAllocFermion(mem_fast); 
  Fermion_t mmp = threadedAllocFermion(mem_fast); 
  Fermion_t r   = threadedAllocFermion(mem_fast); 

  //Initial residual computation & set up
  MprecTW(psi,mp,tmp,DaggerNo,n,L5,0);
  MprecTW(mp,mmp,tmp,DaggerYes,n,L5,0);



  axpy (r, mmp, src,-1.0);
  axpy (p, mmp, src,-1.0);

  cp= norm(r);
  Float ssq =  norm(src);
  Float rsq =  residual* residual*ssq;
  
  //Check if guess is really REALLY good :)
  if ( cp <= rsq ) {
    ThreadBossMessage("CGNE_PV_prec k=0 converged - nice guess %le %le\n",cp,rsq);
    
    threadedFreeFermion(tmp);
    threadedFreeFermion(p);
    threadedFreeFermion(mp);
    threadedFreeFermion(mmp);
    threadedFreeFermion(r);
    return 0;
  }
  //ThreadBossMessage("CGNE_PV_prec k=0 residual %le rsq %le\n",cp,rsq);

  struct timeval start,stop;
  gettimeofday(&start,NULL);


  for (int k=1;k<=max_iter;k++){
    c=cp;

    // Mp
    d = MprecTW(p,mp,tmp,0,n,L5,1);

    // Update solution vector
    a = c/d;
    axpy(psi,p,psi,a);
    
    // MdagMp
    MprecTW(mp,mmp,tmp,1,n,L5,0); 

    //Fused axpy-norm required.    
    //As know "a" in advance could save bandwidth by
    //not writing out mmp, but doing linear combination.
    //in Mprec.
    cp = axpy_norm(r,mmp,r,-a);
    
    // Stopping condition
    if ( cp <= rsq ) { 

      gettimeofday(&stop,NULL);
      struct timeval diff;
      timersub(&stop,&start,&diff);

      ThreadBossMessage("CGNE_PV_prec converged in %d iterations\n",k);
      //ThreadBossMessage("CGNE_PV_prec converged in %d.%6.6d s\n",diff.tv_sec,diff.tv_usec);

      int lx = node_latt[0];
      int ly = node_latt[1];
      int lz = node_latt[2];
      int lt = node_latt[3];
      int sites = (lx*ly*lz*lt*Ls)/2;

      //double flops = 1.0*sites*k*(  2.0*mprecFlopsPerSite() // Two matrix multiples
	//			    + 48*3                  // Three AXPY's
	//			    + 48*2  );              // Two   Norms (one in Mprec)

      //double t = diff.tv_sec*1.0E6 + diff.tv_usec;
      //ThreadBossMessage("CGNE_PV_prec: %d mprec flops/site\n",mprecFlopsPerSite());
      //ThreadBossMessage("CGNE_PV_prec: %le flops\n",flops);
      //ThreadBossMessage("CGNE_PV_prec: %le mflops\n",flops/t);

      MprecTW(psi,mp,tmp,0,n,L5,0);
      MprecTW(mp,mmp,tmp,1,n,L5,0); 
      axpy(tmp,src,mmp,-1.0);
      double true_residual = sqrt(norm(tmp)/norm(src));
      ThreadBossMessage("CGNE_PV_prec: n = %d true residual is %le \n",n,true_residual);

      threadedFreeFermion(tmp);
      threadedFreeFermion(p);
      threadedFreeFermion(mp);
      threadedFreeFermion(mmp);
      threadedFreeFermion(r);
      return k;
    }
    // New (conjugate/M-orthogonal) search direction
    b = cp/c;
    axpy(p,p,r,b);

  }
  ThreadBossMessage("CGNE_PV_prec: CG not converged \n");
  threadedFreeFermion(tmp);
  threadedFreeFermion(p);
  threadedFreeFermion(mp);
  threadedFreeFermion(mmp);
  threadedFreeFermion(r);
  return -1;
}


/*
CG iterates: to express as a polynomial,

r= [1,0,0.... ]
p= [1,0,0.... ]
x= [0,0,0.... ]

Iterate:
  Ap= right_shift(p)
  r = r - a Ap
  x = x + a Ap
  p = r + b p

Poly = psi[0] + psi[1] x^1 etc..
*/
template <class Float>
int bfmbase<Float>::CGNE_single_shift(Fermion_t psi, Fermion_t src,double shift,
				      std::vector<double> & polynomial)
{
  double f;
  double cp,c,a,d,b,qq,b_pred;
  double pnorm;
  int me = thread_barrier();

  Fermion_t p   = threadedAllocFermion(mem_fast); 
  Fermion_t tmp = threadedAllocFermion(mem_fast); 
  Fermion_t mp  = threadedAllocFermion(mem_fast); 
  Fermion_t mmp = threadedAllocFermion(mem_fast); 
  Fermion_t r   = threadedAllocFermion(mem_fast); 

  struct timeval start,stop;
  gettimeofday(&start,NULL);
    
  //Initial residual computation & set up
  axpby(psi,src,src,0,0);
  //fill(psi,0.0);

  axpy (r, src, src,0.0);
  axpy (p, src, src,0.0);
  
  pnorm = axpy_norm(p,p,p,0.0);
  a =pnorm;
  cp=pnorm;
    
  qq=0;
  b=0;
  d=0;
  
  Float ssq    = axpy_norm(src,src,src,0.0);
  Float rsq =  residual* residual*ssq;

  ThreadBossMessage("CGNE_single_shift %le  k=0, cp=%le initial resid %le\n",shift,cp,sqrt(cp/ssq));

  std::vector<double> poly_p(1);
  std::vector<double> poly_r(1);
  std::vector<double> poly_Ap;

  if ( !me ) { 
    poly_p[0]=1.0;
    poly_r[0]=1.0;
    polynomial.resize(0);
  }

  int k;
  for (k=1;k<=max_iter;k++){
      
    this->iter=k;

    c=cp;

    // mmp=[MdagM + shift]p
    //
    //   d=p^dag [MdagM+s] p
    //    =pdag MdagM p + s PdagP
    //
    //   qq= norm( [MdagM+s] p )
    d = Mprec(p,mp,tmp,0,1);// Dag no
        Mprec(mp,mmp,tmp,1);// Dag yes



    qq = axpy_norm(mmp,p,mmp,shift);  
    pnorm=axpy_norm(p,p,p,0.0);//Wasteful. Add norm to cg_psi_p
    d+= shift*pnorm;

    a = c/d;

    b = a*(a*qq-d)/c;
    
    cp = cg_psi_p(psi,p,r,mmp,a,b);

    if(!me){
      //  Ap= right_shift(p)
      poly_Ap.resize(k+1);
      poly_Ap[0]=0.0;
      for(int i=0;i<k;i++){
	poly_Ap[i+1]=poly_p[i];
      }

      //  x = x + a p
      polynomial.resize(k);
      polynomial[k-1]=0.0;
      for(int i=0;i<k;i++){
	polynomial[i] = polynomial[i] + a * poly_p[i];
      }
    
      //  r = r - a Ap
      //  p = r + b p
      poly_r.resize(k+1);
      poly_p.resize(k+1);
      poly_r[k] = poly_p[k] = 0.0;
      for(int i=0;i<k+1;i++){
	poly_r[i] = poly_r[i] - a * poly_Ap[i];
	poly_p[i] = poly_r[i] + b * poly_p[i];
      }
    }

    if ( (k%100==0) ){
      ThreadBossDebug("CGNE_single_shift k=%d, cp=%le %le %le %le %le\n",k,cp,a,b,qq,d); 
    }
    // Stopping condition
    if ( cp <= rsq ) { 

      ThreadBossMessage("CGNE_single_shift converged in %d iterations resid %le\n",k,sqrt(cp/ssq));

      goto cleanup;
	
    }

  }


 cleanup:

  gettimeofday(&stop,NULL);
  struct timeval diff;
  timersub(&stop,&start,&diff);
  
  double flops = mprecFlops()*2.0 + 2.0*axpyNormFlops() + axpyFlops()*2.0;
  flops = flops * (k-1);
  double t = diff.tv_sec*1.0E6 + diff.tv_usec;
  ThreadBossMessage("CGNE_single_shift: shift %le CG terminating after %d iter resid %le : %le mflops %le s\n",
		    shift,k,sqrt(cp/ssq),
		    flops/t,t*1.0e-6,t*1.0e-6);

  for(int i=0;i<polynomial.size();i++){
    //    ThreadBossMessage("CGpolynomial[%d]: %.16le\n",i,polynomial[i]);
  }

  if ( 0 ) { 
    Mprec(psi,mp,tmp,0,1);// Dag no
    Mprec(mp,mmp,tmp,1);// Dag yes
    axpy(mmp,psi,mmp,shift);  
    axpy(mmp,src,mmp,-1.0);
    double true_residual = sqrt(norm(mmp)/norm(src));
    ThreadBossMessage("CGNE_single_shift true residual %le\n",true_residual);
  }

  threadedFreeFermion(tmp);
  threadedFreeFermion(p);
  threadedFreeFermion(mp);
  threadedFreeFermion(mmp);
  threadedFreeFermion(r); 
  return k;

}


/*************************************************************************
**************************************************************************
**************************************************************************/

#include "bfm_cg_multi.h"

#endif

