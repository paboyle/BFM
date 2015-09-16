
#ifndef BFM_CG_UNPREC_H
#define BFM_CG_UNPREC_H

#include <bfm.h>
#include <bfm_qdp.h>

template <class Float>
int bfmbase<Float>::CGNE_Mdag_unprec(Fermion_t sol_guess[2], Fermion_t source[2])
{
  Fermion_t solp[2]; // Both checkerboards
  for(int cb=0;cb<2;cb++){
    solp[cb]   = this->threadedAllocFermion();
  }

  axpy(solp,sol_guess,sol_guess,0.0);

  int me = this->thread_barrier();
  this->ThreadBossLog("CGNE_Mdag_unprec\n");

  int k =CGNE_MdagM_unprec(solp,source);

  G5D_Munprec(solp,sol_guess,DaggerNo); // CGNE,  Sol' = (MdagM)^{-1} Src = M^{1} Mdag^{-1} Src
                                        //        Sol  = M Sol'

  // Verify
  G5D_Munprec(sol_guess,solp,DaggerYes);
  axpy(solp,solp,source,-1.0);

  double nrm = norm(solp);
  double sn  = norm(source);
  if ( (nrm/sn) > 1.0e-4 ) { 
    this->Error("suspicious Mdag_unprec true resid ^%le\n",nrm);
    exit(-1);
  }

  for(int cb=0;cb<2;cb++){
    this->threadedFreeFermion(solp[cb]);
  }
  return k;
}

template <class Float>
int bfmbase<Float>::CGNE_unprec(Fermion_t sol_guess[2], Fermion_t source[2])
{
  return CGNE_M_unprec(sol_guess,source);
}

template <class Float>
int bfmbase<Float>::CGNE_M_unprec(Fermion_t sol_guess[2], Fermion_t source[2])
{
    Fermion_t srcp[2]; // Both checkerboards

    for(int cb=0;cb<2;cb++){
      srcp[cb]   = this->threadedAllocFermion();
    }

    G5D_Munprec(source,srcp,DaggerYes); // CGNE, Src' = Mdag Src, and Sol = (MdagM)^{-1}Mdag Src

    int me = this->thread_barrier();
    this->ThreadBossLog("CGNE_M_unprec\n");

    int k=CGNE_MdagM_unprec(sol_guess,srcp);

  // Verify
  G5D_Munprec(sol_guess,srcp,DaggerNo);
  axpy(srcp,srcp,source,-1.0);
  double sn  = norm(source);
  double nrm = norm(srcp);
  if ( (nrm/sn) > 1.0e-4 ) { 
    this->Error("suspicious M_unprec true resid ^%le\n",nrm);
    exit(-1);
  }

    for(int cb=0;cb<2;cb++){
      this->threadedFreeFermion(srcp[cb]);
    }

    return k;
}

template <class Float>
int bfmbase<Float>::CGNE_MdagM_unprec(Fermion_t sol[2], Fermion_t srcp[2])
  {

    //////////////////////////////////////////////
    //Start the CG for unprec matrix
    //////////////////////////////////////////////
    Fermion_t p[2]; // Both checkerboards
    Fermion_t mp[2];
    Fermion_t mmp[2];
    Fermion_t r[2];
    Fermion_t tmp;

    int me = this->thread_barrier();

    for(int cb=0;cb<2;cb++){
      p[cb]   = this->threadedAllocFermion();
      mp[cb]  = this->threadedAllocFermion();
      mmp[cb] = this->threadedAllocFermion();
      r[cb] = this->threadedAllocFermion();
    }
    tmp=this->threadedAllocFermion();

    double f;
    this->ThreadBossLog("CGNE_unprec: invert solver %s mass %le\n",SolverString(solver),mass);
    f = norm(srcp);
    this->ThreadBossLog("CGNE_unperc: Source norm %le\n",f);
    f = norm(sol);
    this->ThreadBossMessage("CGNE_unprec: Guess  norm %le\n",f);
    
    
    double ssq  =norm(srcp);
    double residual=this->residual;
    double rsq =  residual* residual*ssq;

    // Initial residual r = Src' - MdagM sol
    G5D_Munprec(sol,mp,DaggerNo);
    G5D_Munprec(mp,mmp,DaggerYes);
    f= norm(mmp);
    this->ThreadBossMessage("CGNE_unprec_MdagM: guess=%le\n",f);

    axpy(p,mmp,srcp,-1.0);
    axpy(r,mmp,srcp,-1.0);

    double cp  = norm(r);
    this->ThreadBossMessage("CGNE_unprec_MdagM: guess residual %le\n", cp);
           
    if ( cp <= rsq ) {
      this->ThreadBossMessage("CGNE_unprec converged k=0 - nice guess \n");

      for(int cb=0;cb<2;cb++){
	this->threadedFreeFermion(p[cb]);
	this->threadedFreeFermion(mp[cb]); 
	this->threadedFreeFermion(mmp[cb]);
	this->threadedFreeFermion(r[cb]);
      }
      this->threadedFreeFermion(tmp);
      return 0;
    }  

    double c,a,d,b;

    for (int k=1;k<=max_iter;k++){

      //      this->ThreadBossMessage("CGNE_unprec: iteration %d\n",k);
      c=cp;

      G5D_Munprec(p,mp,DaggerNo);
      d = norm(mp);

      // Update solution vector
      a = c/d;
      axpy(sol,p,sol,a);
    
      // MdagMp
      G5D_Munprec(mp,mmp,DaggerYes); 
      axpy(r,mmp,r,-a);
      cp= norm(r);

      if ( k%100==0 ) {
	this->ThreadBossMessage("CGNE_unprec k=%d r^2=%le\n",k,cp);
      }

      // Stopping condition
      if ( cp <= rsq ) { 
	this->ThreadBossMessage("CGNE_unprec converged k=%d r^2=%le rsq %le\n",k,cp,rsq);
	double nrm = norm(sol);
	this->ThreadBossMessage("CGNE_unprec sol norm  k=%d %le \n",k,nrm);
	nrm = norm(sol,0);
	this->ThreadBossMessage("CGNE_unprec sol[0] norm  k=%d %le \n",k,nrm);
	nrm = norm(sol,this->Ls-1);
	this->ThreadBossMessage("CGNE_unprec sol[Ls-1] norm  k=%d %le \n",k,nrm);

	G5D_Munprec(sol,mp,DaggerNo);
	G5D_Munprec(mp,mmp,DaggerYes); 
	axpy(r,mmp,srcp,-1.0);
	double true_residual = norm(r)/ssq;
	this->ThreadBossMessage("CGNE_unprec true_residual %le\n",sqrt(true_residual));
	this->ThreadBossMessage("CGNE_unprec true_residual^2 %le\n",true_residual);
	/*
	for(int s=0;s<this->cbLs;s++){
	  double nrm=this->norm(r,s);
	  nrm = nrm/ssq;
	  this->ThreadBossMessage("CGNE_unprec s-slice %d = %le\n",s,nrm);
	  
	}
	*/

	for(int cb=0;cb<2;cb++){
	  this->threadedFreeFermion(p[cb]);
	  this->threadedFreeFermion(mp[cb]); 
	  this->threadedFreeFermion(mmp[cb]);
	  this->threadedFreeFermion(r[cb]);
	}
	this->threadedFreeFermion(tmp);

	return k;
      }

      // New (conjugate/M-orthogonal) search direction
      b = cp/c;
      axpy(p,p,r,b); //p= b*p+r;
    }
    this->ThreadBossMessage("CGNE_unprec NOT converged r^2=%le rsq %le\n",cp,rsq);
    for(int cb=0;cb<2;cb++){
      this->threadedFreeFermion(p[cb]);
      this->threadedFreeFermion(mp[cb]); 
      this->threadedFreeFermion(mmp[cb]);
      this->threadedFreeFermion(r[cb]);
    }
    this->threadedFreeFermion(tmp);

    return -1;
    
  }

#endif
