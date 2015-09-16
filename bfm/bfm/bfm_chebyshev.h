#ifndef _BFM_CHEBYSHEV_H_
#define _BFM_CHEBYSHEV_H_

class Chebyshev {
 public:
  int LeftPrecondition;
  int RightPrecondition;
  int order;
  double alpha;
  double hi,lo;
  double *coeffs;

  ~Chebyshev(){ delete[] coeffs; }
  
 Chebyshev(double _lo,double _hi,double _alpha,int _order) : lo(_lo),hi(_hi),order(_order), alpha(_alpha) {

    if( alpha == -0.5 ) {
      RightPrecondition=1;
      LeftPrecondition=1;
    } else if ( alpha == -1 ) {
      RightPrecondition=0;
      LeftPrecondition=1;
    } else if ( alpha == 0 ) {
      RightPrecondition=0;
      LeftPrecondition=0;
    } else { 
      exit(-1);
    }

    if(order < 2) exit(-1);
    coeffs = new double [order];
    for(int j=0;j<order;j++){
      double s=0;
      for(int k=0;k<order;k++){
	double y=cos(M_PI*(k+0.5)/order);
	double x=0.5*(y*(hi-lo)+(hi+lo));
	double f=pow(x,alpha);
	s=s+f*cos( j*M_PI*(k+0.5)/order );
      }
      coeffs[j] = s * 2.0/order;
    }
  }

  template<class Float> void PolyMdagMprecRight(bfmbase<Float> &dop,Fermion_t in,Fermion_t out);
  template<class Float> void PolyMdagMprecLeft (bfmbase<Float> &dop,Fermion_t in,Fermion_t out);
  template<class Float> void PolyMdagMprec(bfmbase<Float> &dop,Fermion_t in,Fermion_t out);
  template<class Float> void SolverMatrix(bfmbase<Float> & dop,Fermion_t in,Fermion_t out);
};

template<class Float>
int bfmbase<Float>::MCR_PolyPrec(Fermion_t psi, Fermion_t src_uprec) // Internal
{
  Chebyshev poly(0.1,150,-0.5,160);

  double f;
  double cp,c,a,d,b,ci;
  complex<double> num;
  int me = thread_barrier();

  if ( this->isBoss() && (!me) ) { 
    this->InverterEnter();
  }

  Fermion_t src = threadedAllocFermion(mem_fast); 
  Fermion_t p   = threadedAllocFermion(mem_fast); 
  Fermion_t tmp = threadedAllocFermion(mem_fast); 
  Fermion_t mp  = threadedAllocFermion(mem_fast); 
  Fermion_t r   = threadedAllocFermion(mem_fast); 
  Fermion_t Ap  = threadedAllocFermion(mem_fast); 
  Fermion_t Ar  = threadedAllocFermion(mem_fast); 

  //Initial residual computation & set up
  double guess = norm(psi);
  double junk;
  poly.SolverMatrix<double>(*this,psi,Ar);

  poly.PolyMdagMprecLeft<double>(*this,src_uprec,src);
  //  Mprec(psi,mp,tmp,DaggerNo);
  //  Mprec(mp,Ar,tmp,DaggerYes);

  axpy (r, Ar, src,-1.0);
  axpy (p, Ar, src,-1.0);

  a =norm(p);
  cp=norm(r);

  Float ssq =  norm(src);
  Float rsq =  residual* residual*ssq;

  ThreadBossMessage("MCR_PolyPrec guess %le \n",guess);
  ThreadBossMessage("MCR_PolyPrec ssq %le rsq %le\n",ssq,rsq);
  ThreadBossMessage("MCR_PolyPrec a %le cp %le\n",a,cp);
  

  //Check if guess is really REALLY good :)
  if ( cp <= rsq ) {
    ThreadBossMessage("MCR_PolyPrec k=0 converged - suspiciously nice guess %le %le\n",cp,rsq);

    threadedFreeFermion(src);
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


  struct timeval start,stop;
  gettimeofday(&start,NULL);


  //  c = Mprec(p,mp,tmp,DaggerNo);
  //  d = Mprec(mp,Ap,tmp,DaggerYes);  //Ap
  poly.SolverMatrix<Float>(*this,p,Ap);
  c = inner_real(p,Ap);
  d = norm(Ap);

  axpy(Ar,Ap,Ap,0);	       //Ar

  //c = real( dot(r,Ar) );
  for (int k=1;k<=max_iter;k++){

    this->iter=k;

    //c = real( dot(r,Ar) );		//c = rAr
    a = c/d;

    axpy(psi,p,psi,a);		//x = x + ap
    cp = axpy_norm(r,Ap,r,-a);		//r = r - aAp
  
    if ( k%10 == 0 ){
      ThreadBossMessage("MCR_PolyPrec: k= %d r^2= %le %le\n",k,cp,sqrt(cp/ssq));
    }

    // Stopping condition
    if ( cp <= rsq ) { 

	gettimeofday(&stop,NULL);
	struct timeval diff;
	timersub(&stop,&start,&diff);
	
	ThreadBossMessage("MCR_PolyPrec converged in %d iterations\n",k);
	ThreadBossMessage("MCR_PolyPrec converged in %d.%6.6d seconds\n",diff.tv_sec,diff.tv_usec);

	double flops = mprecFlops()*2.0 + 2.0*axpyNormFlops() + axpyFlops()*2.0;
	flops = flops * k;

	double t = diff.tv_sec*1.0E6 + diff.tv_usec;
	ThreadBossMessage("MCR_PolyPrec: %d mprec flops/site\n",mprecFlopsPerSite());
	ThreadBossMessage("MCR_PolyPrec: %le flops\n",flops);
	ThreadBossMessage("MCR_PolyPrec: %le mflops per node\n",flops/t);

	poly.PolyMdagMprecRight<double>(*this,psi,tmp);
	axpy(psi,tmp,tmp,0.0);
	Mprec(psi,mp,tmp,0);
	Mprec(mp,Ap,tmp,1); 
	axpy(tmp,src_uprec,Ap,-1.0);
	double true_residual = sqrt(norm(tmp)/norm(src));
	ThreadBossMessage("MCR_PolyPrec: true unprec residual is %le \n",true_residual);

	threadedFreeFermion(src);
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

    //    d = Mprec(r,mp,tmp,DaggerNo); //(r,Ar)
    //        Mprec(mp,Ar,tmp,DaggerYes);
    poly.SolverMatrix<Float>(*this,r,Ar);
    d = inner_real(r,Ar);
    b = d/c;
    c = d;
    axpy(p,p,r,b);
    d = axpy_norm(Ap,Ap,Ar,b);

  }
  
  ThreadBossMessage("MCR_PolyPrec: CG not converged \n");
  threadedFreeFermion(src);
  threadedFreeFermion(tmp);
  threadedFreeFermion(p);
  threadedFreeFermion(mp);
  threadedFreeFermion(r);
  threadedFreeFermion(Ap);
    threadedFreeFermion(Ar);
  if ( this->isBoss() && (!me) ) { 
    this->InverterExit();
  }

  return -1;

}

template<class Float> 
void Chebyshev::PolyMdagMprecRight(bfmbase<Float> &dop,Fermion_t in,Fermion_t out)
{
  int me = dop.thread_barrier();
  if ( RightPrecondition ) { 
    dop.ThreadBossMessage("Right preconditioning matrix x^%f\n",alpha);
    this->PolyMdagMprec(dop,in,out);
  } else { 
    dop.axpy(out,in,in,0.0); 
  }
}
template<class Float> 
void Chebyshev::PolyMdagMprecLeft(bfmbase<Float> &dop,Fermion_t in,Fermion_t out)
{
  int me = dop.thread_barrier();
  if ( LeftPrecondition ) { 
    dop.ThreadBossMessage("Left preconditioning matrix x^%f\n",alpha);
    this->PolyMdagMprec(dop,in,out);
  } else { 
    dop.axpy(out,in,in,0.0); 
  }
}


template<class Float> 
void Chebyshev::PolyMdagMprec(bfmbase<Float> &dop,Fermion_t in,Fermion_t out)
{
    Fermion_t y=    dop.threadedAllocFermion();
    Fermion_t Mtmp= dop.threadedAllocFermion();
    Fermion_t Tnm = dop.threadedAllocFermion();
    Fermion_t Tn  = dop.threadedAllocFermion();
    Fermion_t Tnp = dop.threadedAllocFermion();
    Fermion_t tmp=Tnp;

    int me = dop.thread_barrier();

    dop.ThreadBossMessage("Chebyshev::PolyMdagMprec Evaluating polynomial for x^%f\n",alpha);
    
    double xscale = 2.0/(hi-lo);
    double mscale = -(hi+lo)/(hi-lo);
    
    dop.axpy(Tnm,in,in,0.0);              // Tnm=T0=in

    dop.Mprec(in ,tmp,Mtmp,DaggerNo);  
    dop.Mprec(tmp,y  ,Mtmp,DaggerYes);    
    dop.axpby(Tn ,y,in,xscale,mscale);    // TN=T1 = (xscale MdagM + mscale)in 

    dop.axpby(out,Tnm,Tn,coeffs[0]/2.0,coeffs[1]); // sum = .5c[0]T0 + c[1]T1

    double oo;
    double TTnp;

    double xx=0.4;

    double TTnm= 1.0;

    double yy= (xscale*xx + mscale);
    double TTn = yy;

    oo = coeffs[0]*TTnm/2.0 + coeffs[1] * TTn;
    
    for(int i=2;i<order;i++){
      
      dop.Mprec(Tn,tmp,Mtmp,DaggerNo);  
      dop.Mprec(tmp, y,Mtmp,DaggerYes);    
      dop.axpby(y,y,Tn,xscale,mscale);    // y Tn = [xscale MdagM+mscale] Tn
      dop.axpby(Tnp,y,Tnm,2.0,-1.0);      // Tnp=2yTn - Tnm

      TTnp=2*yy*TTn-TTnm;
      dop.axpy(Tnm,Tn,Tn,0.0); 
      dop.axpy(Tn,Tnp,Tnp,0.0);
      dop.axpy(out,Tn,out,coeffs[i]);//Accumulate

      TTnm=TTn;
      TTn =TTnp;
      oo = oo + coeffs[i]*TTn;
	
    }
    dop.ThreadBossMessage(("Summed chebyshev %le %le\n",pow(xx,alpha),oo);
    dop.threadedFreeFermion(y);
    dop.threadedFreeFermion(Mtmp);
    dop.threadedFreeFermion(Tnm);
    dop.threadedFreeFermion(Tn);
    dop.threadedFreeFermion(Tnp);
};

template<class Float> void Chebyshev::SolverMatrix(bfmbase<Float> & dop,Fermion_t in,Fermion_t out)
{
    
  Fermion_t tmp = dop.threadedAllocFermion();
  Fermion_t Mtmp= dop.threadedAllocFermion();
  
  this->PolyMdagMprecRight(dop,in,tmp);
  
  dop.Mprec(tmp,out,Mtmp,DaggerNo);  
  dop.Mprec(out,tmp,Mtmp,DaggerYes);    

  this->PolyMdagMprecLeft(dop,tmp,out);
  
  // Can turn this into P_L M if inexact deflation
  // P_L M = [ 1  -Mbs Mss^-1] 
  //         [ 0   0         ] 
  // 
  // i) Project tmp to subspace.
  // ii) Promote and subtract to get orthog piece
  // iii) Apply Mss^{-1}
  // iv)  Promote and apply full M
  // v)   Project, Promote,Subtract
  // Add to "1" part.
  
    
  dop.threadedFreeFermion(tmp);
  dop.threadedFreeFermion(Mtmp);
  
};


#endif
