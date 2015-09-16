#ifndef NESTED_4D_H
#define NESTED_4D_H

#include <bfm_qdp_chroma_linop.h>


///////////////////////////////////////
// General 4d effective operator class
///////////////////////////////////////

template <class Float> class bfmCayley4d : public bfmarg {
public:

  bfm_qdp<Float> dop;
  bfm_qdp<Float> dpv;
  multi1d<LatticeColorMatrix> u_qdp;
  g5dParams parms;

  bool isBoss(void ) { return dop.isBoss() ; };

  void init(g5dParams & solver, multi1d<LatticeColorMatrix> &U, double resid, int maxit);
  void end(void);

  double Munprec(const LatticeFermion &src, LatticeFermion &psi, int dag);
  double Munprec(Fermion_t chi[2],Fermion_t psi[2], int dag);

  void fourd_operator(LatticeFermion &src,LatticeFermion &d);
  void sign_approx   (LatticeFermion &src,LatticeFermion &eps);
  void deltaL        (LatticeFermion &src,LatticeFermion &delta);

  void mres(multi1d<Complex> &PPcorr, 
	    multi1d<Complex> &PAcorr, 
	    multi1d<Complex> &PJ5qcorr, 
	    LatticePropagator &prop,
	    LatticePropagator &src,
	    LatticePropagator &midpoint);


  int CGNE_MdagM(LatticeFermion &sol,const LatticeFermion &src);
  int CGNE_MdagM_shift(LatticeFermion &sol,const LatticeFermion &src, Real shift);

  int CGNE_Mdag (LatticeFermion &sol,const LatticeFermion &src);

  int CGNE_M    (LatticeFermion &sol,const LatticeFermion &src);
  int CGNE_M_shift(LatticeFermion &sol,const LatticeFermion &src, Real shift);

};

template<class Float> 
void bfmCayley4d<Float>::end(void)
{
  dop.end();
  dpv.end();
}

template<class Float> 
void bfmCayley4d<Float>::init(g5dParams & solver, multi1d<LatticeColorMatrix> &U, double resid, int maxit)
{
  parms=solver;
  mass=solver.mass;
  Ls = solver.Ls;
  M5 = solver.M5;
  residual=resid;
  max_iter=maxit;

  /**********************************
   * Set up BFM object
   **********************************/

  int lx = QDP::Layout::subgridLattSize()[0];
  int ly = QDP::Layout::subgridLattSize()[1];
  int lz = QDP::Layout::subgridLattSize()[2];
  int lt = QDP::Layout::subgridLattSize()[3];

  multi1d<int> procs = QDP::Layout::logicalSize();

  bfmarg bfma;

  //Physics parameters
  bfma.solver = solver.solver;
  bfma.Ls           = Ls;
  bfma.M5           = M5;
  bfma.mass         = mass;
  bfma.mobius_scale = solver.mobius_scale;
  bfma.zolo_lo      = solver.zolo_lo;
  bfma.zolo_hi      = solver.zolo_hi;
  if ( solver.solver == DWF ) {
    bfma.precon_5d    = 1;
  } else { 
    bfma.precon_5d    = 0;
  }
  bfma.max_iter     = max_iter;
  bfma.residual     = residual;
  //Geometry
  bfma.node_latt[0] = lx;
  bfma.node_latt[1] = ly;
  bfma.node_latt[2] = lz;
  bfma.node_latt[3] = lt;
  QDPIO::cout << "Nested 4d init: mass "<<mass<<endl;

  for(int mu=0;mu<4;mu++){
    if (procs[mu]>1) bfma.local_comm[mu] = 0;
    else             bfma.local_comm[mu] = 1;
  }

  u_qdp=U;
  dop.init(bfma);
  dop.importGauge(U);

  bfma.mass=1.0;
  dpv.init(bfma);
  dpv.importGauge(U);

  this->zolo_delta = dop.zolo_delta;
}

template <class Float> int bfmCayley4d<Float>::CGNE_M(LatticeFermion &psi,const LatticeFermion &src)
{
  LatticeFermion MdagSrc;
  Munprec(src,MdagSrc,DaggerYes);
  return CGNE_MdagM(psi,MdagSrc);
}
template <class Float> int bfmCayley4d<Float>::CGNE_M_shift(LatticeFermion &psi,const LatticeFermion &src, Real shift)
{
  LatticeFermion MdagSrc;
  Munprec(src,MdagSrc,DaggerYes);
  MdagSrc = MdagSrc - shift*src;
  return CGNE_MdagM_shift(psi,MdagSrc, shift);
}
template <class Float> int bfmCayley4d<Float>::CGNE_Mdag(LatticeFermion &psi,const LatticeFermion &src)
{
  LatticeFermion MinvPsi;
  int iter = CGNE_MdagM(MinvPsi,src);
  Munprec(MinvPsi,psi,DaggerNo);
  return iter;
}
template <class Float> int bfmCayley4d<Float>::CGNE_MdagM_shift(LatticeFermion &psi,const LatticeFermion &src, Real shift)
{
  LatticeFermion p;
  LatticeFermion mp;
  LatticeFermion mmp;
  LatticeFermion r;


  Munprec(psi,mp,DaggerNo);
  mp = mp - shift*psi;
  Munprec(mp,mmp,DaggerYes);
  mmp = mmp - shift*mp;

  r = src - mmp;
  p = r;

  Real ssq = norm2(src);
  Real rsq =  residual* residual*ssq;
  Real cp  = norm2(r);

  if ( toDouble(cp) <= toDouble(rsq) ) {
    QDPIO::cout <<"CGNE_unprec converged k=0 - nice guess "<<cp << " " << rsq<<endl;
    return 0;
  }  

  Real c,a,d,b;

  for (int k=1;k<=max_iter;k++){

    c=cp;
    Munprec(p,mp,DaggerNo);
    mp = mp - shift*p;
    d = norm2(mp);
    // Update solution vector
    a = c/d;
    psi = a*p+psi;
    // MdagMp
    Munprec(mp,mmp,DaggerYes); 
    mmp = mmp - shift*mp;

    r = r-a*mmp;
    cp= norm2(r);
    if ( verbose ) 
      QDPIO::cout <<"bfmCayley4d: CGNE_unprec k="<<k<<" r^2 = "<<cp<<endl;

    // Stopping condition
    if ( toDouble(cp) <= toDouble(rsq) ) { 
      QDPIO::cout <<"bfmCayley4d: CGNE_unprec converged k="<<k<<" "<<cp << " " << rsq<<endl;
      Munprec(psi,mp,DaggerNo);
      mp = mp - shift*psi;
      Munprec(mp,mmp,DaggerYes); 
      mmp = mmp - shift*mp;
      r=mmp-src;
      Real true_residual = sqrt(norm2(r)/norm2(src));
      QDPIO::cout <<"bfmCayley4d: CGNE_unprec true residual "<<true_residual<<endl;
      return k;
    }

    // New (conjugate/M-orthogonal) search direction
    b = cp/c;
    p = b*p+r;
  }
  QDPIO::cout <<"bfmCayley4d: CGNE_unprec not converged"<<endl;
  return -1;
}
template <class Float> int bfmCayley4d<Float>::CGNE_MdagM(LatticeFermion &psi,const LatticeFermion &src)
{
  LatticeFermion p;
  LatticeFermion mp;
  LatticeFermion mmp;
  LatticeFermion r;

  Munprec(psi,mp,DaggerNo);
  Munprec(mp,mmp,DaggerYes);

  r = src - mmp;
  p = r;

  Real ssq = norm2(src);
  Real rsq =  residual* residual*ssq;
  Real cp  = norm2(r);

  if ( toDouble(cp) <= toDouble(rsq) ) {
    QDPIO::cout <<"CGNE_unprec converged k=0 - nice guess "<<cp << " " << rsq<<endl;
    return 0;
  }  

  Real c,a,d,b;

  for (int k=1;k<=max_iter;k++){

    c=cp;

    d = Munprec(p,mp,DaggerNo);

    // Update solution vector
    a = c/d;
    psi = a*p+psi;
    
    // MdagMp
    Munprec(mp,mmp,DaggerYes); 

    r = r-a*mmp;
    cp= norm2(r);
    if ( verbose ) 
      QDPIO::cout <<"bfmCayley4d: CGNE_unprec k="<<k<<" r^2 = "<<cp<<endl;

    // Stopping condition
    if ( toDouble(cp) <= toDouble(rsq) ) { 
      QDPIO::cout <<"bfmCayley4d: CGNE_unprec converged k="<<k<<" "<<cp << " " << rsq<<endl;
      Munprec(psi,mp,DaggerNo);
      Munprec(mp,mmp,DaggerYes); 
      r=mmp-src;
      Real true_residual = sqrt(norm2(r)/norm2(src));
      QDPIO::cout <<"bfmCayley4d: CGNE_unprec true residual "<<true_residual<<endl;
      return k;
    }

    // New (conjugate/M-orthogonal) search direction
    b = cp/c;
    p = b*p+r;
  }
  QDPIO::cout <<"bfmCayley4d: CGNE_unprec not converged"<<endl;
  return -1;
}

template<class Float>
double bfmCayley4d<Float>::Munprec(const LatticeFermion &chi,
			       LatticeFermion &psi, 
			       int dag)
{
  multi1d<LatticeFermion> chi_5(Ls), psi_5(Ls);
  /**************************************************************
   * 4d -> 5d system
   **************************************************************/
  // Need chi_5[j] = P_[j][0] chi
  chi_5[0]   = chiralProjectMinus(chi);
  chi_5[Ls-1]= chiralProjectPlus(chi);

  // Zero bulk of chi and result psi
  for(int s=1;s<Ls-1;s++) chi_5[s]=zero;
  for(int s=0;s<Ls;s++)   psi_5[s]=zero;

  /********************
   * Bagel 
   ********************/
  Fermion_t tmp;
  Fermion_t s[2];
  Fermion_t Ms[2];
  Fermion_t MMs[2];

  tmp     = dop.allocFermion();
  for(int cb=0;cb<2;cb++){
    s[cb]   = dop.allocFermion();
    Ms[cb]  = dop.allocFermion();
    MMs[cb] = dop.allocFermion();
    
  }

  for(int cb=0;cb<2;cb++){
    dop.importFermion(chi_5,s[cb],cb);
    dop.importFermion(psi_5,Ms[cb],cb);
    dop.importFermion(psi_5,MMs[cb],cb);
  }
  QDPIO::cout << "bfmCayley4d: required residual = " <<  dop.residual << endl;

  // Test code
  Float f;

  if ( dag ) { 
    /**************************************************************
     * PV^{-d}
     **************************************************************/
    omp_set_num_threads(dop.nthread);
#pragma omp parallel 
    {
#pragma omp for 
      for(int i=0;i<dop.nthread;i++) {
	int mycount = dpv.CGNE_Mdag(Ms,s);
      }     
    } 
    
    /**************************************************************
     * D(m)
     **************************************************************/
#pragma omp parallel 
    {
#pragma omp for 
      for(int i=0;i<dop.nthread;i++) {
	dop.Munprec( Ms, MMs, tmp, DaggerYes) ;
      }
    }
  } else { 
    /**************************************************************
     * D(m)
     **************************************************************/
#pragma omp parallel 
    {
#pragma omp for 
      for(int i=0;i<dop.nthread;i++) {
	dop.Munprec( s, Ms, tmp, DaggerNo) ;
      }
    }

    /**************************************************************
     * PV^{-1}
     **************************************************************/
#pragma omp parallel 
    {
#pragma omp for 
      for(int i=0;i<dop.nthread;i++) {
	int mycount = dpv.CGNE_M(MMs,Ms);
      }
    }
  }

  /********************************************************
   * Merge the chiralities from the walls to reconstruct 4d prop
   ********************************************************/
  dop.freeFermion(tmp);
  for(int cb=0;cb<2;cb++){
    dop.exportFermion(psi_5,MMs[cb],cb);
  }
  for(int cb=0;cb<2;cb++){
    dop.freeFermion(s[cb]);
    dop.freeFermion(Ms[cb]);
    dop.freeFermion(MMs[cb]);
  }
  
  // Need psi = P^{inv}[0][j] psi_5[j]
  psi  = chiralProjectMinus(psi_5[0]);
  psi += chiralProjectPlus (psi_5[Ls-1]);
  
  double nrm =toDouble(norm2(psi));
  return nrm;
}

template<class Float>
double bfmCayley4d<Float>::Munprec(Fermion_t chi[2],
			       Fermion_t psi[2], 
			       int dag)
{
  /********************
   * Bagel 
   ********************/
  Fermion_t tmp;
  Fermion_t Ms[2];

  tmp     = dop.allocFermion();
  for(int cb=0;cb<2;cb++){
    Ms[cb]  = dop.allocFermion();
  }

  omp_set_num_threads(dpv.nthread);
  if ( dag ) { 
    /**************************************************************
     * PV^{-d}
     **************************************************************/
#pragma omp parallel 
    {
#pragma omp for 
      for(int i=0;i<dpv.nthread;i++) {
	int mycount = dpv.CGNE_Mdag(Ms,chi);
      }
    }    
    /**************************************************************
     * D(m)
     **************************************************************/
#pragma omp parallel 
    {
#pragma omp for 
      for(int i=0;i<dop.nthread;i++) {
	dop.Munprec( Ms, psi, tmp, DaggerYes) ;
      }
    }
  } else { 
    /**************************************************************
     * D(m)
     **************************************************************/
#pragma omp parallel 
    {
#pragma omp for 
      for(int i=0;i<dop.nthread;i++) {
	dop.Munprec( chi, Ms, tmp, DaggerNo) ;
      }
    }
    /**************************************************************
     * PV^{-1}
     **************************************************************/
#pragma omp parallel 
    {
#pragma omp for 
      for(int i=0;i<dpv.nthread;i++) {
	int mycount = dpv.CGNE_M(psi,Ms);
      }
    } 
  }
  dop.freeFermion(tmp);
  for(int cb=0;cb<2;cb++){dop.freeFermion(Ms[cb]);}

  return (Float)1.0;
}

////////////////////////////////////////////////////////////////////////
// 4d effective operator for general 5d case
////////////////////////////////////////////////////////////////////////
template<class Float>
void bfmCayley4d<Float>::fourd_operator (LatticeFermion &d,
				     LatticeFermion &src)
{
  int ok=1;
  if ( dop.solver == HwCayleyTanh ) ok=1;
  if ( dop.solver == HtCayleyTanh ) ok=1;
  if ( dop.solver == HmCayleyTanh ) ok=1;
  if ( dop.solver == DWF ) ok=1;
  if ( dop.solver == DWFrb4d ) ok=1;
  if ( !ok ) {
    printf("Bad solver type for bfmCayley4d::fourd_operator\n");
    exit(-1);
  }
  Munprec( src, d, DaggerNo) ;
}

////////////////////////////////////////////////////////////////////////
// delta Ls (Tanh)
////////////////////////////////////////////////////////////////////////
template <class Float>
void bfmCayley4d<Float>::deltaL(LatticeFermion &deltaLs, 
				LatticeFermion &src)
{
  LatticeFermion eps;
  LatticeFermion epseps;

  QDPIO::cout << "Evaluating sign function error"<<endl;
  if ( (this->mass !=0.0)||(this ->dop.mass !=0.0)|| (this ->dop.mass !=mass) ) {
    printf("bfmCayley4d::deltaL -- mass logic bomb %le %le \n",this->mass, this->dop.mass);
    exit(-1);
  }
  
  ////////////////////////////////
  // Here use the fact that dov = (1+m)/2 + (1-m)/2 g5 \epsilon(H)
  // => epsilon_H = g5 [ 2/(1-m) dov - (1+m)/(1-m)]
  ////////////////////////////////
  Gamma G5(15);

  fourd_operator( eps, src);
  eps= G5 * ( (2.0/(1.0-mass))*eps - src*((1.0+mass)/(1.0-mass))) ;


  fourd_operator( epseps, eps);
  epseps = G5 * ( (2.0/(1.0-mass))*epseps - eps*(1.0+mass)/(1.0-mass)) ;

  QDPIO::cout << "Delta_L src norm "<< norm2(src)<<endl;
  QDPIO::cout << "Delta_L eps norm "<< norm2(eps)<<endl;
  QDPIO::cout << "Delta_L eps^2 norm "<< norm2(epseps)<<endl;

#define ZOLO_SCALE
  double scale =1.0;
#ifdef ZOLO_SCALE
  if ( (parms.solver == HtCayleyZolo) || (parms.solver == HwCayleyZolo) ) { 
    scale = 1.0+this->zolo_delta;
    QDPIO::cout <<"Scaling by "<<scale<<" to make sign function error positive definite "<<endl;
  }
#endif
  deltaLs = (src-epseps/(scale*scale))*0.25;
  QDPIO::cout << "Delta_L delta_l norm "<< norm2(deltaLs)<<endl;

#undef REGRESS_DELTA
#ifdef REGRESS_DELTA
  Handle< LinearOperator<LatticeFermion> >  DelQDP = GetDeltaL(u_qdp,parms);
  LatticeFermion check;
  check=zero;
  (*DelQDP)(check,src,PLUS);
  QDPIO::cout<<"************************************************************"<<endl;
  QDPIO::cout<<"Difference from QDP DeltaL "<<norm2(check-deltaLs) << " qdp " <<norm2(check) << " bagel " << norm2(deltaLs) <<endl;
  QDPIO::cout<<"************************************************************"<<endl;
#endif

}

template<class Float>
void bfmCayley4d<Float>::mres(multi1d<Complex> &PPcorr, 
			  multi1d<Complex> &PAcorr, 
			  multi1d<Complex> &PJ5qcorr, 
			  LatticePropagator &prop,
			  LatticePropagator &src,
			  LatticePropagator &midpoint)
{
  LatticeFermion delta_L;  
  LatticeFermion chi;
  LatticeFermion schi;
  LatticeFermion mid;
  LatticeComplex PJ5q    = zero;
  LatticeComplex J5q     = zero;
  LatticeComplex midp    = zero;
  LatticeComplex PP      = zero;
  LatticeComplex PA      = zero;

  for(int color_source=0;color_source<3;color_source++){
  for(int spin_source=0;spin_source<4;spin_source++){

    BfmPropToFerm(prop, chi, color_source, spin_source);
    BfmPropToFerm(src , schi, color_source, spin_source);
    BfmPropToFerm(midpoint , mid, color_source, spin_source);

    Gamma G5(15);
    Gamma Gt(8);
    PP+=localInnerProduct(chi,chi);
    PA+=localInnerProduct(chi,Gt*chi);

    schi=chi*(1.0-mass)+schi; // Add contact term
    deltaL(delta_L,schi);

    J5q  =localInnerProduct(schi,delta_L);
    midp =localInnerProduct(mid,mid);
    QDPIO::cout << "J5q "<<norm2(J5q)<<endl;
    QDPIO::cout << "midp "<<norm2(midp)<<endl;
    QDPIO::cout << "diff "<<norm2(midp-J5q)<<endl;
    PJ5q+=J5q;

  }}

  Set &tslice = GetTimeslice();

  PJ5qcorr = sumMulti(PJ5q,tslice);  
  PPcorr   = sumMulti(PP,tslice);  
  PAcorr   = sumMulti(PA,tslice);  

  return;
}

template <class Float>
int bfm_g5d_DeltaL(LatticeFermion &sol, 
		   LatticeFermion &src,
		   multi1d<LatticeColorMatrix> &U,
		   g5dParams & solver,
		   Real resid,int maxit)
{
  g5dParams mysolver = solver;

  bfmCayley4d<Float> dop;
  
  mysolver.mass=0.0;
  double rr=toDouble(resid);

  dop.init(mysolver,U,rr,maxit);
  dop.deltaL(sol,src);
  dop.end();
}

#endif
