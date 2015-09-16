#ifndef DWF4D_H
#define DWF4D_H

#include <omp.h>
#include "Matrix.h"
#include <cstdlib>
#include <string>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <fstream>
#include <complex>
#include <sys/time.h>

/*
void Gram(LatticeFermion &r, multi1d<LatticeFermion> q, int N){
  for(int i=0;i<N;i++){
    r = r   -  ( innerProduct(q[i],r)/innerProduct(q[i],q[i]) )*q[i];
  }
  r = r*( 1.0/sqrt( norm2(r) ) );
}

void Gram(LatticeFermion &r, multi1d<LatticeFermion> q){
  int N = q.size();
  for(int i=0;i<N;i++){
    r = r   -  ( innerProduct(q[i],r)/innerProduct(q[i],q[i]) )*q[i];
  }
  r = r*( 1.0/sqrt( norm2(r) ) );
}

void Gram(multi1d<LatticeFermion> &q){
  int M = q.size();
  for(int j=0;j<M;j++){
    for(int i=0;i<j;i++){
      q[j] = q[j]   -  ( innerProduct(q[i],q[j])/innerProduct(q[i],q[i]) )*q[i];
    }
    q[j] = q[j]*( 1.0/sqrt( norm2(q[j]) ) );
  }
}

void Gram(vector<LatticeFermion> &q){
  int M = q.size();
  for(int j=0;j<M;j++){
    for(int i=0;i<j;i++){
      q[j] = q[j]   -  ( innerProduct(q[i],q[j])/innerProduct(q[i],q[i]) )*q[i];
    }
    q[j] = q[j]*( 1.0/sqrt( norm2(q[j]) ) );
  }
}


void Gram(multi1d<LatticeFermion> &q, int M){
  for(int j=0;j<M;j++){
    for(int i=0;i<j;i++){
      q[j] = q[j]   -  ( innerProduct(q[i],q[j])/innerProduct(q[i],q[i]) )*q[i];
    }
    q[j] = q[j]*( 1.0/sqrt( norm2(q[j]) ) );
  }
}

void Gram(vector<LatticeFermion> &q, int M){
  for(int j=0;j<M;j++){
    for(int i=0;i<j;i++){
      q[j] = q[j]   -  ( innerProduct(q[i],q[j])/innerProduct(q[i],q[i]) )*q[i];
    }
    q[j] = q[j]*( 1.0/sqrt( norm2(q[j]) ) );
  } 
}

/// q -> q Q
template <class T> void times(multi1d<LatticeFermion> &q, Matrix<T> &Q){
	multi1d<LatticeFermion> S( q.size() );
	for(int j=0;j<q.size();j++){
		S[j] = zero;
		for(int k=0;k<q.size();k++){
			S[j] = S[j] + q[k] * toComplex( Q(k,j) ); 
		}
	}
	for(int j=0;j<q.size();j++){q[j] = S[j];}
}



/// q -> q Q
template <class T> void times(multi1d<LatticeFermion> &q, Matrix<T> &Q, int N){
	multi1d<LatticeFermion> S( N );
	for(int j=0;j<N;j++){
		S[j] = zero;
		for(int k=0;k<N;k++){
			S[j] = S[j] + q[k] * toComplex( Q(k,j) ); 
		}
	}
	for(int j=0;j<N;j++){q[j] = S[j];}
}
*/

/** Print a 4d Lattice Fermion **/
/*
void print_full_vector(LatticeFermion in, int lx, int ly, int lz, int lt){

  multi1d<int> site(4);	
  site[0] = 0;
  site[1] = 0;
  site[2] = 0;
  site[3] = 0;
  int count = 0;  
  for(int i0=0;i0<lx;i0++){
  for(int i1=0;i1<ly;i1++){
  for(int i2=0;i2<lz;i2++){
  for(int i3=0;i3<lt;i3++){
    site[0] = i0;site[1] = i1;site[2] = i2;site[3] = i3;
  
    QDPIO::cout << "(" << i0 << ","<< i1 << ","<< i2 << ","<< i3 << ")" << endl; 
    Fermion F = peekSite(in,site);
    print_fermion2(F);
    count++;	
  }}}}
  QDPIO::cout << endl;
}*/


/*
	4d bfm fermion class
	Removed some functionality in anticipation of bfm changes	
				*/

template <class T> class mybfmDwf4d : public bfmarg {
public:
  int Ls;
  T mass;
  int countit;
  bfm_qdp<T> dop;
  //bfm_qdp<T> dop_tw;
  
  bool isBoss(void ) { return dop.isBoss() ; };

  void end(void);
  void init(T m, T M5, int N5, multi1d<LatticeColorMatrix> U, double resid);
  //void init_HT(T m, T M5, int N5, multi1d<LatticeColorMatrix> U, double resid);

  T Munprec(const LatticeFermion &chi,
		 LatticeFermion &psi, 
		 int dag);

  T Munprec(Fermion_t chi[2],
			 Fermion_t psi[2], 
			 int dag);

  /*T Munprec_TW(Fermion_t chi[2],
			 Fermion_t psi[2], 
			 int dag);
			 
  T Munprec_HT(Fermion_t chi[2],
			Fermion_t psi[2]);
			
  T Munprec_TW_QMR(Fermion_t chi[2],
					  Fermion_t psi[2], 
					  int dag);*/
    
  int CGNE_MdagM(LatticeFermion &sol,const LatticeFermion &src);
  int CGNE_Mdag (LatticeFermion &sol,const LatticeFermion &src);
  int CGNE_M    (LatticeFermion &sol,const LatticeFermion &src);

};
template <class T> void mybfmDwf4d<T>::end(void)
{
  dop.end();
}

template <class T> void mybfmDwf4d<T>::init(T m, T M5, int N5, multi1d<LatticeColorMatrix> U, double resid)
{
  mass=m;
  Ls  =N5;
  countit = 0;
  /**********************************
   * Set up BFM object
   **********************************/

  int lx = QDP::Layout::subgridLattSize()[0];
  int ly = QDP::Layout::subgridLattSize()[1];
  int lz = QDP::Layout::subgridLattSize()[2];
  int lt = QDP::Layout::subgridLattSize()[3];

  bfmarg dwfa;
  //Physics parameters
  dwfa.solver       = DWF;
  dwfa.Ls           = Ls;
  dwfa.M5           = M5;
  dwfa.mass         = mass;
  dwfa.precon_5d    = 1;
  dwfa.max_iter     = 10000;
  #ifdef BLEUGENE
  dwfa.time_report_iter = 10000; 
  #endif
  dwfa.residual     = resid;
  //Geometry
  dwfa.node_latt[0] = lx;
  dwfa.node_latt[1] = ly;
  dwfa.node_latt[2] = lz;
  dwfa.node_latt[3] = lt;

  multi1d<int> procs = QDP::Layout::logicalSize();
  if ( this -> isBoss() ) {
    printf("dwf_4d: %d dim machine\n\t", procs.size());
  }
  
QDPIO::cout << "dwf_4d: required residual = " <<  dwfa.residual << endl;
  
  
  for(int mu=0;mu<4;mu++){
    if ( procs[mu]>1 ) {
      dwfa.local_comm[mu] = 0;
    } else { 
      dwfa.local_comm[mu] = 1;
    }
  }
  if ( this -> isBoss() ) {
    printf("dwf_4d: Local comm = ");
    for(int mu=0;mu<4;mu++){
      printf("%d ", dwfa.local_comm[mu]);
    }
    printf("\n");
  }
  
  multi1d<int> ncoor = QDP::Layout::nodeCoord();

  dop.init(dwfa);
  dop.importGauge(U);
  
  /*dwfa.Ls = 1;
  dwfa.precon_5d = 0;
  dwfa.mass = 1.0;
  dop_tw.init(dwfa);
  dop_tw.importGauge(U);*/

}
/*
template <class T> void mybfmDwf4d<T>::init_HT(T m, T M5, int N5, multi1d<LatticeColorMatrix> U, double resid)
{
  QDPIO::cout << "Init Ht" << endl;
  mass=m;
  Ls  =1;
  countit = 0;

   // Set up BFM object


  int lx = QDP::Layout::subgridLattSize()[0];
  int ly = QDP::Layout::subgridLattSize()[1];
  int lz = QDP::Layout::subgridLattSize()[2];
  int lt = QDP::Layout::subgridLattSize()[3];

  bfmarg dwfa;
  QDPIO::cout << "mass = " << dwfa.mass << " " << mass << endl;
  //Physics parameters
  dwfa.Ls           = Ls;
  dwfa.M5           = 1.0 - mass;
  dwfa.mass         = mass;
  dwfa.precon_5d    = 0;
  dwfa.max_iter     = 10000;
  #ifdef BLEUGENE
  dwfa.time_report_iter = 10000; 
  #endif
  dwfa.residual     = resid;
  //Geometry
  dwfa.node_latt[0] = lx;
  dwfa.node_latt[1] = ly;
  dwfa.node_latt[2] = lz;
  dwfa.node_latt[3] = lt;

  multi1d<int> procs = QDP::Layout::logicalSize();
  if ( this -> isBoss() ) {
    printf("dwf_4d: %d dim machine\n\t", procs.size());
  }
  
QDPIO::cout << "dwf_4d: required residual = " <<  dwfa.residual << endl;
  
  
  for(int mu=0;mu<4;mu++){
    if ( procs[mu]>1 ) {
      dwfa.local_comm[mu] = 0;
    } else { 
      dwfa.local_comm[mu] = 1;
    }
  }
  if ( this -> isBoss() ) {
    printf("dwf_4d: Local comm = ");
    for(int mu=0;mu<4;mu++){
      printf("%d ", dwfa.local_comm[mu]);
    }
    printf("\n");
  }
  
  multi1d<int> ncoor = QDP::Layout::nodeCoord();

  dop.init(dwfa);
  dop.importGauge(U);
  
  dwfa.M5           = 1.0 - (2.0 + mass);
  dop_tw.init(dwfa);
  dop_tw.importGauge(U);

}*/

template <class T> int mybfmDwf4d<T>::CGNE_M(LatticeFermion &psi,const LatticeFermion &src)
{
  LatticeFermion MdagSrc;
  Munprec(src,MdagSrc,DaggerYes);
  return CGNE_MdagM(psi,MdagSrc);
}
template <class T> int mybfmDwf4d<T>::CGNE_Mdag(LatticeFermion &psi,const LatticeFermion &src)
{
  LatticeFermion MinvPsi;
  int iter = CGNE_MdagM(MinvPsi,src);
  Munprec(MinvPsi,psi,DaggerNo);
  return iter;
}
template <class T> int mybfmDwf4d<T>::CGNE_MdagM(LatticeFermion &psi,const LatticeFermion &src)
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
      QDPIO::cout <<"dwf_4d: CGNE_unprec k="<<k<<" r^2 = "<<cp<<endl;

    // Stopping condition
    if ( toDouble(cp) <= toDouble(rsq) ) { 
      QDPIO::cout <<"dwf_4d: CGNE_unprec converged k="<<k<<" "<<cp << " " << rsq<<endl;
      Munprec(psi,mp,DaggerNo);
      Munprec(mp,mmp,DaggerYes); 
      r=mmp-src;
      Real true_residual = sqrt(norm2(r)/norm2(src));
      QDPIO::cout <<"dwf_4d: CGNE_unprec true residual "<<true_residual<<endl;
      return k;
    }

    // New (conjugate/M-orthogonal) search direction
    b = cp/c;
    p = b*p+r;
  }
  QDPIO::cout <<"dwf_4d: CGNE_unprec not converged"<<endl;
  return -1;
}
template <class T> T mybfmDwf4d<T>::Munprec(const LatticeFermion &chi,
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
#pragma omp parallel 
  {
    omp_set_num_threads(dop.nthread);
#pragma omp for 
    for(int i=0;i<dop.nthread;i++) {
  tmp     = dop.allocFermion();
  for(int cb=0;cb<2;cb++){
    s[cb]   = dop.allocFermion();
    Ms[cb]  = dop.allocFermion();
    MMs[cb] = dop.allocFermion();
    
  }
}}
  for(int cb=0;cb<2;cb++){
    dop.importFermion(chi_5,s[cb],cb);
    dop.importFermion(psi_5,Ms[cb],cb);
    dop.importFermion(psi_5,MMs[cb],cb);
  }
  QDPIO::cout << "dwf_4d: required residual = " <<  dop.residual << endl;

  // Test code
  T f;
  int mycountit;

    if ( dag ) { 
      /**************************************************************
       * PV^{-d}
       **************************************************************/
      dop.mass = 1.0;
#pragma omp parallel 
  {
    omp_set_num_threads(dop.nthread);
#pragma omp for 
    for(int i=0;i<dop.nthread;i++) {
	int mycount = dop.CGNE_Mdag(Ms,s);
	mycountit=mycount;
      }    }  
      countit+=mycountit;

      /**************************************************************
       * D(m)
       **************************************************************/
      dop.mass = mass;
#pragma omp parallel 
  {
    omp_set_num_threads(dop.nthread);
#pragma omp for 
    for(int i=0;i<dop.nthread;i++) {
	dop.Munprec( Ms, MMs, tmp, DaggerYes) ;
      }}
      countit++;
    } else { 
      /**************************************************************
       * D(m)
       **************************************************************/
      dop.mass = mass;
#pragma omp parallel 
  {
    omp_set_num_threads(dop.nthread);
#pragma omp for 
    for(int i=0;i<dop.nthread;i++) {
	dop.Munprec( s, Ms, tmp, DaggerNo) ;
      }}
      countit++;
      
      /**************************************************************
       * PV^{-1}
       **************************************************************/
      dop.mass = 1.0;
#pragma omp parallel 
  {
    omp_set_num_threads(dop.nthread);
#pragma omp for 
    for(int i=0;i<dop.nthread;i++) {
	int mycount = dop.CGNE_M(MMs,Ms);
	mycountit=mycount;
      }}
      countit+=mycountit;
    }

    
  /********************************************************
   * Merge the chiralities from the walls to reconstruct 4d prop
   ********************************************************/
#pragma omp parallel 
  {
    omp_set_num_threads(dop.nthread);
#pragma omp for 
    for(int i=0;i<dop.nthread;i++) {
  dop.freeFermion(tmp);
  for(int cb=0;cb<2;cb++){
    dop.exportFermion(psi_5,MMs[cb],cb);
  }
  for(int cb=0;cb<2;cb++){
    dop.freeFermion(s[cb]);
    dop.freeFermion(Ms[cb]);
    dop.freeFermion(MMs[cb]);
  }
}}  
  // Need psi = P^{inv}[0][j] psi_5[j]
  psi = chiralProjectMinus(psi_5[0]);
  psi += chiralProjectPlus (psi_5[Ls-1]);
  
  T nrm =toDouble(norm2(psi));
  return nrm;
}


template <class T> T mybfmDwf4d<T>::Munprec(Fermion_t chi[2],
					  Fermion_t psi[2], 
					  int dag)
{

  /********************
   * Bagel 
   ********************/
  Fermion_t tmp;
  Fermion_t Ms[2];
#pragma omp parallel 
  {
    omp_set_num_threads(dop.nthread);
#pragma omp for 
    for(int i=0;i<dop.nthread;i++) {
  tmp     = dop.allocFermion();
  for(int cb=0;cb<2;cb++){
    Ms[cb]  = dop.allocFermion();
  }
}}
  // Test code
  T f;

  int mycountit;

  if ( dag ) { 
    /**************************************************************
     * PV^{-d}
     **************************************************************/
    if ( dop.isBoss() ) {
      //printf("dwf_4d thread %d\n",omp_get_thread_num());  //TODO BFM VERSION MISMATCH
      fflush(stdout);
    }
    dop.mass = 1.0;
#pragma omp parallel 
  {
    omp_set_num_threads(dop.nthread);

#pragma omp for 
    for(int i=0;i<dop.nthread;i++) {
      int mycount = dop.CGNE_Mdag(Ms,chi);
      mycountit=mycount;
    }    }
    countit+=mycountit;
    /**************************************************************
     * D(m)
     **************************************************************/
    dop.mass = mass;
#pragma omp parallel 
  {
    omp_set_num_threads(dop.nthread);

#pragma omp for 
    for(int i=0;i<dop.nthread;i++) {
      dop.Munprec( Ms, psi, tmp, DaggerYes) ;
    }}
    countit++;
  } else { 
    /**************************************************************
     * D(m)
     **************************************************************/
    dop.mass = mass;
#pragma omp parallel 
  {
    omp_set_num_threads(dop.nthread);

#pragma omp for 
    for(int i=0;i<dop.nthread;i++) {
      dop.Munprec( chi, Ms, tmp, DaggerNo) ;
    }}
    countit++;
    /**************************************************************
     * PV^{-1}
     **************************************************************/
    dop.mass = 1.0;
#pragma omp parallel 
  {
    omp_set_num_threads(dop.nthread);

#pragma omp for 
    for(int i=0;i<dop.nthread;i++) {
      int mycount = dop.CGNE_M(psi,Ms);
      mycountit=mycount;
    } }
    countit+=mycountit;
  }


  /********************************************************
   * Merge the chiralities from the walls to reconstruct 4d prop
   ********************************************************/
#pragma omp parallel 
  {
    omp_set_num_threads(dop.nthread);

#pragma omp for 
    for(int i=0;i<dop.nthread;i++) {
  dop.freeFermion(tmp);
  for(int cb=0;cb<2;cb++){dop.freeFermion(Ms[cb]);}
}}

  return (T)1.0;
}

/*
template <class T> T mybfmDwf4d<T>::Munprec_TW(Fermion_t chi[2],
					  Fermion_t psi[2], 
					  int dag)
{


   // Bagel 

  Fermion_t tmp;
  Fermion_t Ms[2], src[2], sol[2];

  multi1d<LatticeFermion> psi_5(Ls), chi_5(Ls),  tst_5(Ls);
  
  tmp     = dop.allocFermion();
  for(int cb=0;cb<2;cb++){
    Ms[cb]  = dop.allocFermion();
    src[cb] = dop_tw.allocFermion();
    sol[cb] = dop_tw.allocFermion();
  }
  struct timeval start,stop;
  gettimeofday(&start,NULL);  
    
    
  // Test code
  T f;

  QDPIO::cout << "dwf_4d: spawning " <<  NUM_BFM_THREADS << " threads" << endl;
  int mycountit;

  if ( dag ) {

     // PV^{-d}
    if ( dop.isBoss() ) {
      //printf("dwf_4d thread %d\n",omp_get_thread_num());  //TODO BFM VERSION MISMATCH
      fflush(stdout);
    }
    dop.mass = 1.0;
#pragma omp parallel num_threads NUM_BFM_THREADS
        {
      dop.exportFermion(chi_5, chi[0], 0); dop.exportFermion( chi_5, chi[1], 1);
      dop.multU(chi_5, psi_5,  1);      
      int mycount = 0;
      for(int s = 0;s<Ls; s++){
	  dop_tw.importFermion(psi_5[s], src[0], 0); dop_tw.importFermion( psi_5[s], src[1], 1);
	  mycount += dop_tw.CGNE_PV_Mdag(sol, src, s, Ls);
	  dop_tw.exportFermion(psi_5[s], sol[0], 0); dop_tw.exportFermion( psi_5[s], sol[1], 1);
      }
      dop.multU(psi_5,chi_5,  0);
      dop.importFermion(chi_5, Ms[0],0); dop.importFermion(chi_5, Ms[1],1);
          QDPIO::cout << "Converged in " << mycount << " wilson iterations ~~~ " << mycount / Ls << " DWF mults" << endl; 
	  mycountit = mycount / Ls;
    } 
    countit+=mycountit;

     // D(m)

    dop.mass = mass;
#pragma omp parallel num_threads NUM_BFM_THREADS
    {
      dop.Munprec( Ms, psi, tmp, DaggerYes) ;
    }
    countit++;
  } else { 

     // D(m)
    dop.mass = mass;
#pragma omp parallel num_threads NUM_BFM_THREADS
    {
      dop.Munprec( chi, Ms, tmp, DaggerNo) ;
    }
    countit++;

     // PV^{-1}
    dop.mass = 1.0;
#pragma omp parallel num_threads NUM_BFM_THREADS
    {
      dop.exportFermion(chi_5, Ms[0], 0); dop.exportFermion( chi_5, Ms[1], 1);
      QDPIO::cout << "chi_5 " <<  norm2(chi_5) << endl;
      QDPIO::cout << "psi_5 " <<  norm2(psi_5) << endl;
      
      dop.multU(chi_5, psi_5,  1);      
      QDPIO::cout << "chi_5 " <<  norm2(chi_5) << endl;
      QDPIO::cout << "psi_5 " <<  norm2(psi_5) << endl;
      int mycount = 0;
      for(int s = 0;s<Ls; s++){
	  dop_tw.importFermion(psi_5[s], src[0], 0); dop_tw.importFermion( psi_5[s], src[1], 1);
	  mycount += dop_tw.CGNE_PV_M(sol, src, s, Ls);
	  dop_tw.exportFermion(psi_5[s], sol[0], 0); dop_tw.exportFermion( psi_5[s], sol[1], 1);
	  QDPIO::cout << "s = " << s << " chi_5 " <<  norm2(chi_5[s]) << " psi_5 " <<  norm2(psi_5[s]) << endl;
      }
      QDPIO::cout << "chi_5 " <<  norm2(chi_5) << endl;
      QDPIO::cout << "psi_5 " <<  norm2(psi_5) << endl;
      dop.multU(psi_5,chi_5,  0);
      
      QDPIO::cout << "chi_5 " <<  norm2(chi_5) << endl;
      QDPIO::cout << "psi_5 " <<  norm2(psi_5) << endl;
      dop.importFermion(chi_5, psi[0],0); dop.importFermion(chi_5, psi[1],1);
          QDPIO::cout << "Converged in " << mycount << " wilson iterations ~~~ " << mycount / Ls << " DWF mults" << endl; 
	  mycountit = mycount / Ls;
    } 
    countit+=mycountit;
  }



   // Merge the chiralities from the walls to reconstruct 4d prop
  dop.freeFermion(tmp);
  for(int cb=0;cb<2;cb++){
    dop.freeFermion(Ms[cb]);
    dop_tw.freeFermion(src[cb]);
    dop_tw.freeFermion(sol[cb]);
  }
  
  gettimeofday(&stop,NULL);
  struct timeval diff;
  timersub(&stop,&start,&diff);
  QDPIO::cout << "Ls CGNE 4D converged in " << diff.tv_sec << "." << diff.tv_usec << "s" << endl;
  
  return (T)1.0;
}

*/

/*
template <class T> T mybfmDwf4d<T>::Munprec_HT(Fermion_t chi[2],
					  Fermion_t psi[2])
{
  Fermion_t src[2], sol[2], tmp, Mt[2];  
  for(int cb=0;cb<2;cb++){
    src[cb] = dop_tw.allocFermion();
    sol[cb] = dop_tw.allocFermion();
    Mt[cb] = dop_tw.allocFermion();
  }
  tmp = dop_tw.allocFermion();
  
  ///Copy input
  for(int cb=0;cb<2;cb++){
    dop.axpy( src[cb], chi[cb], chi[cb], 0.0);
  }
  struct timeval start,stop;
  gettimeofday(&start,NULL);

      
  // Test code
  QDPIO::cout << "dwf_4d: spawning " <<  NUM_BFM_THREADS << " threads" << endl;

     //				 HT(m)				*
    #pragma omp parallel num_threads NUM_BFM_THREADS 
    //dop.Munprec(src,Mt,tmp,false);  // (Dw) src
    //for(int cb=0;cb<2;cb++) {dop.axpy( Mt[cb], src[cb], Mt[cb], 2.0);}
    countit += dop_tw.CGNE_M(sol, src); // (2 + Dw)^-1 src
    dop.Munprec(sol,Mt,tmp,false);
    dop.chiralproj_axpby(Mt[0], Mt[0], -2.0 , -1,1);
    dop.chiralproj_axpby(Mt[1], Mt[1], -2.0 , -1,1);  
    
    for(int cb=0;cb<2;cb++){
    dop.axpy( psi[cb], Mt[cb], Mt[cb], 0.0);
    }
    countit++;
    

   // Merge the chiralities from the walls to reconstruct 4d prop

  dop_tw.freeFermion(tmp);
  for(int cb=0;cb<2;cb++){
    dop_tw.freeFermion(Mt[cb]);
    dop_tw.freeFermion(src[cb]);
    dop_tw.freeFermion(sol[cb]);
  }
  gettimeofday(&stop,NULL);
  struct timeval diff;
  timersub(&stop,&start,&diff);
  QDPIO::cout << "Ls CGNE 4D converged in " << diff.tv_sec << "." << diff.tv_usec << "s" << endl;
  
  return (T)1.0;
}

*/


///Takes same amount of time really...
/*
template <class T> T mybfmDwf4d<T>::Munprec_TW_QMR(Fermion_t chi[2],
					  Fermion_t psi[2], 
					  int dag)
{

   // Bagel 
  Fermion_t tmp;
  Fermion_t Ms[2], src[2], sol[2];

  multi1d<LatticeFermion> psi_5(Ls), chi_5(Ls),  tst_5(Ls);
  
  tmp     = dop.allocFermion();
  for(int cb=0;cb<2;cb++){
    Ms[cb]  = dop.allocFermion();
    src[cb] = dop_tw.allocFermion();
    sol[cb] = dop_tw.allocFermion();
  }
  struct timeval start,stop;
  gettimeofday(&start,NULL);

      
  // Test code
  T f;

  QDPIO::cout << "dwf_4d: spawning " <<  NUM_BFM_THREADS << " threads" << endl;
  int mycountit;

  if ( dag ) {

     // PV^{-d}
    if ( dop.isBoss() ) {
      //printf("dwf_4d thread %d\n",omp_get_thread_num());  //TODO BFM VERSION MISMATCH
      fflush(stdout);
    }
    dop.mass = 1.0;
#pragma omp parallel num_threads NUM_BFM_THREADS
        {
      dop.exportFermion(chi_5, chi[0], 0); dop.exportFermion( chi_5, chi[1], 1);
      dop.multU(chi_5, psi_5,  1);      
      int mycount = 0;
      for(int s = 0;s<Ls; s++){
	  dop_tw.importFermion(psi_5[s], src[0], 0); dop_tw.importFermion( psi_5[s], src[1], 1);
	  mycount += dop_tw.QMR_PV_Mdag(sol, src, s, Ls);
	  dop_tw.exportFermion(psi_5[s], sol[0], 0); dop_tw.exportFermion( psi_5[s], sol[1], 1);
      }
      dop.multU(psi_5,chi_5,  0);
      dop.importFermion(chi_5, Ms[0],0); dop.importFermion(chi_5, Ms[1],1);
          QDPIO::cout << "Converged in " << mycount << " wilson iterations ~~~ " << mycount / Ls << " DWF mults" << endl; 
	  mycountit = mycount / Ls;
    } 
    countit+=mycountit;

     // D(m)
    dop.mass = mass;
#pragma omp parallel num_threads NUM_BFM_THREADS
    {
      dop.Munprec( Ms, psi, tmp, DaggerYes) ;
    }
    countit++;
  } else { 

     // D(m)
    dop.mass = mass;
#pragma omp parallel num_threads NUM_BFM_THREADS
    {
      dop.Munprec( chi, Ms, tmp, DaggerNo) ;
    }
    countit++;

     // PV^{-1}
    dop.mass = 1.0;
#pragma omp parallel num_threads NUM_BFM_THREADS
    {
      dop.exportFermion(chi_5, Ms[0], 0); dop.exportFermion( chi_5, Ms[1], 1);
      dop.multU(chi_5, psi_5,  1);      
      int mycount = 0;
      for(int s = 0;s<Ls; s++){
	  dop_tw.importFermion(psi_5[s], src[0], 0); dop_tw.importFermion( psi_5[s], src[1], 1);
	  mycount += dop_tw.QMR_PV_M(sol, src, s, Ls);
	  dop_tw.exportFermion(psi_5[s], sol[0], 0); dop_tw.exportFermion( psi_5[s], sol[1], 1);
      }
      dop.multU(psi_5,chi_5,  0);
      dop.importFermion(chi_5, psi[0],0); dop.importFermion(chi_5, psi[1],1);
          QDPIO::cout << "Converged in " << mycount << " wilson iterations ~~~ " << mycount / Ls << " DWF mults" << endl; 
	  mycountit = mycount / Ls;
    } 
    countit+=mycountit;
  }



   // Merge the chiralities from the walls to reconstruct 4d prop
  dop.freeFermion(tmp);
  for(int cb=0;cb<2;cb++){
    dop.freeFermion(Ms[cb]);
    dop_tw.freeFermion(src[cb]);
    dop_tw.freeFermion(sol[cb]);
  }
  
  gettimeofday(&stop,NULL);
  struct timeval diff;
  timersub(&stop,&start,&diff);
  QDPIO::cout << "Ls CGNE 4D converged in " << diff.tv_sec << "." << diff.tv_usec << "s" << endl;
  
  return (T)1.0;
}

*/






#endif
