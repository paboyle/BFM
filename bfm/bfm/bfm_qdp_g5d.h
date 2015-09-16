#ifndef _BFM_QDP_G5D_H_
#define _BFM_QDP_G5D_H_

#include <bfm_cg_unprec.h>

#include <omp.h>

//////////////////////////////////////////////////////
// Prototype Generalised 5d overlap solver
//////////////////////////////////////////////////////

template <class Float>
int bfm_g5d_CG_unprec(LatticeFermion &sol, 
		      LatticeFermion &src,
		      multi1d<LatticeColorMatrix> &U,
		      bfmActionParams & solver,
		      double residual,int max_iter);


template <class Float>
int bfm_g5d_CG_unprec(LatticeFermion &sol, 
		      LatticeFermion &src,
		      multi1d<LatticeColorMatrix> &U,		     
		      LatticeComplex &PA,
		      LatticeComplex &PP,
		      LatticeComplex &PAconsv,
		      LatticeComplex &PJ5q,
		      bfmActionParams & solver,
		      double residual,int max_iter);

template <class Float>
int bfm_g5d_CG_unprec(LatticeFermion &sol, 
		      LatticeFermion &src,
		      multi1d<LatticeColorMatrix> &U,		     
		      bfmActionParams & solver,
		      double residual,int max_iter);

template <class Float>
int bfm_g5d_DeltaL(LatticeFermion &sol, 
		   LatticeFermion &src,
		   multi1d<LatticeColorMatrix> &U,
		   bfmActionParams & solver,
		   Real resid,int max_iter);


////////////////////////////////////////////////////////////////////////
// Implementation Conjugate gradient (unprec) for GENERAL 5d overlap
////////////////////////////////////////////////////////////////////////
typedef multi1d<LatticeColorMatrix> U;
typedef LatticeFermion T4;
typedef multi1d<LatticeFermion> T5;

template <class Float>
int bfm_g5d_CG_unprec(LatticeFermion &sol, 
		      LatticeFermion &src,
		      multi1d<LatticeColorMatrix> &U,
		      bfmActionParams & solver,
		      double residual,int max_iter)
{
  LatticeComplex PJ5q;
  LatticeComplex PA;
  LatticeComplex PP;
  LatticeComplex PAconsv;

  return bfm_g5d_CG_unprec<Float>(sol, src, U, PA,PP,PAconsv,PJ5q,solver,residual,max_iter);
}


template <class Float>
int bfm_g5d_CG_unprec(LatticeFermion &sol, 
		      LatticeFermion &src,
		      multi1d<LatticeColorMatrix> &U,
		      LatticeComplex &PA,
		      LatticeComplex &PP,
		      LatticeComplex &PAconsv,
		      LatticeComplex &PJ5q,
		      bfmActionParams & solver,
		      double residual,int max_iter)
{
  // Set up BAGEL object
  int lx = QDP::Layout::subgridLattSize()[0];
  int ly = QDP::Layout::subgridLattSize()[1];
  int lz = QDP::Layout::subgridLattSize()[2];
  int lt = QDP::Layout::subgridLattSize()[3];


  bfmarg bfma;

  //Physics parameters
  bfma.Ls           = solver.Ls;
  bfma.M5           = solver.M5;
  bfma.mass         = solver.mass;
  bfma.precon_5d    = 0;
  bfma.zolo_lo      = solver.zolo_lo;
  bfma.zolo_hi      = solver.zolo_hi;
  bfma.mobius_scale = solver.mobius_scale;

  bfma.solver       = solver.solver;
  bfma.max_iter     = max_iter;
  bfma.residual     = toDouble(residual);

  //Geometry
  bfma.node_latt[0] = lx;
  bfma.node_latt[1] = ly;
  bfma.node_latt[2] = lz;
  bfma.node_latt[3] = lt;

  multi1d<int> procs = QDP::Layout::logicalSize();
  for(int mu=0;mu<4;mu++){
    if (procs[mu]>1) bfma.local_comm[mu] = 0;
    else             bfma.local_comm[mu] = 1;
  }

  QDPIO::cout << "bfm_g5d_cg_unprec : Initialising BAGEL-2 G5D solver "<<endl;

  bfm_qdp<Float> bfm; 
  bfm.init(bfma);

  Fermion_t sol_t[2];
  Fermion_t src_t[2];
  Fermion_t check_t[2];
  Fermion_t tmp;

  bfm.importGauge(U);
  tmp=bfm.allocFermion();
  

  // Possible G5 in matrix convention
  LatticeFermion src_tmp = src;
  if ( ( solver.solver == HwPartFracZolo)
     ||( solver.solver == HwPartFracTanh) 
     ||( solver.solver == HwContFracZolo) 
     ||( solver.solver == HwContFracTanh) 
       ){ 
    Gamma G5(15);
    src_tmp = G5*src;
  }
  
  // Import to Bagel
  for(int cb=0;cb<2;cb++){
    sol_t[cb] = bfm.allocFermion();
    src_t[cb] = bfm.allocFermion();
    check_t[cb] = bfm.allocFermion();
    bfm.importPhysicalFermion(src_tmp,src_t[cb],cb);
    bfm.importPhysicalFermion(sol,sol_t[cb],cb);
  }


  // Some source preparation for various cases
  omp_set_num_threads(bfm.nthread);

  if ( ( solver.solver == HmCayleyTanh ) 
     ||( solver.solver == HwCayleyTanh ) 
     ||( solver.solver == HtCayleyTanh ) 
     ||( solver.solver == HtCayleyZolo ) 
     ||( solver.solver == HwCayleyZolo ) 
       ){ 


#pragma omp parallel 
  {
    omp_set_num_threads(bfm.nthread);

#pragma omp for 
    for(int i=0;i<bfm.nthread;i++) {
      bfm.G5D_Dminus(src_t,check_t,DaggerNo);

      for(int cb=0;cb<2;cb++){
	bfm.copy(src_t[cb],check_t[cb]);
      }
    }
  }}


  // Run the inverter
#pragma omp parallel 
  {
#pragma omp for 
    for(int i=0;i<bfm.nthread;i++) {
 
      double n = bfm.norm(src_t);

      bfm.CGNE_unprec(sol_t, src_t);

    }
  }

  // snarf J5q
  PJ5q=zero;
  PAconsv=zero;
  PA=zero;
  PP=zero;
  if ( ( solver.solver == HmCayleyTanh ) 
     ||( solver.solver == HwCayleyTanh ) 
     ||( solver.solver == HtCayleyTanh ) 
     ||( solver.solver == HtCayleyZolo ) 
     ||( solver.solver == HwCayleyZolo ) 
     ||( solver.solver == DWF)
       ){ 

    LatticeFermion  Midp;
    LatticeFermion  Midm;
    LatticeFermion  Mid;

    for(int cb=0;cb<2;cb++){
      int sm = solver.Ls/2-1;
      int sp = solver.Ls/2;
      bfm.exportFermion(Midm,sol_t[cb],cb,sm);
      bfm.exportFermion(Midp,sol_t[cb],cb,sp);
    }
    Mid  = chiralProjectPlus(Midm);
    Mid += chiralProjectMinus(Midp);
    PJ5q = localInnerProduct(Mid,Mid);
    QDPIO::cout << "DWF midpoint vol sum = " <<norm2(Mid) <<endl;
    QDPIO::cout << "DWF midpoint vol sum = " <<norm2(Mid) <<endl;

  }

  if(0) { // Conserved current not supported yet in general case
  if ( ( solver.solver == HmCayleyTanh ) 
     ||( solver.solver == HwCayleyTanh ) 
     ||( solver.solver == HtCayleyTanh ) 
     ||( solver.solver == HtCayleyZolo ) 
     ||( solver.solver == HwCayleyZolo ) 
       ){ 

#pragma omp parallel 
    {
#pragma omp for 
    for(int i=0;i<bfm.nthread;i++) {
 
      bfm.G5D_Dplus(sol_t,check_t,DaggerNo);

      bfm.solveMobiusDminus=1;

      if(i==0) QDPIO::cout<<"Calling Dminus inverse"<<endl;
      //Right hand term for general conserved current
      // D+ -> check_t
      // D-Inv D+->src_t
      // (1+b D-Inv D+) -> check_t
      bfm.CGNE_unprec(src_t,check_t);
      for(int cb=0;cb<2;cb++){
	for(int s=0;s<solver.Ls;s++){
	  bfm.axpby_ssp(check_t[cb],bfm.bs[s],sol_t[cb],bfm.cs[s],src_t[cb],s,s);
	}
      }

      //Left hand term for general conserved current
      if(i==0) QDPIO::cout<<"Calling Dminus inverse"<<endl;
      bfm.CGNE_unprec(src_t,sol_t);
      bfm.solveMobiusDminus=0;

    }}
  
    LatticeComplex C;
    LatticeFermion    LH;
    LatticeFermion    RH;
    LatticeFermion    p5d;
    LatticeFermion us_p5d;
    for(int s=0;s<solver.Ls;s++){

      
      int mu=Nd-1;
      int g5=Ns*Ns-1;
      int gt=8;

      // src_t == LH
      // check_t==RH
      for(int cb=0;cb<2;cb++){
	bfm.exportFermion(LH,src_t[cb],cb,solver.Ls-1-s);
	bfm.exportFermion(RH,check_t[cb],cb,s);
      }
	  p5d= LH;
      us_p5d = U[mu]*shift(RH,FORWARD,mu);

      C = 0.5*localInnerProduct(Gamma(g5)*p5d,Gamma(gt)*us_p5d);
      C-= 0.5*localInnerProduct(Gamma(g5)*p5d,us_p5d);

      p5d = RH;
      us_p5d = U[mu]*shift(LH,FORWARD,mu);

      C+= 0.5*localInnerProduct(Gamma(g5)*us_p5d,Gamma(gt)*p5d);
      C+= 0.5*localInnerProduct(Gamma(g5)*us_p5d,p5d);

      if (s < solver.Ls/2) PAconsv -= C;
      else          PAconsv += C;
      
    }
  }  
  }  

  // Export the solutoin
  sol=zero;
  for(int cb=0;cb<2;cb++){
    LatticeFermion  tmp=zero;
    bfm.exportPhysicalFermion(tmp,sol_t[cb],cb);
    sol+=tmp;
    bfm.freeFermion(sol_t[cb]);
    bfm.freeFermion(src_t[cb]);
    bfm.freeFermion(check_t[cb]);
  }
  bfm.freeFermion(tmp);
  bfm.end();

  
  // Continued and part frac subtract contact term
  if ( ( solver.solver == HwPartFracZolo)
     ||( solver.solver == HwPartFracTanh) ){
    sol = sol*2;
  }
  if ( ( solver.solver == HwPartFracZolo)
     ||( solver.solver == HwPartFracTanh) 
     ||( solver.solver == HwContFracZolo) 
     ||( solver.solver == HwContFracTanh) 
       ){ 

    sol = sol * 2.0/(1.0-solver.mass); // CF/PF Matrix norm factor
                                // Matrix is [ 1+m/1-m + g5 eps]
                                // Want [(1+m)/2+(1-m)/2 g5 eps] ^{-1}
                                // so 2/(1-m) factor on inverse.

    sol = sol -src;             // Contact term subtraction
    sol = sol / (1.0 - solver.mass);   
  }

  PA = localInnerProduct(sol,Gamma(1<<Nd-1)*sol);
  PP = localInnerProduct(sol,sol);
  QDPIO::cout << "DWF PP = "<<norm2(sol)<<endl;

}  

template <class Float>
int bfm_g5d_CG_unprec(LatticePropagator &sol, 
		      LatticePropagator &src,
		      multi1d<LatticeColorMatrix> &U,
		      int Ls,
		      int DodwfConservedVector,
		      Real mass,
		      Real M5,int Ncg, Real resid[],int max_iter[])
{
  return -1;
}





#endif
