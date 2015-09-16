#ifndef _BFM_QDP_CAYLEY_PREC_H_
#define _BFM_QDP_CAYLEY_PREC_H_

#include <bfm_cg_unprec.h>
#include <bfm_cg_mixed_prec.h>
#include <omp.h>

//////////////////////////////////////////////////////
// Prototype Generalised 5d overlap solver
//////////////////////////////////////////////////////

template <class Float>
int bfm_Cayley_CG(LatticeFermion &sol, 
		  LatticeFermion &src,
		  multi1d<LatticeColorMatrix> &U,		     
		  LatticeComplex &PA,
		  LatticeComplex &PP,
		  LatticeComplex &PAconsv,
		  LatticeComplex &PJ5q,
		  bfmActionParams parms,
		  double residual,int max_iter);


//////////////////////////////////////////////////////
// Inverts the odd parity of MatPCdagMatPC.
// Maps to FMatEvlInv in CPS
//////////////////////////////////////////////////////
template <class Float>
int bfm_Cayley_CG_precMdagM_oo(multi1d<LatticeFermion> &sol, 
			       multi1d<LatticeFermion> &src,
			       multi1d<LatticeColorMatrix> &U,
			       std::string output_stem,
			       bfmActionParams parms, Real residual, int max_iter);


////////////////////////////////////////////////////////////////////////
// Implementation Conjugate gradient for Cayley form 5d preconditioned
////////////////////////////////////////////////////////////////////////
typedef multi1d<LatticeColorMatrix> U;

template <class Float>
int bfm_Cayley_CG(LatticeFermion &sol, 
		  LatticeFermion &src,
		  multi1d<LatticeColorMatrix> &U,
		  LatticeComplex &PA,
		  LatticeComplex &PP,
		  LatticeComplex &PAconsv,
		  LatticeComplex &PJ5q,
		  bfmActionParams parms,double residual,int max_iter)
{
  // Set up BAGEL object
  int lx = QDP::Layout::subgridLattSize()[0];
  int ly = QDP::Layout::subgridLattSize()[1];
  int lz = QDP::Layout::subgridLattSize()[2];
  int lt = QDP::Layout::subgridLattSize()[3];


  bfmarg bfma;
  //Physics parameters
  bfma.Ls           = parms.Ls;
  bfma.M5           = toDouble(parms.M5);
  bfma.mass         = toDouble(parms.mass);
  bfma.precon_5d    = 0;
  bfma.time_report_iter=50;

  // Following three need to be propagated to interface
  bfma.zolo_lo      = parms.zolo_lo;
  bfma.zolo_hi      = parms.zolo_hi;
  bfma.mobius_scale = parms.mobius_scale;

  bfma.solver       = parms.solver;
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
  bfma.solver= parms.solver;

  QDPIO::cout << "bfm_Cayley_CG : Initialising BAGEL-2 G5D solver "<<endl;

  bfm_qdp<Float> bfm; 
  bfm.init(bfma);

  Fermion_t sol_t[2];
  Fermion_t src_t[2];
  Fermion_t check_t[2];
  Fermion_t tmp;

  bfm.importGauge(U);
  tmp=bfm.allocFermion();

  // Import to Bagel
  for(int cb=0;cb<2;cb++){
    sol_t[cb] = bfm.allocFermion();
    src_t[cb] = bfm.allocFermion();
    check_t[cb] = bfm.allocFermion();
    bfm.importPhysicalFermion(src,src_t[cb],cb);
    bfm.importPhysicalFermion(sol,sol_t[cb],cb);
  }

  // Some source preparation for various cases
  QDPIO::cout << "setting thread count to "<<bfm.nthread<<endl;
  omp_set_num_threads(bfm.nthread);

  if ( ( parms.solver == HmCayleyTanh ) 
     ||( parms.solver == HwCayleyTanh ) 
     ||( parms.solver == HtCayleyTanh ) 
     ||( parms.solver == HtCayleyZolo ) 
     ||( parms.solver == HwCayleyZolo ) 
       ){ 

#pragma omp parallel 
  {

#pragma omp for 
    for(int i=0;i<bfm.nthread;i++) {
      bfm.G5D_Dminus(src_t,check_t,DaggerNo);

      for(int cb=0;cb<2;cb++){
	bfm.copy(src_t[cb],check_t[cb]);
      }
    }
  }
  } else {
    QDPIO::cout << "Unexpected solver type "<< endl;
    exit(-1);
  }
  // Run the inverter
#pragma omp parallel 
  {
#pragma omp for 
    for(int i=0;i<bfm.nthread;i++) {
      bfm.CGNE_M(sol_t, src_t);
    }
  }

  // snarf J5q
  PJ5q=zero;
  PAconsv=zero;
  PA=zero;
  PP=zero;

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
  if ( ( parms.solver == HwPartFracZolo)
     ||( parms.solver == HwPartFracTanh) ){
    sol = sol*2;
  }
  if ( ( parms.solver == HwPartFracZolo)
     ||( parms.solver == HwPartFracTanh) 
     ||( parms.solver == HwContFracZolo) 
     ||( parms.solver == HwContFracTanh) 
       ){ 

    sol = sol * 2.0/(1.0-parms.mass); // CF/PF Matrix norm factor
                                // Matrix is [ 1+m/1-m + g5 eps]
                                // Want [(1+m)/2+(1-m)/2 g5 eps] ^{-1}
                                // so 2/(1-m) factor on inverse.

    sol = sol -src;             // Contact term subtraction
    sol = sol / (1.0 - parms.mass);   
  }

  PA = localInnerProduct(sol,Gamma(1<<Nd-1)*sol);
  PP = localInnerProduct(sol,sol);

}  

int bfm_Cayley_mixed_precision_CG(LatticeFermion &sol, 
				  LatticeFermion &src,
				  multi1d<LatticeColorMatrix> &U,		     
				  LatticeComplex &PA,
				  LatticeComplex &PP,
				  LatticeComplex &PAconsv,
				  LatticeComplex &PJ5q,
				  bfmActionParams parms,
				  double residual,int max_iter);

#endif


