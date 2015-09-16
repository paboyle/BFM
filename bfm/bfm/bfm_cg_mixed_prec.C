#include <bfm.h>
#include <bfm_qdp.h>
#include <bfm_cg_mixed_prec.h>

int bfm_mixed_precision_CG(  bfm_qdp<double> &bfm_double,
			     bfm_qdp<float>  &bfm_single,
			     double residual,int max_iter,
			     Fermion_t sol_d[2],
			     Fermion_t src_d[2]
			     )
{


  int max_outer = 10;
  int iters = 0;

  Fermion_t tmp_t[2];
  Fermion_t sol_s[2];
  Fermion_t defect_s[2];
  Fermion_t mtmp;
  // Import to Bagel
  for(int cb=0;cb<2;cb++){
    tmp_t[cb]    = bfm_double.allocFermion();
    defect_s[cb] = bfm_single.allocFermion();
    sol_s[cb]    = bfm_single.allocFermion();
  }
  mtmp = bfm_double.allocFermion();
  int converged = 0;
#pragma omp parallel 
  {
#pragma omp for 
    for(int i=0;i<bfm_double.nthread;i++) {

      int me = bfm_double.thread_barrier();

      double src_norm = bfm_double.norm(src_d);
      double defect_norm;
      double true_residual;
      bfm_double.ThreadBossDebug("BfmMixedPrecision: source norm is %le\n",src_norm);

      for(int outer = 0 ; outer < max_outer; outer ++) {
	
	int inner_iters;
	       
	// Iterate:
	// 
	//
	// defect_k    = eta - M sol_d_k 
	//
	// sol_d_{k+1}    = sol_d_k + M^{-1} defect_k
	//
	//
	// This implies
	//      => M sol_d_{k+1} - eta = M sol_d_k + defect_k  - eta
       	//                             = M sol_d_k - eta + eta - M sol_d_k
	//                             -> 0
	
	//////////////////////////////////////////////////////////
	// Double precision compute double precision residual vec
	//////////////////////////////////////////////////////////

	// M sol_d -> tmp_t

	if ( bfm_double.SPIcomms() && !me ) bfm_double.comm_init();  // Start double comms
	bfm_double.thread_barrier();

	bfm_double.Munprec(sol_d,tmp_t,mtmp,DaggerNo);

	defect_norm=0.0;
	// defect_k    = eta - M sol_d_k  = eta - tmp_t
	for(int cb=0;cb<2;cb++){
	  bfm_double.axpy(tmp_t[cb],tmp_t[cb],src_d[cb],-1.0);
	  defect_norm += bfm_double.norm(tmp_t[cb]);
	  bfm_double.precisionChange(tmp_t[cb],defect_s[cb],DoubleToSingle,cb);
	}
	true_residual = sqrt(defect_norm/src_norm);

	bfm_double.ThreadBossMessage("BfmMixedPrecision[%d]: defect norm is %le\n",outer,defect_norm);
	bfm_double.ThreadBossMessage("BfmMixedPrecision[%d]: defect norm is %le\n",outer,defect_norm);
	bfm_double.ThreadBossMessage("BfmMixedPrecision[%d]: true residual is %le\n",outer,true_residual);
	if ( true_residual < residual ) {
	  bfm_double.ThreadBossMessage("BfmMixedPrecision[%d]: converged after %d iterations and residual is %le\n",outer,iters,true_residual);
	  outer = max_outer;
	  converged=1;
	}

	if ( bfm_double.SPIcomms() && !me ) bfm_double.comm_end(); // End double comms
	bfm_double.thread_barrier();

	if ( !converged ) { 

	  if ( bfm_single.SPIcomms() && !me ) bfm_single.comm_init(); // Start single comms
	  bfm_single.thread_barrier();

	  ////////////////////////////////////////////////////
	  // Single precision inner CG
	  ////////////////////////////////////////////////////
	  double target_residual = residual/true_residual / 10;
	  if ( target_residual < 1.0e-6 ) target_residual = 1.0e-6;
	  bfm_single.ThreadBossMessage("BfmMixedPrecision[%d]: setting target residual to %le\n",outer,target_residual);
	  
	  for(int cb=0;cb<2;cb++){
	    bfm_single.copy(sol_s[cb]   ,defect_s[cb]); // Guess = source
	  }
	  
	  bfm_single.residual = target_residual;
	  inner_iters = bfm_single.CGNE_M(sol_s, defect_s); // Do the solve
	  if( me == 0 ) {
	    iters += inner_iters;
	  }
	  if ( bfm_single.SPIcomms() && !me ) bfm_single.comm_end(); // End single comms
	  bfm_single.thread_barrier();

	  ////////////////////////////////////////////////////
	  // Convert to double and subtract from sol_d
	  ////////////////////////////////////////////////////
	  // sol_d_{k+1}    = sol_d_k + M^{-1} defect_k
	  for(int cb=0;cb<2;cb++){
	    bfm_double.precisionChange(sol_s[cb],tmp_t[cb],SingleToDouble,cb);
	    bfm_double.axpy(sol_d[cb],sol_d[cb],tmp_t[cb],+1.0);
	  }
	}

      }
      
    }
  }
  if( converged ) return iters;
  bfm_double.Error("CG not converged\n");
  exit(-1);
}

int dwf_mixed_precision_CG(multi1d<LatticeFermion> &sol, 
			   multi1d<LatticeFermion> &src,
			   multi1d<LatticeColorMatrix> &U,		     
			   int Ls, Real mass, Real M5,double residual,int max_iter)
{
  // Set up BAGEL object
  int lx = QDP::Layout::subgridLattSize()[0];
  int ly = QDP::Layout::subgridLattSize()[1];
  int lz = QDP::Layout::subgridLattSize()[2];
  int lt = QDP::Layout::subgridLattSize()[3];

  bfmarg bfma;

  //Physics parameters
  bfma.Ls           = Ls;
  bfma.M5           = toDouble(M5);
  bfma.mass         = toDouble(mass);
  bfma.precon_5d    = 1;

  // Following three need to be propagated to interface
  bfma.solver       = DWF;
  bfma.max_iter     = max_iter;
  bfma.residual     = residual;

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

  QDPIO::cout << "dwf_mixedprecision_CG : Initialising BAGEL-2 G5D solver "<<endl;

  QDPIO::cout << "dwf_mixedprecision_CG : single precision matrix "<<endl;
  bfm_qdp<float> bfm_single; 
  bfm_single.init(bfma);
  bfm_single.importGauge(U);
  if (bfm_single.SPIcomms() )  bfm_single.comm_end();

  QDPIO::cout << "dwf_mixedprecision_CG : double precision matrix "<<endl;
  bfm_qdp<double> bfm_double; 
  bfm_double.init(bfma);
  bfm_double.importGauge(U);
  if (bfm_double.SPIcomms() )  bfm_double.comm_end(); // No overlapping

  Fermion_t sol_d[2];
  Fermion_t src_d[2];

  // Import to Bagel
  for(int cb=0;cb<2;cb++){
    sol_d[cb] = bfm_double.allocFermion();
    src_d[cb] = bfm_double.allocFermion();
  }

  // Some source preparation for various cases
  QDPIO::cout << "dwf_mixedprecision_CG: setting OpenMP thread count to "<<bfm_double.nthread<<endl;
  omp_set_num_threads(bfm_double.nthread);

  // Start with initial source, zero guess
  for(int cb=0;cb<2;cb++){
    bfm_double.importFermion(src,src_d[cb],cb);
  }


#pragma omp parallel 
  {
#pragma omp for 
    for(int i=0;i<bfm_double.nthread;i++) {
      for(int cb=0;cb<2;cb++){
	bfm_double.set_zero(sol_d[cb]);
      }
    }
  }

  int iter = bfm_mixed_precision_CG(bfm_double,bfm_single,residual, max_iter,sol_d,src_d);
  for(int cb=0;cb<2;cb++){
    bfm_double.freeFermion(sol_d[cb]);
    bfm_double.freeFermion(src_d[cb]);
  }
  return iter;
}

int dwf_mixed_precision_CG(LatticeFermion &sol, 
			   LatticeFermion &src,
			   multi1d<LatticeColorMatrix> &U,		     
			   int Ls, Real mass, Real M5,double residual,int max_iter)
{
  // Set up BAGEL object
  int lx = QDP::Layout::subgridLattSize()[0];
  int ly = QDP::Layout::subgridLattSize()[1];
  int lz = QDP::Layout::subgridLattSize()[2];
  int lt = QDP::Layout::subgridLattSize()[3];

  bfmarg bfma;

  //Physics parameters
  bfma.Ls           = Ls;
  bfma.M5           = toDouble(M5);
  bfma.mass         = toDouble(mass);
  bfma.precon_5d    = 1;

  // Following three need to be propagated to interface
  bfma.solver       = DWF;
  bfma.max_iter     = max_iter;
  bfma.residual     = residual;

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

  QDPIO::cout << "dwf_mixedprecision_CG : Initialising BAGEL-2 G5D solver "<<endl;

  QDPIO::cout << "dwf_mixedprecision_CG : single precision matrix "<<endl;
  bfm_qdp<float> bfm_single; 
  bfm_single.init(bfma);
  bfm_single.importGauge(U);
  if (bfm_single.SPIcomms() )  bfm_single.comm_end(); // No overlapping

  QDPIO::cout << "dwf_mixedprecision_CG : double precision matrix "<<endl;
  bfm_qdp<double> bfm_double; 
  bfm_double.init(bfma);
  bfm_double.importGauge(U);
  if (bfm_double.SPIcomms() )  bfm_double.comm_end(); // No overlapping

  Fermion_t sol_d[2];
  Fermion_t src_d[2];

  // Import to Bagel
  for(int cb=0;cb<2;cb++){
    sol_d[cb] = bfm_double.allocFermion();
    src_d[cb] = bfm_double.allocFermion();
  }

  // Some source preparation for various cases
  QDPIO::cout << "dwf_mixedprecision_CG: setting OpenMP thread count to "<<bfm_double.nthread<<endl;
  omp_set_num_threads(bfm_double.nthread);

  // Start with initial source, zero guess
  for(int cb=0;cb<2;cb++){
    bfm_double.importPhysicalFermion(src,src_d[cb],cb);
  }


#pragma omp parallel 
  {
#pragma omp for 
    for(int i=0;i<bfm_double.nthread;i++) {
      for(int cb=0;cb<2;cb++){
	bfm_double.set_zero(sol_d[cb]);
      }
    }
  }

  int iter = bfm_mixed_precision_CG(bfm_double,bfm_single,residual, max_iter,sol_d,src_d);
  for(int cb=0;cb<2;cb++){
    bfm_double.freeFermion(sol_d[cb]);
    bfm_double.freeFermion(src_d[cb]);
  }
  return iter;
}

int bfm_Cayley_mixed_precision_CG(LatticeFermion &sol, 
				  LatticeFermion &src,
				  multi1d<LatticeColorMatrix> &U,		     
				  LatticeComplex &PA,
				  LatticeComplex &PP,
				  LatticeComplex &PAconsv,
				  LatticeComplex &PJ5q,
				  bfmActionParams parms,
				  double residual,int max_iter)
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

  // Following three need to be propagated to interface
  bfma.zolo_lo      = parms.zolo_lo;
  bfma.zolo_hi      = parms.zolo_hi;
  bfma.mobius_scale = parms.mobius_scale;

  bfma.solver       = parms.solver;
  bfma.max_iter     = max_iter;
  bfma.residual     = residual;

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

  QDPIO::cout << "bfm_Cayley_mixedprecision_CG : Initialising BAGEL-2 G5D solver "<<endl;

  QDPIO::cout << "bfm_Cayley_mixedprecision_CG : single precision matrix "<<endl;
  bfm_qdp<float> bfm_single; 
  bfm_single.init(bfma);
  bfm_single.importGauge(U);

  QDPIO::cout << "bfm_Cayley_mixedprecision_CG : double precision matrix "<<endl;
  bfm_qdp<double> bfm_double; 
  bfm_double.init(bfma);
  bfm_double.importGauge(U);

  Fermion_t sol_d[2];
  Fermion_t src_d[2];
  Fermion_t tmp_t[2];

  // Import to Bagel
  for(int cb=0;cb<2;cb++){
    sol_d[cb] = bfm_double.allocFermion();
    src_d[cb] = bfm_double.allocFermion();
    tmp_t[cb] = bfm_double.allocFermion();
  }

  // Some source preparation for various cases
  QDPIO::cout << "bfm_cayley_mixedprecision_CG: setting OpenMP thread count to "<<bfm_double.nthread<<endl;
  omp_set_num_threads(bfm_double.nthread);

  if ( ( parms.solver != HmCayleyTanh ) 
     &&( parms.solver != HwCayleyTanh ) 
     &&( parms.solver != HtCayleyTanh ) 
     &&( parms.solver != HtCayleyZolo ) 
     &&( parms.solver != HwCayleyZolo ) 
       ){ 
    QDPIO::cout << "Unexpected solver type "<< endl;
    exit(-1);
  }
  
  // Start with initial source, zero guess
  for(int cb=0;cb<2;cb++){
    bfm_double.importPhysicalFermion(src,src_d[cb],cb);
  }


#pragma omp parallel 
  {
#pragma omp for 
    for(int i=0;i<bfm_double.nthread;i++) {

      // Basis change for Cayley source
      // Zero initial guess
      bfm_double.G5D_Dminus(src_d,tmp_t,DaggerNo); 

      for(int cb=0;cb<2;cb++){
	bfm_double.copy(src_d[cb],tmp_t[cb]);  
	bfm_double.set_zero(sol_d[cb]);
      }

    }
  }

  int iter =  bfm_mixed_precision_CG(bfm_double,bfm_single,residual, max_iter,sol_d,src_d);
  for(int cb=0;cb<2;cb++){
    bfm_double.freeFermion(sol_d[cb]);
    bfm_double.freeFermion(src_d[cb]);
    bfm_double.freeFermion(tmp_t[cb]);
  }
  return iter;
}
