#include "bfm.h"
#include "bfm_qdp.h"
#include "bfm_wrapper.h"
#include <math.h>
#if defined(QDP_USE_OMP_THREADS)
#warning QDP using OpenMP threading
#include <omp.h>
#endif

template<class Float> 
void BfmWrapper::SpectralRange(bfm_qdp<Float>  &bfm)
{
  int Ls = bfm.Ls;
  multi1d<T4> gaus(Ls);
  for(int s=0;s<Ls;s++)
    gaussian(gaus[s]);


  Fermion_t tmp1,tmp2;
  Fermion_t b;
  Fermion_t Ab;
  double mu;
  tmp1   = bfm.allocFermion();
  tmp2   = bfm.allocFermion();
  b      = bfm.allocFermion();
  Ab     = bfm.allocFermion();

  bfm.importFermion(gaus,b,0);
  
  int dagyes= 1;
  int dagno = 0;
  int donrm = 1;
  double abnorm;
  double absq;
  double n;
  ///////////////////////////////////////
  // Get the highest evalue of MdagM
  ///////////////////////////////////////
#pragma omp parallel for  
for(int i=0;i<bfm.nthread;i++) {
  int me = bfm.thread_barrier();

  for(int k=0;k<100;k++){

    n = bfm.norm(b);

    absq = bfm.Mprec(b,tmp2,tmp1,dagno,donrm);
           bfm.Mprec(tmp2,Ab,tmp1,dagyes); //MdagM
    
    mu = absq/n;
    if(bfm.isBoss() && (!me) ) { 
      printf("PowerMethod %d mu=%le absq=%le n=%le \n",k,mu,absq,n);
    }

    abnorm = sqrt(absq);

    // b = Ab / |Ab|
    bfm.axpby(b,Ab,Ab,0.0,1.0/abnorm);
  }
 }

 /////////////////////////////////////////
 // Get the lowest eigenvalue of mdagm
 // Apply power method to [lambda_max - MdagM]
 /////////////////////////////////////////
 double lambda_max = mu;

 bfm.importFermion(gaus,b,0);

#pragma omp parallel for  
for(int i=0;i<bfm.nthread;i++) {
  int me = bfm.thread_barrier();

  for(int k=0;k<100;k++){

    n = bfm.norm(b);

    absq = bfm.Mprec(b,tmp2,tmp1,dagno,donrm);
    bfm.Mprec(tmp2,Ab,tmp1,dagyes); 
    bfm.axpby(Ab,Ab,b,-1,lambda_max);

    absq = lambda_max * n - absq;

    mu = absq/n;

    if(bfm.isBoss() && (!me) ) { 
      printf("PowerMethod for low EV %d mu=%le absq=%le n=%le \n",k,mu,absq,n);
      printf("PowerMethod  %d lambda_min  =%le\n",k,lambda_max - mu);
    }

    abnorm = sqrt(absq);

    // b = Ab / |Ab|
    bfm.axpby(b,Ab,Ab,0.0,1.0/abnorm);

  }
}

  bfm.freeFermion(tmp1);
  bfm.freeFermion(tmp2);
  bfm.freeFermion(b);
  bfm.freeFermion(Ab);


}

  int BfmWrapper::bfmInvert4d_mixed(T4 &sol,const T4 &chi)
  {
    QDP_error_exit("bfmInvert4d_mixed not implemented");
  }

template<class Float>
int BfmWrapper::bfmInvert_mixed(bfm_qdp<Float> &bfm_dp,
				Fermion_t sol_d[2],
				Fermion_t src_d[2])
{
  bfmarg bfma = *( (bfmarg *) & bfm_dp);
  bfm_qdp<float> bfm_sp;
  bfm_sp.init(bfma);
  int iter = bfmInvert_mixed(bfm_dp,bfm_sp,sol_d,src_d);
  bfm_sp.end();
  return iter;
}


template<class Float>
int BfmWrapper::bfmInvert_mixed(bfm_qdp<Float>  &bfm_dp,
				bfm_qdp<float> &bfm_sp,
				Fermion_t sol_d[2],
				Fermion_t src_d[2])
{
    Fermion_t tmp_t[2];
    Fermion_t sol_s[2];
    Fermion_t mtmp;
    Fermion_t defect_s[2];

    double residual   = toDouble(invParam.RsdTarget[0]);
    bfm_dp.residual = bfm_sp.residual = residual;

    for(int cb=0;cb<2; cb++){
      tmp_t[cb] = bfm_dp.allocFermion();
      sol_s[cb] = bfm_sp.allocFermion();
      defect_s[cb] = bfm_sp.allocFermion();
    }
    mtmp =  bfm_dp.allocFermion();

    // If preconditioned MdagM inversion set the other parity to zero
    if ( invParam.BfmMatrix == BfmMat_MdagM_oo ) { 
      bfm_dp.master_fill(src_d[0],0.0);
      bfm_dp.master_fill(sol_d[0],0.0);
      bfm_sp.master_fill(defect_s[0],0.0);
      bfm_dp.master_fill(tmp_t[0],0.0);
    }
    int converged = 0;
    int iter = 0;
    int max_outer = 10;

    QDPIO::cerr << "bfmInvert_mixed: starting thread loop for "<<bfm_dp.nthread << " threads and target residual "<< residual<<endl;

    // Run the inverter
#pragma omp parallel for
    for(int i=0;i<bfm_dp.nthread;i++) {

      int me = bfm_dp.thread_barrier();

      double src_norm = bfm_dp.norm(src_d);
      double defect_norm;
      double true_residual;

      if ( bfm_dp.isBoss() && (!me) ) {
	printf("bfmInvert_mixed: source norm is %le\n",src_norm);
      }


      for(int outer = 0 ; outer < max_outer; outer ++) {

	int inner_iters;

	if ( bfm_dp.SPIcomms() && (!me) ) {
	  bfm_dp.comm_init();  // Start double comms
	}
	bfm_dp.thread_barrier();


	defect_norm=0.0;
	if ( invParam.BfmMatrix == BfmMat_MdagM_oo ) { 

	  bfm_dp.Mprec(sol_d[1],tmp_t[0],mtmp,DaggerNo);
	  bfm_dp.Mprec(tmp_t[0],tmp_t[1],mtmp,DaggerYes);
	  bfm_dp.axpy(tmp_t[1],tmp_t[1],src_d[1],-1.0);
	  defect_norm = bfm_dp.norm(tmp_t[1]);
	  bfm_dp.precisionChange(tmp_t[1],defect_s[1],DoubleToSingle,1);

	} else if ( invParam.BfmMatrix == BfmMat_M ) { 

	  bfm_dp.Munprec(sol_d,tmp_t,mtmp,DaggerNo);

	  for(int cb=0;cb<2;cb++){
	    bfm_dp.axpy(tmp_t[cb],tmp_t[cb],src_d[cb],-1.0);
	    defect_norm += bfm_dp.norm(tmp_t[cb]);
	    bfm_dp.precisionChange(tmp_t[cb],defect_s[cb],DoubleToSingle,cb);
	  }
	} 
	true_residual = sqrt(defect_norm/src_norm);

	if ( bfm_dp.isBoss() && (!me) ) printf("bfmInvert_mixed[%d]: true residual is %le\n",outer,true_residual); fflush(stdout);
	if ( true_residual < residual ) {
	  if ( bfm_dp.isBoss() && !me ) 
	    printf("bfmInvert_mixed[%d]: converged after %d iterations and residual is %le\n",outer,iter,true_residual);
	  outer = max_outer;
	  converged=1;
	}

	// Hack to reinit comms buffers that are shared between the two
	if ( bfm_dp.SPIcomms() && !me ) bfm_dp.comm_end(); // End double comms
	bfm_dp.thread_barrier();

	if ( !converged ) { 

	  if ( bfm_sp.SPIcomms() && !me ) bfm_sp.comm_init(); // Start single comms
	  bfm_sp.thread_barrier();

	  ////////////////////////////////////////////////////
	  // Single precision inner CG
	  ////////////////////////////////////////////////////
	  double target_residual = residual/true_residual * 0.75;
	  double delta = toDouble(invParam.Delta);
	  if ( target_residual < delta ) target_residual = delta;
	  if ( bfm_sp.isBoss() && (!me) )
	    printf("bfmInvert_mixed[%d]: setting target residual to %le\n",outer,target_residual);
	  
	  for(int cb=0;cb<2;cb++){
	    bfm_sp.copy(sol_s[cb]   ,defect_s[cb]); // Guess = source
	  }
	  bfm_sp.residual = target_residual;

	  if ( invParam.BfmMatrix == BfmMat_MdagM_oo ) { 
	    if ( bfm_sp.isBoss() && (!me) )
	      printf("bfmInvert_mixed: CGNE_MdagM prec solver \n");
	    inner_iters = bfm_sp.CGNE_prec(sol_s[1], defect_s[1]); // Do the solve
	  } else if ( invParam.BfmMatrix == BfmMat_M ) { 
	    if ( bfm_sp.isBoss() && (!me) )
	      printf("bfmInvert_mixed: CGNE_M solver \n");
	  
	    inner_iters = bfm_sp.CGNE_M(sol_s, defect_s); // Do the solve
	    double n_e = bfm_sp.norm(sol_s[0]);
	    double n_o = bfm_sp.norm(sol_s[1]);
	    if ( bfm_sp.isBoss() && (!me) )
	      printf("bfmInvert_mixed: CGNE_M solver solution  %le %le \n",n_e,n_o);
	    
	  } else if ( invParam.BfmMatrix == BfmMat_Mdag ) { 
	    if ( bfm_sp.isBoss() && (!me) )
	      printf("bfmInvert_mixed: CGNE_Mdag solver \n");
	  
	    inner_iters = bfm_sp.CGNE_Mdag(sol_s, defect_s); // Do the solve
	    double n_e = bfm_sp.norm(sol_s[0]);
	    double n_o = bfm_sp.norm(sol_s[1]);
	    if ( bfm_sp.isBoss() && (!me) )
	      printf("bfmInvert_mixed: CGNE_Mdag solver solution  %le %le \n",n_e,n_o);
	    
	  } else { 
	    exit(-1);
	  }

	  if( me == 0 ) {
	    iter += inner_iters;
	  }
	  if ( bfm_sp.SPIcomms() && !me ) bfm_sp.comm_end(); // End single comms
	  bfm_sp.thread_barrier();

	  ////////////////////////////////////////////////////
	  // Convert to double and subtract from sol_d
	  ////////////////////////////////////////////////////
	  // sol_d_{k+1}    = sol_d_k + M^{-1} defect_k
	  for(int cb=0;cb<2;cb++){
	    bfm_dp.precisionChange(sol_s[cb],tmp_t[cb],SingleToDouble,cb);
	    bfm_dp.axpy(sol_d[cb],sol_d[cb],tmp_t[cb],+1.0);
	  }
	}

      }

    }
    for(int cb=0;cb<2;cb++) {
      bfm_dp.freeFermion(tmp_t[cb]);
      bfm_sp.freeFermion(defect_s[cb]);
      bfm_sp.freeFermion(sol_s[cb]);
    }
    bfm_dp.freeFermion(mtmp);
    return iter;
}


  int BfmWrapper::bfmInvert5d_mixed(multi1d<T4> &sol,const multi1d<T4> &chi)
  {
    int res;
    multi1d<T4> src = chi;

    // Set up BAGEL object
    int lx = QDP::Layout::subgridLattSize()[0];
    int ly = QDP::Layout::subgridLattSize()[1];
    int lz = QDP::Layout::subgridLattSize()[2];
    int lt = QDP::Layout::subgridLattSize()[3];
    
    bfmarg bfma;
#if defined(QDP_USE_OMP_THREADS)
    bfma.Threads(omp_get_max_threads());
#else
    bfma.Threads(1);
#endif
    //    bfma.Verbose(0);

    //Physics parameters
    bfmActionParams *bfmap = (bfmActionParams *) &bfma;
    *bfmap = invParam.BAP;
    
    // Algorithm & code control
    bfma.time_report_iter=-100;
    bfma.max_iter     = invParam.MaxIter;
    int max_iter      = bfma.max_iter;

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
    
    // Bfm object
    bfm_qdp<float> bfm_sp; 
    bfm_qdp<double> bfm_dp; 
    bfm_sp.init(bfma);
    bfm_dp.init(bfma);
    bfm_sp.importGauge(links);
    bfm_dp.importGauge(links);
    
    //Fermion import
    Fermion_t sol_d[2];
    Fermion_t src_d[2];
    
    for(int cb=0;cb<2; cb++){
      sol_d[cb] = bfm_dp.allocFermion();
      src_d[cb] = bfm_dp.allocFermion();
      bfm_dp.importFermion(src,src_d[cb],cb);
      bfm_dp.importFermion(sol,sol_d[cb],cb);
    }

    res=bfmInvert_mixed<double>(bfm_dp,bfm_sp,sol_d,src_d);

    // MdagM solution
    for(int cb=0;cb<2;cb++) {
      bfm_dp.exportFermion(sol,sol_d[cb],cb);
    }
    for(int cb=0;cb<2;cb++) {
      bfm_dp.freeFermion(sol_d[cb]);
      bfm_dp.freeFermion(src_d[cb]);
    }
    bfm_dp.end();
    bfm_sp.end();
    return res;
  }

template<class Float> 
int BfmWrapper::bfmCayleyDminus(multi1d<T4> &inout)
{
    // Set up BAGEL object
    int lx = QDP::Layout::subgridLattSize()[0];
    int ly = QDP::Layout::subgridLattSize()[1];
    int lz = QDP::Layout::subgridLattSize()[2];
    int lt = QDP::Layout::subgridLattSize()[3];
    
    bfmarg bfma;

#if defined(QDP_USE_OMP_THREADS)
    bfma.Threads(omp_get_max_threads());
#else
    bfma.Threads(1);
#endif
    //    bfma.Verbose(0);

    //Physics parameters
    bfmActionParams *bfmap = (bfmActionParams *) &bfma;
    *bfmap = invParam.BAP;
    
    // Algorithm & code control
    bfma.time_report_iter=-100;
    bfma.max_iter     = invParam.MaxIter;
    bfma.residual     = toDouble(invParam.RsdTarget[0]);

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
    
    // Bfm object
    bfm_qdp<Float> bfm; 
    bfm.init(bfma);

    //Gauge field import
    bfm.importGauge(links);
    
    //Fermion import
    Fermion_t sol_t[2];
    Fermion_t src_t[2];
    
    for(int cb=0;cb<2; cb++){
      sol_t[cb] = bfm.allocFermion();
      src_t[cb] = bfm.allocFermion();
      bfm.importFermion(inout,src_t[cb],cb);
    }

#pragma omp parallel for
    for(int i=0;i<bfm.nthread;i++) {
      bfm.G5D_Dminus(src_t,sol_t,0);
    }

    for(int cb=0;cb<2;cb++) {
      bfm.exportFermion(inout,sol_t[cb],cb);
    }

    for(int cb=0;cb<2;cb++) {
      bfm.freeFermion(sol_t[cb]);
      bfm.freeFermion(src_t[cb]);
    }
    bfm.end();
    return 0;
}


  template<class Float>
  int BfmWrapper::bfmInvert5d(multi1d<T4> &sol,const multi1d<T4> &chi)
  {
    int res;
    multi1d<T4> src = chi;

    // Set up BAGEL object
    int lx = QDP::Layout::subgridLattSize()[0];
    int ly = QDP::Layout::subgridLattSize()[1];
    int lz = QDP::Layout::subgridLattSize()[2];
    int lt = QDP::Layout::subgridLattSize()[3];
    
    bfmarg bfma;

#if defined(QDP_USE_OMP_THREADS)
    bfma.Threads(omp_get_max_threads());
#else
    bfma.Threads(1);
#endif
    //    bfma.Verbose(0);

    //Physics parameters
    bfmActionParams *bfmap = (bfmActionParams *) &bfma;
    *bfmap = invParam.BAP;
    
    // Algorithm & code control
    bfma.time_report_iter=-100;
    bfma.max_iter     = invParam.MaxIter;
    bfma.residual     = toDouble(invParam.RsdTarget[0]);

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
    
    // Bfm object
    bfm_qdp<Float> bfm; 
    bfm.init(bfma);

    //Gauge field import
    bfm.importGauge(links);
    
    //Fermion import
    Fermion_t sol_t[2];
    Fermion_t src_t[2];
    
    for(int cb=0;cb<2; cb++){
      sol_t[cb] = bfm.allocFermion();
      src_t[cb] = bfm.allocFermion();
      bfm.importFermion(src,src_t[cb],cb);
      bfm.importFermion(sol,sol_t[cb],cb);
    }
    
    QDPIO::cout << "BfmWrapper guess norm is "<< norm2(sol)<<endl;
    QDPIO::cerr << "bfmInvert5d: solver beginning"<<endl;
    // Run the inverter
    sol=zero;
    int iter;

    // MdagM solution
    if ( invParam.BfmMatrix == BfmMat_MdagM_oo ) { 
#pragma omp parallel for 
      for(int i=0;i<bfm.nthread;i++) {
	int iter_thr = bfm.CGNE_prec_MdagM(sol_t[1], src_t[1]);
	if ( i==0 ) iter = iter_thr;
      }
      res = iter;
      bfm.exportFermion(sol,sol_t[1],1);

    // M solution 
    // MdagM solution

    } else if ( invParam.BfmMatrix == BfmMat_CayleyDminusUnprec ) { 

      bfm.solveMobiusDminus=1;
#pragma omp parallel for 
      for(int i=0;i<bfm.nthread;i++) {
	int iter_thr = bfm.CGNE_unprec(sol_t,src_t);
	if ( i==0 ) iter = iter_thr;
      }
      bfm.solveMobiusDminus=0;
      res = iter;
      bfm.exportFermion(sol,sol_t[0],0);
      bfm.exportFermion(sol,sol_t[1],1);

    // M solution 
    } else if ( invParam.BfmMatrix == BfmMat_M ) { 
#pragma omp parallel for 
      for(int i=0;i<bfm.nthread;i++) {
	int iter_thr = bfm.CGNE_M(sol_t, src_t);
	if ( i==0 ) iter = iter_thr;
      }
      for(int cb=0;cb<2;cb++) {
	bfm.exportFermion(sol,sol_t[cb],cb);
      }
    } else if ( invParam.BfmMatrix == BfmMat_Mdag ) { 

#pragma omp parallel for 
      for(int i=0;i<bfm.nthread;i++) {
	int iter_thr = bfm.CGNE_Mdag(sol_t, src_t);
	if ( i==0 ) iter = iter_thr;
      }
      for(int cb=0;cb<2;cb++) {
	bfm.exportFermion(sol,sol_t[cb],cb);
      }
    } else { 
      printf("BfmMatrix %d \n",invParam.BfmMatrix); fflush(stdout);
      QDP_error_exit("BfmWrapper :: Invert5d don't know the matrix type given");
    }
    for(int cb=0;cb<2;cb++) {
      bfm.freeFermion(sol_t[cb]);
      bfm.freeFermion(src_t[cb]);
    }
    res=iter;
    bfm.end();
    return res;
  }
  
  template<class Float> int BfmWrapper::bfmInvert4d(T4 &sol,const T4 &chi)
  {
    int res;
    T4 src = chi;
    
    // Set up BAGEL object
    int lx = QDP::Layout::subgridLattSize()[0];
    int ly = QDP::Layout::subgridLattSize()[1];
    int lz = QDP::Layout::subgridLattSize()[2];
    int lt = QDP::Layout::subgridLattSize()[3];
    
    bfmarg bfma;
#if defined(QDP_USE_OMP_THREADS)
    bfma.Threads( omp_get_max_threads());
#else
    bfma.Threads(1);
#endif
    //    bfma.Verbose(0);
    
    //Physics parameters
    bfmActionParams *bfmap = (bfmActionParams *) &bfma;
    *bfmap = invParam.BAP;
    
    // Algorithm & code control
    bfma.time_report_iter=-100;
    bfma.max_iter     = invParam.MaxIter;
    bfma.residual     = toDouble(invParam.RsdTarget[0]);
    
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

    QDPIO::cout << "Initialising BAGEL-2 solver "<<endl;

    // Bfm object
    bfm_qdp<Float> bfm; 
    bfm.init(bfma);
    
    //Gauge field import
    bfm.importGauge(links);
    
    //begin karthee clover
    //Missing clover term
    if ( invParam.BAP.solver  == CloverFermion) {
      // Import the ClovDiag, ClovOffDiag
      // Import the ClovInvDiag, ClovInvOffDiag
      bfm.importClover(CloverDiag,CloverOffDiag, bfm.A);
      bfm.importClover(CloverInvDiag,CloverInvOffDiag, bfm.Ainv);
    }
    //end karthee clover

    //Fermion import
    Fermion_t sol_t[2];
    Fermion_t src_t[2];

    for(int cb=0;cb<2; cb++){
      sol_t[cb] = bfm.allocFermion();
      src_t[cb] = bfm.allocFermion();
      bfm.importFermion(src,src_t[cb],cb);
      bfm.importFermion(sol,sol_t[cb],cb);
    }

    // Run the inverter
    sol=zero;
    int iter;

    // MdagM solution
    if ( invParam.BfmMatrix == BfmMat_MdagM_oo ) { 
      QDPIO::cout << "Calling BAGEL BfmMat_MdagM_oo eo CG inverter"<<endl;
#pragma omp parallel for 
      for(int i=0;i<bfm.nthread;i++) {
	int iter_thr = bfm.CGNE_prec_MdagM(sol_t[1], src_t[1]);
	if ( i==0 ) iter = iter_thr;
      }
      res = iter;
      bfm.exportFermion(sol,sol_t[1],1);

    // M solution 
    } else if ( invParam.BfmMatrix == BfmMat_M ) { 

      QDPIO::cout << "Calling BAGEL BfmMat_M eo CG inverter"<<endl;
#pragma omp parallel for 
      for(int i=0;i<bfm.nthread;i++) {
	int iter_thr = bfm.CGNE_M(sol_t, src_t);
	if ( i==0 ) iter = iter_thr;
      }
      for(int cb=0;cb<2;cb++) {
	bfm.exportFermion(sol,sol_t[cb],cb);
      }
    } else if ( invParam.BfmMatrix == BfmMat_Mdag ) { 


      QDPIO::cout << "Calling BAGEL BfmMat_Mdag eo CG inverter"<<endl;
#pragma omp parallel for 
      for(int i=0;i<bfm.nthread;i++) {
	int iter_thr = bfm.CGNE_Mdag(sol_t, src_t);
	if ( i==0 ) iter = iter_thr;
      }
      for(int cb=0;cb<2;cb++) {
	bfm.exportFermion(sol,sol_t[cb],cb);
      }
    } else { 
      printf("BfmMatrix %d \n",invParam.BfmMatrix); fflush(stdout);
      QDP_error_exit("BfmWrapper :: Invert4d don't know the matrix type given");
    }
    
    for(int cb=0;cb<2;cb++) {
      bfm.freeFermion(sol_t[cb]);
      bfm.freeFermion(src_t[cb]);
    }

    res = iter;
   
    bfm.end();
    return res;

  }


  template<class Float>
  int BfmWrapper::bfmMultiInvert5d(multi1d< multi1d<T4> >& psi, 
						     const multi1d<Real>& shifts, 
						     const multi1d<Real>& Residuals, 
						     const multi1d<T4>& chi)
  {
    int nshift = shifts.size();
    Fermion_t sol_t[nshift];
    Fermion_t src_t;
    double    dshifts[nshift];
    double    alpha[nshift];
    double    mresidual[nshift];
    int       dontsum=0;

    int Ls = invParam.BAP.Ls;
    
    if ( chi.size() != Ls ) { 
      QDP_error_exit("Ls mismatch in bfmMultiInvert5d");
    }

    psi.resize(nshift);
    for(int shift=0;shift<nshift;shift++){
      dshifts[shift] =toDouble(shifts[shift]);
      alpha[shift]   = 1.0;
      mresidual[shift]= toDouble(Residuals[shift]);

      psi[shift].resize(Ls);
      for(int s=0;s<Ls;s++){
	psi[shift][s] = zero;
      }
    }

    int res;
    multi1d<T4> src = chi;
    
    // Set up BAGEL object
    int lx = QDP::Layout::subgridLattSize()[0];
    int ly = QDP::Layout::subgridLattSize()[1];
    int lz = QDP::Layout::subgridLattSize()[2];
    int lt = QDP::Layout::subgridLattSize()[3];
    
    bfmarg bfma;

#if defined(QDP_USE_OMP_THREADS)
    bfma.Threads(omp_get_max_threads());
#else
    bfma.Threads(1);
#endif
    //    bfma.Verbose(0);

    //Physics parameters
    bfmActionParams *bfmap = (bfmActionParams *) &bfma;
    *bfmap = invParam.BAP;
    
    // Algorithm & code control
    bfma.time_report_iter=-100;
    bfma.max_iter     = invParam.MaxIter;
  
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
    
    QDPIO::cout << "Initialising BAGEL-2 solver "<<endl;
    // Bfm object
    bfm_qdp<Float> bfm; 
    bfm.init(bfma);

    //Gauge field import
    bfm.importGauge(links);
    
    //Fermion import
    int cb = 1;
    for(int shift=0;shift<nshift;shift++){
      sol_t[shift] = bfm.allocFermion();
      if (sol_t[shift] == NULL ) { 
	QDP_error_exit("Allocate failed\n");
      }
      bfm.importFermion(psi[shift],sol_t[shift],cb);
    }
    src_t = bfm.allocFermion();
    if (src_t == NULL ) { 
      QDP_error_exit("Allocate failed\n");
    }
    bfm.importFermion(src,src_t,cb);
    
    // Run the inverter
    int iter;

#pragma omp parallel for 
    for(int i=0;i<bfm.nthread;i++) {
      int iter_thr = bfm.CGNE_prec_MdagM_multi_shift(sol_t,src_t,
						     dshifts,
						     alpha,
						     nshift,
						     mresidual,
						     dontsum);
      if ( i==0 ) iter = iter_thr;
    }
    res = iter;
    
    for(int shift=0;shift<nshift;shift++) {
      bfm.exportFermion(psi[shift],sol_t[shift],cb);
      bfm.freeFermion(sol_t[shift]);
    }
    bfm.freeFermion(src_t);
    
    res=iter;

    bfm.end();
    return res;
  }

  template<class Float>  int BfmWrapper::bfmMultiInvert4d(multi1d<T4> & psi, 
						     const multi1d<Real>& shifts, 
						     const multi1d<Real>& Residuals, 
						     const T4& chi)
  {
    int nshift = shifts.size();
    Fermion_t sol_t[nshift];
    Fermion_t src_t;
    double    dshifts[nshift];
    double    alpha[nshift];
    double    mresidual[nshift];
    int       dontsum=0;

    for(int shift=0;shift<nshift;shift++){
      dshifts[shift] =toDouble(shifts[shift]);
      alpha[shift]   = 1.0;
      mresidual[shift]= toDouble(Residuals[shift]);
      psi[shift] = zero;
    }

    int res;
    T4 src = chi;
    
    // Set up BAGEL object
    int lx = QDP::Layout::subgridLattSize()[0];
    int ly = QDP::Layout::subgridLattSize()[1];
    int lz = QDP::Layout::subgridLattSize()[2];
    int lt = QDP::Layout::subgridLattSize()[3];
    
    bfmarg bfma;

#if defined(QDP_USE_OMP_THREADS)
    bfma.Threads(omp_get_max_threads());
#else
    bfma.Threads(1);
#endif
    //    bfma.Verbose(0);

    //Physics parameters
    bfmActionParams *bfmap = (bfmActionParams *) &bfma;
    *bfmap = invParam.BAP;
    
    // Algorithm & code control
    bfma.time_report_iter=-100;
    bfma.max_iter     = invParam.MaxIter;
  
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
    
    QDPIO::cout << "Initialising BAGEL-2 solver "<<endl;
    
    // Bfm object
    bfm_qdp<Float> bfm; 
    bfm.init(bfma);

    //Gauge field import
    bfm.importGauge(links);
    
    //begin karthee clover
    //Missing clover term
    if ( invParam.BAP.solver  == CloverFermion) {
      // Import the ClovDiag, ClovOffDiag
      // Import the ClovInvDiag, ClovInvOffDiag
      bfm.importClover(CloverDiag,CloverOffDiag, bfm.A);
      bfm.importClover(CloverInvDiag,CloverInvOffDiag, bfm.Ainv);
    }
    //end karthee clover

    int cb = 1;
    for(int shift=0;shift<nshift;shift++){
      sol_t[shift] = bfm.allocFermion();
      bfm.importFermion(psi[shift],sol_t[shift],cb);
    }
    src_t = bfm.allocFermion();
    bfm.importFermion(src,src_t,cb);
    
    // Run the inverter
    int iter;

    QDPIO::cout << "Calling BAGEL multishift inverter"<<endl;
#pragma omp parallel for  
    for(int i=0;i<bfm.nthread;i++) {
      int iter_thr = bfm.CGNE_prec_MdagM_multi_shift(sol_t,src_t,
						     dshifts,
						     alpha,
						     nshift,
						     mresidual,
						     dontsum);
      if ( i==0 ) iter = iter_thr;
    }
    res = iter;
    
    
    for(int shift=0;shift<nshift;shift++) {
      bfm.exportFermion(psi[shift],sol_t[shift],cb);
      bfm.freeFermion(sol_t[shift]);
    }
    bfm.freeFermion(src_t);
    bfm.end();
    return res;
  }

  void QuarkProp(LatticeFermion &result,const LatticeFermion &source)
  {
    
  }

template<class Float>  void BfmWrapper::bfmTwoFlavorForce5d(multi1d<T4> &phi,multi1d<LatticeColorMatrix> &force)
{
  exit(-1);
}
template<class Float>  void BfmWrapper::bfmTwoFlavorRatioForce5d(multi1d<T4> &phi,multi1d<LatticeColorMatrix> &force,
								 Real  M_pv, Real  M_f)
{
  double m_pv = toDouble(M_pv);
  double m_f  = toDouble(M_f );
    // Set up BAGEL object
    int lx = QDP::Layout::subgridLattSize()[0];
    int ly = QDP::Layout::subgridLattSize()[1];
    int lz = QDP::Layout::subgridLattSize()[2];
    int lt = QDP::Layout::subgridLattSize()[3];
    
    bfmarg bfma;

#if defined(QDP_USE_OMP_THREADS)
    bfma.Threads(omp_get_max_threads());
#else
    bfma.Threads(1);
#endif
    //    bfma.Verbose(0);

    //Physics parameters
    bfmActionParams *bfmap = (bfmActionParams *) &bfma;
    *bfmap = invParam.BAP;
    
    // Algorithm & code control
    bfma.time_report_iter=-100;
    bfma.max_iter     = invParam.MaxIter;
    bfma.residual     = toDouble(invParam.RsdTarget[0]); 
    QDPIO::cout << "BfmWrapper: residual "<< bfma.residual<<endl;
  
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
    
    // Bfm object
    bfm_qdp<Float> M; 
    bfm_qdp<Float> PV; 
    bfma.mass = m_f;
    M.init(bfma);

    bfma.mass = m_pv;
    PV.init(bfma);


    //Gauge field import
    M.importGauge(links);
    PV.importGauge(links);

    Matrix_t force_t[2];
    Fermion_t phi_t;
    phi_t = M.allocFermion();
    force_t[0] = M.allocMatrix();
    force_t[1] = M.allocMatrix();

#pragma omp parallel for  
    for(int i=0;i<M.nthread;i++) {
      int cb=1;
      M.importFermion(phi,phi_t,cb); // Odd odd force term

      this->bfmTwoFlavorRatioForce<Float>(M,PV,phi_t, force_t);

      M.exportForce(force_t[0],force,0);
      M.exportForce(force_t[1],force,1);
    }

    M.freeFermion(phi_t);
    M.freeMatrix(force_t[0]);
    M.freeMatrix(force_t[1]);
    M.end();
    PV.end();

}

    // Possible single prec
 
template<class Float>  void BfmWrapper::bfmTwoFlavorRatioForce(bfm_qdp<Float> &M, bfm_qdp<Float> &PV, 
							       Fermion_t phi, Matrix_t force[2])
{
  Fermion_t X  = M.threadedAllocFermion();
  Fermion_t Y  = M.threadedAllocFermion();
  Fermion_t tmp1= M.threadedAllocFermion();
  Fermion_t tmp2= M.threadedAllocFermion();

  // S = phi^dag V (Mdag M)^-1 V^dag  phi
  // dS/du = phi^dag dV (Mdag M)^-1 V^dag  phi
  //       - phi^dag V (Mdag M)^-1 [ Mdag dM + dMdag M ]  (Mdag M)^-1 V^dag  phi
  //       + phi^dag V (Mdag M)^-1 dV^dag  phi

  int dagno=0;
  int dagyes=1;

  int me = M.thread_barrier();

  M.zeroMatrix(force[0]);
  M.zeroMatrix(force[1]);
  
  if ( M.isBoss() && (!me) ) {
    printf("TwoFlavorRatioForce target residual %le sizeoffloat %d\n",M.residual,sizeof(Float));
  }
  
  double n;
#define PRINT_VEC(A) n= PV.norm(A) ; if ( M.isBoss() && ( !me ) ) { printf("%s %le \n",#A,n);};
  //Get X = (Mdag M)^-1_oo V^dag_oo  phi
  PV.Mprec(phi,tmp2,tmp1,dagyes,0);

  M.fill(X,0.0);
  if ( invParam.BfmPrecision == Bfm32bitDefectCorrect ) { 
    Fermion_t src[2];
    src[0] = src[1] = tmp2;
    Fermion_t sol[2];
    sol[1] = X;
    sol[0] = tmp1;
    bfmInvert_mixed<Float>(M,sol,src);
  } else { 
    M.CGNE_prec(X, tmp2); 
  }
  //Y = (Mdag)^-1 V^dag  phi
  M.Mprec(X,Y,tmp2,dagno,0);

  // phi^dag V (Mdag M)^-1 dV^dag  phi
  PV.MprecDeriv(X,phi,force,dagyes,1.0);
  
  // phi^dag dV (Mdag M)^-1 V^dag  phi
  PV.MprecDeriv(phi,X,force,dagno,1.0);

  //        phi^dag V (Mdag M)^-1 Mdag dM   (Mdag M)^-1 V^dag  phi
  //        phi^dag V (Mdag M)^-1 dMdag M   (Mdag M)^-1 V^dag  phi
  // ... need a relative sign here. Force it into X.
  M.MprecDeriv(Y,X,force,dagno,-1.0);
  M.MprecDeriv(X,Y,force,dagyes,-1.0);

  M.threadedFreeFermion(X);
  M.threadedFreeFermion(Y);
  M.threadedFreeFermion(tmp2);
  M.threadedFreeFermion(tmp1);

}
// S_f = chi^dag* P(V^dag*V)/Q(V^dag*V)* N(M^dag*M)/D(M^dag*M)* P(V^dag*V)/Q(V^dag*V)* chi
//
// Here, M is some 5D operator and V is the Pauli-Villars field
// N and D makeup the rat. poly of the M term and P and & makeup the rat.poly of the denom term
//
// Need
// dS_f/dU =  chi^dag d[P/Q]  N/D   P/Q  chi
//         +  chi^dag   P/Q d[N/D]  P/Q  chi
//         +  chi^dag   P/Q   N/D d[P/Q] chi
//
// Here P/Q \sim R_{1/4}  ~ (V^dagV)^{1/4}
// Here N/D \sim R_{-1/2} ~ (M^dagM)^{-1/2}
//
// P/Q is expressed as partial fraction expansion:
//
//           a0 + \sum_k ak/(V^dagV + bk)
//
// d[P/Q] is then
//
//          \sum_k -ak [V^dagV+bk]^{-1}  [ dV^dag V + V^dag dV ] [V^dag V + bk]^{-1}
//
// and similar for N/D.
// 
// Need
//       MpvPhi_k   = [Vdag V + bk]^{-1} chi
//
//       MpvPhi     = {a0 +  \sum_k ak [Vdag V + bk]^{-1} }chi
//
//       MfMpvPhi_k = [MdagM+bk]^{-1} MpvPhi
//      
//       MfMpvPhi   = {a0 +  \sum_k ak [Mdag M + bk]^{-1} } MpvPhi
//
//       MpvMfMpvPhi_k = [Vdag V + bk]^{-1} MfMpvchi
//
// With these building blocks
//
//       dS/dU =  
//                 \sum_k -ak MpvPhi_k^dag        [ dV^dag V + V^dag dV ] MpvMfMpvPhi_k           <- deriv on P left
//             +   \sum_k -ak MpvMfMpvPhi_k^\dag  [ dV^dag V + V^dag dV ] MpvPhi_k
//             +   \sum_k -ak MfMpvPhi_k^dag      [ dM^dag M + M^dag dM ] MfMpvPhi_k
//

template<class Float>  void BfmWrapper::bfmOneFlavorRatioRationalForce(bfm_qdp<Float> &M, bfm_qdp<Float> &PV, 
								       Fermion_t phi, 
								       Matrix_t force[2],
								       int n_pv, double a0_pv, double ak_pv[],double bk_pv[], double r_pv[],
								       int n_f , double a0_f , double ak_f[] ,double bk_f[], double r_f[])
{
  Fermion_t      MpvPhi_k[n_pv];
  Fermion_t    MfMpvPhi_k[n_f];
  Fermion_t MpvMfMpvPhi_k[n_pv];

  Fermion_t      MpvPhi;
  Fermion_t    MfMpvPhi;
  Fermion_t MpvMfMpvPhi;
  Fermion_t           Y;
  Fermion_t         tmp;

  double one_pv[n_pv];
  double one_f[n_f];

  for(int k=0;k<n_pv;k++) one_pv[k] = 1.0;
  for(int k=0;k<n_f ;k++) one_f[k]  = 1.0;

  //Allocate
  for(int k=0;k<n_pv;k++){
         MpvPhi_k[k] = M.threadedAllocFermion();
    MpvMfMpvPhi_k[k] = M.threadedAllocFermion();
  }
  for(int k=0;k<n_f;k++){
    MfMpvPhi_k[k] = M.threadedAllocFermion();
  }

  MpvPhi       = M.threadedAllocFermion();
  MfMpvPhi     = M.threadedAllocFermion();
  MpvMfMpvPhi  = M.threadedAllocFermion();
  Y            = M.threadedAllocFermion();
  tmp          = M.threadedAllocFermion();

  int me = PV.thread_barrier();


  // The multishift inversions
  PV.CGNE_prec_MdagM_multi_shift(MpvPhi_k,
				 phi,
				 bk_pv,
				 one_pv,
				 n_pv,
				 r_pv,
				 0);
  PV.axpby(MpvPhi,phi,phi,0.0,a0_pv);
  for(int k=0; k < n_pv; k++) { 
    PV.axpby(MpvPhi,MpvPhi,MpvPhi_k[k],1.0,ak_pv[k]);
  }
#define RESID_CHECK
#ifdef RESID_CHECK
  double n_sum;
  double n_k;
  double e_k;
  
  n_sum=PV.norm(MpvPhi);
  n_sum=sqrt(n_sum);
  for(int k=0; k < n_pv; k++) { 
    n_k = PV.norm(MpvPhi_k[k]);
    n_k = sqrt(n_k);
    e_k = r_pv[k] * ak_pv[k] * n_k;
    if(PV.isBoss() && (!me)){
      printf("PV inverse error [%d] = %le  / %le (%le,%le)\n",k,e_k,n_sum,n_k,ak_pv[k]);
    }
  }
#endif
  
  M.CGNE_prec_MdagM_multi_shift(MfMpvPhi_k,
			      MpvPhi,
			      bk_f,
			      one_f,
			      n_f,
			      r_f,
			      0);
  M.axpby(MfMpvPhi,MpvPhi,MpvPhi,0.0,a0_f);
  for(int k=0; k < n_f; k++) { 
    M.axpby(MfMpvPhi,MfMpvPhi,MfMpvPhi_k[k],1.0,ak_f[k]);
  }
#ifdef RESID_CHECK
  n_sum=PV.norm(MpvPhi);
  n_sum=sqrt(n_sum);
  for(int k=0; k < n_f; k++) { 
    n_k = PV.norm(MpvPhi_k[k]);
    n_k = sqrt(n_k);
    e_k = r_f[k] * ak_f[k] * n_k;
    if(PV.isBoss() && (!me)){
      printf("Mf inverse error [%d] = %le  / %le (%le,%le)\n",k,e_k,n_sum,n_k,ak_f[k]);
    }
  }
#endif

  PV.CGNE_prec_MdagM_multi_shift(MpvMfMpvPhi_k,
				 MfMpvPhi,
				 bk_pv,
				 one_pv,
				 n_pv,
				 r_pv,
				 0);
  PV.axpby(MpvMfMpvPhi,MfMpvPhi,MfMpvPhi,0.0,a0_pv);
  for(int k=0; k < n_pv; k++) { 
    PV.axpby(MpvMfMpvPhi,MpvMfMpvPhi,MpvMfMpvPhi_k[k],1.0,ak_pv[k]);
  }

  int dagno=0;
  int dagyes=1;
  // Loop over poles in PV sector
  for(int k=0; k < n_pv; k++) {

    //                 \sum_k -ak MpvPhi_k^dag        [ dV^dag V + V^dag dV ] MpvMfMpvPhi_k           <- deriv on P left
    //             +   \sum_k -ak MpvMfMpvPhi_k^\dag  [ dV^dag V + V^dag dV ] MpvPhi_k

    double coeff = -ak_pv[k];
    PV.Mprec(MpvPhi_k[k],Y,tmp,dagno,0);
    PV.MprecDeriv(MpvMfMpvPhi_k[k],Y,force,dagyes, coeff);

    PV.Mprec(MpvMfMpvPhi_k[k],Y,tmp,dagno);   // V as we take Ydag 
    PV.MprecDeriv(Y, MpvPhi_k[k],force, dagno, coeff);        // dV

    PV.Mprec(MpvMfMpvPhi_k[k],Y,tmp,dagno,0);
    PV.MprecDeriv(MpvPhi_k[k], Y,force,dagyes, coeff);
      
    PV.Mprec(MpvPhi_k[k], Y, tmp, dagno,0);   // V as we take Ydag 
    PV.MprecDeriv(Y, MpvMfMpvPhi_k[k],force, dagno, coeff);        // dV

  }

  //             +   \sum_k -ak MfMpvPhi_k^dag      [ dM^dag M + M^dag dM ] MfMpvPhi_k
  for(int k=0; k < n_f; k++){

    double coeff = -ak_f[k];
    // The  d(M^dag)*M  term
    M.Mprec(MfMpvPhi_k[k],Y,tmp,dagno);               
    M.MprecDeriv(MfMpvPhi_k[k],Y,force,dagyes, coeff);

    // The  M^dag*d(M)  term
    M.Mprec(MfMpvPhi_k[k],Y,tmp,dagno);               
    M.MprecDeriv(Y,MfMpvPhi_k[k],force,dagno, coeff);

  }

  //Free vectors
  for(int k=0;k<n_pv;k++){
    M.threadedFreeFermion(MpvPhi_k[k]);
    M.threadedFreeFermion(MpvMfMpvPhi_k[k]);
  }
  for(int k=0;k<n_f;k++){
    M.threadedFreeFermion(MfMpvPhi_k[k]);
  }

  M.threadedFreeFermion(MpvPhi);
  M.threadedFreeFermion(MfMpvPhi);
  M.threadedFreeFermion(MpvMfMpvPhi);
  M.threadedFreeFermion(Y);
  M.threadedFreeFermion(tmp);

}

void BfmWrapper::RationalSanityCheck(double power, int n_p, double a0, double ak[], double bk[])
{
  double x=0.4134;

  double partfrac = a0;
  for(int i=0; i < n_p; i++){
    partfrac+= ak[i]/(x+bk[i]);
  }
  double exact=pow(x,power);
  QDPIO::cout << "Rational Sanity Check for power "<<power << " " << exact << " " << partfrac << endl;
  if ( fabs(partfrac-exact) > 1.0e-5 ) {
    printf("Insane Rational pow(%le,%le) = %le  : partfrac  %le\n",x,power,pow(x,power),partfrac);
    printf("a0 %le\n",a0);
    for(int i=0; i < n_p; i++){
      printf("a%d %le b%d %le\n",i+1,ak[i],i+1,bk[i]);
    }
    exit(-1);
  }
}

 
template<class Float>  
void BfmWrapper::bfmOneFlavorRatioRationalForce(multi1d<T4> &phi,
						multi1d<LatticeColorMatrix> &force,
						Real M_pv, Real M_f,
						Real nrm_pv,multi1d<Real> &shifts_pv, multi1d<Real> &residue_pv,
						Real nrm_f,multi1d<Real> &shifts_f, multi1d<Real> &residue_f
						)
{
  double m_pv = toDouble(M_pv);
  double m_f  = toDouble(M_f );
  int n_pv = shifts_pv.size();
  int n_f  = shifts_f.size();

  double  ak_f[n_f];
  double  bk_f[n_f];

  double  ak_pv[n_pv];
  double  bk_pv[n_pv];

  double  rk_pv[n_pv];
  double  rk_f[n_f];

  double a0_f  = toDouble(nrm_f);
  double a0_pv = toDouble(nrm_pv);
  double residual = toDouble(invParam.RsdTarget[0]); 

 for(int pole=0;pole<n_f;pole++){
    ak_f[pole]=toDouble(residue_f[pole]);
    bk_f[pole]=toDouble(shifts_f[pole]);
    rk_f[pole] =  residual;
  }
 for(int pole=0;pole<n_pv;pole++){
    ak_pv[pole]=toDouble(residue_pv[pole]);
    bk_pv[pole]=toDouble(shifts_pv[pole]);
    rk_pv[pole] = residual;
  }
 for(int k=0;k<4;k++){
   rk_f[k] *= 10;
   rk_pv[k] *= 10;
 }
 rk_f[0] *= 2;

  RationalSanityCheck(-0.5,n_f,a0_f,ak_f,bk_f);
  RationalSanityCheck(0.25,n_pv,a0_pv,ak_pv,bk_pv);
  QDPIO::cout << "BfmWrapper: Checked the rational functions look rational" <<endl;

    // Set up BAGEL object
    int lx = QDP::Layout::subgridLattSize()[0];
    int ly = QDP::Layout::subgridLattSize()[1];
    int lz = QDP::Layout::subgridLattSize()[2];
    int lt = QDP::Layout::subgridLattSize()[3];
    
    bfmarg bfma;

#if defined(QDP_USE_OMP_THREADS)
    bfma.Threads(omp_get_max_threads());
#else
    bfma.Threads(1);
#endif
    //    bfma.Verbose(0);

    //Physics parameters
    bfmActionParams *bfmap = (bfmActionParams *) &bfma;
    *bfmap = invParam.BAP;
    
    // Algorithm & code control
    bfma.time_report_iter=-100;
    bfma.max_iter     = invParam.MaxIter;
    bfma.residual     = residual;
    QDPIO::cout << "BfmWrapper: residual "<< bfma.residual<<endl;
  
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
    
    // Bfm object
    bfm_qdp<Float> M; 
    bfm_qdp<Float> PV; 
    bfma.mass = m_f;
    M.init(bfma);
    bfma.mass = m_pv;
    PV.init(bfma);

    //Gauge field import
    M.importGauge(links);
    PV.importGauge(links);

    Matrix_t force_t[2];
    force_t[0] = M.allocMatrix();   
    force_t[1] = M.allocMatrix();    

    Fermion_t phi_t;
    phi_t = M.allocFermion();

    static int range_check;
    if ( range_check == 0 ) { 
      QDPIO::cout << "BfmWrapper: checking spectral range "<<endl;
      this->SpectralRange<Float>(M);
      QDPIO::cout << "BfmWrapper: checking spectral range "<<endl;
      this->SpectralRange<Float>(PV);
      range_check=1;
    }

#pragma omp parallel for  
    for(int i=0;i<M.nthread;i++) {
      int cb=1;
      M.importFermion(phi,phi_t,cb); // Odd odd force term
      M.zeroMatrix(force_t[0]);
      M.zeroMatrix(force_t[1]);



      this->bfmOneFlavorRatioRationalForce<Float>(M, PV, 
						  phi_t,force_t,
						  n_pv,a0_pv,ak_pv,bk_pv,rk_pv,
						  n_f ,a0_f ,ak_f ,bk_f,rk_f);

      M.exportForce(force_t[0],force,0);
      M.exportForce(force_t[1],force,1);
    }

    M.freeFermion(phi_t);
    M.freeMatrix(force_t[0]);
    M.freeMatrix(force_t[1]);
    M.end();
    PV.end();
}


template void BfmWrapper::bfmTwoFlavorForce5d<float>(multi1d<T4> &phi,multi1d<LatticeColorMatrix> &force);
template void BfmWrapper::bfmTwoFlavorRatioForce5d<float>(multi1d<T4> &phi,multi1d<LatticeColorMatrix> &force,Real  M_pv, Real  M_f);

template void BfmWrapper::bfmTwoFlavorForce5d<double>(multi1d<T4> &phi,multi1d<LatticeColorMatrix> &force);
template void BfmWrapper::bfmTwoFlavorRatioForce5d<double>(multi1d<T4> &phi,multi1d<LatticeColorMatrix> &force,Real  M_pv, Real  M_f);


template void BfmWrapper::bfmOneFlavorRatioRationalForce<float>(multi1d<T4> &phi,
								multi1d<LatticeColorMatrix> &force,
								Real M_pv, Real M_f,
								Real,multi1d<Real> &shifts_pv, multi1d<Real> &residue_pv,
								Real,multi1d<Real> &shifts_f, multi1d<Real> &residue_f
								);

template void BfmWrapper::bfmOneFlavorRatioRationalForce<double>(multi1d<T4> &phi,
								multi1d<LatticeColorMatrix> &force,
								Real M_pv, Real M_f,
								 Real,multi1d<Real> &shifts_pv, multi1d<Real> &residue_pv,
								 Real,multi1d<Real> &shifts_f, multi1d<Real> &residue_f
								);

  template  int BfmWrapper::bfmInvert5d<float> (multi1d<T4> &sol,const multi1d<T4> &chi);
  template  int BfmWrapper::bfmInvert5d<double>(multi1d<T4> &sol,const multi1d<T4> &chi);
  template  int BfmWrapper::bfmInvert4d<float>(T4 &sol,const T4 &chi);
  template  int BfmWrapper::bfmInvert4d<double>(T4 &sol,const T4 &chi);
template  int BfmWrapper::bfmCayleyDminus<float> (multi1d<T4> &sol);
template  int BfmWrapper::bfmCayleyDminus<double>(multi1d<T4> &sol);

  template  int BfmWrapper::bfmMultiInvert5d<float>(multi1d< multi1d<T4> >& psi, 
						     const multi1d<Real>& shifts, 
						     const multi1d<Real>& Residuals, 
						     const multi1d<T4>& chi);
  template  int BfmWrapper::bfmMultiInvert5d<double>(multi1d< multi1d<T4> >& psi, 
						     const multi1d<Real>& shifts, 
						     const multi1d<Real>& Residuals, 
						     const multi1d<T4>& chi);

  template  int BfmWrapper::bfmMultiInvert4d<float>(multi1d<T4> & psi, 
						     const multi1d<Real>& shifts, 
						     const multi1d<Real>& Residuals, 
						     const T4& chi);
  template  int BfmWrapper::bfmMultiInvert4d<double>(multi1d<T4> & psi, 
						     const multi1d<Real>& shifts, 
						     const multi1d<Real>& Residuals, 
						     const T4& chi);

