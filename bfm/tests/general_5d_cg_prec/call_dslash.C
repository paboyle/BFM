#include <bfm.h>
#include <bfm_qdp.h>
#include <bfm_qdp_cayley_prec.h>
#include <bfm_qdp_dwf.h>

#define Printf if ( QMP_is_primary_node() ) printf

#include <bfm_qdp_chroma_linop.h>

#include <bfmcommqmp.h>
#include <bfmcommspi.h>

typedef bfmcommQMP<double> bfm_qmp;

void test_solver(BfmSolver solver, LatticeFermion &src);


int main (int argc,char **argv )
{
  Chroma::initialize(&argc,&argv);

  WilsonTypeFermActs4DEnv::registerAll(); 
  LinOpSysSolverArrayEnv::registerAll();

  bfmarg::UseCGdiagonalMee(1);

  /********************************************************
   * Command line parsing
   ********************************************************
   */
#define COMMANDLINE
#ifdef COMMANDLINE
  if ( argc != 6 ) { 
   Printf("Usage: %s lx ly lz lt Ls\n All must be even\n",argv[0]);
   Printf("argc is %d\n",argc);
   for ( int i=0;i<argc;i++){
      Printf("%d %s\n",i,argv[i]);
   }
      exit(-1);
  }
#endif

  /********************************************************
   * Setup QDP
   ********************************************************
   */
  multi1d<int> nrow(Nd);
#ifdef COMMANDLINE
  nrow[0] = atoi(argv[1]);
  nrow[1] = atoi(argv[2]);
  nrow[2] = atoi(argv[3]);
  nrow[3] = atoi(argv[4]);
#else
  nrow[0] = 4;
  nrow[1] = 2;
  nrow[2] = 2;
  nrow[3] = 4;
#endif

  bfmarg::Threads(1);
  bfmarg::Verbose(BfmMessage|BfmPerformance);
  Layout::setLattSize(nrow);
  Layout::create();


  LatticeFermion src;
  gaussian(src);

  for(int i=0;i<1;i++){
    //    test_solver(HtCayleyZolo,src);
    //    test_solver(HwCayleyZolo,src);
    test_solver(HtCayleyTanh,src);
    //    test_solver(HmCayleyTanh,src);
    //test_solver(HwCayleyTanh,src);
    //    test_solver(HwPartFracZolo,src);
    //    test_solver(HwContFracZolo,src);
  }
}

void test_solver(BfmSolver solver,LatticeFermion &src)
{

  g5dParams parms;

  int Ls=8;
  double M5=1.8;
  double mq=0.01;
  double wilson_lo = 0.05;
  double wilson_hi = 6.8;
  double shamir_lo = 0.025;
  double shamir_hi = 1.7;
  double ht_scale=1.7;
  double hw_scale=1.0;

  if ( solver == HtCayleyTanh ) { 
    Printf("Testing HtCayleyTanh aka DWF\n");
    parms.ShamirCayleyTanh(mq,M5,Ls);
  } else if ( solver == HmCayleyTanh ) { 
    parms.ScaledShamirCayleyTanh(mq,M5,Ls,ht_scale);
    Printf("Testing HmCayleyTanh Moebius\n");
  } else if ( solver == HwCayleyTanh ) { 
    parms.WilsonCayleyTanh(mq,M5,Ls,hw_scale);
    Printf("Testing HwCayleyTanh aka Borici\n");
  } else if ( solver == HtCayleyZolo ) { 
    parms.ShamirCayleyZolo(mq,M5,Ls,shamir_lo,shamir_hi);
    Printf("Testing HtCayleyZolo\n");
  } else if ( solver == HwCayleyZolo ) { 
    parms.WilsonCayleyZolo(mq,M5,Ls,wilson_lo,wilson_hi);
    Printf("Testing HwCayleyZolo aka Chiu\n");
  } else if ( solver == HwPartFracZolo ) { 
    Printf("Testing HwPartFracZolo aka KEK\n");
    Ls = 9;
    parms.WilsonPartFracZolo(mq,M5,Ls,wilson_lo,wilson_hi);
    Printf("Partial fraction: setting Ls =%d\n",Ls);
  } else if ( solver == HwContFracZolo ) { 
    Printf("Testing HwContFracZolo \n");
    Ls = 9;
    parms.WilsonContFracZolo(mq,M5,Ls,wilson_lo,wilson_hi);
  } else { 
    exit(-1);
  }

  multi1d<LatticeColorMatrix> u(4);
  HotSt(u);

  LatticeFermion sol;
  LatticeFermion sol_sp;
  LatticeFermion check;
  check=zero;

  int ChromaRegress = 0;

  ////////////////////////////////
  // check the results with Chroma
  ////////////////////////////////
  if( ChromaRegress ) {
    QDPIO::cout << "===================== " <<endl;
    Printf("Getting Chroma Linop\n"); fflush(stdout);
    Handle< SystemSolver<LatticeFermion> > sysolv =    GetSolver(u, parms);

    Printf("Solving with CHROMA to create regression results ...\n ...and impress you with how fast and wonderful BAGEL is\n"); fflush(stdout);
    
    SystemSolverResults_t r=(*sysolv)(check,src);

  QDPIO::cout << "===================== "<<r.resid <<endl;
  }

  ////////////////////////////////
  // solve with BAGEL
  ////////////////////////////////
  LatticeComplex PP,PA,PAc,PJ5;


  QDPIO::cout << "===================== "<<endl;
  QDPIO::cout << "DOUBLE Source normalisation is "<< norm2(src) <<endl;
  QDPIO::cout << "===================== " <<endl;
  sol=zero;
  bfm_Cayley_CG<double>(sol,src,u,PP,PA,PAc,PJ5,parms,1.0e-18,80000);
  QDPIO::cout << "Solved with new  DWF sol    ="<<norm2(sol)<<endl;


  QDPIO::cout << "Solved with DOUBLE  bfm solution norm    ="<<norm2(sol)<<endl;
  Real n2diff;
  if ( ChromaRegress ) {
    n2diff = norm2(check-sol);
    QDPIO::cout << "Comparison: sol - check = " <<n2diff<<endl;
    QDPIO::cout << "Comparison: sol = "<<norm2(sol)<< " check = " <<norm2(check)<<endl;
    if( toDouble(n2diff)/toDouble(norm2(sol)) > 1.0e-8 ) { 
      QDP_error_exit("Suspicous difference");
    } else { 
      QDPIO::cout << "looks ok"<<endl;
    }
  }

  QDPIO::cout << "===================== "<<endl;
  QDPIO::cout << "SINGLE Source normalisation is "<< norm2(src) <<endl;
  QDPIO::cout << "===================== "<<endl;
  sol_sp=zero;
  bfm_Cayley_CG<float>(sol_sp,src,u,PP,PA,PAc,PJ5,parms,1.0e-8,80000);

  if ( ChromaRegress ) {
    QDPIO::cout << "Solved with SINGLE  bfm solution norm    ="<<norm2(sol_sp)<<endl;
    n2diff = norm2(check-sol_sp);
    QDPIO::cout << "Comparison: sol - check = " <<n2diff<<endl;
    QDPIO::cout << "Comparison: sol = "<<norm2(sol_sp)<< " check = " <<norm2(check)<<endl;
    if( toDouble(n2diff)/toDouble(norm2(sol_sp)) > 1.0e-4 ) { 
    QDP_error_exit("Suspicous difference");
    } else { 
      QDPIO::cout << "looks ok"<<endl;
      
    }
  }
  QDPIO::cout << "===================== " <<endl;
  QDPIO::cout << "The end " <<endl;
  QDPIO::cout << "===================== " <<endl;
  
}






