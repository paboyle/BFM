#include <bfm.h>
#include <bfm_qdp.h>
#include <bfm_qdp_cayley_prec.h>
#include <bfm_qdp_dwf.h>

#define Printf if ( QMP_is_primary_node() ) printf

#include <bfm_qdp_chroma_linop.h>

#include <bfmcommqmp.h>
#include <bfmcommspi.h>

typedef LatticeFermion T;
typedef LatticeFermion Phi;
typedef multi1d<LatticeColorMatrix> U;
typedef multi1d<LatticeColorMatrix> P;
typedef multi1d<LatticeColorMatrix> Q;
typedef bfm bfm_t;

typedef bfmcommQMP<double> bfm_qmp;

void test_solver(BfmSolver solver, LatticeFermion &src);


int main (int argc,char **argv )
{
  Chroma::initialize(&argc,&argv);

  WilsonTypeFermActs4DEnv::registerAll(); 
  LinOpSysSolverArrayEnv::registerAll();

  bfmarg::Verbose(1);
  bfmarg::UseCGdiagonalMee(0);

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
  bfmarg::Verbose(0);
  Layout::setLattSize(nrow);
  Layout::create();


  LatticeFermion src;
  gaussian(src);

  for(int i=0;i<1;i++){
    test_solver(HmCayleyTanh,src);
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
  double ht_scale=2.0;
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
  } else { 
    exit(-1);
  }

  multi1d<LatticeColorMatrix> u(4);
  HotSt(u);

  multi1d<LatticeFermion> X(Ls);
  multi1d<LatticeFermion> Y(Ls);
  for(int s=0;s<Ls;s++){
    gaussian(X[s]);
    gaussian(Y[s]);
  }

  multi1d<LatticeColorMatrix> Fb(4);
  multi1d<LatticeColorMatrix> F (4);

  int lx = QDP::Layout::subgridLattSize()[0];
  int ly = QDP::Layout::subgridLattSize()[1];
  int lz = QDP::Layout::subgridLattSize()[2];
  int lt = QDP::Layout::subgridLattSize()[3];

  multi1d<int> procs = QDP::Layout::logicalSize();
  /********************************************************
   * Setup DWF operator
   ********************************************************
   */
  bfmarg  dwfa;
  dwfa.node_latt[0]  = lx;
  dwfa.node_latt[1]  = ly;
  dwfa.node_latt[2]  = lz;
  dwfa.node_latt[3]  = lt;
  for(int mu=0;mu<4;mu++){
    if ( procs[mu]>1 ) {
      dwfa.local_comm[mu] = 0;
    } else { 
      dwfa.local_comm[mu] = 1;
    }
  }

  bfm_qmp dwf_qmp;

  dwfa.mass = parms.mass;
  dwfa.M5   = parms.M5;
  dwfa.Csw = 0.0;
  dwfa.precon_5d=0;
  dwfa.residual=1.0e-5;
  dwfa.max_iter=5000;
  dwfa.Ls = parms.Ls;
  dwfa.solver= parms.solver;
  dwfa.zolo_lo   = parms.zolo_lo;
  dwfa.zolo_hi   = parms.zolo_hi;
  dwfa.mobius_scale = parms.mobius_scale;

  dwf_qmp.init(dwfa);
  dwf_qmp.importGauge(u);

  Fermion_t X_t;
  Fermion_t Y_t;
  Matrix_t  Force_t[2];

  X_t     = dwf_qmp.allocFermion();
  Y_t     = dwf_qmp.allocFermion();
  for(int cb=0;cb<2;cb++){
    Force_t[cb] = dwf_qmp.allocMatrix();
    printf("Force_t pointer %lx\n",Force_t[cb]);
  }


  multi1d<int> bcs(Nd); bcs[0] = bcs[1] = bcs[2] = bcs[3] = 1;
  Handle<FermBC<T,U,U> > fbc(new SimpleFermBC< T, U, U >(bcs));
  Handle<CreateFermState<T,U,U> > cfs( new CreateSimpleFermState<T,U,U>(fbc));
  Handle<FermState<T,U,U> > fs ((*cfs)(u));


  //! Params for NEFF
  EvenOddPrecKNOFermActArrayParams kparams;

  Real scale = dwf_qmp.mobius_scale;
  kparams.OverMass = dwf_qmp.M5;
  kparams.Mass     = dwf_qmp.mass;
  kparams.a5       = 1.0;
  kparams.coefs.resize(Ls);
  for(int s=0;s<Ls;s++)
    kparams.coefs[s] = scale;
  kparams.N5       = dwf_qmp.Ls;

  EvenOddPrecKNOFermActArray FA(cfs,kparams);
  Handle< DiffLinearOperatorArray<Phi,P,Q> > M(FA.linOp(fs));
  
  for( int dag=0;dag<2;dag++){

    Fb=zero;
    F=zero;

    dwf_qmp.importFermion(X,X_t,1);
    dwf_qmp.importFermion(Y,Y_t,1);

    dwf_qmp.zeroMatrix(Force_t[0]);
    dwf_qmp.zeroMatrix(Force_t[1]);


    dwf_qmp.MprecDeriv(X_t,Y_t,Force_t,dag);

    for(int cb=0;cb<2;cb++)
      dwf_qmp.exportForce(Force_t[cb],Fb,cb);
    
    Printf("Calling S_f deriv\n");
    if( dag==1 ) 
      M->deriv(F, X, Y, MINUS);
    else
      M->deriv(F, X, Y, PLUS);
      
    for(int mu=0;mu<4;mu++){
      QDPIO::cout << "dag "<< dag<<"mu "<<mu<<" n2diff " << norm2(Fb[mu]-F[mu])<<endl;
    }

#if 0
    QDPIO::cout << "Calling force term" << endl;
    dwf_qmp.zeroMatrix(Force_t[0]);
    dwf_qmp.zeroMatrix(Force_t[1]);
    dwf_qmp.TwoFlavorRatioForce(X_t,Force_t,1.0,dwf_qmp.mass);
    for(int cb=0;cb<2;cb++)
      dwf_qmp.exportForce(Force_t[cb],Fb,cb);
    for(int mu=0;mu<4;mu++){
      QDPIO::cout << "dag "<< dag<<"mu "<<mu<<" TwoFlavorForce " << norm2(Fb[mu])<<endl;
    }
#endif

    //    { 
    // EvenOddPrecConstDetTwoFlavorRatioConvConvWilsonTypeFermMonomial5D Monomial();
    //    }

  }

    dwf_qmp.freeMatrix(Force_t[0]);
    dwf_qmp.freeMatrix(Force_t[1]);
    exit(0);


}






