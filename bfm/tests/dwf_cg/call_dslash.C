#include <chroma.h>
#include <actions/ferm/invert/syssolver_linop_cg_array.h>
#include <bfm.h>

typedef bfm bfm_t;
//typedef commQMPbagel1<double> bfm_t;

using namespace Chroma;

typedef LatticeFermion T;
typedef multi1d<LatticeFermion> T5;
typedef multi1d<LatticeColorMatrix> U;

//#define Printf if ( QMP_is_primary_node() ) printf
#define Printf printf

int main (int argc,char **argv )
{
  Chroma::initialize(&argc,&argv);


  /********************************************************
   * Command line parsing
   ********************************************************
   */
  if ( argc != 7 ) { 
   Printf("Usage: %s lx ly lz lt Ls name\n All must be even\n",argv[0]);
   Printf("argc is %d\n",argc);
   for ( int i=0;i<argc;i++){
      Printf("%d %s\n",i,argv[i]);
   }
      exit
(-1);

  }
  /********************************************************
   * Setup QDP
   ********************************************************
   */
  multi1d<int> nrow(Nd);
  nrow[0] = atoi(argv[1]);
  nrow[1] = atoi(argv[2]);
  nrow[2] = atoi(argv[3]);
  nrow[3] = atoi(argv[4]);
  int Ls = atoi(argv[5]);



  Layout::setLattSize(nrow);
  Layout::create();

  int lx = QDP::Layout::subgridLattSize()[0];
  int ly = QDP::Layout::subgridLattSize()[1];
  int lz = QDP::Layout::subgridLattSize()[2];
  int lt = QDP::Layout::subgridLattSize()[3];

  /********************************************************
   * Setup DWF operator
   ********************************************************
   */
  bfmarg dwfa;
  bfm_t  dwf;
  dwfa.solver = DWFrb4d;
  dwfa.node_latt[0]  = lx;
  dwfa.node_latt[1]  = ly;
  dwfa.node_latt[2]  = lz;
  dwfa.node_latt[3]  = lt;
  dwfa.verbose=1;

  multi1d<int> procs = QDP::Layout::logicalSize();
  Printf("%d dim machine\n\t", procs.size());
  for(int mu=0;mu<4;mu++){
    Printf("%d ", procs[mu]);
    if ( procs[mu]>1 ) {
      dwfa.local_comm[mu] = 0;
    } else { 
      dwfa.local_comm[mu] = 1;
    }
  }
  Printf("\nLocal comm = ");
  for(int mu=0;mu<4;mu++){
    Printf("%d ", dwfa.local_comm[mu]);
  }
  Printf("\n");
  
  multi1d<int> ncoor = QDP::Layout::nodeCoord();

  Real M5(1.8);
  Real mq(0.1);

  dwfa.precon_5d = 0;
  dwfa.Ls   = Ls;
  dwfa.M5   = toDouble(M5);
  dwfa.mass = toDouble(mq);
  dwfa.Csw  = 0.0;
  dwfa.max_iter = 1000;
  dwfa.residual = 1.e-8;
  Printf("Initialising bfm operator\n");
  dwf.init(dwfa);

  /********************************************************
   * Gaussian gauge field
   ********************************************************
   */

  multi1d<LatticeColorMatrix>  u(Nd);
  //HotSt(u);
  ///Read Gauge field
  {
  	ArchivGauge_t Header ;
  	readArchiv(Header,u,argv[6]);
  }
  multi1d<LatticeColorMatrixF> uf(Nd);
  for(int m=0; m < Nd; ++m){
    uf[m] = u [m];
    u [m] = uf[m];
  }

  /********************************************************
   * Gaussian source and result vectors
   ********************************************************
   */
  multi1d<LatticeFermion> psi(Ls);
  multi1d<LatticeFermion> chi(Ls);
  multi1d<LatticeFermion> chi_qdp(Ls);

  for(int s=0;s<Ls;s++) gaussian(psi[s]);
  for(int s=0;s<Ls;s++) gaussian(chi[s]);
  for(int s=0;s<Ls;s++) gaussian(chi_qdp[s]);



  /********************************************************
   * Bagel internal single checkerboard vectors
   ********************************************************
   */
   Fermion_t psi_h[2];
   psi_h[0] = dwf.allocFermion();
   psi_h[1] = dwf.allocFermion();
   Fermion_t chi_h[2];
   chi_h[0] = dwf.allocFermion();
   chi_h[1] = dwf.allocFermion();

  /********************************************************
   * Import gauge field to BAGEL
   ********************************************************
   */
  dwf.importGauge(u);
  dwf.inner(psi_h[0], chi_h[0]);
#undef DEBUG

  /*Import this checkerboard of source field to bagel*/
  Printf("Importing psi field cb %d\n",1);
  dwf.importFermion(psi,psi_h[0],0);
  dwf.importFermion(psi,psi_h[1],1);
 
  // Fill the other checkerboard of result with noise
  for(int s=0;s<Ls;s++) gaussian(chi[s]);
  for(int s=0;s<Ls;s++) chi_qdp[s] = chi[s];
  dwf.importFermion(chi,chi_h[0],0);
  dwf.importFermion(chi,chi_h[1],1);

#define NITER 1
  for(int i=0;i<NITER;i++) {
    //    dwf.axpy(chi_h[0],psi_h[0],psi_h[0],0.0)  ;
    //    dwf.axpy(chi_h[1],psi_h[0],psi_h[0],0.0)  ;
    dwf.CGNE(chi_h,psi_h);
  }

  dwf.exportFermion(chi,chi_h[0],0);
  dwf.exportFermion(chi,chi_h[1],1);
 
#define CHECKIT
#ifdef  CHECKIT 
  /********************************************************
   * QDP Linop
   ********************************************************
   */
   multi1d<int> bcs(Nd);
   bcs[0] = bcs[1] = bcs[2] = bcs[3] = 1;

   Handle< FermBC<T,U,U> > fbc(new SimpleFermBC< T, U, U >(bcs));
   Handle<CreateFermState<T,U,U> > cfs( new CreateSimpleFermState<T,U,U>(fbc));
   EvenOddPrecDWFermActArray  S_f(cfs, M5, mq, Ls);
   Handle< FermState<T,U,U> > fs( S_f.createState(u) );

   Handle< EvenOddPrecLinearOperatorArray<T,U,U> > M(S_f.precLinOp(fs,mq));
   Handle< LinearOperatorArray<T> > HM(S_f.precLinOp(fs,mq));
  // Check the result
  SysSolverCGParams invParam;
  invParam.RsdCG        = 1.0e-8;
  invParam.RsdCGRestart = 1.0e-8;
  invParam.MaxCG        = dwfa.max_iter ;
  invParam.MaxCGRestart = dwfa.max_iter ;

  Handle<LinOpSystemSolverArray<T> > Hsolver = 
    new LinOpSysSolverCGArray<T>(HM,invParam);
  Handle<PrecFermAct5DQprop<T,U,U> > solver =
    new PrecFermAct5DQprop<T,U,U>(M,Hsolver);
  SystemSolverResults_t  stats =  (*solver)(chi_qdp,psi);
  int n_count = stats.n_count;

  Double n2 = 0.0;
  for(int s=0;s<Ls;s++) { 
    Double sn2 = norm2(chi[s]-chi_qdp[s]);
    n2+=sn2;
#ifdef DEBUG
    QDPIO::cout << "|| Bagel - QDP || ["<<s<<"] = "<< sn2 << endl;
#endif
  }
  QDPIO::cout << "|| Bagel - QDP || = "<< n2 << endl;
  
  n2 = 0.0;
  for(int s=0;s<Ls;s++) { 
    Double sn2 = norm2(chi[s]);
    n2+=sn2;
#ifdef DEBUG
    QDPIO::cout << "|| Bagel || ["<<s<<"] = "<< sn2 << endl;
#endif
  }
  QDPIO::cout << "|| Bagel || = "<< n2 << endl;
  
  n2 = 0.0;
  for(int s=0;s<Ls;s++) {
    Double sn2 = norm2(chi_qdp[s]);
    n2+=sn2;
#ifdef DEBUG
    QDPIO::cout << "|| QDP || ["<<s<<"] = "<< sn2 << endl;
#endif
  }
  QDPIO::cout << "|| QDP || = "<< n2 << endl;
#endif



}








