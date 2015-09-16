#include <chroma.h>
#include <actions/ferm/invert/syssolver_linop_cg_array.h>
#include <bfm.h>
#include <bfm_qdp.h>

typedef bfm_qdp<double> bfm_t;

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
#undef COMMANDLINE
#ifdef COMMANDLINE
  if ( argc != 6 ) { 
   Printf("Usage: %s lx ly lz lt Ls\n All must be even\n",argv[0]);
   Printf("argc is %d\n",argc);
   for ( int i=0;i<argc;i++){
      Printf("%d %s\n",i,argv[i]);
   }
      exit
(-1);

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
  int Ls = atoi(argv[5]);
#else
  nrow[0] = 4;
  nrow[1] = 4;
  nrow[2] = 4;
  nrow[3] = 4;
  int Ls  = 4;
#endif

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
  dwfa.solver = DWF;
  dwfa.node_latt[0]  = lx;
  dwfa.node_latt[1]  = ly;
  dwfa.node_latt[2]  = lz;
  dwfa.node_latt[3]  = lt;
  dwfa.verbose=1;
  dwfa.reproduce=1;
  bfmarg::Threads(64);
  bfmarg::Reproduce(0);
  bfmarg::ReproduceChecksum(0);
  bfmarg::ReproduceMasterCheck(0);
  bfmarg::Verbose(1);

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

  dwfa.precon_5d = 1;
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
  HotSt(u);
  multi1d<LatticeColorMatrixF> uf(Nd);
  for(int m=0; m < Nd; ++m){
    uf[m] = u [m];
    u [m] = uf[m];
  }
  printf("Setup gauge field\n");
  fflush(stdout);
  /********************************************************
   * Gaussian source and result vectors
   ********************************************************
   */
  multi1d<LatticeFermion> source(Ls);
  for(int s=0;s<Ls;s++) gaussian(source[s]);
  Fermion_t src = dwf.allocFermion();
  Fermion_t test= dwf.allocFermion();
  printf("Filling with zeroes field\n");
  fflush(stdout);
  dwf.master_fill(test,0.0);

#define NMULTI (3)
   Fermion_t psi[NMULTI];
   double masses[NMULTI]   = {0.0,0.1,0.2};
   double alpha [NMULTI]    = {1.0,1.0,1.0};
   double residuals[NMULTI] = {1.0e-8,1.0e-8,1.0e-8};
   psi[0] = dwf.allocFermion();
   psi[1] = dwf.allocFermion();
   psi[2] = dwf.allocFermion();
   dwf.master_fill(psi[0],0.0);
   dwf.master_fill(psi[1],0.0);
   dwf.master_fill(psi[2],0.0);

  /********************************************************
   * Import gauge field to BAGEL
   ********************************************************
   */
  dwf.importGauge(u);
  dwf.importFermion(source,src,0);

  printf("Calling half cb inverter\n"); fflush(stdout);
  dwf.inv_type=CG_PREC_MDAGM;
  dwf.qdp_chi_h[0]=psi[0];
  dwf.qdp_chi_h[1]=psi[0];
  dwf.qdp_psi_h[0]=src;
  dwf.qdp_psi_h[1]=src;
  bfm_spawn_cg(dwf);

  printf("Calling multi-shift inverter\n");fflush(stdout);
  dwf.inv_type=CG_PREC_MDAGM_MULTI;
  dwf.qdp_chi_multi_h=psi;
  dwf.shifts=masses;
  dwf.alpha =alpha;
  dwf.nshift=NMULTI;
  dwf.mresidual=residuals;
  dwf.single=0;
  bfm_spawn_cg(dwf);
  //  dwf.CGNE_prec_MdagM_multi_shift(psi,
  //				  src,
  //				  masses,
  //				  alpha,
  //				  NMULTI,
  //				  residuals,
  //				  0);

  Printf("Done\n"); 

}








