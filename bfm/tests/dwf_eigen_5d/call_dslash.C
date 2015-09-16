#include <chroma.h>
#include <actions/ferm/invert/syssolver_linop_cg_array.h>
#include <bfm.h>
#include <bfm_qdp.h>
#include <eigen/Krylov_5d.h>

typedef bfm bfm_t;
//typedef commQMPbagel1<double> bfm_t;

using namespace Chroma;

typedef LatticeFermion T;
typedef multi1d<LatticeFermion> T5;
typedef multi1d<LatticeColorMatrix> U;

typedef double number;

//#define Printf if ( QMP_is_primary_node() ) printf
#define Printf printf

int main (int argc,char **argv )
{
  Chroma::initialize(&argc,&argv);


  /********************************************************
   * Command line parsing
   ********************************************************
   */
  if ( argc < 10 ) { 
   Printf("Usage: %s lx ly lz lt Ls Nwanted Tn off rad\n All must be even\n",argv[0]);
   Printf("Good test: D_oo 4 4 4 4 4 10 30 1.0 6.0");
   Printf("argc is %d\n",argc);
   for ( int i=0;i<argc;i++){
      Printf("%d %s\n",i,argv[i]);
   }
      exit(-1);
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

  int threads = 16;
  bfmarg::Threads(threads);
  bfmarg::Reproduce(0);
  bfmarg::ReproduceChecksum(0);
  bfmarg::ReproduceMasterCheck(0);
  bfmarg::Verbose(0);

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
  bfmarg::Threads(threads);

  multi1d<int> ncoor = QDP::Layout::nodeCoord();

  Real M5(1.8);
  Real mq(0.01);
  dwfa.solver       = DWF; //DWFrb4d = 4d preconditioning. DWF = 5d preconditioning
  dwfa.precon_5d = 1;
  dwfa.Ls   = Ls;
  dwfa.M5   = toDouble(M5);
  dwfa.mass = toDouble(mq);
  dwfa.Csw  = 0.0;
  dwfa.max_iter = 100000;

  dwfa.residual = 1.e-10;
  Printf("Initialising bfm operator\n"); 

  bfm_qdp<number> dwf_qdp;
  dwf_qdp.init(dwfa);

  /********************************************************
   * Gaussian gauge field
   ********************************************************
   */

  multi1d<LatticeColorMatrix>  u(Nd);
  for(int i=0;i<Nd; i++){
   u[i] = 1.0; 
  }
  HotSt(u);
  ArchivGauge_t Header ; readArchiv(Header,u,argv[10]);

  /********************************************************
   * Anti-Periodic Boundary conditions
   ********************************************************
   */

    multi1d<int> bcond (Nd);
    bcond[0] = bcond[1] = bcond[2] = 1; bcond[3] = -1;
    Handle<FermBC<T,U,U> > fbca (new SimpleFermBC<T,U,U>(bcond));
    fbca->modify(u); 
    dwf_qdp.importGauge(u);
    struct timeval start,stop,diff ;

  /********************************************************
   * EigenCalculation: Calculate K eigenvectors
   ********************************************************
   */
   int K = atoi(argv[6])+10;
   int M = K + 10;
   int Tn = atoi(argv[7]);

   multi1d<multi1d<LatticeFermion> > Eigs(0); 
   multi1d<Complex> vals(0);
   
	QDPIO::cout << "Set up Lanczos 5d..." << endl;
  	Lanczos_5d<number> eig(dwf_qdp);		///Constructor
	///Eigensolver parameters	
	//eig.u = u;				///Read gauge field
	eig.M = M;				///Nsearch > Nwanted
	eig.K = K;				///Nwanted
	eig.conv = 1e-8;			///allowed residual
	eig.D = DDAGD;	 			///Type of multiplication;
	eig.Tn = Tn;				///Order of Chebyshev
	eig.ch_off = atof(argv[8]);		///> max wanted eigenvalue of G5D 
	eig.ch_rad = atof(argv[9]);		///> spectral radius of G5D 
	eig.sh = false; eig.ch_shift = 0.0;	///shift to get internal eigenvalues
	eig.prec = 0;				///preconditioned : 1 or not : 0
	eig.lock = false;
	eig.get = atoi(argv[6]);
	eig.maxits = 10000;
	QDPIO::cout << "Call eigensolver..." << endl;

  gettimeofday(&start,NULL);

	eig.Run();

  gettimeofday(&stop,NULL);
  timersub(&stop,&start,&diff);
  QDPIO::cout <<  "eig ran in " << diff.tv_sec << "." << diff.tv_usec << " " << M << " " << Tn << endl;



 
}








