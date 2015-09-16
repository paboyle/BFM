#include <chroma.h>
#include <actions/ferm/invert/syssolver_linop_cg_array.h>
#include <bfm_qdp.h>
#include <eigen/Krylov_4d.h>

//typedef bfm bfm_t;
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

  if ( argc < 10 ) { 
   Printf("Usage: %s lx ly lz lt Ls Nwanted Tn off rad\n All must be even\n",argv[0]);
   Printf("Good test: D_ov 4 4 4 4 4 10 0 0 0");
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

  /********************************************************
   * Setup DWF operator
   ********************************************************
   */
   bfmarg dwfa;
  bfm_qdp<double> dwf;

  double M5 = 1.8;
  double mq = 0.0;

  //Physics parameters
  dwfa.solver       = HtCayleyTanh;
  dwfa.Ls           = Ls;
  dwfa.mass         = mq;
  dwfa.M5           = M5;
  dwfa.precon_5d    = 0;
  dwfa.max_iter     = 10000 ;
  dwfa.residual     = 1e-4;
  //Geometry
  dwfa.node_latt[0] = QDP::Layout::subgridLattSize()[0];
  dwfa.node_latt[1] = QDP::Layout::subgridLattSize()[1];
  dwfa.node_latt[2] = QDP::Layout::subgridLattSize()[2];
  dwfa.node_latt[3] = QDP::Layout::subgridLattSize()[3];

  multi1d<int> procs = QDP::Layout::logicalSize();
  QDPIO::cout << procs.size() << " dim machine\n\t" << endl;
  for(int mu=0;mu<4;mu++){
    QDPIO::cout << procs[mu] << " ";
    if ( procs[mu]>1 ) {
      dwfa.local_comm[mu] = 0;
    } else { 
      dwfa.local_comm[mu] = 1;
    }
  }
  QDPIO::cout << "\nLocal comm = ";
  for(int mu=0;mu<4;mu++){
    QDPIO::cout << dwfa.local_comm[mu] << " ";
  }
  QDPIO::cout << endl; 
  
  multi1d<int> ncoor = QDP::Layout::nodeCoord();

  dwf.init(dwfa);


  
 Seed seed;
 QDP::RNG::savern(seed);
 QDPIO::cout<<"HADRON: seed declare "<< seed <<endl;
 QDP::RNG::setrn(seed);

  /********************************************************
   * Gaussian gauge field
   ********************************************************
   */

  multi1d<LatticeColorMatrix>  u(Nd);
  for(int i=0;i<Nd; i++){
   u[i] = 1.0; 
  }
  HotSt(u);
  if(argc == 11) {ArchivGauge_t Header ; readArchiv(Header,u,argv[10]);}

  /********************************************************
   * Anti-Periodic Boundary conditions
   ********************************************************
   */

    multi1d<int> bcond (Nd);
    bcond[0] = bcond[1] = bcond[2] = 1; bcond[3] = -1;
    Handle<FermBC<T,U,U> > fbc (new SimpleFermBC<T,U,U>(bcond));
    fbc->modify(u);
    dwf.importGauge(u);

  /********************************************************
   * EigenCalculation: Calculate K eigenvectors
   ********************************************************
   */
   int K = atoi(argv[6]);
   int M = K + 20;
   int Tn = atoi(argv[7]);

   multi1d<LatticeFermion> Eigs(0); 
   multi1d<Complex> vals(0);
   
	QDPIO::cout << "Set up Lanczos 4d..." << endl;
  	Lanczos_4d<double> eig(dwf);
	///Eigensolver parameters
	eig.u = u;				
	eig.M = M;
	eig.K = K;
	eig.conv = 1e-4;
	eig.D = G5D;
	eig.Tn = Tn;
	eig.ch_off = atof(argv[8]);
	eig.ch_rad = atof(argv[9]);
	eig.sh = false; eig.ch_shift = 0.0;
	QDPIO::cout << "Call eigensolver..." << endl;
	eig.Run();
	//eig.eigs_to_qdp(Eigs,vals);


  /********************************************************
   * Calculate 10 extra just for fun
   ********************************************************
   */
  /* K += 10;
   M += 10;

   
	QDPIO::cout << "Set up Lanczos 4d again..." << endl;
  	Lanczos_4d<double> eig2(dwf);
	///Eigensolver parameters
        eig2.u = u;		
	eig2.M = M;
	eig2.K = K;
	eig2.conv = 1e-4;
	eig2.D = G5D;
	eig2.Tn = Tn;
	eig2.ch_off = atof(argv[8]);
	eig2.ch_rad = atof(argv[9]);
	eig2.sh = false; eig.ch_shift = 0.0;
	///Set continue parameter = true
	eig2.cont = true;
	QDPIO::cout << "Call eigensolver..." << endl;
	eig2.Initial_State(Eigs, vals);
	eig2.Run();

	Printf("Done\n");*/

}








