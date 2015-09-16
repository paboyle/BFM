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

  multi1d<LatticeFermion> psi(Ls);
  multi1d<LatticeFermion> chi(Ls);
  multi1d<LatticeFermion> chi_deflated(Ls);

  for(int s=0;s<Ls;s++) gaussian(psi[s]);
  for(int s=0;s<Ls;s++) gaussian(chi[s]);
  for(int s=0;s<Ls;s++) gaussian(chi_deflated[s]);

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

dwf.init(dwfa);
dwf.importGauge(u);

   Fermion_t psi_h[2];
   psi_h[0] = dwf.allocFermion();
   psi_h[1] = dwf.allocFermion();
   Fermion_t chi_h[2];
   chi_h[0] = dwf.allocFermion();
   chi_h[1] = dwf.allocFermion();

#undef DEBUG

//Point Source
LatticePropagator Src = zero;

  multi1d<int> site(4);
  site[0] = 0;
  site[1] = 0;
  site[2] = 0;
  site[3] = 0;

  Real one = 1.0;
  Propagator F = one;
  pokeSite(Src,F,site);

      LatticeFermion ferm;
      PropToFerm(Src, ferm, 0,0);
  for(int s=0;s<Ls;s++) psi[s] = zero;
	psi[0]   = chiralProjectPlus(ferm);
	psi[Ls-1]= chiralProjectMinus(ferm);


/***************************/

 // Printf("Importing psi field cb %d\n",1);
  dwf.importFermion(psi,psi_h[0],0);
  dwf.importFermion(psi,psi_h[1],1);
  dwf.importFermion(chi,chi_h[0],0);
  dwf.importFermion(chi,chi_h[1],1); 

int Niter = 1;
for(int n=0; n<Niter; n++){

#pragma omp parallel 
  {
    omp_set_num_threads(threads);
#pragma omp for 
    for(int i=0;i<threads;i++) {
      //dwf.CGNE_M(chi_h,psi_h);
      dwf.MCR(chi_h,psi_h);
    }
  }

}

  	dwf.exportFermion(chi,chi_h[0],0);
  	dwf.exportFermion(chi,chi_h[1],1);


dwf.freeFermion(psi_h[0]);
dwf.freeFermion(psi_h[1]);
dwf.freeFermion(chi_h[0]);
dwf.freeFermion(chi_h[1]);

dwf.end();


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
	eig.M = M;				///Nsearch > Nwanted
	eig.K = K;				///Nwanted
	eig.conv = 1e-10;			///allowed residual
	eig.D = DDAGD;	 			///Type of multiplication;
	eig.Tn = Tn;				///Order of Chebyshev
	eig.ch_off = atof(argv[8]);		///> max wanted eigenvalue of G5D 
	eig.ch_rad = atof(argv[9]);		///> spectral radius of G5D 
	eig.sh = false; eig.ch_shift = 0.0;	///shift to get internal eigenvalues
	eig.prec = 1;				///preconditioned : 1 or not : 0
	eig.lock = false;
	eig.get = atoi(argv[6]);
	eig.maxits = 10000;
	QDPIO::cout << "Call eigensolver..." << endl;

  gettimeofday(&start,NULL);

	eig.Run();

  gettimeofday(&stop,NULL);
  timersub(&stop,&start,&diff);
  QDPIO::cout <<  "eig ran in " << diff.tv_sec << "." << diff.tv_usec << " " << M << " " << Tn << endl;



  /********************************************************
   * Deflate with calculated evecs
   ********************************************************
   */

   Fermion_t psi_hh[2];
   psi_hh[0] = dwf_qdp.allocFermion();
   psi_hh[1] = dwf_qdp.allocFermion();
   Fermion_t chi_hh[2];
   chi_hh[0] = dwf_qdp.allocFermion();
   chi_hh[1] = dwf_qdp.allocFermion();

/***************************/

 // Printf("Importing psi field cb %d\n",1);
  dwf_qdp.importFermion(psi,psi_hh[0],0);
  dwf_qdp.importFermion(psi,psi_hh[1],1);
  dwf_qdp.importFermion(chi_deflated,chi_hh[0],0);
  dwf_qdp.importFermion(chi_deflated,chi_hh[1],1); 

#pragma omp parallel 
  {
    omp_set_num_threads(threads);
#pragma omp for 
    for(int i=0;i<threads;i++) {
  	//dwf_qdp.CGNE(chi_hh,psi_hh,eig.bq,eig.bl);
  	dwf_qdp.MCR(chi_hh,psi_hh,eig.bq,eig.bl);
	}
}

  	dwf_qdp.exportFermion(chi_deflated,chi_hh[0],0);
  	dwf_qdp.exportFermion(chi_deflated,chi_hh[1],1);

	Real ndiff = 0.0;
	for(int s=0;s<Ls;s++) { ndiff = norm2(chi[s] - chi_deflated[s]); }
	QDPIO::cout << "|chi - chi_deflated| = " << sqrt(ndiff) << endl;
	Printf("Done\n");
 
dwf_qdp.freeFermion(psi_hh[0]);
dwf_qdp.freeFermion(psi_hh[1]);
dwf_qdp.freeFermion(chi_hh[0]);
dwf_qdp.freeFermion(chi_hh[1]);

dwf_qdp.end();


}








