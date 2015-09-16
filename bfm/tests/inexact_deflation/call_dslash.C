#include <chroma.h>
#include <actions/ferm/invert/syssolver_linop_cg_array.h>
#include <bfm.h>
#include <bfm_qdp.h>
#include <eigen/Krylov_5d.h>
#include "littleDiracOp.h"
#include <bfm_mcr.h>

#include <bfm_chebyshev.h>

#include <chroma.h>
#include <actions/ferm/invert/syssolver_linop_cg_array.h>
#include <actions/ferm/invert/syssolver_linop_aggregate.h>
#include <bfm_qdp_chroma_linop.h>

typedef bfm bfm_t;
//typedef commQMPbagel1<double> bfm_t;
#define SOLVER MCR_PolyPrec
#define UNPREC_SOLVER MCR
using namespace Chroma;

typedef LatticeFermion T;
typedef multi1d<LatticeFermion> T5;
typedef multi1d<LatticeColorMatrix> U;

typedef double number;

#define Printf if ( QMP_is_primary_node() ) printf
//#define Printf printf
void LuscherSubspace (int Ls,Fermion_t f, bfm_qdp<double> &dwf);
void ResidualSubspace (int Ls,Fermion_t f, bfm_qdp<double> &dwf);
void NestedInverse (
		    BfmLittleDiracOperator &ldop,
		    Fermion_t src,
		    Fermion_t sol,
		    Fermion_t resid,
		    bfm_qdp<double> &dwf_qdp);

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
  dwfa.solver = DWF;
  dwfa.node_latt[0]  = lx;
  dwfa.node_latt[1]  = ly;
  dwfa.node_latt[2]  = lz;
  dwfa.node_latt[3]  = lt;
  dwfa.verbose=0;

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
  //  dwfa.solver    = DWFrb4d; //DWFrb4d = 4d preconditioning. DWF = 5d preconditioning
  //  dwfa.precon_5d = 0;
  dwfa.Ls   = Ls;
  dwfa.M5   = toDouble(M5);
  dwfa.mass = toDouble(mq);
  dwfa.Csw  = 0.0;
  dwfa.max_iter = 100000;

  dwfa.residual = 1.e-8;
  Printf("Initialising bfm operator\n"); 

  bfm_qdp<number> dwf_qdp;
  dwf_qdp.init(dwfa);

  /********************************************************
   * Gaussian gauge field
   ********************************************************
   */

  multi1d<LatticeColorMatrix>  u(Nd);
  ArchivGauge_t Header ; readArchiv(Header,u,argv[10]);

  multi1d<LatticeFermion> psi(Ls);
  multi1d<LatticeFermion> chi(Ls);
  multi1d<LatticeFermion> chi_deflated(Ls);

  for(int s=0;s<Ls;s++) gaussian(psi[s]);
  for(int s=0;s<Ls;s++) chi[s]=zero;
  for(int s=0;s<Ls;s++) chi_deflated[s]=zero;

  /********************************************************
   * Anti-Periodic Boundary conditions
   ********************************************************
   */

    multi1d<int> bcond (Nd);
    bcond[0] = bcond[1] = bcond[2] = 1; bcond[3] = 1;
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

#ifdef POINT_SOURCE
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
#else
  LatticeFermion ferm;
  gaussian(ferm);
#endif

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
  Printf("Calling undeflated unprec solver");
for(int n=0; n<Niter; n++){

#pragma omp parallel 
  {
    omp_set_num_threads(threads);
#pragma omp for 
    for(int i=0;i<threads;i++) {
      dwf.UNPREC_SOLVER(chi_h,psi_h);
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
   int NumberSubspace=atoi(argv[6]);
   int K = NumberSubspace+10;
   int M = K + 10;
   int Tn = atoi(argv[7]);

   multi1d<Complex> vals(0);
   
	QDPIO::cout << "Set up Lanczos 5d..." << endl;
  	Lanczos_5d<number> eig(dwf_qdp);		///Constructor
	///Eigensolver parameters	
	eig.M = M;				///Nsearch > Nwanted
	eig.K = K;				///Nwanted
	eig.conv = 1e-8;			///allowed residual
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

	multi1d< multi1d<LatticeFermion> > subspace(NumberSubspace);
	multi1d<Complex> evals;
	eig.eigs_to_qdp(subspace,evals);
	
  gettimeofday(&stop,NULL);
  timersub(&stop,&start,&diff);
  QDPIO::cout <<  "eig ran in " << diff.tv_sec << "." << diff.tv_usec << " " << M << " " << Tn << endl;
  for(int i=0;i<subspace.size();i++){
    QDPIO::cout << "Evec "<<i<<" norm " << norm2(subspace[i])<< " lambda " << evals[i]<<endl;
  }


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

  Printf("Calling deflated unprec solver \n");
#pragma omp parallel 
  {
    omp_set_num_threads(threads);
#pragma omp for 
    for(int i=0;i<threads;i++) {
  	dwf_qdp.UNPREC_SOLVER(chi_hh,psi_hh,eig.bq,eig.bl);
    }
  }

  dwf_qdp.exportFermion(chi_deflated,chi_hh[0],0);
  dwf_qdp.exportFermion(chi_deflated,chi_hh[1],1);

  Real ndiff = 0.0;
  for(int s=0;s<Ls;s++) { ndiff = norm2(chi[s] - chi_deflated[s]); }
  QDPIO::cout << "|chi - chi_deflated| = " << sqrt(ndiff) << endl;
  Printf("Done\n");


// Try a little dirac operator
   Handle<CreateFermState<T4,U,U> > cfs( new CreateSimpleFermState<T4,U,U>(fbca));
   EvenOddPrecDWFermActArray  S_f(cfs, M5, mq, Ls);
   Handle< FermState<T4,U,U> > fs( S_f.createState(u) );
   Handle< LinearOperatorArray<T4> > linop(S_f.precLinOp(fs,mq));
   int block[5];
   block[0]=4;
   block[1]=4;
   block[2]=4;
   block[3]=4;
   block[4]=4;

   Fermion_t src = dwf_qdp.allocFermion(); 
   Fermion_t tmp = dwf_qdp.allocFermion(); 
   Fermion_t Mtmp= dwf_qdp.allocFermion(); 
   Fermion_t resid= dwf_qdp.allocFermion(); 

#pragma omp parallel 
  {
    omp_set_num_threads(threads);
#pragma omp for 
    for(int i=0;i<threads;i++) {
      // Schur complement / decompose entry to get same source
      dwf_qdp.MooeeInv(psi_hh[Even],tmp,DaggerNo);
      dwf_qdp.Meo     (tmp,src,Odd,DaggerNo);
      dwf_qdp.axpy    (tmp,src,psi_hh[Odd],-1.0);
      dwf_qdp.Mprec(tmp,src,Mtmp,DaggerYes);  
    }
  }

  QDPIO::cout << "Creating little Dirac op for "<< NumberSubspace<< " vectors"<<endl;

  int NumberSubspaceTrunc=NumberSubspace;
  Fermion_t subspaceTrunc[NumberSubspaceTrunc];

  enum { SubspaceEigen, SubspaceLuscher, SubspaceResiduals};
  int subspaceType = SubspaceLuscher;

  for(int i=0;i<NumberSubspaceTrunc;i++)    subspaceTrunc[i] = dwf_qdp.allocFermion();

  if( subspaceType == SubspaceEigen ) {
    QDPIO::cout << "Eigen subspace"<<endl;
    for(int i=0;i<NumberSubspaceTrunc;i++)      
      dwf_qdp.importFermion(subspace[i],subspaceTrunc[i],1);
  } else if ( subspaceType==SubspaceLuscher ) { 
    QDPIO::cout << "Luscher subspace"<<endl;
    for(int i=0;i<NumberSubspaceTrunc;i++)      
      LuscherSubspace(Ls,subspaceTrunc[i],dwf_qdp);
  } else if ( subspaceType==SubspaceResiduals ) { 
    QDPIO::cout << "Residual subspace"<<endl;
    for(int i=0;i<NumberSubspaceTrunc;i++)      
      ResidualSubspace(Ls,subspaceTrunc[i],dwf_qdp);
  }

  BfmLittleDiracOperator ldop(Ls,NumberSubspaceTrunc,block,subspaceTrunc,&dwf_qdp);
  int ns = ldop.SubspaceLocalDimension();
  QDPIO::cout << "local  subspace dimension is "<< ns<<endl;
  QDPIO::cout << "global subspace dimension is "<< ldop.SubspaceDimension()<<endl;
  QMP_barrier();

  std::vector<std::complex<double>  > res(ns);
  std::vector<std::complex<double>  > Ares(ns);

  QDPIO::cout << "Completeness of eigenbasis test"<<endl;
  for(int e=0;e<subspace.size();e++){
    Fermion_t tmpO = dwf_qdp.allocFermion();
    Fermion_t tmpI = dwf_qdp.allocFermion();
    multi1d<LatticeFermion> tmpF(Ls);
    for(int s=0;s<Ls;s++) tmpF[s]=zero;
    dwf_qdp.importFermion(subspace[e],tmpI,1);
    ldop.ProjectToSubspace(tmpI,res);
    ldop.PromoteFromSubspace(res,tmpO);
    dwf_qdp.exportFermion(tmpF,tmpO,1);
    QDPIO::cout << e<<" completeness "<< norm2(tmpF-subspace[e]) << " " <<norm2(tmpF)<< " " << norm2(subspace[e]) << endl;
    dwf_qdp.freeFermion(tmpO);
    dwf_qdp.freeFermion(tmpI);
  }

    FlightLog("Computing little operator");
    //    ldop.ComputeLittleMatrix();
  //  FlightLog("Computing little using coloring operator");
    ldop.ComputeLittleMatrixColored();
    FlightLog("Projecting to subspace");
    ldop.ProjectToSubspace(src,res);
    FlightLog("Applying inverse");
    ldop.ApplyInverse(res,Ares);
    FlightLog("Promoting");
    ldop.PromoteFromSubspace(Ares,chi_hh[Odd]);
    FlightLog("Done");
 
#pragma omp parallel 
  {
#pragma omp for 
    for(int i=0;i<threads;i++) {
      dwf_qdp.SOLVER(chi_hh[Odd],src);
      dwf_qdp.Mprec(chi_hh[Odd],tmp,Mtmp,DaggerNo);  
      dwf_qdp.Mprec(tmp,resid,Mtmp,DaggerYes);  
      dwf_qdp.axpy (resid,resid,src,-1.0);
    }
}

  QDPIO::cout << "Completeness of residual test"<<endl;
  FILE * fp;
  if ( dwf_qdp.isBoss() ) fp=fopen("resid.dat","w");
  for(double rr=1.0e-4;rr>0.9e-12;rr=rr/10.0){

    dwf_qdp.residual=rr;
#pragma omp parallel 
  {
#pragma omp for 
    for(int i=0;i<threads;i++) {
      //zero the solution
      if ( 1 ) dwf_qdp.axpby(chi_hh[Odd],chi_hh[Odd],chi_hh[Odd],0.0,0.0);
      dwf_qdp.SOLVER(chi_hh[Odd],src);
      dwf_qdp.Mprec(chi_hh[Odd],tmp,Mtmp,DaggerNo);  
      dwf_qdp.Mprec(tmp,resid,Mtmp,DaggerYes);  
      dwf_qdp.axpy (resid,resid,src,-1.0);

    }
  }

  
  {
    ldop.ProjectToSubspace(resid,res);
    ldop.PromoteFromSubspace(res,chi_hh[0]);
    double n2diff;
    double n2r;
#pragma omp parallel 
    {
#pragma omp for 
      for(int i=0;i<threads;i++) {
	dwf_qdp.axpy(chi_hh[0],chi_hh[0],resid,-1.0);
	n2diff = dwf_qdp.norm(chi_hh[0]);
	n2r    = dwf_qdp.norm(resid);
      }
    }
    if ( dwf_qdp.isBoss() ) fprintf(fp,"%le %le %le %le %le\n",rr,n2diff,n2r,n2diff/n2r,1.0-n2diff/n2r);
    QDPIO::cout << "diff/nrm = "<<n2diff<<"/"<<n2r <<endl;
  }
  }
  if ( dwf_qdp.isBoss() ) fclose(fp);
  QDPIO::cout << " Nested inverse test" <<endl;
#pragma omp parallel 
    {
#pragma omp for 
      for(int i=0;i<threads;i++) {
	dwf_qdp.set_zero(chi_hh[Odd]);
      }
    }
    NestedInverse (ldop,
		   src,
		   chi_hh[Odd],resid,dwf_qdp);

  QDPIO::cout << "Undeflated inversion to 10^-12\n ";

  dwf_qdp.residual = 1.0e-12;
#pragma omp parallel 
    {
#pragma omp for 
      for(int i=0;i<threads;i++) {
	dwf_qdp.set_zero(chi_hh[Odd]);
	int iter = dwf_qdp.SOLVER(chi_hh[Odd],src);
      }
    }

  dwf_qdp.freeFermion(psi_hh[0]);
  dwf_qdp.freeFermion(psi_hh[1]);
  dwf_qdp.freeFermion(chi_hh[0]);
  dwf_qdp.freeFermion(chi_hh[1]);
  dwf_qdp.end();

}

//In MCR and CG the sequence of estimated residuals
//is a function of current residual, search direction and iteration number.
//
//Wish to return the current search direction and true residual from the solver
//and provide a routine to step this. Use this as a basis for the inexact deflation.
//

void NestedInverse (BfmLittleDiracOperator &ldop,
		    Fermion_t src,
		    Fermion_t sol,
		    Fermion_t resid,
		    bfm_qdp<double> &dwf_qdp)
{

  Fermion_t tmp  = dwf_qdp.allocFermion();
  Fermion_t xc   = dwf_qdp.allocFermion();
  Fermion_t dxc   = dwf_qdp.allocFermion();
  Fermion_t Mtmp = dwf_qdp.allocFermion();

  double ns,nr,nrr;
  int threads = dwf_qdp.threads;
#pragma omp parallel 
    {
#pragma omp for 
      for(int i=0;i<threads;i++) {
	ns =dwf_qdp.norm(src);
	dwf_qdp.set_zero(sol);
      }
    }

    QDPIO::cout << "Source norm " << sqrt(ns) << endl;
  for(double r=1.0e-6;r>0.9e-12;r/=1000) {


    // Compute true residual
#pragma omp parallel 
    {
#pragma omp for 
      for(int i=0;i<threads;i++) {
	dwf_qdp.Mprec(sol,tmp,Mtmp,DaggerNo);  
	dwf_qdp.Mprec(tmp,resid,Mtmp,DaggerYes);  
	dwf_qdp.axpy (resid,src,resid,-1.0); // r=M sol-src
	nr =dwf_qdp.norm(resid);
      }
    }
    QDPIO::cout << "True residual " << sqrt(nr/ns) << endl;
    
    // Evaluate Dss^-1 resid
    {
      int ns = ldop.SubspaceLocalDimension();
      std::vector<std::complex<double>  > res(ns);
      std::vector<std::complex<double>  > Ares(ns);
      ldop.ProjectToSubspace(resid,res);
      ldop.ApplyInverse(res,Ares);
      ldop.PromoteFromSubspace(Ares,xc);  // xc=P Ainv R (M sol-src)  xcoarse
    }

    double dxcn;
    std::complex<double> alpha;
#pragma omp parallel 
    {
#pragma omp for 
      for(int i=0;i<threads;i++) {    

// sol' = sol - Ainv [Msol-src] 
//    => Msol' - src = [Msol-src] - M Ainv[Msol - src] 
	dwf_qdp.axpy (sol,xc,sol,-1.0);
	// recompute resid
	dwf_qdp.Mprec(sol,tmp,Mtmp,DaggerNo);  
	dwf_qdp.Mprec(tmp,resid,Mtmp,DaggerYes);  
	dwf_qdp.axpy (resid,resid,src,-1.0);
	nrr =dwf_qdp.norm(resid);
      }
    }
    QDPIO::cout << "Deflation residual " << sqrt(nrr/ns) <<endl;
     
    dwf_qdp.residual = r;  // Solve to present residual
#pragma omp parallel 
    {
#pragma omp for 
      for(int i=0;i<threads;i++) {    
	int iter = dwf_qdp.SOLVER(sol,src);
      }
    }

  }


  dwf_qdp.freeFermion(dxc);
  dwf_qdp.freeFermion(xc);
  dwf_qdp.freeFermion(tmp);
  dwf_qdp.freeFermion(Mtmp);
}
void ResidualSubspace (int Ls,Fermion_t vec, bfm_qdp<double> &dwf)
{
  multi1d<LatticeFermion> gauss(Ls);
  for(int s=0;s<Ls;s++) gaussian(gauss[s]);
  Fermion_t src = dwf.allocFermion();
  Fermion_t sol = dwf.allocFermion();
  dwf.importFermion(gauss,src,1);
  double restore_residual = dwf.residual;
  dwf.residual = 1.0e-6;
#pragma omp parallel 
  {
    int threads = omp_get_num_threads();
#pragma omp for 
    for(int i=0;i<threads;i++) {
      dwf.MCR_prec(vec,src,&sol,1);
    }
  }
  dwf.freeFermion(sol);
  dwf.freeFermion(src);
  dwf.residual= restore_residual;
}

void LuscherSubspace (int Ls,Fermion_t sol, bfm_qdp<double> &dwf)
{
  multi1d<LatticeFermion> gauss(Ls);
  for(int s=0;s<Ls;s++){
    gaussian(gauss[s]);
  }
  Fermion_t src = dwf.allocFermion();
  dwf.importFermion(gauss,src,1);
  double restore_residual = dwf.residual;
  dwf.residual = 1.0e-5;
#pragma omp parallel 
  {
    int threads = omp_get_num_threads();
#pragma omp for 
    for(int t=0;t<threads;t++) {
      for(int i=0;i<5;i++){
	dwf.SOLVER(sol,src);
	dwf.axpy(src,sol,sol,0.0);
      }
    }
  }
  dwf.freeFermion(src);
  dwf.residual= restore_residual;
}



