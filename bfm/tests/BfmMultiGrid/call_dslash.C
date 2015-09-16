#include <chroma.h>
#include <actions/ferm/invert/syssolver_linop_cg_array.h>
#include <bfm.h>
#include <bfm_qdp.h>
#include <eigen/Krylov_5d.h>
#include <bfm_mcr.h>
#include "BfmMultiGrid.h"


#include <chroma.h>
#include <actions/ferm/invert/syssolver_linop_cg_array.h>
#include <actions/ferm/invert/syssolver_linop_aggregate.h>
#include <bfm_qdp_chroma_linop.h>

//define SOLVER MCR_prec
#define SOLVER CGNE_prec
#define UNPREC_SOLVER CGNE_M
using namespace Chroma;

typedef LatticeFermion T;
typedef multi1d<LatticeFermion> T5;
typedef multi1d<LatticeColorMatrix> U;

typedef double number;
typedef float  relax_number;

extern double LittleDopSolverResidual;
extern int    LittleDopSolver;
extern int    SubspaceSurfaceDepth;
extern int    InnerKrylovIterMax;
extern double InnerKrylovShift;


#define Printf if ( QMP_is_primary_node() ) printf
//#define Printf printf
void LuscherSubspace (int Ls,Fermion_t f, bfm_qdp<double> &dwf);
void ResidualSubspace (int Ls,Fermion_t f, bfm_qdp<double> &dwf);
void ResidualSequenceSubspace (int Ls,int nv,Fermion_t *f, bfm_qdp<double> &dwf);
void NestedInverse (
		    BfmMultiGrid &ldop,
		    Fermion_t src,
		    Fermion_t sol,
		    Fermion_t resid,
		    bfm_qdp<double> &dop);

int main (int argc,char **argv )
{
  Chroma::initialize(&argc,&argv);


  /********************************************************
   * Command line parsing
   ********************************************************
   */
  if ( argc < 11 ) { 
   Printf("Usage: %s lx ly lz lt Ls Nwanted Tn off rad\n All must be even\n",argv[0]);
   Printf("Good test: D_oo 4 4 4 4 4 10 30 1.0 6.0 BLCOK");
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
  int NumberSubspace=atoi(argv[6]);
  SubspaceSurfaceDepth=atoi(argv[7]);
  Printf("SurfaceDepth %d\n",SubspaceSurfaceDepth);
   int blocksize = atoi(argv[10]);
   int block[5];
   if(nrow[0]==48){
     block[0]=blocksize;
     block[1]=blocksize;
     block[2]=blocksize;
     block[3]=blocksize;
   } else { 
     block[0]=blocksize;
     block[1]=blocksize;
     block[2]=blocksize;
     block[3]=blocksize;
   }
   block[4]=Ls;

  Layout::setLattSize(nrow);
  Layout::create();

  int lx = QDP::Layout::subgridLattSize()[0];
  int ly = QDP::Layout::subgridLattSize()[1];
  int lz = QDP::Layout::subgridLattSize()[2];
  int lt = QDP::Layout::subgridLattSize()[3];

  int threads =   omp_get_max_threads();
  bfmarg::Threads(threads);
  bfmarg::Reproduce(0);
  bfmarg::ReproduceChecksum(0);
  bfmarg::ReproduceMasterCheck(0);
  bfmarg::Verbose(11);
  bfmarg::UseCGdiagonalMee(0);

  /********************************************************
   * Setup DWF operator
   ********************************************************
   */
  bfmarg dwfa;
  bfm_t  dwf;
  bfm_qdp<number> dop;
  dwfa.solver = DWF;
  dwfa.node_latt[0]  = lx;
  dwfa.node_latt[1]  = ly;
  dwfa.node_latt[2]  = lz;
  dwfa.node_latt[3]  = lt;
  dwfa.verbose=11;

  multi1d<int> procs = QDP::Layout::logicalSize();
  Printf("%d dim machine\n\t", procs.size());
  for(int mu=0;mu<4;mu++){
    if ( procs[mu]>1 ) {
      dwfa.local_comm[mu] = 0;
    } else { 
      dwfa.local_comm[mu] = 1;
    }
  }

  bfmarg::Threads(threads);

  multi1d<int> ncoor = QDP::Layout::nodeCoord();

  double M5=1.8;
  double mq=0.00078;
  dwfa.max_iter = 100000;
  dwfa.residual = 1.e-8;
  dwfa.ScaledShamirCayleyTanh(mq,M5,Ls,2.0);

  dop.init(dwfa);

  bfmarg dwfa_prec = dwfa;
  //  bfm_qdp<relax_number> dop_sp;
  //  dop_sp.init(dwfa_prec);

  /********************************************************
   * gauge field
   ********************************************************
   */

  multi1d<LatticeColorMatrix>  u(Nd);
  ArchivGauge_t Header ; readArchiv(Header,u,argv[11]);


  /********************************************************
   * Anti-Periodic Boundary conditions
   ********************************************************
   */

    multi1d<int> bcond (Nd);
    bcond[0] = bcond[1] = bcond[2] = 1; bcond[3] = -1;
    Handle<FermBC<T,U,U> > fbca (new SimpleFermBC<T,U,U>(bcond));
    fbca->modify(u); 
    dop.importGauge(u);
    //    dop_sp.importGauge(u);
    struct timeval start,stop,diff ;

    ///////////////////////////
    // Fermion fields gaussian  source
    ///////////////////////////

  multi1d<LatticeFermion> psi(Ls);
  multi1d<LatticeFermion> chi(Ls);

  for(int s=0;s<Ls;s++) psi[s]=zero;
  for(int s=0;s<Ls;s++) chi[s]=zero;

  Fermion_t psi_h[2];
  psi_h[0] = dop.allocFermion();
  psi_h[1] = dop.allocFermion();
  Fermion_t chi_h[2];
  chi_h[0] = dop.allocFermion();
  chi_h[1] = dop.allocFermion();
    
  LatticeFermion ferm;
  gaussian(ferm);
  ferm=where(Layout::latticeCoordinate(3) == 0, ferm,LatticeFermion(zero));

  for(int s=0;s<Ls;s++) psi[s] = zero;
  psi[0]   = chiralProjectPlus(ferm);
  psi[Ls-1]= chiralProjectMinus(ferm);
  
/***************************/

 // Printf("Importing psi field cb %d\n",1);
   dop.importFermion(psi,psi_h[0],0);
   dop.importFermion(psi,psi_h[1],1);
   dop.importFermion(chi,chi_h[0],0);
   dop.importFermion(chi,chi_h[1],1); 

   Fermion_t src = dop.allocFermion(); 
   Fermion_t tmp = dop.allocFermion(); 
   Fermion_t Mtmp= dop.allocFermion(); 
   Fermion_t resid= dop.allocFermion(); 

#pragma omp parallel 
  {
#pragma omp for 
    for(int i=0;i<threads;i++) {
      // Schur complement / decompose entry to get same source
      dop.MooeeInv(psi_h[Even],tmp,DaggerNo);
      dop.Meo     (tmp,src,Odd,DaggerNo);
      dop.axpy    (tmp,src,psi_h[Odd],-1.0);
      dop.Mprec(tmp,src,Mtmp,DaggerYes);  
      dop.fill(Mtmp,0.0);
    }
  }

  std::string file(argv[8]);
   
#if 0
   int Niter = 1;
   Printf("Calling undeflated unprec solver");
   for(int n=0; n<Niter; n++){

   std::string filecgne=file+".CGNE";
   dop.InverterLoggingBegin(filecgne);
#pragma omp parallel 
     {
#pragma omp for 
       for(int i=0;i<threads;i++) {
	 dop.SOLVER(Mtmp,src);
       }
     }
   }
   dop.InverterLoggingEnd();
   dop.InverterRegisterExactSolution(Mtmp);
#endif

/***************************/
 // Printf("Importing psi field cb %d\n",1);

  ///////////////////
  // Little dirac op
  ///////////////////
  QDPIO::cout << "Creating little Dirac op for "<< NumberSubspace<< " vectors"<<endl;

  int lq = atoi(argv[9]);
#if 1
  int quadrant[4] = { block[0],block[1],block[2],block[3] };
#else 
  int quadrant[4] = { lq,lq,block[2],block[3] };
#endif
  BfmMultiGrid ldop(Ls,NumberSubspace,block,quadrant,&dop);
  ldop.RelaxSubspace<number>(&dop);
  ldop.UsePcgPrecLinop(&dop);

  int ns = ldop.SubspaceLocalDimension();
  std::vector<std::complex<double>  > res(ns);
  std::vector<std::complex<double>  > Ares(ns);
  
  FlightLog("Computing little  operator using coloring");
#pragma omp parallel 
  {
#pragma omp for 
    for(int i=0;i<threads;i++) {
      ldop.ComputeLittleMatrixColored();
    }
  }

  int nv=128;
  QDPIO::cout << "Initialising second level deflation for Nvec="<< nv<<endl;
  ldop.LdopDeflationBasisInit(nv);

  FlightLog("SOLVING WITH POLYPREC");

  std::string filecg=file+".CG_Poly";
  std::string filemcr=file+".MCR";
  std::string filemcraug=file+".MCR_aug";

  LittleDopSolverResidual = 1.0e-11;
  InnerKrylovIterMax  = 8;
  InnerKrylovShift    = 1.0;
  filecg=file+".CG_PolyInnerBS6_8shift1.0";
  if(0){
  dop.InverterLoggingBegin(filecg);
#pragma omp parallel 
  {
#pragma omp for 
    for(int i=0;i<threads;i++) {
      dop.axpy(chi_h[Odd],src,src,0.0);
      ldop.CG_PolyPrec(chi_h[Odd],src,chi_h[Even]);
    }
  }
  dop.InverterLoggingEnd();
  }

  // 4^4 block
  // 8,1.0  @ 24 Vec => 566
  // 8,1.0  @ 36 Vec => 384
  // 8,3.0  @ 36 Vec => 568 <---800s @1k
  // 8,10.0 @ 36 Vec => 553
  // 16,10.0 @ 36 Vec => 550
  // 16,0.1  @ 36 Vec => ~ 400 but no converge
  // 
  // 6^4 
  // 8,1.0 Deep     @ 36 Vec => 478     [ est 670s @1k => ~2x over CGNE, better in single )
  //  stdout[0]: BfmMultiGrid::CG_PolyPrec MssInv 550.638561  s
  //  stdout[0]: BfmMultiGrid::CG_PolyPrec Proj   331.069631  s <--- 10% proj 25% inner => 48 vecs @1k
  //  stdout[0]: BfmMultiGrid::CG_PolyPrec other  2380.156044 s
  //  stdout[0]: 478
  //
  // 8,0.5 Deep     @ 36 Vec => 504
  //
  //
  // Compared to alpha SA -- i) Ac -> [1-alpha A] A [1-alpha A]  -- smooth edges of "chunks"
  //                        ii) Determine half my vectors using deflated linop iteration/residual.
  //                            -> alternative to RelaxSubspace
  // 
  // Try 6^4, Nv=128 [suppress MssInv O/H] Surface [suppress Proj] vecs=36
  // 8,1.0 36 Deep Relaxed subspace => 533
  // 6,3.0 36 Deep Relaxed subspace => 585
  //
  // Deep again, but source restricted to surfaces
  // Try 6^4, Nv=128 [suppress MssInv O/H] Surface [suppress Proj] vecs=36
  // 8,1.0 36 Deep Relaxed subspace => 486
  // 6,3.0 36 Deep Relaxed subspace => 535
  //
  //
  // 4,4,6,6 Deep 1.0: => 399
  // Updating subspace as we run is required => variable preconditioner


  InnerKrylovIterMax  = 8;
  InnerKrylovShift    = 1.0;
  ldop.PcgShift       = 1.0;
  ldop.PcgType        = PcgAdef2;

  filecg=file+".PcgADef2_e-4";
  dop.InverterLoggingBegin(filecg);
#pragma omp parallel 
  {
#pragma omp for 
    for(int i=0;i<threads;i++) {
      dop.axpy(chi_h[Odd],src,src,0.0);
      ldop.Pcg(chi_h[Odd],src,chi_h[Even],1.0e-12,1.0e-4);
    }
  }
  dop.InverterLoggingEnd();


  InnerKrylovIterMax  = 8;
  InnerKrylovShift    = 1.0;
  ldop.PcgShift       = 1.0;
  ldop.PcgType        = PcgAdef2;

  filecg=file+".PcgADef2_e-3";
  dop.InverterLoggingBegin(filecg);
#pragma omp parallel 
  {
#pragma omp for 
    for(int i=0;i<threads;i++) {
      dop.axpy(chi_h[Odd],src,src,0.0);
      ldop.Pcg(chi_h[Odd],src,chi_h[Even],1.0e-12,1.0e-3);
    }
  }
  dop.InverterLoggingEnd();


  InnerKrylovIterMax  = 8;
  InnerKrylovShift    = 1.0;
  ldop.PcgShift       = 1.0;
  ldop.PcgType        = PcgMssDef;

  filecg=file+".PcgMssDef";
  dop.InverterLoggingBegin(filecg);
#pragma omp parallel 
  {
#pragma omp for 
    for(int i=0;i<threads;i++) {
      dop.axpy(chi_h[Odd],src,src,0.0);
      ldop.Pcg(chi_h[Odd],src,chi_h[Even],1.0e-12,1.0e-4);
    }
  }
  dop.InverterLoggingEnd();



  dop.InverterLoggingBegin(filemcr);
#pragma omp parallel 
  {
#pragma omp for 
    for(int i=0;i<threads;i++) {
      dop.axpy(chi_h[Odd],src,src,0.0);
      ldop.MCR_PolyPrec(chi_h[Odd],src,resid);

    }
  }
  dop.InverterLoggingEnd();


  exit(0);

  InnerKrylovIterMax  = 8;
  InnerKrylovShift    =0.5;
  filecg=file+".CG_PolyInner8shift0.5";
  dop.InverterLoggingBegin(filecg);
#pragma omp parallel 
  {
#pragma omp for 
    for(int i=0;i<threads;i++) {
      dop.axpy(chi_h[Odd],src,src,0.0);
      ldop.CG_PolyPrec(chi_h[Odd],src,chi_h[Even]);
    }
  }
  dop.InverterLoggingEnd();


  InnerKrylovIterMax  = 16;
  InnerKrylovShift    =0.3;
  filecg=file+".CG_PolyInner16Shift0.3";
  dop.InverterLoggingBegin(filecg);
#pragma omp parallel 
  {
#pragma omp for 
    for(int i=0;i<threads;i++) {
      dop.axpy(chi_h[Odd],src,src,0.0);
      ldop.CG_PolyPrec(chi_h[Odd],src,chi_h[Even]);
    }
  }
  dop.InverterLoggingEnd();



  ldop.AddToSubspace(resid);

  QDPIO::cout << "Reinitialising second level deflation for Nvec="<< nv<<endl;
  ldop.LdopDeflationBasisInit(nv);
  
  dop.InverterLoggingBegin(filemcraug);
#pragma omp parallel 
  {
#pragma omp for 
    for(int i=0;i<threads;i++) {
      dop.axpy(chi_h[Odd],src,src,0.0);
      ldop.MCR_PolyPrec(chi_h[Odd],src,resid);
    }
  }
  dop.InverterLoggingEnd();


  FlightLog("DONE");
  exit(0);

}

//In MCR and CG the sequence of estimated residuals
//is a function of current residual, search direction and iteration number.
//
//Wish to return the current search direction and true residual from the solver
//and provide a routine to step this. Use this as a basis for the inexact deflation.
//

void NestedInverse (BfmMultiGrid &ldop,
		    Fermion_t src,
		    Fermion_t sol,
		    Fermion_t resid,
		    bfm_qdp<double> &dwf_qdp)
{

  Fermion_t tmp  = dwf_qdp.allocFermion();
  Fermion_t xc   = dwf_qdp.allocFermion();
  Fermion_t dxc   = dwf_qdp.allocFermion();
  Fermion_t Mtmp = dwf_qdp.allocFermion();

  double nsrc,nr,nrr,nrrs;
  int threads = dwf_qdp.threads;
#pragma omp parallel 
  {
#pragma omp for 
    for(int i=0;i<threads;i++) {
      nsrc =dwf_qdp.norm(src);
      dwf_qdp.set_zero(sol);
    }
  }
  
  int ns = ldop.SubspaceLocalDimension();
  std::vector<std::complex<double>  > res(ns);
  std::vector<std::complex<double>  > Ares(ns);
  QDPIO::cout << "Source norm " << sqrt(ns) << endl;
  
  LittleDopSolver=LittleDopSolverMCR;
  LittleDopSolverResidual=1.0e-5;
  int iters=0;
  while(1){

    dwf_qdp.max_iter = 150;
    // Compute true residual
#pragma omp parallel 
    {
#pragma omp for 
      for(int i=0;i<threads;i++) {
	dwf_qdp.Mprec(sol,tmp,Mtmp,DaggerNo);  
	dwf_qdp.Mprec(tmp,resid,Mtmp,DaggerYes);  
	dwf_qdp.axpy (resid,src,resid,-1.0);                  // r = M sol - src
	nr =dwf_qdp.norm(resid);
      }
    }
    QDPIO::cout << "True residual " << sqrt(nr/nsrc) << endl;
    if ( sqrt(nr/nsrc) < 1.0e-8 ) { 
      QDPIO::cout << "converged with residual " << sqrt(nr/nsrc) << " after "<<iters <<" iterations"<<endl;
      return;
    }
    // Evaluate Dss^-1 resid
    {
#pragma omp parallel 
    {
#pragma omp for 
      for(int i=0;i<threads;i++) {
	ldop.ProjectToSubspace(resid,res);
	ldop.ApplyInverse(res,Ares);
	ldop.PromoteFromSubspace(Ares,xc);        // xc=P Ainv R (M sol-src)  xcoarse

	// sol' =  sol - P Ainv R [Msol-src] 
	//      <=> M sol' - src = [M sol - src] - M Ainv [M sol - src]  ~ 0
	dwf_qdp.axpy (sol,xc,sol,-1.0);
      }
    }
    }

    double dxcn;
#pragma omp parallel 
    {
#pragma omp for 
      for(int i=0;i<threads;i++) {    

	// recompute resid
	dwf_qdp.Mprec(sol,tmp,Mtmp,DaggerNo);  
	dwf_qdp.Mprec(tmp,resid,Mtmp,DaggerYes);  
	dwf_qdp.axpy (resid,src,resid,-1.0);       //  M sol - src

	nrr =dwf_qdp.norm(resid);
	ldop.ProjectToSubspace(resid,res);
	nrrs=ldop.norm_vec(res);
      }
    }
    QDPIO::cout << "Deflation residual   " << sqrt(nrr/nsrc) <<endl;
    QDPIO::cout << "Deflation residual_s " << sqrt(nrrs/nsrc) <<endl;
     
    dwf_qdp.residual = 1.0e-8;  // Solve to present residual
#pragma omp parallel 
    {
#pragma omp for 
      for(int i=0;i<threads;i++) {    
	int iter = dwf_qdp.SOLVER(sol,src);
	iters+=iter;
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
      dwf.MCR_prec(sol,src,&vec,1);
    }
  }
  dwf.freeFermion(sol);
  dwf.freeFermion(src);
  dwf.residual= restore_residual;
}

void ResidualSequenceSubspace (int Ls,int nvec,Fermion_t *vec, bfm_qdp<double> &dwf)
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
      dwf.MCR_prec(sol,src,vec,nvec);
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
  dwf.residual = 1.0e-6;
#pragma omp parallel 
  {
    int threads = omp_get_num_threads();
#pragma omp for 
    for(int t=0;t<threads;t++) {
      for(int i=0;i<1;i++){
	dwf.SOLVER(sol,src);
	dwf.axpy(src,sol,sol,0.0);
      }
    }
  }
  dwf.freeFermion(src);
  dwf.residual= restore_residual;
}



