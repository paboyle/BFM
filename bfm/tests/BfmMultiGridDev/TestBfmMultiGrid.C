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
#include <spi/include/kernel/spec.h>

#define SOLVER CGNE_prec
#define UNPREC_SOLVER CGNE_M
using namespace Chroma;

typedef LatticeFermion T;
typedef multi1d<LatticeFermion> T5;
typedef multi1d<LatticeColorMatrix> U;

typedef double number;
typedef float relax_number;


#define Printf if ( QMP_is_primary_node() ) printf

int main (int argc,char **argv )
{

  Chroma::initialize(&argc,&argv);

  srand48(QMP_get_node_number()*7-1);

  /********************************************************
   * Command line parsing
   ********************************************************
   */
  if ( argc < 11 ) { 
   Printf("Usage: %s Config Outname lx ly lz lt Ls Nvec Depth Bx By Bz Bt \n All must be even\n",argv[0]);
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
  char * confname = argv[1];
  char * outname  = argv[2];

  multi1d<int> nrow(Nd);
  nrow[0] = atoi(argv[3]);
  nrow[1] = atoi(argv[4]);
  nrow[2] = atoi(argv[5]);
  nrow[3] = atoi(argv[6]);
  int Ls  = atoi(argv[7]);


  int NumberSubspace=atoi(argv[8]);
  double SubspaceSurfaceDepth=atoi(argv[9]);
  int block[5];
  block[0]=atoi(argv[10]);
  block[1]=atoi(argv[11]);
  block[2]=atoi(argv[12]);
  block[3]=atoi(argv[13]);
  block[4]=Ls;

  QDPIO::cout << "Vol "
	      << nrow[0] << " "
	      << nrow[1] << " "
	      << nrow[2] << " "
	      << nrow[3] << " "
	      << Ls   
	      <<endl;

  QDPIO::cout << "Block "
	      << block[0] << " "
	      << block[1] << " "
	      << block[2] << " "
	      << block[3] << " "
	      << block[4]
	      <<endl;

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
  bfm_qdp<relax_number> dop_sp;
  dop_sp.init(dwfa_prec);



  /********************************************************
   * gauge field
   ********************************************************
   */

  multi1d<LatticeColorMatrix>  u(Nd);
  if ( !strcmp(confname,"HotSt") ) { 
    HotSt(u);
  } else { 
    ArchivGauge_t Header ; 
    readArchiv(Header,u,confname);
  }
  /********************************************************
   * Anti-Periodic Boundary conditions
   ********************************************************
   */

  multi1d<int> bcond (Nd);
  bcond[0] = bcond[1] = bcond[2] = 1; bcond[3] = -1;
  Handle<FermBC<T,U,U> > fbca (new SimpleFermBC<T,U,U>(bcond));
  fbca->modify(u); 
  dop.importGauge(u);
  dop_sp.importGauge(u);



    ///////////////////////////
    // Fermion fields gaussian  source
    ///////////////////////////

  struct timeval start,stop,diff ;
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
  dop.comm_init();
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

  std::string file(outname);
  std::string filecg=file+".CGNE";
#if 0
   Printf("Calling undeflated unprec solver");
   dop.InverterLoggingBegin(filecg);
#pragma omp parallel 
     {
#pragma omp for 
       for(int i=0;i<threads;i++) {
	 dop.SOLVER(Mtmp,src);
       }
     }
   dop.InverterLoggingEnd();
   dop.InverterRegisterExactSolution(Mtmp);
#endif


  ///////////////////
  // Little dirac op
  ///////////////////

  // Setup Parms.
  BfmMultiGridParams bmgp;
  bmgp.NumberSubspace=NumberSubspace;

  for(int mu=0;mu<4;mu++) bmgp.Block[mu] = block[mu];
  bmgp.Block[4]=Ls;

  bmgp.Ls=Ls;
  bmgp.SubspaceRationalLo      = 0.001;
  bmgp.SubspaceRationalResidual= 3.0e-7;
  bmgp.SubspaceSurfaceDepth    = SubspaceSurfaceDepth ;

  bmgp.PreconditionerKrylovResidual     = 1.0e-6; // Run to non-converge for fixed depth
  bmgp.PreconditionerKrylovIterMax      = 8;
  bmgp.PreconditionerKrylovShift        = 1.0;

  bmgp.LittleDopSolverResidualInner = 1.0e-2;
  bmgp.LittleDopSolverResidualVstart= 1.0e-6;
  bmgp.LittleDopSolverResidualSubspace=1.0e-7;

  bmgp.PcgType                 = PcgAdef2f;


  QDPIO::cout << "Creating little Dirac op for "<< NumberSubspace<< " vectors"<<endl;
  //  int quadrant[4] = { block[0],block[1],block[2],block[3] };
  BfmMultiGrid<double> ldop(Ls,bmgp.NumberSubspace,bmgp.Block,bmgp.Block,&dop,&dop_sp);
  ldop.SetParams(bmgp);

  // Complicated init sequence
  ldop.RelaxSubspace<number>(&dop);
  ldop.ComputeLittleMatrixColored();

  // Deprecated params.

  ldop.LdopDeflationBasisInit(8);
  ldop.LdopDeflationBasisInit(16);
  ldop.LdopDeflationBasisInit(24);
  ldop.LdopDeflationBasisInit(32);
  ldop.LdopDeflationBasisInit(48);
  ldop.LdopDeflationBasisInit(64);
  ldop.LdopDeflationBasisInit(128);
  ldop.LdopDeflationBasisDiagonalise(128);
  ldop.SinglePrecSubspace();


  BfmMultiGrid<float>  ldop_f(Ls,bmgp.NumberSubspace,bmgp.Block,bmgp.Block,&dop,&dop_sp);
  ldop_f.CloneSubspace<double>(ldop);
  bmgp.LittleDopSolverResidualVstart=1.0e-5;
  bmgp.LittleDopSolverResidualInner =5.0e-3;
  ldop_f.SetParams(bmgp);

#pragma omp parallel 
  {
#pragma omp for 
    for(int i=0;i<threads;i++) {
      dop.axpy(chi_h[Odd],src,src,0.0);
      ldop_f.Pcg(chi_h[Odd],src,chi_h[Even]);
    }
  }


  ldop_f.SetParams(bmgp);
#pragma omp parallel 
  {
#pragma omp for 
    for(int i=0;i<threads;i++) {
      dop.axpy(chi_h[Odd],src,src,0.0);
      dop.residual=3.0e-5;
      ldop_f.Pcg(chi_h[Odd],src,chi_h[Even]);
      dop.residual=3.0e-7;
      ldop_f.Pcg(chi_h[Odd],src,chi_h[Even]);
      dop.residual=1.0e-8;
      ldop_f.Pcg(chi_h[Odd],src,chi_h[Even]);
    }
  }


  bmgp.LittleDopSolverResidualVstart=1.0e-5;
  bmgp.LittleDopSolverResidualInner =1.0e-4;
  bmgp.PreconditionerKrylovIterMax      = 35;
  bmgp.PreconditionerKrylovShift        = 0.03;
  ldop_f.SetParams(bmgp);
#pragma omp parallel 
  {
#pragma omp for 
    for(int i=0;i<threads;i++) {
      dop.axpy(chi_h[Odd],src,src,0.0);
      ldop_f.Pcg(chi_h[Odd],src,chi_h[Even]);
    }
  }

  bmgp.LittleDopSolverResidualVstart=1.0e-4;
  bmgp.LittleDopSolverResidualInner =1.0e-4;
  bmgp.PreconditionerKrylovIterMax      = 50;
  bmgp.PreconditionerKrylovShift        = 0.03;
  ldop_f.SetParams(bmgp);
#pragma omp parallel 
  {
#pragma omp for 
    for(int i=0;i<threads;i++) {
      dop.axpy(chi_h[Odd],src,src,0.0);
      ldop_f.Pcg(chi_h[Odd],src,chi_h[Even]);
    }
  }


  exit(0);

}


