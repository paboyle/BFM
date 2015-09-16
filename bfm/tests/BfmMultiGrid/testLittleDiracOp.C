#include <qdp.h>
#include <bfm.h>
#include "BfmMultiGrid.h"

#include <bfm_qdp_chroma_linop.h>

#include <bfm_qdp.h>
#include <bfm_qdp_g5d.h>
#include <bfm_qdp_dwf.h>
#include <bfm_qdp_chroma_linop.h>
#include <bfm_wrapper.h>
#include <eigen/Krylov_5d.h>

#define Printf if ( QMP_is_primary_node() ) printf

#include <bfmcommqmp.h>
#include <bfmcommspi.h>
typedef bfmcommQMP<double> bfm_qmp;

void test_solver(BfmSolver solver);

int main (int argc,char **argv )
{
  Chroma::initialize(&argc,&argv);

  WilsonTypeFermActs4DEnv::registerAll(); 
  LinOpSysSolverArrayEnv::registerAll();

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


  bfmarg::Threads(16);
  Layout::setLattSize(nrow);
  Layout::create();



  test_solver(DWF);

}

void test_solver(BfmSolver solver)
{

  g5dParams parms;

  int Ls=16;
  double M5=1.8;
  double mq=0.0001;
  double wilson_lo = 0.05;
  double wilson_hi = 6.8;
  double shamir_lo = 0.025;
  double shamir_hi = 1.7;
  double ht_scale=1.7;
  double hw_scale=1.0;

  if ( solver != DWF ) { 
    exit(0);
    Printf("Should be testing HtCayleyTanh aka DWF\n");
  }
  parms.pDWF(mq,M5,Ls);

  multi1d<LatticeColorMatrix> u(4);
  HotSt(u);
  //  ArchivGauge_t Header ; readArchiv(Header,u,"ckpoint_lat.3000");  

  multi1d<LatticeFermion> src(Ls);

/* Rudy calculate some eigenvectors */


  BfmWrapperParams BWP;
  BWP.BfmInverter = BfmInv_CG; 
  BWP.BfmMatrix   = BfmMat_M;
  BWP.BfmPrecision= Bfm64bit;
  BWP.MaxIter     = 10000;
  BWP.RsdTarget.resize(1);
  BWP.RsdTarget[0]= 1.0e-9;
  BWP.Delta = 1.0e-4;
  BWP.BAP = parms;
  BfmWrapper bfm(BWP);

    bfmarg bfma;
#if defined(QDP_USE_OMP_THREADS)
    bfma.Threads(omp_get_max_threads());
#else
    bfma.Threads(16);
#endif
    bfma.Verbose(0);

    //Physics parameters
    bfmActionParams *bfmap = (bfmActionParams *) &bfma;
    *bfmap = bfm.invParam.BAP;
    
    // Algorithm & code control
    bfma.time_report_iter=-100;
    bfma.max_iter     = bfm.invParam.MaxIter;
    bfma.residual     = toDouble(bfm.invParam.RsdTarget[0]);

  int lx = QDP::Layout::subgridLattSize()[0];
  int ly = QDP::Layout::subgridLattSize()[1];
  int lz = QDP::Layout::subgridLattSize()[2];
  int lt = QDP::Layout::subgridLattSize()[3];
    //Geometry
    bfma.node_latt[0] = lx;
    bfma.node_latt[1] = ly;
    bfma.node_latt[2] = lz;
    bfma.node_latt[3] = lt;
    
    multi1d<int> procs = QDP::Layout::logicalSize();
    for(int mu=0;mu<4;mu++){
      if (procs[mu]>1) bfma.local_comm[mu] = 0;
      else             bfma.local_comm[mu] = 1;
    }
    
    // Bfm object
    bfm_qdp<double> bfm_eig; 
    bfm_eig.init(bfma);

    //Gauge field import
    bfm_eig.importGauge(u);

    //Subspace
#define NumberGaussian (1)
  Fermion_t subspace[NumberGaussian];
  Fermion_t check;
  Fermion_t mp;
  Fermion_t mmp;
  Fermion_t tmp_t;
  check = bfm_eig.allocFermion();
     mp = bfm_eig.allocFermion();
    mmp = bfm_eig.allocFermion();
  tmp_t = bfm_eig.allocFermion();
  bfm_eig.importFermion(src,check,1);

  QDPIO::cout << "Ls = "<<Ls<<endl;
  for(int g=0;g<NumberGaussian;g++){
    for(int s=0;s<Ls;s++){
      gaussian(src[s]);
    }
    subspace[g]=bfm_eig.allocFermion();
    bfm_eig.importFermion(src,subspace[g],1); // Half parity gaussian
    if ( g==0) {
      bfm_eig.importFermion(src,check,1);
    }
    for(int s=0;s<Ls;s++){
      src[s]=zero;
    }
    bfm_eig.exportFermion(src,subspace[g],1);
    QDPIO::cout << "Subspace norm " << norm2(src)<<endl;
  }
  for(int s=0;s<Ls;s++){
    gaussian(src[s]);
  }
  QDPIO::cout << "Got here " << endl;

  //  Handle< LinearOperatorArray<T> > linop =GetLinOp(u, parms);
  int block[5];
  for(int i=0;i<5;i++) block[i]=4;

  QDPIO::cout << "Initialised dirac op"<<endl;
  BfmMultiGrid ldop(Ls,NumberGaussian,block,subspace,&bfm_eig);

  int ns = ldop.SubspaceDimension();
  QDPIO::cout << "subspace dimension is "<< ns<<endl;
  ns = ldop.SubspaceLocalDimension();
  QDPIO::cout << "subspace dimension per node is "<< ns<<endl;

  std::vector<std::complex<double> > decomp(ns);
  ldop.ProjectToSubspace(check,decomp);
  if (QMP_is_primary_node()){
    FILE * fp = fopen("coeff.dat","w");
    for(int s=0;s<ns;s++){
      fprintf(fp,"coeff %d %le %le\n",s,real(decomp[s]),imag(decomp[s]));
    }
    fclose(fp);
  }
  for(int s=0;s<ns;s++){
    QDPIO::cout << "coeff "<<s<<" " << real(decomp[s]) << " " << imag(decomp[s])<<endl;
  }
  ldop.PromoteFromSubspace(decomp,mp);
  double n;
#pragma omp parallel 
  {
    omp_set_num_threads(bfm_eig.nthread);
#pragma omp for 
    for(int t=0;t<bfm_eig.nthread;t++) {
      bfm_eig.axpy(check,mp,check,-1);
      n = bfm_eig.norm(check);
    }
  }
  QDPIO::cout << "project/promote n2diff "<< n<<endl;
  QMP_barrier();

QDPIO::cout << "Computing little dirac matrix"<<endl;
  ldop.ComputeLittleMatrixColored();

  QDPIO::cout << "Done"<<endl;

  std::vector<std::complex<double> > Aphi(ns);
  //        phi^dag DdagD phi = |Dphi|^2 with phi a subspace vector
  //        should be equal to Project/Apply/Promote + inner product

#pragma omp parallel 
  {
#pragma omp for 
    for(int t=0;t<bfm_eig.nthread;t++) {
      bfm_eig.Mprec(subspace[0],mp,tmp_t,0);
    }
  }

  QDPIO::cout << "Applied BFM matrix "<<endl;

  double n2;
#pragma omp parallel 
  {
    omp_set_num_threads(bfm_eig.nthread);
#pragma omp for 
    for(int t=0;t<bfm_eig.nthread;t++) {
      n2 = bfm_eig.norm(mp);
    }
  }

  QDPIO::cout << "Applied BFM matrix "<<n2<<endl;

  ldop.ProjectToSubspace(subspace[0],decomp);
  QDPIO::cout << "Projected to subspace "<<endl;
  ldop.Apply(decomp,Aphi);
  QDPIO::cout << "Applied A "<<endl;
  ldop.PromoteFromSubspace(Aphi,check);
  QDPIO::cout << "Promoted "<<endl;

  complex<double> inn;
#pragma omp parallel 
  {
#pragma omp for 
    for(int t=0;t<bfm_eig.nthread;t++) {
      inn = bfm_eig.inner(subspace[0],check);
    }
  }

  QDPIO::cout << "phi^dag Ddag D phi check " << n2 << " " <<real(inn) << imag(inn) <<endl;

  std::vector<std::complex<double> > AinvAphi(ns);
  ldop.ProjectToSubspace(subspace[0],decomp);
  ldop.Apply(decomp,Aphi);
  for(int s=0;s<ns;s++){
    QDPIO::cout << "Aphi "<<s<<" " << real(Aphi[s]) <<" " << imag(Aphi[s])<<endl;
  }
  ldop.PromoteFromSubspace(Aphi,check);

#pragma omp parallel 
  {
#pragma omp for 
    for(int t=0;t<bfm_eig.nthread;t++) {
      bfm_eig.Mprec(subspace[0],mp,tmp_t,0);
      bfm_eig.Mprec(mp,mmp,tmp_t,1);
    }
  }
  ldop.ProjectToSubspace(mmp,decomp);
  ldop.PromoteFromSubspace(decomp,mmp);
#pragma omp parallel 
  {
#pragma omp for 
    for(int t=0;t<bfm_eig.nthread;t++) {
      bfm_eig.axpy(check,mmp,check,-1.0);
      n2 = bfm_eig.norm(check);
    }
  }
  QDPIO::cout << "PMdagMP check n2diff "<< n2<<endl;


  QMP_barrier();
  QDPIO::cout << "Applying inverse"<<endl;
  ldop.ApplyInverse(Aphi,AinvAphi);
  QMP_barrier();
  for(int s=0;s<ns;s++){
    QDPIO::cout << "AinvAphi "<<s<<" " << real(AinvAphi[s]) << " " << imag(AinvAphi[s])<<endl;
  }
  ldop.PromoteFromSubspace(AinvAphi,check);

#pragma omp parallel 
  {
#pragma omp for 
    for(int t=0;t<bfm_eig.nthread;t++) {
      bfm_eig.axpy(check,subspace[0],check,-1.0);
      n2 = bfm_eig.norm(check);
    }
  }
  QDPIO::cout << "AinvA check n2diff "<< n2<<endl;
  

}
