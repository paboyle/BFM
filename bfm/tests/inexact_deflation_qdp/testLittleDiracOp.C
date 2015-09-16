#include <qdp.h>
#include <bfm.h>
#include "littleDiracOp.h"

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


  bfmarg::Threads(1);
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
  //ArchivGauge_t Header ; readArchiv(Header,u,"ckpoint_lat.3000");  

  multi1d<LatticeFermion> src(Ls);
  multi1d<LatticeFermion> check(Ls);
  multi1d<LatticeFermion> checkD(Ls);

/* Rudy
   Calculate some eigenvectors */
#define NumberGaussian (1)


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
    bfma.Threads(1);
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


  multi1d< multi1d<LatticeFermion> > subspace(NumberGaussian);
  QDPIO::cout << "Ls = "<<Ls<<endl;
  for(int g=0;g<NumberGaussian;g++){
    subspace[g].resize(Ls);
    for(int s=0;s<Ls;s++){
      gaussian(subspace[g][s]);
      int par = s&0x1;
      subspace[g][s][rb[par]] = zero;
    }
    QDPIO::cout << "Subspace norm " << norm2(subspace[g])<<endl;
  }

  for(int s=0;s<Ls;s++){
    gaussian(src[s]);
  }


QDPIO::cout << "Got here " << endl;

//  Handle< LinearOperatorArray<T> > linop =GetLinOp(u, parms);
  check=subspace[0];
  multi1d<int> block(5);
  for(int i=0;i<5;i++) block[i]=4;

#define TESTS
#ifdef TESTS

  QDPIO::cout << "Initialised dirac op"<<endl;
  LittleDiracOperator ldop(block,subspace);

  int ns = ldop.SubspaceDimension();
  QDPIO::cout << "subspace dimension is "<< ns<<endl;

  std::vector<std::complex<double> > decomp(ns);

  ldop.ProjectToSubspace(subspace[0],decomp);
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
  ldop.PromoteFromSubspace(decomp,check);
  QDPIO::cout << "project/promote n2diff "<< norm2(check-subspace[0])<<endl;

  QDPIO::cout << "Computing little dirac matrix"<<endl;
  ldop.ComputeLittleMatrix(bfm_eig);

#if 1
  Fermion_t tmp_t = bfm_eig.allocFermion();
  Fermion_t   p_t = bfm_eig.allocFermion();
  Fermion_t  mp_t = bfm_eig.allocFermion();
  Fermion_t mmp_t = bfm_eig.allocFermion();
  bfm_eig.importFermion(subspace[0],p_t,1);
#pragma omp parallel 
  {
    omp_set_num_threads(bfm_eig.nthread);
#pragma omp for 
    for(int i=0;i<bfm_eig.nthread;i++) {
      bfm_eig.Mprec(p_t,mp_t,tmp_t,DaggerNo);  
      bfm_eig.Mprec(mp_t,mmp_t,tmp_t,DaggerYes);  
    }
  }

  for(int s=0;s<Ls;s++){
    check[s][rb[0]]=zero;
    checkD[s][rb[0]]=zero;
  }
  bfm_eig.exportFermion(check,mp_t,1);
  bfm_eig.exportFermion(checkD,mmp_t,1);
  bfm_eig.freeFermion( tmp_t );
  bfm_eig.freeFermion(   p_t );
  bfm_eig.freeFermion(  mp_t );
  bfm_eig.freeFermion( mmp_t );
#else
  (*linop)(check,subspace[0],PLUS);
#endif
  Real n2 = norm2(check);

  std::vector<std::complex<double> > Aphi(ns);
  //        phi^dag DdagD phi = |Dphi|^2 with phi a subspace vector
  //        should be equal to Project/Apply/Promote + inner product
  ldop.ProjectToSubspace(subspace[0],decomp);
  ldop.Apply(decomp,Aphi);
  ldop.PromoteFromSubspace(Aphi,check);
  Complex inn = innerProduct(subspace[0],check);
  QDPIO::cout << "phi^dag Ddag D phi check " << n2 << " " <<inn <<endl;

  ldop.ProjectToSubspace(checkD,decomp);
  ldop.PromoteFromSubspace(decomp,checkD);
  n2 = norm2(check - checkD);
  QDPIO::cout << "P Ddag D phi - A phi diff" << n2 <<endl;

  //
  // Subspace Inverse x A applied on vector already lying in subspace
  // should be the same as the orig vector
  //

  std::vector<std::complex<double> > AinvAphi(ns);
  ldop.ProjectToSubspace(subspace[0],decomp);
  ldop.Apply(decomp,Aphi);
  for(int s=0;s<ns;s++){
    QDPIO::cout << "Aphi "<<s<<" " << real(Aphi[s]) <<" " << imag(Aphi[s])<<endl;
  }
  ldop.ApplyInverse(Aphi,AinvAphi);
  for(int s=0;s<ns;s++){
    QDPIO::cout << "AinvAphi "<<s<<" " << real(AinvAphi[s]) << " " << imag(AinvAphi[s])<<endl;
  }
  ldop.PromoteFromSubspace(AinvAphi,check);
  QDPIO::cout << "AinvA check n2diff "<< norm2(check-subspace[0])<<endl;

#else


  multi1d<LatticeFermion> subspacetmp(Ls);
  
  multi1d<LatticeFermion> guess(Ls);
  multi1d<LatticeFermion> eta(Ls);
  multi1d<LatticeFermion> eta_s(Ls);
  multi1d<LatticeFermion> sol(Ls);
  multi1d<LatticeFermion> defect(Ls);
  multi1d<LatticeFermion> tmp(Ls);
  multi1d<LatticeFermion> residual(Ls);
  LatticeFermion ferm;
  gaussian(ferm);
  eta=zero;
  eta[0]   = chiralProjectPlus (ferm);
  eta[Ls-1]= chiralProjectMinus(ferm);


  bfm.links = u;

  //  subspace[0] = residual;

  // Prepare subspace Luscher's way


  for(int g=0;g<NumberGaussian;g++){
    for(int i=0;i<5;i++){
      bfm.invParam.BfmMatrix = BfmMat_Mdag;
      bfm.bfmInvert5d<double>(subspacetmp,subspace[g]);
      bfm.invParam.BfmMatrix = BfmMat_M;
      bfm.bfmInvert5d<double>(subspacetmp,subspace[g]);
      subspace[g]=subspacetmp;
    }
  }  

  LittleDiracOperator ldop(block,subspace);
  int ns = ldop.SubspaceDimension();

  QDPIO::cout << "subspace dimension is "<< ns<<endl;
  ldop.ComputeLittleMatrix(linop);

  QDPIO::cout << "UNDEFLATED"<<endl;
  bfm.bfmInvert5d<double>(sol,eta);
  (*linop)(tmp,sol,PLUS);
  (*linop)(residual,tmp,MINUS);
  residual=residual-eta;
  QDPIO::cout << "Residual norm "<<norm2(residual)<<endl;

  multi1d<DComplex>  res(ns);
  multi1d<DComplex> Ares(ns);

  ldop.ProjectToSubspace(eta,res);
  ldop.PromoteFromSubspace(res,eta_s);
  ldop.ApplyInverse(res,Ares);
  ldop.PromoteFromSubspace(Ares,guess);
  QDPIO::cout << "test Guess norm "<<norm2(guess)<<"  Source norm "<< norm2(eta) << " eta_s "<< norm2(eta_s)<< endl;

  defect=zero;
  (*linop)(tmp,guess,PLUS);
  (*linop)(defect,tmp,MINUS);
  
  Real dn =norm2(defect-eta);
  Real rn =norm2(eta) ;
  QDPIO::cout << "n2diff = "<< dn/rn << " eta  norm " << rn << " resid "<<dn <<endl;

  ldop.ProjectToSubspace(defect,res);
  ldop.PromoteFromSubspace(res,defect);
  dn =norm2(defect-eta_s);
  rn =norm2(eta_s) ;
  QDPIO::cout << "n2diff = "<< dn/rn << " eta_s norm " << rn << " resid "<<dn <<endl;

  
  QDPIO::cout << "DEFLATED"<<endl;
  bfm.bfmInvert5d<double>(guess,eta);
  
  sol = sol + defect;
  (*linop)(tmp,sol,PLUS);
  (*linop)(residual,tmp,MINUS);
  residual=residual-eta;
  
  QDPIO::cout << "Residual norm "<<norm2(residual)<<endl;

#endif

}
