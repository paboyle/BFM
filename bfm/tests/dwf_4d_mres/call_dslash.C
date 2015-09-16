#include <chroma.h>
#include <bfm.h>
#include <bfm_qdp.h>
#include <bfm_qdp_g5d.h>

using namespace Chroma;
typedef multi1d<LatticeColorMatrix> U;

#define Printf if ( QMP_is_primary_node() ) printf

typedef LatticeFermion T;
typedef multi1d<LatticeFermion> T5;
typedef multi1d<LatticeColorMatrix> U;

#include <bfmcommqmp.h>
#include <bfmcommspi.h>
typedef bfmcommQMP<double> bfm_qmp;

void Tinv(LatticeFermion &sol, 
	  LatticeFermion &src, 
	  multi1d<LatticeColorMatrix> u,
	  double mass,
	  double M5,int inv=1);

void test_solver(BfmSolver solver);
Handle< UnprecLinearOperatorArray<T,U,U> > GetLinOp(U u, int solver,int Ls, double M5, double mq,double eps);

int main (int argc,char **argv )
{
  Chroma::initialize(&argc,&argv);
  WilsonTypeFermActs4DEnv::registerAll(); 

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
  nrow[0] = 32;
  nrow[1] = 32;
  nrow[2] = 32;
  nrow[3] = 32;
  int Ls = 16;
#endif

  Layout::setLattSize(nrow);
  Layout::create();


  test_solver(HtCayleyTanh);

}

void test_solver(BfmSolver solver)
{
  int Ls=8;

  multi1d<LatticeColorMatrix> u(4);
  HotSt(u);

  double M5=1.0;
  double mq=0.0;
  LatticeFermion src;
  LatticeFermion sol;
  LatticePropagator dwfprop;
  LatticePropagator prop;
  LatticePropagator source_prop;
  src=zero;
  sol=zero;

  double residual=1.0e-10;
  int    max_iter=10000;

  for(int c=0;c<3;c++){
  for(int s=0;s<4;s++){

    ColorVector cv=zero;
    Fermion f     =zero;
    LatticeFermion src = zero;

    Complex cone=cmplx(Real(1.0),Real(0.0));
    pokeColor(cv,cone,c);
    pokeSpin(f,cv,s);

    multi1d<int> coords(4);
    QDPIO::cout << "source: "<< f << endl;
    coords[0] = coords[1] = coords[2] = coords[3] = 0;

    pokeSite(src,f,coords);

    bfm_g5d_CG_unprec<double>(sol,src,u,solver,Ls,mq, M5,residual,max_iter);

    BfmFermToProp(src,source_prop, c,s);
    BfmFermToProp(sol,prop, c,s);

  }}

  QDPIO::cout<<endl;
  QDPIO::cout<<"================================="<<endl;
  QDPIO::cout<<"Mres using DWF"<<endl;
  QDPIO::cout<<"================================="<<endl;
  QDPIO::cout<<endl;
  // Do things the mres way
  int Ncg=1;
  Real resid[1];  resid[0]=1.0e-9;
  int maxit[1];   maxit[0]=10000;
  LatticePropagator midpoint;
  //dwf_CG<double,double>(dwfprop,source_prop,u,"/dwf5d",Ls,0,mq,M5,Ncg,resid,maxit,midpoint);

  QDPIO::cout<<endl;
  QDPIO::cout<<"================================="<<endl;
  QDPIO::cout<<"Mres using midpoint to restore surface"<<endl;
  QDPIO::cout<<"================================="<<endl;
  QDPIO::cout<<endl;
  
  
  for(int c=0;c<3;c++){
  for(int sp=0;sp<4;sp++){
  {
    double mass=mq;
    LatticeFermion ovlap_surface;
    LatticeFermion dov_reconstruct;

    LatticeFermion tmp;
    LatticeFermion tmp1;

    BfmPropToFerm(dwfprop,ovlap_surface,c,sp);
    BfmPropToFerm(source_prop,src,c,sp);
    ovlap_surface = ovlap_surface* (1.0-mass) + src;

    BfmPropToFerm(midpoint,tmp,c,sp);
    QDPIO::cout << "Reconstruction midpoint J5q norm "<<norm2(tmp)<<endl;

    for(int s=0;s<Ls/2;s++){
      Tinv(tmp1,tmp,u,mass,M5,1);
      tmp=tmp1;
    }
    dov_reconstruct = tmp;

    BfmPropToFerm(midpoint,tmp,c,sp);
    for(int s=0;s<Ls/2;s++){
      Tinv(tmp1,tmp,u,mass,M5,0);
      tmp=tmp1;
    }
    dov_reconstruct += tmp;
  
    tmp=dov_reconstruct - ovlap_surface;
    QDPIO::cout << "Reconstruction of the surface prop error"<<norm2(tmp)<<endl;
    QDPIO::cout << "Reconstruction of the surface prop dov_reconstruc"<<
      norm2(dov_reconstruct)<<endl;
    QDPIO::cout << "Reconstruction of the surface prop ovlap_surface"<<norm2(ovlap_surface)<<endl;

  }}}

  // Do things the deltaL way

  QDPIO::cout<<endl;
  QDPIO::cout<<"================================="<<endl;
  QDPIO::cout<<"Mres using deltaL"<<endl;
  QDPIO::cout<<"================================="<<endl;
  QDPIO::cout<<endl;

  {
    multi1d<Complex> PPcorr;
    multi1d<Complex> PAcorr; 
    multi1d<Complex> PJ5qcorr;

    bfmDwf4d Cayley;

    Cayley.init(DWF, mq, M5, Ls, u,residual,max_iter);
    Cayley.mres(PPcorr,PAcorr,PJ5qcorr,prop,source_prop,midpoint);
    Cayley.end();
    
    int length=PJ5qcorr.size();
    multi1d<Real>    rcorr(length);

    for(int t=0;t<length;t++) {
      rcorr[t]=real(PPcorr[t]);
      QDPIO::cout << "PP["<< t<<"] = " << rcorr[t]<<endl;
    }
    for(int t=0;t<length;t++) {
      rcorr[t]=real(PAcorr[t]);
      QDPIO::cout << "PA["<< t<<"] = " << rcorr[t]<<endl;
    }
    for(int t=0;t<length;t++) {
      rcorr[t]=real(PJ5qcorr[t]);
      QDPIO::cout << "PJ5q["<< t<<"] = " << rcorr[t]<<endl;
    }
    
  }

  Printf("Done!!\n");

}

  





void Tinv(LatticeFermion &sol, 
	  LatticeFermion &src, 
	  multi1d<LatticeColorMatrix> u,
	  double mass,
	  double M5, 
	  int inv)
{
  // Set up BAGEL object
  double residual=1.0e-9;
  int max_iter   = 10000;

  BfmSolver solver = DWFTransfer;
  if ( inv )solver = DWFTransferInv;

  int lx = QDP::Layout::subgridLattSize()[0];
  int ly = QDP::Layout::subgridLattSize()[1];
  int lz = QDP::Layout::subgridLattSize()[2];
  int lt = QDP::Layout::subgridLattSize()[3];

  multi1d<int> procs = QDP::Layout::logicalSize();

  bfmarg bfma;

  //Physics parameters
  bfma.Ls           = 1;
  bfma.M5           = toDouble(M5);
  bfma.mass         = toDouble(mass);
  bfma.precon_5d    = 0;
  bfma.solver       = solver;
  bfma.list_engine  = 0;
  bfma.max_iter     = max_iter;
  bfma.residual     = toDouble(residual);

  //Geometry
  bfma.node_latt[0] = lx;
  bfma.node_latt[1] = ly;
  bfma.node_latt[2] = lz;
  bfma.node_latt[3] = lt;
  for(int mu=0;mu<4;mu++){
    if (procs[mu]>1) bfma.local_comm[mu] = 0;
    else             bfma.local_comm[mu] = 1;
  }

  QDPIO::cout << "bfm_qdp:: Initialising BAGEL-2 G5D solver for transfer matrix "<<endl;

  bfm_qdp<double> bfm; 

  bfm.init(bfma);

  bfm.importGauge(u);

  Fermion_t sol_t[2];
  Fermion_t src_t[2];
  for(int cb=0;cb<2;cb++){
    sol_t[cb] = bfm.allocFermion();
    src_t[cb] = bfm.allocFermion();
    bfm.importPhysicalFermion(src,src_t[cb],cb);
    bfm.importPhysicalFermion(sol,sol_t[cb],cb);
  }
  printf("src_norm = %le\n",bfm.norm(src_t));
  bfm.G5D_Transfer(src_t,sol_t);
  printf("sol_norm = %le\n",bfm.norm(sol_t));
  
  for(int cb=0;cb<2;cb++){
    bfm.exportPhysicalFermion(sol,sol_t[cb],cb);
    bfm.freeFermion(sol_t[cb]);
    bfm.freeFermion(src_t[cb]);
  }
  bfm.end();

  // Could get CHROMA linop and verify now. but need 5d solution.
}



