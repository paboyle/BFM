#include <chroma.h>
#include <actions/ferm/invert/syssolver_linop_cg_array.h>
#include <actions/ferm/invert/syssolver_linop_aggregate.h>
#include <bfm.h>
#include <bfm_qdp.h>
#include <bfm_qdp_g5d.h>
#include <bfm_qdp_dwf.h>
#include <nested_4d.h>

using namespace Chroma;

#define Printf if ( QMP_is_primary_node() ) printf

typedef multi1d<LatticeColorMatrix> U;
typedef LatticeFermion T;
typedef multi1d<LatticeFermion> T5;

#include <bfmcommqmp.h>
#include <bfmcommspi.h>
typedef bfmcommQMP<double> bfm_qmp;

void test_solver(multi1d<LatticeColorMatrix> & u,BfmSolver solver,LatticeFermion & src, std::string configuration,int Ls);

  double M5=1.8;
  double mq=0.005;
 double tol=1.0e-10;

int main (int argc,char **argv )
{
  Chroma::initialize(&argc,&argv);

  WilsonTypeFermActs4DEnv::registerAll(); 
  LinOpSysSolverArrayEnv::registerAll();

  /********************************************************
   * Command line parsing
   ********************************************************
   */
  if ( argc != 6 ) { 
   Printf("Usage: %s lx ly lz lt <conf>\n All must be even\n",argv[0]);
   Printf("argc is %d\n",argc);
   for ( int i=0;i<argc;i++){
      Printf("%d %s\n",i,argv[i]);
   }
      exit(-1);
  }
  std::string configuration(argv[5]);

  /********************************************************
   * Setup QDP
   ********************************************************
   */
  multi1d<int> nrow(Nd);
  nrow[0] = atoi(argv[1]);
  nrow[1] = atoi(argv[2]);
  nrow[2] = atoi(argv[3]);
  nrow[3] = atoi(argv[4]);
  int threads = 64;

  Layout::setLattSize(nrow);
  Layout::create();

  bfmarg::Threads(threads);
  bfmarg::Reproduce(0);
  bfmarg::ReproduceChecksum(0);
  bfmarg::ReproduceMasterCheck(0);
  bfmarg::Verbose(0);

  multi1d<LatticeColorMatrix> u(4);
  ArchivGauge_t Header;
  readArchiv(Header,u,configuration);
  //HotSt(u);

  std::string output;

  LatticeFermion src=zero;
  //Z2
  Real nrm = 1./sqrt(2.0);
  Real nnrm = -1./sqrt(2.0);
  for(int cc=0;cc<3;cc++){
    LatticeSpinVector sv;
    for(int ss=0;ss<4;ss++){

      LatticeComplex cx;
      LatticeReal re,im;
      LatticeReal rnd;
      
      random(rnd);
      re = where(rnd>0.5,LatticeReal(nrm),LatticeReal(nnrm));
      random(rnd);
      im = where(rnd>0.5,LatticeReal(nrm),LatticeReal(nnrm));
      cx = cmplx(re,im);
      
      pokeSpin(sv,cx,ss);
    }
    pokeColor(src,sv,cc);
  }

  int Lss[] = { 16, 12, 8};

  for (int i=0;i<3;i++){
    QDPIO::cout <<"*************************" <<endl;
    QDPIO::cout <<"* Ls = "<<Lss[i] <<endl;
    QDPIO::cout <<"*************************" <<endl;
    int Ls = Lss[i];
    test_solver(u,HtCayleyTanh,src,configuration,Ls);
    test_solver(u,HmCayleyTanh,src,configuration,Ls);
    test_solver(u,HtCayleyZolo,src,configuration,Ls);
    test_solver(u,HwCayleyZolo,src,configuration,Ls);
  }

  exit(0);
}

void test_solver(multi1d<LatticeColorMatrix> &u,
		 BfmSolver solver,
		 LatticeFermion & src,
		 std::string configuration,
		 int Ls)

{  
  g5dParams parms;

  double M5=1.8;
  double mq=0.1;
  double wilson_lo = 0.05;
  double wilson_hi = 6.8;
  double shamir_lo = 0.025;
  double shamir_hi = 1.7;
  double ht_scale=2.0;
  double hw_scale=1.0;

  if ( solver == HtCayleyTanh ) { 
    Printf("Testing HtCayleyTanh aka DWF\n");
    parms.ShamirCayleyTanh(mq,M5,Ls);
  } else if ( solver == HmCayleyTanh ) { 
    Printf("Testing HmCayleyTanh Moebius\n");
    parms.ScaledShamirCayleyTanh(mq,M5,Ls,ht_scale);
  } else if ( solver == HwCayleyTanh ) { 
    Printf("Testing HwCayleyTanh aka Borici\n");
    parms.WilsonCayleyTanh(mq,M5,Ls,hw_scale);
  } else if ( solver == HtCayleyZolo ) { 
    Printf("Testing HtCayleyZolo\n");
    parms.ShamirCayleyZolo(mq,M5,Ls,shamir_lo,shamir_hi);
  } else if ( solver == HwCayleyZolo ) { 
    Printf("Testing HwCayleyZolo aka Chiu\n");
    parms.WilsonCayleyZolo(mq,M5,Ls,wilson_lo,wilson_hi);
  } else { 
    Printf("Unknown testcase \n");
    exit(-1);
  }

  LatticeFermion sol=zero;
  LatticeFermion delta=zero;
  Complex PJ5q=zero;

  ////////////////////////////////
  // solve with BAGEL
  ////////////////////////////////
  LatticeComplex PA;
  LatticeComplex PP;
  LatticeComplex PAc;
  LatticeComplex PJ5;
      
  if ( solver == DWF ) { 

    LatticePropagator solp=zero;
    LatticePropagator srcp=zero;
    BfmFermToProp(src, srcp, 0,0);

    int its[1] = { 40000 };
    Real ttol[1]= {tol};
    Real mmq(mq);
    Real MM5(M5);
    dwf_CG<double,double>(solp,srcp,u,"pooh.xml",Ls,0,mmq,MM5,1,ttol,its);
    BfmPropToFerm(solp, sol, 0,0);

  } else { 
    bfm_Cayley_CG<double>(sol,src,u,PP,PA,PAc,PJ5,parms,1.0e-12,80000);
  }

  QDPIO::cout << "Adding contact term and normalising"<<endl;
  LatticeFermion      ovlap_surface = sol * (1.0-mq) + src;
  
  bfm_g5d_DeltaL<double>(delta,ovlap_surface,u,parms,tol,40000);

  Real PPsum   = norm2(sol);
  PJ5q = innerProduct(ovlap_surface,delta);

  QDPIO::cout << "DeltaL " << PJ5q << " PP " << PPsum <<endl;
  QDPIO::cout << "Mres   " << PJ5q / PPsum <<endl;

}

