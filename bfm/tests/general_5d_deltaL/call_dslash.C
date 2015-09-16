#include <chroma.h>
#include <actions/ferm/invert/syssolver_linop_cg_array.h>
#include <actions/ferm/invert/syssolver_linop_aggregate.h>
#include <bfm.h>
#include <bfm_qdp.h>
#include <bfm_qdp_g5d.h>
#include <bfm_qdp_dwf.h>

using namespace Chroma;

#define Printf if ( QMP_is_primary_node() ) printf

typedef multi1d<LatticeColorMatrix> U;
typedef LatticeFermion T;
typedef multi1d<LatticeFermion> T5;

#include <bfmcommqmp.h>
#include <bfmcommspi.h>
typedef bfmcommQMP<double> bfm_qmp;

void test_solver(multi1d<LatticeColorMatrix> & u,BfmSolver solver,std::string output, std::string configuration,int Ls);

  double M5=1.8;
  double mq=0.01;
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
#define COMMANDLINE
#ifdef COMMANDLINE
  if ( argc != 8 ) { 
   Printf("Usage: %s lx ly lz lt Ls <conf>\n All must be even\n",argv[0]);
   Printf("argc is %d\n",argc);
   for ( int i=0;i<argc;i++){
      Printf("%d %s\n",i,argv[i]);
   }
      exit(-1);
  }
  std::string configuration(argv[6]);
  std::string outputdir(argv[7]);
#else
  std::string configuration("/pooh");
  std::string outputdir("pooh");
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

  int threads = 64;

  Layout::setLattSize(nrow);
  Layout::create();

  bfmarg::Threads(threads);
  bfmarg::Reproduce(0);
  bfmarg::ReproduceChecksum(0);
  bfmarg::ReproduceMasterCheck(0);
  bfmarg::Verbose(1);

  multi1d<LatticeColorMatrix> u(4);
  if ( 1 ) { 
    ArchivGauge_t Header;
    readArchiv(Header,u,configuration);
  } else {
    HotSt(u);
  }


  int Ls=  atoi(argv[5]);

  std::string output;

  int Lss[] = { 16 , 12, 8};

  for (int i=0;i<1;i++){

    Ls = Lss[i];

    /*
    output=outputdir+"/HtCayleyTanh.xml";
    test_solver(u,DWF,output,configuration,Ls);

    output=outputdir+"/HtCayleyTanh.xml";
    test_solver(u,HtCayleyTanh,output,configuration,Ls);

    output=outputdir+"/HmCayleyTanh.xml";
    test_solver(u,HmCayleyTanh,output,configuration,Ls);
    */

    Ls++;
    output=outputdir+"/HwContFracZolo.xml";
    test_solver(u,HwContFracZolo,output,configuration,Ls);

    /* ZOLO MIN MAX not right
    output=outputdir+"/HtCayleyZolo.xml";
    test_solver(u,HtCayleyZolo,output,configuration,Ls);
    */

  }

  exit(0);
}

void test_solver(multi1d<LatticeColorMatrix> &u,
		 BfmSolver solver,
		 std::string output, std::string configuration,int Ls)

{  
  if ( solver == HtCayleyTanh ) { 
    Printf("Testing HtCayleyTanh aka DWF\n");
  } else if ( solver == DWF ) { 
    Printf("Testing DWF\n");
  } else if ( solver == HmCayleyTanh ) { 
    Printf("Testing HmCayleyTanh Moebius\n");
  } else if ( solver == HwCayleyTanh ) { 
    Printf("Testing HwCayleyTanh aka Borici\n");
  } else if ( solver == HtCayleyZolo ) { 
    Printf("Testing HtCayleyZolo\n");
  } else if ( solver == HwCayleyZolo ) { 
    Printf("Testing HwCayleyZolo aka Chiu\n");
  } else if ( solver == HwPartFracZolo ) { 
    Printf("Testing HwPartFracZolo aka KEK\n");
    if (Ls&0x1 == 0 ) {
      Printf("Partial fraction: Ls should be odd\n",Ls);
      exit(0);
    }
  } else if ( solver == HwContFracZolo ) { 
    Printf("Testing HwContFracZolo \n");
    if (Ls&0x1 == 0 ) {
      Printf("Partial fraction: Ls should be odd\n",Ls);
      exit(0);
    }
  }


  LatticeFermion sol;
  LatticeFermion delta;
  Real      PP=zero;
  Complex PJ5q=zero;

  for(int c=0;c<1;c++){
    for(int s=0;s<1;s++){
    
      LatticeFermion src=zero;

#undef POINT_SOURCE
#ifdef POINT_SOURCE
      multi1d<int> site(4);
      site[0]=site[1]=site[2]=site[3]=0;
      ColorVector cv=zero;
      Fermion f=zero;
      pokeColor(cv,Real(1.0),c);
      pokeSpin(f,cv,s);
      pokeSite(src,f,site);
#else 
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
#endif

      sol=zero;

      ////////////////////////////////
      // solve with BAGEL
      ////////////////////////////////
      LatticeComplex PA_t;
      LatticeComplex PP_t;
      LatticeComplex PAconsv_t;
      LatticeComplex PJ5q_t;
      
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
      bfm_g5d_CG_unprec<double>(sol,src,u,
				PA_t,PP_t,PAconsv_t,PJ5q_t,
				solver,Ls,mq, M5,tol,40000);
      }

      QDPIO::cout << "Adding contact term and normalising"<<endl;
      LatticeFermion      ovlap_surface = sol * (1.0-mq) + src;

      bfm_g5d_DeltaL<double>(delta,ovlap_surface,u,
			     solver,Ls,mq,M5,tol,40000);

      PP   = norm2(sol);
      PJ5q = innerProduct(ovlap_surface,delta);

      QDPIO::cout << "DeltaL " << PJ5q << " PP " << PP <<endl;
      QDPIO::cout << "Mres   " << PJ5q / PP <<endl;

      
    }
  }

}

