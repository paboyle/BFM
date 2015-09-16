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
 double tol=1.0e-7;

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

    std::string output=outputdir+"/HtCayleyZolo.xml";
    test_solver(u,HtCayleyZolo,output,configuration,Ls);

    //    output=outputdir+"/HwCayleyZolo.xml";
    //    test_solver(u,HwCayleyZolo,output,configuration,Ls);
    //output=outputdir+"/HtCayleyTanh.xml";
    //test_solver(u,HtCayleyTanh,output,configuration,Ls);
    //output=outputdir+"/HmCayleyTanh.xml";
    //test_solver(u,HmCayleyTanh,output,configuration,Ls);
    //output=outputdir+"/HwCayleyTanh.xml";
    //test_solver(u,HwCayleyTanh,output,configuration,Ls);
    return;

    {

      LatticePropagator sol=zero;
      LatticePropagator src=zero;

      // Point source at origin
      multi1d<int> site(4);
      site[0]=site[1]=site[2]=site[3]=0;
      Real one = 1.0;
      Propagator F = one;
      pokeSite(src,F,site);

      output=outputdir+"/DWF.xml";

      Real resid[1] = { Real(tol) };
      int  max_iter[1] = { 40000 };
      dwf_CG<double,double>(sol,src,u,output,Ls,false,
			Real (mq),
			Real (M5),1,resid,max_iter);
    
    }

}

void test_solver(multi1d<LatticeColorMatrix> &u,
		 BfmSolver solver,
		 std::string output, std::string configuration,int Ls)
{
  
  if ( solver == HtCayleyTanh ) { 
    Printf("Testing HtCayleyTanh aka DWF\n");
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
    Ls = 17;
    Printf("Partial fraction: setting Ls =%d\n",Ls);
  } else if ( solver == HwContFracZolo ) { 
    Printf("Testing HwContFracZolo \n");
    Ls = 17;
  }



  LatticeFermion sol;
  LatticeFermion check;
  LatticeComplex PA=zero;
  LatticeComplex PP=zero;
  LatticeComplex PAconsv=zero;
  LatticeComplex PJ5q=zero;

  for(int c=0;c<3;c++){
    for(int s=0;s<4;s++){
    
      LatticeFermion src=zero;
      multi1d<int> site(4);
      site[0]=site[1]=site[2]=site[3]=0;

      ColorVector cv=zero;
      Fermion f=zero;
      pokeColor(cv,Real(1.0),c);
      pokeSpin(f,cv,s);
      pokeSite(src,f,site);

      sol=zero;

      ////////////////////////////////
      // solve with BAGEL
      ////////////////////////////////

      {
	LatticeComplex PA_t=zero;
	LatticeComplex PP_t=zero;
	LatticeComplex PAconsv_t=zero;
	LatticeComplex PJ5q_t=zero;
      
	bfm_g5d_CG_unprec<double>(sol,src,u,
				  PA_t,PP_t,PAconsv_t,PJ5q_t,
				  solver,Ls,mq, M5,tol,40000);

	PA+=PA_t;
	PP+=PP_t;
	PAconsv+=PAconsv_t;
	PJ5q+=PJ5q_t;
      }
    }
  }

  {
    Set &tslice = GetTimeslice();
    int length = tslice.numSubsets();
    multi1d<DComplex> corr(length);
    multi1d<Real>    rcorr(length);
	
    XMLFileWriter xml(output);
    push(xml,"DWF_prop");
    push(xml,"DWF_observables");

    corr = sumMulti(PP,tslice);  
    for(int t=0;t<length;t++)rcorr[t]=real(corr[t]);
    write(xml,"PP",rcorr);

    corr = sumMulti(PA,tslice);  
    for(int t=0;t<length;t++)rcorr[t]=real(corr[t]);
    write(xml,"PA",rcorr);
      
    corr = sumMulti(PAconsv,tslice);  
    for(int t=0;t<length;t++)rcorr[t]=real(corr[t]);
    write(xml,"PAconsv",rcorr);

    corr = sumMulti(PJ5q,tslice);  
    for(int t=0;t<length;t++)rcorr[t]=real(corr[t]);
    write(xml,"PJ5q",rcorr);
    
    pop(xml);     
    pop(xml);      
  }

}


