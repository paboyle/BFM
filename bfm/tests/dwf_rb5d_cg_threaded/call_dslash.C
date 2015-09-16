 
#include <chroma.h>
#include <actions/ferm/invert/syssolver_linop_cg_array.h>
#include <bfm.h>
#include <bfm_qdp.h>

#include "/bgsys/drivers/ppcfloor/hwi/include/bqc/nd_rese_dcr.h"


typedef bfm_dp bfm_t;

using namespace Chroma;

typedef LatticeFermion T;
typedef multi1d<LatticeFermion> T5;
typedef multi1d<LatticeColorMatrix> U;

#define Printf if ( QMP_is_primary_node() ) printf
//#define Printf printf

int main (int argc,char **argv )
{
  Chroma::initialize(&argc,&argv);


  /********************************************************
   * Command line parsing
   ********************************************************
   */
#define COMMANDLINE
#ifdef COMMANDLINE
  if ( argc != 7 ) { 
   Printf("Usage: %s lx ly lz lt Ls threads\n All must be even\n",argv[0]);
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
  int Ls = atoi(argv[5]);
  int threads = atoi(argv[6]);
#else
#if 0
  nrow[0] = 16;
  nrow[1] = 16;
  nrow[2] = 16;
  nrow[3] = 16;
  int Ls = 16;
#else
  nrow[0] = 32;
  nrow[1] = 32;
  nrow[2] = 64;
  nrow[3] = 8;
  int Ls = 8;
#endif
  int threads = 64;
#endif

  Layout::setLattSize(nrow);
  Layout::create();

  bfmarg::Threads(threads);
  bfmarg::Reproduce(1);
  bfmarg::ReproduceChecksum(0);
  bfmarg::ReproduceMasterCheck(0);
  bfmarg::Verbose(0);
  
  /********************************************************
   * Gaussian gauge field
   ********************************************************
   */

  multi1d<LatticeColorMatrix>  u(Nd); HotSt(u);
  multi1d<LatticeColorMatrixF> uf(Nd);
  for(int m=0; m < Nd; ++m){
    uf[m] = u [m];
    u [m] = uf[m];
  }

  /********************************************************
   * Gaussian source and result vectors
   ********************************************************
   */
  multi1d<LatticeFermion> psi(Ls);
  multi1d<LatticeFermion> chi(Ls);
  multi1d<LatticeFermion> chi_qdp(Ls);

  int howmany=0;
  while (howmany++ < 20000) {
    //while (howmany++ < 3) {

    for(int s=0;s<Ls;s++) gaussian(psi[s]);
    for(int s=0;s<Ls;s++) gaussian(chi_qdp[s]);

    /********************************************************
     * Import gauge field to BAGEL
     ********************************************************
     */

    Real M5 = 1.8;
    Real mass = 0.1;
    Real mq = mass;

#define Ncg 1
    Real residuals[Ncg] ; 
    int  max_iter[Ncg];

    

    for(int i=0;i<Ncg;i++) residuals[i] = 1.0e-12;
    for(int i=0;i<Ncg;i++) max_iter[i]  = 1000;

    QDPIO::cout << "*******************************************"<<endl;
    QDPIO::cout << "* Double precision                         "<<endl;
    QDPIO::cout << "*******************************************"<<endl;
    for(int s=0;s<Ls;s++) chi[s]=zero;
    dwf_restarted_invert<float,double>(CG_PREC_M,
				       chi, 
				       psi,
				       u,
				       "pooh.xml",
				       Ls,
				       mass,
				       M5,Ncg,residuals,max_iter);


    QDPIO::cout << "*******************************************"<<endl;
    QDPIO::cout << "* Single precision                         "<<endl;
    QDPIO::cout << "*******************************************"<<endl;
    for(int s=0;s<Ls;s++) chi[s]=zero;
    dwf_restarted_invert<float,float>(CG_PREC_M,
				       chi, 
				       psi,
				       u,
				       "pooh.xml",
				       Ls,
				       mass,
				       M5,Ncg,residuals,max_iter);


    QDPIO::cout << "*******************************************"<<endl;
    QDPIO::cout << "* Mixed precision                         "<<endl;
    QDPIO::cout << "*******************************************"<<endl;
    for(int s=0;s<Ls;s++) chi[s]=zero;

    double residual = 1.0e-12;
    int maxit=10000;
    dwf_mixed_precision_CG(chi,psi,u,Ls,mass,M5,residual,maxit);



    char link_name[ND_RESE_DCR_num][10] = { "A-", "A+", "B-", "B+", "C-", "C+", "D-", "D+", "E-", "E+", "IO" };
    uint32_t i;
    for (i = 0; i < ND_RESE_DCR_num; i++)
      {
	uint64_t val_re = DCRReadUser(ND_RESE_DCR(i, RE_LINK_ERR_CNT));
	uint64_t val_retran = DCRReadUser(ND_RESE_DCR(i, SE_RETRANS_CNT));
	if (val_re || val_retran)
	  printf("LINK Errors on RESE %s Recv Count = %ld Retran = %ld\n",
		 link_name[i],val_re,val_retran);
      }
 
#undef CHECKIT
#ifdef  CHECKIT 
    /********************************************************
     * QDP Linop
     ********************************************************
     */
    multi1d<int> bcs(Nd);
    bcs[0] = bcs[1] = bcs[2] = bcs[3] = 1;

    Handle< FermBC<T,U,U> > fbc(new SimpleFermBC< T, U, U >(bcs));
    Handle<CreateFermState<T,U,U> > cfs( new CreateSimpleFermState<T,U,U>(fbc));
    EvenOddPrecDWFermActArray  S_f(cfs, M5, mq, Ls);
    Handle< FermState<T,U,U> > fs( S_f.createState(u) );

    Handle< EvenOddPrecLinearOperatorArray<T,U,U> > M(S_f.precLinOp(fs,mq));
    Handle< LinearOperatorArray<T> > HM(S_f.precLinOp(fs,mq));

  // Check the result
    SysSolverCGParams invParam;
    invParam.RsdCG        = 1.e-7;
    invParam.RsdCGRestart = 1.e-7;
    invParam.MaxCG        = 1000;
    invParam.MaxCGRestart = 1000;

    Handle<LinOpSystemSolverArray<T> > Hsolver = new LinOpSysSolverCGArray<T>(HM,invParam);
    Handle<PrecFermAct5DQprop<T,U,U> > solver =  new PrecFermAct5DQprop<T,U,U>(M,Hsolver);
    SystemSolverResults_t  stats =  (*solver)(chi_qdp,psi);
    int n_count = stats.n_count;

    Double n2 = 0.0;
    for(int s=0;s<Ls;s++) { 
      Double sn2 = norm2(chi[s]-chi_qdp[s]);
      n2+=sn2;
#ifdef DEBUG
      QDPIO::cout << "|| Bagel - QDP || ["<<s<<"] = "<< sn2 << endl;
#endif
    }
    QDPIO::cout << "|| Bagel - QDP || = "<< n2 << endl;
    if ( toDouble(n2) >  1.0e-5 ) { 
      QDPIO::cout << "Difference is suspicious" << endl;
      exit(-1);
    }
    n2 = 0.0;
    for(int s=0;s<Ls;s++) { 
      Double sn2 = norm2(chi[s]);
      n2+=sn2;
#ifdef DEBUG
      QDPIO::cout << "|| Bagel || ["<<s<<"] = "<< sn2 << endl;
#endif
    }
    QDPIO::cout << "|| Bagel || = "<< n2 << endl;
  
    n2 = 0.0;
    for(int s=0;s<Ls;s++) {
      Double sn2 = norm2(chi_qdp[s]);
      n2+=sn2;
#ifdef DEBUG
      QDPIO::cout << "|| QDP || ["<<s<<"] = "<< sn2 << endl;
#endif
    }
    QDPIO::cout << "|| QDP || = "<< n2 << endl;
#endif
  }
  Printf("Done\n"); 

}








