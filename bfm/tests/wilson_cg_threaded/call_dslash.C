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

//#define Printf if ( QMP_is_primary_node() ) printf
#define Printf printf

int main (int argc,char **argv )
{
  Chroma::initialize(&argc,&argv);


  /********************************************************
   * Command line parsing
   ********************************************************
   */
#define COMMANDLINE
#ifdef COMMANDLINE
  if ( argc != 6 ) { 
   Printf("Usage: %s lx ly lz lt threads\n All must be even\n",argv[0]);
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
  int threads = atoi(argv[5]);
#else
#if 0
  nrow[0] = 16;
  nrow[1] = 16;
  nrow[2] = 16;
  nrow[3] = 16;
#else
  nrow[0] = 32;
  nrow[1] = 32;
  nrow[2] = 64;
  nrow[3] = 8;
#endif
  int threads = 64;
#endif

  Layout::setLattSize(nrow);
  Layout::create();

  bfmarg::Threads(threads);
  bfmarg::Reproduce(1);
  bfmarg::ReproduceChecksum(0);
  bfmarg::ReproduceMasterCheck(0);
  bfmarg::Verbose(1);
  
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
  //  u[3] = zero;

  /********************************************************
   * Gaussian source and result vectors
   ********************************************************
   */
  LatticeFermion psi;
  LatticeFermion chi;
  LatticeFermion chi_qdp;

  int howmany=0;
  while (howmany++ < 3000) {

    gaussian(psi);
    gaussian(chi);
    gaussian(chi_qdp);

    /********************************************************
     * Import gauge field to BAGEL
     ********************************************************
     */

    Real M5 = 0;
    Real mass = -0.1;
    Real mq = mass;

#define Ncg 1
    Real residuals[Ncg] ; 
    int  max_iter[Ncg];

    for(int i=0;i<Ncg;i++) residuals[i] = 1.0e-12;
    for(int i=0;i<Ncg;i++) max_iter[i]  = 1000;

    wilson_CG<float,double>(chi,psi,u,mass,Ncg, residuals,max_iter);

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
 
  }
  Printf("Done\n"); 

}








