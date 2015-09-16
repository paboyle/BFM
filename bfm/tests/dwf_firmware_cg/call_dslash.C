#include <bfm.h>
#include <stdio.h>
#include <stdlib.h>

#define _BGQ_TESTINT_DCR__THREAD_ACTIVE0_RB  (0x0C805A) 
#define _BGQ_TESTINT_DCR__THREAD_ACTIVE1_RB  (0x0C805B) 
#define _BGQ_PRIV_DCR_BASE 0x3ffe0000000LL

inline uint64_t  BgDcrReadPriv(uint64_t dcr)
{
    register volatile uint64_t *p = (uint64_t *)(_BGQ_PRIV_DCR_BASE);
    return(p[dcr]);
    
}

int NTHREAD;

void * thr_main(void *p);

Fermion_t psi_h;
Fermion_t chi_h;

extern double src[];
extern double gauge[];


typedef double Float;
class singletest : public bfm {
public:

  void staticGauge(void)
  {
    this->u = gauge;
  }


};


singletest    *dwf;
singletest    dwf_private;
double sum_array[1024];

void setup(void);
void compute(int me);

ThreadModelFirmware TMF;

extern "C" {
void test_exit(int);
int test_main(int argc, char **argv);

int test_main(int argc, char **argv)
{
  char buffer[1024];

  uint64_t A0 = BgDcrReadPriv(_BGQ_TESTINT_DCR__THREAD_ACTIVE0_RB);
  uint64_t A1 = BgDcrReadPriv(_BGQ_TESTINT_DCR__THREAD_ACTIVE1_RB);
  
  int active=0;
  for(int i=0;i<64;i++){
    if ( A0 & (0x1UL<<(63-i)) ) active++;
    if ( A1 & (0x1UL<<(63-i)) ) active++;
  }

  NTHREAD = active ;

  TMF.thread_init_noalloc(NTHREAD);
  TMF.sum_array= sum_array;

  int me  = TMF.thread_barrier();

  if ( me == 0 ) {

    printf("Hello, world\n"); fflush(stdout);
    printf("%d active threads\n",active);

    malloc(0x1000000);
    dwf = new (& dwf_private) singletest();
    
    printf("Calling setup\n");
    setup();
  }

  TMF.thread_barrier();

  compute(me);

  test_exit(0);
}
}
void setup(void)
{
  int lx = 4;
  int ly = 4;
  int lz = 4;
  int lt = 4;

  bfmarg dwfa;
 dwfa.solver = DWFrb4d;
  dwfa.threads = NTHREAD;

  dwfa.node_latt[0]  = lx;
  dwfa.node_latt[1]  = ly;
  dwfa.node_latt[2]  = lz;
  dwfa.node_latt[3]  = lt;

  dwfa.local_comm[0]  = 1;
  dwfa.local_comm[1]  = 1;
  dwfa.local_comm[2]  = 1;
  dwfa.local_comm[3]  = 1;

  dwfa.Ls = 4;
  dwfa.mass = 0.0;
  dwfa.Csw  = 0.0;

  printf("Initialising bfm operator\n");
  dwf->init(dwfa);

  printf("Initialising arrays\n");

  psi_h = src;
  chi_h = dwf->allocFermion();
  
  dwf->staticGauge(); 

  printf("Setup complete\n");

}

void compute(int me)
{
  double delta;
  double n1;
  double n2;
  int idx=0;

  if ( me == 0 ) {
    printf("Calling CG\n");
  }
      
  dwf->thread_barrier();

  dwf->CGNE_prec(chi_h,psi_h);

  uint64_t csum = this->checksum(chi_h);
  if ( me == 0 ) printf("Result checksum is %lx\n",csum);

  test_exit(0);

}







