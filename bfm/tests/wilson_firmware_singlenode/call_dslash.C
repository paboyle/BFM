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
Fermion_t check;
Fermion_t diff ;

extern double src[];
extern double cb0dag0[];
extern double cb0dag1[];
extern double cb1dag0[];
extern double cb1dag1[];
extern double gauge[];

double * arrays[4] = { 
  cb0dag0,
  cb0dag1,
  cb1dag0,
  cb1dag1
};

typedef double Float;
class singletest : public bfm {
public:

  void staticGauge(void)
  {
    this->u = gauge;
  }

  void unitGauge(void)
  {
    int vol = node_latt[0]
            * node_latt[1]
            * node_latt[2]
            * node_latt[3];
    int words = 18 * vol;

    int bytes = sizeof(QDPdouble)*words;

    QDPdouble *gauge = (QDPdouble *)malloc(bytes);
    bzero(gauge,bytes);

    // No easy way to support bagel-2      
    // without CHROMA/QDP as I used QDP to do the cshift & adjoint
    // for double stored links. just use unit gauge.
    for(int v=0;v<vol;v++ ) {
      gauge[v*18]    = 1.0;
      gauge[v*18+8]  = 1.0;
      gauge[v*18+16] = 1.0;
    }
    int ndir=4;
    if ( doubleStored() ) { 
      ndir = 8;
    }
    for(int mu=0;mu<ndir;mu++) {
      importGauge(gauge, mu); 
    }
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

  int nrow[4];
  nrow[0] = lx;
  nrow[1] = ly;
  nrow[2] = lz;
  nrow[3] = lt;

  bfmarg dwfa;
  dwfa.solver = WilsonFermion;
  dwfa.threads = NTHREAD;

  dwfa.node_latt[0]  = lx;
  dwfa.node_latt[1]  = ly;
  dwfa.node_latt[2]  = lz;
  dwfa.node_latt[3]  = lt;

  dwfa.local_comm[0]  = 1;
  dwfa.local_comm[1]  = 1;
  dwfa.local_comm[2]  = 1;
  dwfa.local_comm[3]  = 1;

  dwfa.Ls = 1;
  dwfa.mass = 0.0;
  dwfa.Csw  = 0.0;

  printf("Initialising bfm operator\n");
  dwf->init(dwfa);

  printf("Initialising arrays\n");

  psi_h = src;
  chi_h = dwf->allocFermion();
  check = dwf->allocFermion();
  diff  = dwf->allocFermion();
  
  dwf->staticGauge(); 

  printf("Setup complete\n");

}

void compute(int me)
{
  // cb is cb of result, 1-cb is cb of input field

  double delta;
  double n1;
  double n2;
  int idx=0;

  for(int cb=0;cb<2;cb++){
    for(int dag=0;dag<2;dag++){

      if ( me == 0 ) {
        printf("Checking cb=%d dag=%d %lx \n",cb,dag,(unsigned long)arrays[idx]);
        // dwf->importFermion(arrays[idx],check,0);
      }
      check = arrays[idx];
      
      dwf->thread_barrier();

      dwf->dslash(psi_h,
		  chi_h, 
		  cb,dag);
      //      dwf->CGNE_prec(chi_h,psi_h);
      uint64_t csum = this->checksum(chi_h);

      if(me == 0) printf("Result checksum is %lx\n",csum);

      delta = dwf->axpy_norm(diff,check,chi_h,-1.0);
      
      if ( delta > 1.0e-8 ) {
        printf("Norm of difference is bad \n",delta);
	test_exit(-1);
      } else {
	if ( me == 0 ) printf("Norm of difference is good \n",delta);
      }
      idx++;
    }
  }

}







