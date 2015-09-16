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


extern integer shift_table_0[];
extern integer shift_table_1[];
extern integer cb_table_0[];
extern integer cb_table_1[];
extern integer face_table_0_0_0[];
extern integer face_table_0_0_1[];
extern integer face_table_0_0_2[];
extern integer face_table_0_0_3[];
extern integer face_table_0_1_0[];
extern integer face_table_0_1_1[];
extern integer face_table_0_1_2[];
extern integer face_table_0_1_3[];
extern integer face_table_1_0_0[];
extern integer face_table_1_0_1[];
extern integer face_table_1_0_2[];
extern integer face_table_1_0_3[];
extern integer face_table_1_1_0[];
extern integer face_table_1_1_1[];
extern integer face_table_1_1_2[];
extern integer face_table_1_1_3[];

integer *linked_shift_table[2] = {
  shift_table_0,
  shift_table_1
};
integer *linked_cb_table[2] = {
  cb_table_0,
  cb_table_1
};
integer *linked_cb_table[2] = {
  face_table_0_0_0,
  face_table_0_0_1,
  face_table_0_0_2,
  face_table_0_0_3,
  face_table_0_1_0,
  face_table_0_1_1,
  face_table_0_1_2,
  face_table_0_1_3,
  face_table_1_0_0,
  face_table_1_0_1,
  face_table_1_0_2,
  face_table_1_0_3,
  face_table_1_1_0,
  face_table_1_1_1,
  face_table_1_1_2,
  face_table_1_1_3
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

virtual void pointers_init(void)
{

  int mu;
  int cb, pm;
  int shift;
  int x,y,z,t;
  int lx,ly,lz,lt;
 
  int bound_index[Ncb][NMinusPlus][ND];
  int local_p;
  int local_addr[ND];
  int shift_addr[ND];

  int cbsite,shift_cbsite;
  int offset,tab,table_size;

  lx = this->simd_latt[0];
  ly = this->simd_latt[1];
  lz = this->simd_latt[2];
  lt = this->simd_latt[3];

  for ( cb = 0 ; cb<2 ; cb++) {
    table_size = NMinusPlus*(this->simd_cbvol+1)*Nmu;
    this->shift_table[cb]= linked_shift_table[cb];
  }
  for ( cb = 0 ; cb<2 ; cb++) {
    table_size = (this->simd_cbvol+8);
    this->cb_table[cb]= linked_cb_table[cb];
  }

  for(cb=0;cb<2;cb++){
    for(pm=0;pm<2;pm++){
      for(mu=0; mu<4;mu++){

	this->face_table[cb][pm][mu]= linked_face_table[cb][pm][mu];

      }
    }
  }
  
  /*
   * This bit points to the first element after the end of the
   * body for comm_offset[Minus][0] 
   */
  this->comm_offset[Minus][0]  = this->simd_cbvol; 
  /* Now we do the rest of the comms offsets packing densely.
     nbound[mu-1] sites after the previous */
  this->comm_offset[Minus][1]  = this->comm_offset[Minus][0] + this->simd_nbound[0];
  this->comm_offset[Minus][2]  = this->comm_offset[Minus][1] + this->simd_nbound[1];
  this->comm_offset[Minus][3]  = this->comm_offset[Minus][2] + this->simd_nbound[2];

  /* Now the PLUS directions. They are the same as the minus directions 
     pushed by allbound */
  this->comm_offset[Plus][ 0 ] = this->comm_offset[Minus][0] + this->simd_allbound;
  this->comm_offset[Plus][ 1 ] = this->comm_offset[Minus][1] + this->simd_allbound;
  this->comm_offset[Plus][ 2 ] = this->comm_offset[Minus][2] + this->simd_allbound;
  this->comm_offset[Plus][ 3 ] = this->comm_offset[Minus][3] + this->simd_allbound;



  return;
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
  dwfa.solver=DWFrb4d;
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
		  cb,dag,0);
      //      dwf->CGNE_prec(chi_h,psi_h);
      
      uint64_t csum = this->checksum(chi_h);
      if ( me == 0 ) printf("Result checksum is %lx\n",csum);

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







