#include <unistd.h>
#include <stdio.h>
#include <bfm.h>


//
//
class britney : public bfm {
public:

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
    free(gauge);
  }
  void randGauge(void)
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
      for(int i=0;i<18;i++) gauge[v*18+i]    = drand48();
    }
    int ndir=4;
    if ( doubleStored() ) { 
      ndir = 8;
    }
    for(int mu=0;mu<ndir;mu++) {
      importGauge(gauge, mu); 
    }
    free(gauge);
  }
  void randFermion(Fermion_t x_t){

    double *x = (double *) x_t;
    int nspinco = 12;
    double *xx;
    double nrm = 0.0;

    for(int site=0;site<simd_cbvol;site++){
      for(int s=0;s<cbLs;s++) {
     
	xx = &x[s*nspinco*nsimd*2   + site*nspinco*nsimd*2*Ls];

	for(int spinco=0;spinco<nspinco;spinco++){
	  for(int simd_el=0;simd_el<nsimd;simd_el++){
	    for(int reim=0;reim<2;reim++){
	      int idx = reim + simd_el*2 + spinco*nsimd*2;
	      xx[idx] = (drand48()-0.5)*2; //[-1,1] uniform
	    }}}
      }
    }  

  };
};

int rand_dim(void)
{
  int dim = lrand48() &0xE;
  if ( dim > 12 ) dim=12;
  if ( dim < 4 ) dim=4;
  return dim;
}
#define Printf printf

void random_test(void);

int main (int argc,char **argv )
{
  for(int i=0;i<10;i++){
    random_test();
  }
}
void random_test(void)
{

  /********************************************************
   * Setup
   ********************************************************
   */

  int lx = rand_dim();
  int ly = rand_dim();
  int lz = rand_dim();
  int lt = rand_dim()&0xC;
  int Ls = rand_dim();
  /*
  int lx = 4;
  int ly = 4;
  int lz = 4;
  int lt = 4;
  int Ls = 4;
  */
  double M5 = 1.8;
  double mass = 0.01;
  double mq = mass;


  /********************************************************
   * Setup DWF operator
   ********************************************************
   */
  bfmarg dwfa;
  britney  dwf;

  dwfa.node_latt[0]  = lx;
  dwfa.node_latt[1]  = ly;
  dwfa.node_latt[2]  = lz;
  dwfa.node_latt[3]  = lt;

  for(int mu=0;mu<4;mu++){
      dwfa.local_comm[mu] = 1;
  }
  
  dwfa.precon_5d = 0;
  dwfa.Ls   = Ls;
  dwfa.M5   = M5;
  dwfa.mass = mq;
  dwfa.Csw  = 0.0;
  dwfa.max_iter = 1000;
  dwfa.residual = 1.e-8;
  Printf("Initialising bfm operator for %dx%dx%dx%dx%d DWF\n",lx,ly,lz,lt,Ls);
  dwf.init(dwfa);

  /********************************************************
   * Bagel internal single checkerboard vectors
   ********************************************************
   */
   Fermion_t psi_h[2];
   psi_h[0] = dwf.allocFermion();
   psi_h[1] = dwf.allocFermion();

   Fermion_t chi_h[2];
   chi_h[0] = dwf.allocFermion();
   chi_h[1] = dwf.allocFermion();

  dwf.unitGauge();

  /*Import this checkerboard of source field to bagel*/
  Printf("Importing psi field cb %d\n",1);
  dwf.randFermion(psi_h[0]);
  dwf.randFermion(psi_h[1]);
 
  // Fill the other checkerboard of result with noise
  dwf.randFermion(chi_h[0]);
  dwf.randFermion(chi_h[1]);

  Printf("Checksums chi %16.16lx %16.16lx\n",dwf.checksum(chi_h[0]),dwf.checksum(chi_h[1]));
  Printf("Checksums psi %16.16lx %16.16lx\n",dwf.checksum(psi_h[0]),dwf.checksum(psi_h[1]));

#define NITER 1
  for(int i=0;i<NITER;i++) {
    dwf.CGNE(chi_h,psi_h);
  }
  Printf("Checksums chi %16.16lx %16.16lx\n",dwf.checksum(chi_h[0]),dwf.checksum(chi_h[1]));
  Printf("Checksums psi %16.16lx %16.16lx\n",dwf.checksum(psi_h[0]),dwf.checksum(psi_h[1]));
  Printf("Done\n"); 
  dwf.freeFermion(chi_h[0]);
  dwf.freeFermion(chi_h[1]);
  dwf.freeFermion(psi_h[0]);
  dwf.freeFermion(psi_h[1]);
  dwf.end();

}








