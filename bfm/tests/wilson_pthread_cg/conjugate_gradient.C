#include <unistd.h>
#include <stdio.h>
#include <bfm.h>

#define NTHREAD (1)

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
  }
  void randFermion(Fermion_t x_t){

    int vol = node_latt[0]
            * node_latt[1]
            * node_latt[2]
            * node_latt[3];
    int words = 24 * vol;
    int bytes = sizeof(QDPdouble)*words;

    QDPdouble *psi = (QDPdouble *)malloc(bytes);

    for(int w=0;w<words;w++ ) {
      psi[w] = drand48();
    }

    importFermion(psi, x_t, 0);

  };
};

// Lazy - just make global
void *thr_cg ( void * ptr);
bfmarg   dwfa;
britney  dwf;
Fermion_t psi_h[2];
Fermion_t chi_h[2];

#define Printf printf

int main (int argc,char **argv )
{

  /********************************************************
   * Setup
   ********************************************************
   */

  int lx = 4;
  int ly = 4;
  int lz = 4;
  int lt = 4;
  int Ls = 1;
  
  double M5 = -0.01;
  double mass = 0.01;
  double mq = mass;


  /********************************************************
   * Setup DWF operator
   ********************************************************
   */
  dwfa.threads = NTHREAD;

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

  Printf("Initialising bfm operator\n");
  dwf.init(dwfa);

  /********************************************************
   * Bagel internal single checkerboard vectors
   ********************************************************
   */
   psi_h[0] = dwf.allocFermion();
   psi_h[1] = dwf.allocFermion();

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

  pthread_t threads[NTHREAD];

  printf("Creating threads\n");
  for(int t=0;t<NTHREAD;t++){
    pthread_create(&threads[t],NULL,thr_cg,NULL);
  }
  printf("Joining threads\n");
  for(int t=0;t<NTHREAD;t++){
    pthread_join(threads[t],NULL);
  }

  Printf("Done\n"); 

}

void *thr_cg ( void * ptr)
{
#define NITER 1
  printf("CG thread\n");
  for(int i=0;i<NITER;i++) {
    printf("CG thread norms %le %le \n",dwf.norm(chi_h[0]),
	                                dwf.norm(psi_h[0])); fflush(stdout);
    dwf.axpy(chi_h[0],psi_h[0],psi_h[0],0.0)  ;
    dwf.axpy(chi_h[1],psi_h[0],psi_h[0],0.0)  ;
    dwf.CGNE(chi_h,psi_h);
    printf("CG thread done\n"); fflush(stdout);
  }
  return NULL;
}

