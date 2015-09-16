#include <unistd.h>
#include <stdio.h>
#include <bfm.h>
#include <mpi.h>

//
// Britney incarnation overrides the global sum  
// to instead broadcast result on node zero and do bit identical compare
//
// Any nodes that dribble a bit will print "Oops, I did it again!".  
//
// And, no, I prophetically made this up prior to Britney's current issues, back
// in 2004. So there.
//
#define NUM_BFM_THREADS (64)
class britney : public bfm {
public:

  void randGauge(void)
  {
  for(int dir=0;dir<8;dir++){
  int x[4];
  int Ndircoco=72; /*8 directions stored only*/
  int Ncoco=9;

  for ( x[3]=0;x[3]<this->node_latt[3];x[3]++ ) { 
  for ( x[2]=0;x[2]<this->node_latt[2];x[2]++ ) { 
  for ( x[1]=0;x[1]<this->node_latt[1];x[1]++ ) { 
  for ( x[0]=0;x[0]<this->node_latt[0];x[0]++ ) { 
     
    int bbase = dir*9;
    int bidx;
    double * bagel = this->u;
      for ( int coco=0;coco<9;coco++ ) { 
      for ( int reim=0;reim<2;reim++ ) { 

        bidx = this->bagel_idx(x,reim,coco+bbase,Ndircoco,0);
	bagel[bidx] = drand48() * 0.001;

      }}
      bidx = this->bagel_idx(x,0,0+bbase,Ndircoco,0);      bagel[bidx]+=1.0;
      bidx = this->bagel_idx(x,0,4+bbase,Ndircoco,0);      bagel[bidx]+=1.0;
      bidx = this->bagel_idx(x,0,8+bbase,Ndircoco,0);      bagel[bidx]+=1.0;
  }}}}
  }}
  void randFermion(Fermion_t x_t){

    double *x = (double *) x_t;
    int nspinco = 12;
    double *xx;
    double nrm = 0.0;

    for(int site=0;site<simd_cbvol;site++){
      for(int s=0;s<Ls;s++) {
     
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
  virtual bool isBoss()
  {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    if ( rank == 0 ) return true;
    return false;
  }

  virtual void comm_gsum(double &val) {  

    int me = this->thread_barrier();
    if ( me == 0 ) { 

      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD,&rank);

      double bcast=0.0;
      double dest;
      if ( isBoss() ) bcast=val;

      MPI_Allreduce((void *)&bcast, (void *)&dest, 1,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      if ( val != dest ) {
	printf("Mismatch between rank %d and rank 0 %le %le\n",rank,val,dest);
	exit(-1);
      }

    }
    this->thread_barrier();
  };
};

bool isBoss()
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    if ( rank == 0 ) return true;
    return false;
}

#define Printf if ( isBoss() ) printf

int main (int argc,char **argv )
{

  MPI_Init(&argc,&argv);
  
  srand48(0x12348765);
  /********************************************************
   * Setup
   ********************************************************
   */
  int lx = 8;
  int ly = 8;
  int lz = 8;
  int lt = 8;
  int Ls = 16;
  
  double M5 = 1.8;
  double mass = 0.001;
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
  
  dwfa.precon_5d = 1;
  dwfa.Ls   = Ls;
  dwfa.M5   = M5;
  dwfa.mass = mq;
  dwfa.Csw  = 0.0;
  dwfa.max_iter = 30000;
  dwfa.residual = 1.e-24;

  int threads = NUM_BFM_THREADS;
  Printf("Initialising bfm operator for %d threads\n", threads);

  SPI_PhysicalThreadMapping(1);

  bfmarg::Threads(threads);
  bfmarg::Reproduce(1);
  bfmarg::ReproduceChecksum(1);
  bfmarg::ReproduceMasterCheck(0);
  bfmarg::Verbose(1);

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

  dwf.randGauge();

  /*Import this checkerboard of source field to bagel*/
  Printf("Importing psi field cb %d\n",0);
  dwf.randFermion(psi_h[0]);
  Printf("Importing psi field cb %d\n",1);
  dwf.randFermion(psi_h[1]);
 
  Printf("Filling result with noise\n",1);
  // Fill the other checkerboard of result with noise
  dwf.randFermion(chi_h[0]);
  dwf.randFermion(chi_h[1]);

#define NITER 200

  for(int i=0;i<NITER;i++) {

#pragma omp parallel num_threads NUM_BFM_THREADS
    {
      dwf.axpy(chi_h[0],psi_h[0],psi_h[0],0.0)  ;
      dwf.axpy(chi_h[1],psi_h[0],psi_h[0],0.0)  ;
      dwf.CGNE(chi_h,psi_h);
    }

  }
  double d=1.0;

  MPI_Barrier(MPI_COMM_WORLD);

  Printf("PASSED\n"); 
  MPI_Finalize();
}








