#include <chroma.h>
#include <bfm.h>
#include <sys/time.h>

using namespace Chroma;
typedef multi1d<LatticeColorMatrix> U;
#define Printf if ( QMP_is_primary_node() ) printf
//#define Printf printf

#include <bfmcommqmpbagel2.h>
#include <bfmcommspibagel2.h>
typedef commSPIbagel2<double> bfm_spi;

#undef MPI_TEST
#ifdef MPI_TEST
void mpi_test(bfm & dwf);
#endif
#define Float double
int main (int argc,char **argv )
{
  Chroma::initialize(&argc,&argv);

  /********************************************************
   * Command line parsing
   ********************************************************
   */
  if ( argc != 6 ) { 
   Printf("Usage: %s lx ly lz lt Ls\n All must be even\n",argv[0]);
   Printf("argc is %d\n",argc);
   for ( int i=0;i<argc;i++){
      Printf("%d %s\n",i,argv[i]);
   }
      exit(-1);
  }

  /********************************************************
   * Setup QDP
   ********************************************************
   */
  multi1d<int> nrow(Nd);
  nrow[0] = atoi(argv[1]);
  nrow[1] = atoi(argv[2]);
  nrow[2] = atoi(argv[3]);
  nrow[3] = atoi(argv[4]);

  Layout::setLattSize(nrow);
  Layout::create();

  int lx = QDP::Layout::subgridLattSize()[0];
  int ly = QDP::Layout::subgridLattSize()[1];
  int lz = QDP::Layout::subgridLattSize()[2];
  int lt = QDP::Layout::subgridLattSize()[3];
  int Ls = atoi(argv[5]);

  /********************************************************
   * Setup DWF operator
   ********************************************************
   */
  bfmarg dwfa;
  bfm_spi dwf;
  dwfa.solver = DWF;
  dwfa.node_latt[0]  = lx;
  dwfa.node_latt[1]  = ly;
  dwfa.node_latt[2]  = lz;
  dwfa.node_latt[3]  = lt;

  multi1d<int> procs = QDP::Layout::logicalSize();
  Printf("%d dim machine\n\t", procs.size());
  for(int mu=0;mu<4;mu++){
    if (procs[mu]>1) dwfa.local_comm[mu] = 0;
    else             dwfa.local_comm[mu] = 1;
    Printf("%d ", procs[mu]);
  }

  Real M5(1.8);
  Real mq(0.1);

  dwfa.Ls   = Ls;
  dwfa.M5   = toDouble(M5);
  dwfa.mass = toDouble(mq);
  dwfa.Csw  = 0.0;
  dwfa.precon_5d  = 1;

  Printf("Initialising bfm operator\n");
  dwf.init(dwfa);

  /********************************************************
   * Gaussian gauge field
   ********************************************************
   */

  multi1d<LatticeColorMatrix> u(Nd); HotSt(u);

  /********************************************************
   * Gaussian source and result vectors
   ********************************************************
   */
  multi1d<LatticeFermion> psi(Ls);
  multi1d<LatticeFermion> chi(Ls);
  multi1d<LatticeFermion> chi_qdp(Ls);

  for(int s=0;s<Ls;s++) gaussian(psi[s]);


  /********************************************************
   * Bagel internal single checkerboard vectors
   ********************************************************
   */
  Fermion_t psi_h = dwf.allocFermion();
  Fermion_t chi_h = dwf.allocFermion();

  /********************************************************
   * Import gauge field to BAGEL
   ********************************************************
   */
  dwf.importGauge(u);

  LatticeFermion noise; gaussian(noise);


  //Calculate bytes sent and received
  int allbound = 0;
  int maxbound = 0;

  for (int mu=0; mu<4;mu++) if ( !dwf.local_comm[mu] ) allbound += dwf.simd_nbound[mu];
  for (int mu=0; mu<4;mu++) if ( (!dwf.local_comm[mu]) && dwf.simd_nbound[mu]>maxbound ) maxbound = dwf.simd_nbound[mu];
  
  double XBytes = 2*allbound * dwf.simd() * 12* sizeof(Float)*dwf.cbLs;
  double MaxFaceBytes = maxbound * dwf.simd() * 12* sizeof(Float)*dwf.cbLs;
  Printf("off node comms send %lf bytes\n",XBytes);
  Printf("off node comms receive %lf bytes\n",XBytes);
  for(int mu=0;mu<4;mu++) { 
    if ( !dwf.local_comm[mu] ) 
      Printf("comms[%d] %d sites, %d bytes\n",mu,dwf.simd_nbound[mu],dwf.simd_nbound[mu]*12*sizeof(Float)*dwf.cbLs*dwf.simd());
  }
  Printf("Biggest face is %lf bytes\n",MaxFaceBytes);

  for(int s=0;s<Ls;s++) chi[s] = zero;
  for(int s=0;s<Ls;s++) chi_qdp[s] = chi[s];

#if 0
  int dag = 0;
  {
    /*Import this checkerboard of source field to bagel*/
    Printf("************************************\n");
    Printf(" Full communications \n");
    Printf("************************************\n");
    QMP_barrier();
    struct timeval start,stop;
    gettimeofday(&start,NULL);
    QMP_barrier();
#define NITER 2000
    for(int i=0;i<NITER;i++){
      dwf.comm(0,psi_h,dag);
    }
    QMP_barrier();
    gettimeofday(&stop,NULL);
    struct timeval diff;
    timersub(&stop,&start,&diff);
    double t = diff.tv_sec + 1.0e-6*diff.tv_usec;
    Printf("%d comms in %f s\n",NITER,t);
    Printf("%f MB/s\n",XBytes*NITER*2./t * 1.e-6);
    Printf("largest face %f MB/s\n",MaxFaceBytes*NITER/t *1.e-6);
  }
#endif

  {
    /*Import this checkerboard of source field to bagel*/
    Printf("************************************\n");
    Printf(" QMP calls only \n");
    Printf("************************************\n");
    QMP_barrier();
    struct timeval start,stop;
    gettimeofday(&start,NULL);
#define NITER 1000
    for(int i=0;i<NITER;i++){
      dwf.comm_spi_start();
      dwf.comm_spi_complete();
    }
    QMP_barrier();
    gettimeofday(&stop,NULL);
    struct timeval diff;
    timersub(&stop,&start,&diff);
    double t = diff.tv_sec + 1.0e-6*diff.tv_usec;
    Printf("%d comms in %f s\n",NITER,t);
    Printf("%f MB/s\n",XBytes*NITER*2./t * 1.e-6);
    Printf("largest face %f MB/s\n",MaxFaceBytes*NITER/t *1.e-6);
  }
  QMP_barrier();
#ifdef MPI_TEST
  mpi_test(dwf);
#endif
  QMP_barrier();
  Printf("Done\n"); 
  QMP_barrier();

}

#ifdef MPI_TEST
#include <mpi.h>
#include <rts.h>

/*
 * Some trials of MPI bandwidth as I cannot believe how bad QMP is
 */
void mpi_test(bfm & dwf)
{
  BGLPersonality pers;

  int bgl_lex, bgl_dims[4], bgl_coor[4];

  rts_get_personality(&pers, sizeof(pers));

  // Need TXYZ MPI node mapping in BGLMPI_MAPPING environment variable
  // Not easy to find this out at runtime
  bgl_dims[1] = pers.xSize;
  bgl_dims[2] = pers.ySize;
  bgl_dims[3] = pers.zSize;
  bgl_dims[0] = (pers.opFlags & BGLPERSONALITY_OPFLAGS_VIRTUALNM) ? 2 : 1;


  bgl_coor[1] = pers.xCoord;
  bgl_coor[2] = pers.yCoord;
  bgl_coor[3] = pers.zCoord;
  bgl_coor[0] = rts_get_processor_id();

  bgl_lex = bgl_coor[0] +
    bgl_coor[1] * bgl_dims[0] +
    bgl_coor[2] * bgl_dims[0] * bgl_dims[1] +
    bgl_coor[3] * bgl_dims[0] * bgl_dims[1] * bgl_dims[2] ;

  int mpi_id ;
  MPI_Comm_rank(MPI_COMM_WORLD,&mpi_id);

  if ( mpi_id != bgl_lex ) { 
    printf("MPI mapping failure : %d %d %d %d \n",bgl_coor[0],bgl_coor[1],bgl_coor[2],bgl_coor[3]);
    printf("MPI mapping failure : MPI %d vs. %d \n",mpi_id,bgl_lex);
    printf("MPI mapping failure : Use TXYZ\n");
    exit(0);
  }

  Float * send_bufs[8];
  Float * recv_bufs[8];

  MPI_Request reqs[16];
  MPI_Status  stats[16];

  int allbytes = 0;
  int nreq = 0;
  for(int mu=0;mu<4;mu++){
    for(int pm=0;pm<2;pm++){

      int nbrcoor[4] = { bgl_coor[0],bgl_coor[1],bgl_coor[2],bgl_coor[3] };

      int delta = pm*2-1;
      int rtag = 17+pm;
      int stag = 18-pm;

      nbrcoor[mu] = (bgl_coor[mu]+delta+bgl_dims[mu])%bgl_dims[mu];
	
      int nbr_id  = nbrcoor[0] +
	nbrcoor[1] * bgl_dims[0] +
	nbrcoor[2] * bgl_dims[0] * bgl_dims[1] +
	nbrcoor[3] * bgl_dims[0] * bgl_dims[1] * bgl_dims[2] ;

      int bytes = dwf.simd_nbound[mu]*12*sizeof(Float)*dwf.cbLs;

      send_bufs[2*mu+pm]   = (Float *)memalign(256,bytes);
      recv_bufs[2*mu+pm]   = (Float *)memalign(256,bytes*2);
#define RECV_STRIDED
#ifdef RECV_STRIDED   
      int block = 12*sizeof(Float);
      int blocks= dwf.simd_nbound[mu] * dwf.cbLs;
      int stride = 2*block;
      MPI_Datatype mpi_type;
      MPI_Type_vector(blocks,block,stride,MPI_BYTE,&mpi_type);
      MPI_Type_commit(&mpi_type);

      if ( MPI_Recv_init(recv_bufs[2*mu+pm],1,mpi_type,nbr_id,rtag,MPI_COMM_WORLD,&reqs[nreq++]) != MPI_SUCCESS ) {
	printf("Oops recv init\n");
      }
#else
      if ( MPI_Recv_init(recv_bufs[2*mu+pm],bytes,MPI_BYTE,nbr_id,rtag,MPI_COMM_WORLD,&reqs[nreq++]) != MPI_SUCCESS ) {
	printf("Oops recv init\n");
      }
#endif
      if ( MPI_Send_init(send_bufs[2*mu+pm],bytes,MPI_BYTE,nbr_id,stag,MPI_COMM_WORLD,&reqs[nreq++]) != MPI_SUCCESS ){
	printf("Oops send init\n");
      }
      allbytes+=2*bytes;
      Printf("MPI init mu=%d pm%d %d requests %d \n",mu,pm,nreq,nbr_id);

    }
  }

  struct timeval start,stop;
  gettimeofday(&start,NULL);
  for(int i=0;i<NITER;i++){
    MPI_Startall(nreq,reqs);
    if ( MPI_Waitall(nreq,reqs,stats) != MPI_SUCCESS ) { 
      printf("MPI barfed\n");
      exit(0);
    }
  }  
  gettimeofday(&stop,NULL);
  struct timeval diff;
  timersub(&stop,&start,&diff);
  double t = diff.tv_sec + 1.0e-6*diff.tv_usec;
  Printf("%d comms %d bytes in %f s\n",NITER,allbytes,t);
  Printf("%f MB/s\n",((double)NITER*allbytes)/t*1.e-6,t);
  fflush(stdout);
  QMP_barrier();
  QMP_abort(0);

}
#endif 




