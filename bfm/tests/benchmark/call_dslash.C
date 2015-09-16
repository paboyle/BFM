#include <unistd.h>
#include <stdio.h>
#include <chroma.h>
#include <bfm.h>
#include <bfm_qdp.h>
#include <sys/time.h>

template<class Float> void benchmark(int lx,int ly, int lz, int lt, int Ls,int local, int th);

double usecond(void);
double usecond(void)
{
  struct timeval tv;
  gettimeofday(&tv,NULL);
  double now = 1.0e6*tv.tv_sec + 1.0*tv.tv_usec;
  return now;
}


template<class Float>
class britney : public bfm_internal<Float> {
public:

  void randGauge(void)
  {
    Float gauge[18] = { 
      1,0,0,0,0,0,
      0,0,1,0,0,0,
      0,0,0,0,1,0
    };

    for(int dir=0;dir<8;dir++){
    int Ndircoco=72; /*8 directions stored*/
    int Ncoco=9;
    for (int site=0;site<this->node_latt[0]*this->node_latt[1]*this->node_latt[2]*this->node_latt[3];site++ ) { 
    
      int x[4] ;
      int s=site;
      x[0]=s%this->node_latt[0];    s=s/this->node_latt[0];
      x[1]=s%this->node_latt[1];    s=s/this->node_latt[1];
      x[2]=s%this->node_latt[2];    s=s/this->node_latt[2];
      x[3]=s%this->node_latt[3];
      
      for ( int coco=0;coco<9;coco++ ) { 
	for ( int reim=0;reim<2;reim++ ) { 

	  Float * bagel = this->u;
	  int bbase = dir*9;
	  int cidx = reim+coco*2;
	  int bidx = this->bagel_idx(x,reim,coco+bbase,Ndircoco,0);
	  bagel[bidx] = gauge[cidx];
	  

	}}
    }
    }
  }
  void randFermion(Fermion_t x_t){

    Float *x = (Float *) x_t;
    int nspinco = 12;
    Float *xx;
    double nrm = 0.0;

    for(int site=0;site<this->simd_cbvol;site++){
      for(int s=0;s<this->cbLs;s++) {
     
	xx = &x[s*nspinco*this->nsimd*2   + site*nspinco*this->nsimd*2*this->cbLs];

	for(int spinco=0;spinco<nspinco;spinco++){
	  for(int simd_el=0;simd_el<this->nsimd;simd_el++){
	    for(int reim=0;reim<2;reim++){
	      int idx = reim + simd_el*2 + spinco*this->nsimd*2;
	      xx[idx] = (drand48()-0.5)*2; //[-1,1] uniform
	    }}}
      }
    }  

  };
};


//#define Printf printf
#define Printf if ( QMP_is_primary_node() ) printf


int main (int argc,char **argv )
{
  Chroma::initialize(&argc,&argv);

  int lens[]  = { 24, 12, 8, 4 };
  int slens[] = { 16,4,2,1 };
  int nlen = 1;
  int nslen = 3;

  int local = atoi(argv[1]);
  int fsize = atoi(argv[2]);
  int nthr  = atoi(argv[3]);
#if 0
  for(int ix=0;ix<nlen;ix++){
    //  for(int iy=0;iy<nlen;iy++){
    //  for(int iz=0;iz<nlen;iz++){
    //  for(int it=0;it<nlen;it++){
    int iy = ix;
    int iz = ix;
    int it = ix;
    {{{
    for(int is=0;is<nslen;is++){
      int lx = lens[ix];
      int ly = lens[iy];
      int lz = lens[iz];
      int lt = lens[it];
      int Ls = slens[is];
#else
      {{{{{
      int lx = 24;
      int ly = 24;
      int lz = 24;
      int lt = 24;
      int Ls = 24;
#endif

      if ( fsize ) {
	benchmark<float>(lx,ly,lz,lt,Ls,local,nthr);
      } else { 
	benchmark<double>(lx,ly,lz,lt,Ls,local,nthr);
      }
	      }abort();
  }}}}

}

template<class Float> void benchmark (int lx,int ly, int lz, int lt, int Ls,int local,int nthr)
{
  /********************************************************
   * Setup
   ********************************************************
   */
  
  double M5 = 1.8;
  double mass = 0.01;
  double mq = mass;

  /********************************************************
   * Setup DWF operator
   ********************************************************
   */
  bfmarg::Threads(nthr);
  bfmarg::Verbose(0);
  bfmarg dwfa;
  britney<Float>  dwf;

  dwfa.node_latt[0]  = lx;
  dwfa.node_latt[1]  = ly;
  dwfa.node_latt[2]  = lz;
  dwfa.node_latt[3]  = lt;

  multi1d<int> procs = QDP::Layout::logicalSize();
  for(int mu=0;mu<4;mu++){
      dwfa.local_comm[mu] = local;
      if ( procs[mu] == 1 ) {
	dwfa.local_comm[mu] = 1;
      }
  }
  
  g5dParams parms;
  parms.ScaledShamirCayleyTanh(mq,M5,Ls,2.0);
  dwfa.mass = parms.mass;
  dwfa.M5   = parms.M5;
  dwfa.Csw = 0.0;
  dwfa.precon_5d=0;
  dwfa.Ls = parms.Ls;
  dwfa.solver= parms.solver;
  dwfa.zolo_lo   = parms.zolo_lo;
  dwfa.zolo_hi   = parms.zolo_hi;
  dwfa.mobius_scale = parms.mobius_scale;

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
  dwf.randFermion(psi_h[0]);
  dwf.copy(psi_h[1],psi_h[0]);
  dwf.copy(chi_h[0],psi_h[0]);
  dwf.copy(chi_h[1],psi_h[0]);
 

  
  int Nloop = 1000;

  Printf("*********************************************\n"); 
  Printf("Benchmarking Dwilson<%d> for volume (%d %d %d %d) Ls=%d : local_comm=%d\n",
	 sizeof(Float),lx,ly,lz,lt,Ls,local);
  Printf("*********************************************\n"); 
  double flops = 1344.0*lx*ly*lz*lt*Ls*0.5 * Nloop;
  double bytes = (24.0*8.0+24.0+18.0*8.0)*lx*ly*lz*lt*Ls*0.5 * Nloop*sizeof(Float);

  double t1,t2;
  int addscale=1;
#pragma omp parallel 
  {
#pragma omp for 
    for(int t=0;t<dwf.nthread;t++) {
      for(int i=0;i<100;i++) {
	dwf.dslash_generic_serial(psi_h[0],chi_h[0],chi_h[0],0,0,0,addscale,0.0,1.0);
      }
      dwf.thread_barrier();
      if (t==0) t1 = usecond();
      for(int i=0;i<Nloop;i++) {
	dwf.dslash_generic_serial(psi_h[0],chi_h[0],chi_h[0],0,0,0,addscale,0.0,1.0);
      }
    }
  }
  t2 = usecond();
  Printf("Dwilson:         %f MB/s %le bytes %f Mflop/s %le usec\n",bytes/(t2-t1),bytes/Nloop,flops/(t2-t1),t2-t1);
  
  if(0){
  // Axpy
  bytes = 100.*24.*sizeof(Float)*3.0*dwf.node_cbvol*dwf.cbLs;
  flops = 100.*24.*2.0*dwf.node_cbvol*dwf.cbLs;
#pragma omp parallel 
  {
#pragma omp for 
    for(int t=0;t<dwf.nthread;t++) {
      t1=usecond();
      dwf.thread_barrier();
      for(int i=0;i<100;i++) {
	dwf.axpy(chi_h[0],psi_h[0],psi_h[1],0.0)  ;
      }
      dwf.thread_barrier();
      t2=usecond();
    }
  }
  Printf("axpy   :         %f MB/s %le bytes %f Mflop/s %le usec\n",bytes/(t2-t1),bytes/100.0,flops/(t2-t1),t2-t1);

  bytes = 100.*24.*sizeof(Float)*3.0*dwf.node_cbvol*dwf.cbLs;
  flops = 100.*24.*4.0*dwf.node_cbvol*dwf.cbLs;
#pragma omp parallel 
  {
#pragma omp for 
    for(int t=0;t<dwf.nthread;t++) {
      t1=usecond();
      for(int i=0;i<100;i++) {
	dwf.axpy_norm(chi_h[0],psi_h[0],psi_h[1],0.0)  ;
      }
      t2=usecond();
    }
  }
  Printf("axpy_norm        %f MB/s %le bytes %f Mflop/s %le usec\n",bytes/(t2-t1),bytes/100.0,flops/(t2-t1),t2-t1);

  bytes = 100.*24.*sizeof(Float)*7.0*dwf.node_cbvol*dwf.cbLs;
  flops = 100.*24.*8.0*dwf.node_cbvol*dwf.cbLs;
#pragma omp parallel 
  {
#pragma omp for 
    for(int t=0;t<dwf.nthread;t++) {
      t1=usecond();
      for(int i=0;i<100;i++) {
	dwf.cg_psi_p(psi_h[0],psi_h[1],chi_h[0],chi_h[1],0.0,0.0)  ;
      }
      t2=usecond();
    }
  }
  Printf("cg_psi_p_update  %f MB/s %le bytes %f Mflop/s %le usec\n",bytes/(t2-t1),bytes/100.0,flops/(t2-t1),t2-t1);

  // Benchmark mobius performance
  dwf.randGauge();

  // MooeeInv
  bytes = 100.*24.*sizeof(Float)*2.0*dwf.node_cbvol*dwf.cbLs;
  flops = 100.*(24+24+48+24)*dwf.node_cbvol*(dwf.cbLs-1); //U Um D Lm L

#pragma omp parallel 
  {
#pragma omp for 
    for(int t=0;t<dwf.nthread;t++) {
      t1=usecond();
      for(int i=0;i<100;i++){
	dwf.G5D_MooeeInv(psi_h[0],psi_h[1],0);
      }
    }
  }

  t2=usecond();
  Printf("MooeeInv         %f MB/s %le bytes  %f Mflop/s %le usec\n",bytes*1.0/(t2-t1),bytes/100.0,flops*1.0/(t2-t1),t2-t1);
  fflush(stdout);
  // Mooee
  bytes = 100.*24.*sizeof(Float)*2.0*dwf.node_cbvol*dwf.cbLs;
  flops = 100.*(48.+24.)*dwf.node_cbvol*dwf.cbLs; //U Um D Lm L
  fflush(stdout);
#pragma omp parallel 
  {
#pragma omp for 
    for(int t=0;t<dwf.nthread;t++) {
      t1=usecond();
      for(int i=0;i<100;i++){
	dwf.G5D_Mooee(psi_h[0],psi_h[1],0);
      }
    }
  }
  t2=usecond();
  Printf("Mooee            %f MB/s %le bytes %le Mflop/s %le usec\n",bytes*1.0/(t2-t1),bytes/100.0,flops*1.0/(t2-t1),t2-t1);

  Printf("*********************************************\n"); 
  }
   dwf.freeFermion(psi_h[0]);
   dwf.freeFermion(psi_h[1]);
   dwf.freeFermion(chi_h[0]);
   dwf.freeFermion(chi_h[1]);

  dwf.end();

  return;

  Printf("Benchmarking CG\n");
#define NITER 1
#pragma omp parallel 
  {
#pragma omp for 
    for(int t=0;t<dwf.nthread;t++) {
      for(int i=0;i<NITER;i++) {
	dwf.axpy(chi_h[0],psi_h[0],psi_h[0],0.0)  ;
	dwf.axpy(chi_h[1],psi_h[0],psi_h[0],0.0)  ;
	dwf.CGNE(chi_h,psi_h);
      }
    }
  }



}








