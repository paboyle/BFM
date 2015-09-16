
#include <chroma.h>
#include <bfm.h>
#include <omp.h>
using namespace Chroma;

void compare_result (LatticeFermion A, LatticeFermion B, bfm & dwf);


int main (int argc,char **argv )
{
  Chroma::initialize(&argc,&argv);


  /********************************************************
   * Command line parsing
   ********************************************************
   */
  if ( argc != 6 ) { 
   printf("Usage: %s lx ly lz lt Ls\n All must be even\n",argv[0]);
   printf("argc is %d\n",argc);
    for ( int i=0;i<argc;i++)
      printf("%d %s\n",i,argv[i]);
    exit(-1);
  }

  int lx = atoi(argv[1]);
  int ly = atoi(argv[2]);
  int lz = atoi(argv[3]);
  int lt = atoi(argv[4]);
  int Ls = atoi(argv[5]);

  QDPIO::cout << "cbvol = "<<lx*ly*lz*lt*Ls/2 <<endl;
  /********************************************************
   * Setup QDP
   ********************************************************
   */
  multi1d<int> nrow(Nd);
  nrow[0] = lx;
  nrow[1] = ly;
  nrow[2] = lz;
  nrow[3] = lt;

  Layout::setLattSize(nrow);
  Layout::create();

  /********************************************************
   * Setup DWF operator
   ********************************************************
   */
  bfmarg dwfa;
  dwfa.solver = DWF;
  bfm    dwf;

  dwfa.node_latt[0]  = lx;
  dwfa.node_latt[1]  = ly;
  dwfa.node_latt[2]  = lz;
  dwfa.node_latt[3]  = lt;

  dwfa.local_comm[0]  = 1;
  dwfa.local_comm[1]  = 1;
  dwfa.local_comm[2]  = 1;
  dwfa.local_comm[3]  = 0;

  dwfa.Ls = Ls;
  dwfa.mass = .1;
  dwfa.M5   = 1.8;
  dwfa.Csw  = 0.0;
  dwfa.precon_5d = 1;

  printf("Initialising bfm operator\n");fflush(stdout);
  dwf.init(dwfa);
  printf("Operator is initialised\n");fflush(stdout);

  /********************************************************
   * Gaussian source and result vectors
   ********************************************************
   */
  multi1d<LatticeFermion> psi(Ls);
  multi1d<LatticeFermion> psi_qdp(Ls);
  multi1d<LatticeFermion> x(Ls);
  multi1d<LatticeFermion> y(Ls);

  for(int s=0;s<Ls;s++) {
#if 0
        gaussian(x[s]);
	gaussian(y[s]);
#else
	x[s]=zero;
	y[s]=zero;
#endif
    psi[s]=zero;
  }

  printf("Allocating Bagel fermions\n");fflush(stdout);

  /********************************************************
   * Bagel internal single checkerboard vectors
   ********************************************************
   */
  Fermion_t psi_h = dwf.allocFermion();
  Fermion_t x_h = dwf.allocFermion();
  Fermion_t y_h = dwf.allocFermion();

  int cb = 0;
  printf("Importing Bagel fermions\n");fflush(stdout);

  dwf.importFermion(x,x_h,cb);
  dwf.importFermion(y,y_h,cb);
  double a = 2.0;
  double b = 3.0;
  int Nloop = 1000;
  int cbsites = (lx*ly*lz*lt*Ls)/2;
  double flops,t;
  struct timeval diff; 
  struct timeval start,stop;
  Real n2diff;
  Real n2qdp;
  Real n2check;
  double n2dwf;
  double bytes = 1.*Nloop*cbsites*24*8*3; // RRW
  double foot  = 1.*cbsites*24*8*3; // RRW
  
  printf("Starting Linear Combination tests with footprint %le bytes\n",foot);fflush(stdout);

  ////////////////////////////////////////
  // a x + i b gamma_5 x
  ////////////////////////////////////////
  dwf.axpibg5x(psi_h,x_h,a,b);
  dwf.exportFermion(psi,psi_h,cb);
  {
  n2diff = 0.0;
  int scb = cb;
  for(int s=0;s<Ls;s++) {
    Real bb(b);
    Real z(0.0);
    Gamma G5(15);
    LatticeFermion tmp = G5 * x[s];
    psi_qdp[s] = a*x[s] + cmplx(z,bb)*tmp;
    psi_qdp[s][rb[1-scb]] = zero;
    n2diff += norm2(psi[s]-psi_qdp[s],rb[scb]);
    //    printf("axpibg5x %d %le\n",s,toDouble(norm2(psi[s]-psi_qdp[s],rb[scb])));
    //    printf("axpibg5x %d %le %le\n",s,toDouble(norm2(psi[s],rb[scb])),toDouble(norm2(psi_qdp[s],rb[scb])));
    scb = 1-scb;
  }  
  }  
  printf("n2diff = %le\n",toDouble(n2diff));

  /////////////////////////////////////////
  // AXPY
  /////////////////////////////////////////
  printf("x_h %le %le\n",dwf.norm(x_h),
	 toDouble(norm2(x[0],rb[0])));

  //Benchmark
  double *pp = (double *)psi_h;
  double *xp = (double *)x_h;
  double *yp = (double *)y_h;
#define THREADS (128)
  omp_set_num_threads(THREADS);
  int work = (cbsites*24)/THREADS;
  t=0;
  for(int i=0;i<Nloop;i++){
#pragma omp parallel for 
    for(int tid=0;tid<THREADS;tid++){
      if (i==0) {
	printf("Preloading cache with %d/%d threads\n",
	     omp_get_thread_num(),
	     omp_get_max_threads());
      }
      for(int nn=tid*work;nn<(tid+1)*work;nn++){
	pp[nn]+=0;
	xp[nn]+=1;
	yp[nn]+=2;
      }
    }

    //Benchmark
    gettimeofday(&start,NULL);
    dwf.axpy(psi_h,x_h,y_h,a);
    gettimeofday(&stop,NULL);
    timersub(&stop,&start,&diff); 
    t= t+diff.tv_sec + 1.0e-6 * diff.tv_usec;
  }
  flops = 1.*Nloop * cbsites * 24*2;
  printf("axpy : %le s : %le Gflop/s : %le GB/s\n",t,flops/t*1.0e-9,bytes/t*1.0e-9);

  //Check
  dwf.exportFermion(psi,psi_h,cb);
  n2diff = 0.0;
  int scb = cb;
  for(int s=0;s<Ls;s++) {
    psi_qdp[s] = a*x[s] + y[s];
    psi_qdp[s][rb[1-scb]] = zero;
    n2diff += norm2(psi[s]-psi_qdp[s],rb[scb]);
    //    printf("%d %le\n",s,toDouble(norm2(psi[s]-psi_qdp[s],rb[scb])));
    //    printf("%d %le %le\n",s,toDouble(norm2(psi[s],rb[scb])),toDouble(norm2(psi_qdp[s],rb[scb])));
    scb = 1-scb;
  }  
  printf("n2diff = %le\n",toDouble(n2diff));

  /////////////////////////////////////////
  // AXPBY
  /////////////////////////////////////////


  gettimeofday(&start,NULL);
  for(int i=0;i<Nloop;i++){
    dwf.axpby(psi_h,x_h,y_h,a,b);
  }
  gettimeofday(&stop,NULL);
  timersub(&stop,&start,&diff); t = diff.tv_sec + 1.0e-6 * diff.tv_usec;
  flops = 1.*Nloop * cbsites * 24*3;
  printf("axpby : %le s : %le Gflop/s : %le GB/s\n",t,flops/t*1.0e-9,bytes/t*1.0e-9);

  //Check
  dwf.exportFermion(psi,psi_h,cb);
  n2diff = 0.0;
  scb = 0;
  for(int s=0;s<Ls;s++) {
    psi_qdp[s] = x[s]*a + y[s]*b;
    psi_qdp[s][rb[1-scb]] = zero;
    n2diff += norm2(psi[s]-psi_qdp[s],rb[scb]);
    scb = 1-scb;
  }  
  printf("n2diff = %le\n",toDouble(n2diff));

  /////////////////////////////////////////
  // AXPY NORM
  /////////////////////////////////////////

  // Benchmark
  gettimeofday(&start,NULL);
  for(int i=0;i<Nloop;i++){
    n2dwf = dwf.axpy_norm(psi_h,x_h,y_h,a);
  }
  gettimeofday(&stop,NULL);
  timersub(&stop,&start,&diff); 
  t = diff.tv_sec + 1.0e-6 * diff.tv_usec;
  flops = 1.*Nloop * cbsites * 24*4;
  printf("axpy_norm : %le s : %le Gflop/s : %le GB/s\n",t,flops/t*1.0e-9,bytes/t*1.0e-9);

  //Check
  dwf.exportFermion(psi,psi_h,cb);
  n2diff = 0.0;
  n2qdp  = 0.0;
  n2check  = 0.0;
  scb = 0;
  for(int s=0;s<Ls;s++) {
    psi_qdp[s] = x[s]*a + y[s];
    psi_qdp[s][rb[1-scb]] = zero;
    n2check += norm2(psi[s],rb[scb]);
    n2diff += norm2(psi[s]-psi_qdp[s],rb[scb]);
    n2qdp+=norm2(psi_qdp[s],rb[scb]);
    scb=1-scb;
  }  
  printf("n2diff = %le\n",toDouble(n2diff));
  printf("n2qdp  = %le\n",toDouble(n2qdp));
  printf("n2dwf  = %le\n",n2dwf);
  printf("n2check  = %le\n",toDouble(n2check));

  /////////////////////////////////////////
  // AXPBY NORM
  /////////////////////////////////////////

  // Benchmark
  gettimeofday(&start,NULL);
  for(int i=0;i<Nloop;i++){
    n2dwf = dwf.axpby_norm(psi_h,x_h,y_h,a,b);
  }
  gettimeofday(&stop,NULL);
  timersub(&stop,&start,&diff); t = diff.tv_sec + 1.0e-6 * diff.tv_usec;
  flops = 1.*Nloop * cbsites * 24*5;
  printf("axpby_norm : %le s : %le Gflop/s : %le GB/s\n",t,flops/t*1.0e-9,bytes/t*1.0e-9);
  
  //Check
  dwf.exportFermion(psi,psi_h,cb);
  n2diff = 0.0;
  n2qdp  = 0.0;
  n2check= 0.0;
  scb = 0;
  for(int s=0;s<Ls;s++) {
    psi_qdp[s] = x[s]*a + y[s]*b;
    psi_qdp[s][rb[1-scb]] = zero;
    n2check += norm2(psi[s],rb[scb]);
    n2diff+= norm2(psi[s]-psi_qdp[s],rb[scb]);
    n2qdp +=norm2(psi_qdp[s],rb[scb]);
    scb = 1-scb;
  }  
  printf("n2diff = %le\n",toDouble(n2diff));
  printf("n2qdp  = %le\n",toDouble(n2qdp));
  printf("n2dwf  = %le\n",n2dwf);
  printf("n2check  = %le\n",toDouble(n2check));

  /////////////////////////////////////////
  // CAXPY 
  /////////////////////////////////////////
  a=2.0;
  b=0.2;
  dwf.caxpy(psi_h,x_h,y_h,a,b);

  printf("caxpy:\n");
  // Check
  dwf.exportFermion(psi,psi_h,cb);
  n2diff = 0.0;
  n2qdp  = 0.0;
  n2check= 0.0;
  scb = 0;
  for(int s=0;s<Ls;s++) {
    Real rra(a);
    Real rrb(b);
    psi_qdp[s] = cmplx(rra,rrb)*x[s] + y[s];
    psi_qdp[s][rb[1-scb]] = zero;
    n2diff+= norm2(psi[s]-psi_qdp[s],rb[scb]);
    scb = 1-scb;
  }
  printf("n2diff = %le\n",toDouble(n2diff));


  printf("Done\n"); 
}



void compare_result (LatticeFermion A, LatticeFermion B, bfm & dwf)
{

  multi1d<int> x(4);

  int Nspinco=12;

  printf("Result vectors \n");
  for ( x[3]=0; x[3]<dwf.node_latt[3];x[3]++ ) { 
  for ( x[2]=0; x[2]<dwf.node_latt[2];x[2]++ ) { 
  for ( x[1]=0; x[1]<dwf.node_latt[1];x[1]++ ) { 
  for ( x[0]=0; x[0]<dwf.node_latt[0];x[0]++ ) { 

    int xx[4];
    xx[0]=x[0];
    xx[1]=x[1];
    xx[2]=x[2];
    xx[3]=x[3];
    printf("site %d %d %d %d\n",x[0],x[1],x[2],x[3]);
    for ( int sp=0;sp<4;sp++ ) { 
    for ( int co=0;co<3;co++ ) { 

        int spco = co+3*sp;

	int reim=0;
        Fermion ferm = peekSite(A,x); 
        ColorVector cv = peekSpin(ferm,sp);
        Complex ca = peekColor(cv,co);

        ferm = peekSite(B,x); 
        cv = peekSpin(ferm,sp);
        Complex cb = peekColor(cv,co);

        printf("%le %le\n",toDouble(real(ca)),toDouble(real(cb)));
        printf("%le %le\n",toDouble(imag(ca)),toDouble(imag(cb)));

    }}

  }}}}

}






