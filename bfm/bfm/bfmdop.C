/****************************************************************************/
/* PAB Dec 2007                                                             */
/*
  // Methods implemented here
  virtual void init(bfmarg &arg) ; 
  virtual void end (void) ;

  virtual void dslash(Fermion_t chi, 
		      Fermion_t psi, 
		      int cb,
		      int dag,int dperp=1) ;


  virtual void Dperp(Fermion_t psi, 
		     Fermion_t chi, 
		     int cb,
		     int dag);

  virtual void pointers_init(void) ;
  virtual void pointers_end(void) ; 
*/
/****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h> 
#include <fcntl.h>
#include <complex>

#include "bfm.h"
#include <sys/time.h>

#include <signal.h>



void sa_action(int sig,siginfo_t *si,void * ptr)
{

         ucontext_t * uc= (ucontext_t *)ptr;
  struct sigcontext *sc = (struct sigcontext *)&uc->uc_mcontext;

  printf("Caught signal %d\n",si->si_signo);
  printf("  mem address %lx\n",si->si_addr);
  printf("         code %d\n",si->si_code);

  printf("  instruction %lx\n",sc->rip);

  if ( si->si_signo == SIGSEGV ) {
    printf("Oops... a sigsegv\n");
    fflush(stdout);
    exit(-1);
  }
#define REG(A)  printf("  %s %lx\n",#A, sc-> A);
  REG(rdi);
  REG(rsi);
  REG(rbp);
  REG(rbx);
  REG(rdx);
  REG(rax);
  REG(rcx);
  REG(rsp);
  REG(rip);


  REG(r8);
  REG(r9);
  REG(r10);
  REG(r11);
  REG(r12);
  REG(r13);
  REG(r14);
  REG(r15);
  fflush(stdout);
  return;
};


template <class Float>
void bfmbase<Float>::init(bfmarg &arg)  
{
  int words;             /* size of the spinor field on the         */
  	  	  	 /* sublattice checkerboard                 */

  int slx;                      /* x-direction size of node sublattice */
  int sly;                      /* y-direction size of node sublattice */
  int slz;                      /* z-direction size of node sublattice */
  int slt;                      /* t-direction size of node sublattice */
  int i;
  int mu;

/*--------------------------------------------------------------------------*/
/* Set sublattice direction sizes                                           */
/*--------------------------------------------------------------------------*/
  * ((bfmarg *)this) = arg;
  precision_test=0;

  /*-----------------------------------------------------------------
   * Advise how many threads will call numerical operations
   *-----------------------------------------------------------------
   */
  this->BossLog("Initialising BFM solver %s for %d threads\n",
		    this->SolverString(this->solver),this->threads); 
  thread_init(this->threads);  
  
  this->halo_depth=1;
  if ( this->solver == WilsonNN ) { 
    this->halo_depth =2;
  }

  this->BossLog("Trapping SIGSEGV interrupts\n");
  stack_t mystack;
  mystack.ss_sp = bfm_alloc(SIGSTKSZ);
  mystack.ss_size=SIGSTKSZ;
  mystack.ss_flags=0;
  sigaltstack(&mystack, NULL);

  struct sigaction sa,osa;
  sigemptyset (&sa.sa_mask);
  sa.sa_sigaction= sa_action;
  sa.sa_flags    = SA_ONSTACK|SA_SIGINFO;

  sigaction(SIGSEGV,&sa,NULL);
  sigaction(SIGTRAP,&sa,NULL);

  /*-------------------------------------------------------------
   * Introduce 5d precon ; use ncb to refer to 4d checkerboards
   *-------------------------------------------------------------
   */
  if ( this->precon_5d && (this->Ls==1) ) { 
    this->Error("bfm::init: precon_5d & Ls=1 logic bomb\n");
    exit(-1);
  }
  if ( this->precon_5d ) this->ncb=1;
  else this->ncb=2;
  if ( this->precon_5d ) this->cbLs=this->Ls/2;
  else this->cbLs=this->Ls;

  this->node_cbvol = this->node_latt[0]
                               *this->node_latt[1]
                               *this->node_latt[2]
                               *this->node_latt[3]
                               /this->ncb;

  /***************************************************/
  /*This should be selected depending on architecture*/
  this->BossLog("bfm::init: nsimd=%d integer %d short_integer %d Float %d\n",simd(),
		sizeof(integer), sizeof(short_integer), sizeof(Float));
  this->nsimd = this->simd();

  if (this->nsimd == 1){
    this->simd_grid[0] = 1; // Query. Is it possible to just move this around?
    this->simd_grid[1] = 1;
    this->simd_grid[2] = 1;
    this->simd_grid[3] = 1;
  } else if (this->nsimd == 2 ) {
    this->simd_grid[0] = 1; // Query. Is it possible to just move this around?
    this->simd_grid[1] = 1;
    this->simd_grid[2] = 1;
    this->simd_grid[3] = 2;
  } else if (this->nsimd == 4 ) {
    this->simd_grid[0] = 1; // Query. Is it possible to just move this around?
    this->simd_grid[1] = 1;
    this->simd_grid[2] = 2;
    this->simd_grid[3] = 2;
  } else if (this->nsimd == 8 ) {
    this->simd_grid[0] = 1; // Query. Is it possible to just move this around?
    this->simd_grid[1] = 2;
    this->simd_grid[2] = 2;
    this->simd_grid[3] = 2;
  } else if (this->nsimd == 16 ) {
    this->simd_grid[0] = 2; // Query. Is it possible to just move this around?
    this->simd_grid[1] = 2;
    this->simd_grid[2] = 2;
    this->simd_grid[3] = 2;
  } else {
    exit(0);
  }

  complex_i_simd = bfm_alloc(2*8*sizeof(double)*this->nsimd);
  int ns=this->nsimd;
  for(int i=0;i<ns;i++){

      complex_i_simd[ 0*ns*2 + i*2+0]=0;
      complex_i_simd[ 0*ns*2 + i*2+1]=1;

      complex_i_simd[ 1*ns*2 + i*2+0]=0;
      complex_i_simd[ 1*ns*2 + i*2+1]=1;

      complex_i_simd[ 2*ns*2 + i*2+0]=0;
      complex_i_simd[ 2*ns*2 + i*2+1]=0;

      complex_i_simd[ 3*ns*2 + i*2+0]=2.0;
      complex_i_simd[ 3*ns*2 + i*2+1]=2.0;

  }
  /***************************************************/

  this->simd_latt[0] = slx = this->node_latt[0]/this->simd_grid[0];
  this->simd_latt[1] = sly = this->node_latt[1]/this->simd_grid[1];
  this->simd_latt[2] = slz = this->node_latt[2]/this->simd_grid[2];
  this->simd_latt[3] = slt = this->node_latt[3]/this->simd_grid[3];

  this->base_parity =(this->ncoor[0]*this->node_latt[0] 
              + this->ncoor[1]*this->node_latt[1]
              + this->ncoor[2]*this->node_latt[2]
              + this->ncoor[3]*this->node_latt[3])&0x1;

/*--------------------------------------------------------------------------*/
/* Set periodic wrap back or not                                            */
/*--------------------------------------------------------------------------*/
  this->local_comm[0] = arg.local_comm[0];
  this->local_comm[1] = arg.local_comm[1];
  this->local_comm[2] = arg.local_comm[2];
  this->local_comm[3] = arg.local_comm[3];

/*-----------------------------------------------------------------------*/
/* compute the subgrd volume of each chkbd ... at least two local dims   */
/* must be even for this code to be correct.                             */
/*-----------------------------------------------------------------------*/
  this->simd_cbvol = (slx * sly * slz * slt)/this->ncb;
  
  this->simd_nbound[0] = this->halo_depth*(sly * slz * slt)/this->ncb; 
  this->simd_nbound[1] = this->halo_depth*(slx * slz * slt)/this->ncb;
  this->simd_nbound[2] = this->halo_depth*(slx * sly * slt)/this->ncb;
  this->simd_nbound[3] = this->halo_depth*(slx * sly * slz)/this->ncb;

  this->simd_allbound  = this->simd_nbound[0]
                       + this->simd_nbound[1]
                       + this->simd_nbound[2]
                       + this->simd_nbound[3];

  if ( this->simd_nbound[0] * slx * this->ncb != (slx*sly*slz*slt)*this->halo_depth ) {
    this->Error("bfm::init Even x logic bomb\n");
    exit(-1);
  }
  if ( this->simd_nbound[1] * sly * this->ncb != (slx*sly*slz*slt)*this->halo_depth ) {
    this->Error("bfm::init Even y logic bomb\n");
    exit(-1);
  }
  if ( this->simd_nbound[2] * slz * this->ncb != (slx*sly*slz*slt)*this->halo_depth ) {
    this->Error("bfm::init Even z logic bomb\n");
    exit(-1);
  }
  if ( this->simd_nbound[3] * slt * this->ncb != (slx*sly*slz*slt)*this->halo_depth ) {
    this->Error("bfm::init Even t logic bomb\n");
    exit(-1);
  }

  /*------------------------------------------------------------------------*/
  /* Check shape                                                            */
  /*------------------------------------------------------------------------*/
  if ( (slx&1) || (slz&0x1) || (slt&0x1) ) {
    this->Error("bfm<Float>::init Bagel is refusing to run as a sub latt is odd\n");
    exit(-1);
  }



/*--------------------------------------------------------------------------*/
/* Reserve memory for gauge field                                           */
/*--------------------------------------------------------------------------*/
  // double stored, two checkerboards
  words = GAUGE_SIZE * this->nsimd * this->simd_cbvol *2*2;
  
  this->BossDebug("bfm::u allocate from heap\n"); 

  this->u = (Float *)bfm_alloc(words*sizeof(Float));
  if(this->u == 0){
    this->Error("bfm::u allocate\n");
    exit(-1);
  }

/*--------------------------------------------------------------------------*/
/* Reserve memory for Clover terms                                          */
/*--------------------------------------------------------------------------*/
    if(this->solver == CloverFermion) {
      words = CLOVER_MAT_SIZE * this->nsimd * this->simd_cbvol *2 *2;
 
      this->BossDebug("bfm::A & bfm::Ainv allocate from heap : %d\n", words);
      this->A = (Float *)bfm_alloc(words*sizeof(Float));
      if(this->A == 0) {
	this->Error("bfm::A allocate\n");
	exit(-1);
      }
      this->Ainv = (Float *)bfm_alloc(words*sizeof(Float));
      if(this->Ainv == 0) {
	this->Error("bfm::Ainv allocate\n");
	exit(-1);
      }
    }

/*----------------------------------------------------------------------*/
/* Build the pointer table                                              */
/*----------------------------------------------------------------------*/
  this->pointers_init();
  this->lebesgue_order();

/*----------------------------------------------------------------------*/
/* Initialise the comms                                                 */
/*----------------------------------------------------------------------*/
  this->comm_init();
  this->BossMessage("comms initialised\n");
  // Possible Overlap implementation init
  this->GeneralisedFiveDimInit();
  this->BossLog("Initialisation complete\n");
}

template <class Float>
void bfmbase<Float>::end (void)
{
  this->comm_end();
  this->pointers_end();
  this->GeneralisedFiveDimEnd();
  bfm_free(this->u);
  bfm_free(complex_i_simd);
}


/*
 * Implement all the platform specific routines for BG/L,P,Q
 * These are:
 *  bagel_idx -- for mapping the layout in import/export routines
 *  dslash    -- calls the BG assembler
 *  Mooee     -- preconditioning is aware of the data layout 
 *  MooeeInv  -- preconditioning is aware of the data layout 
 *  comm_init,comm_end,comm
 */

static double two_mass_simd1[4] __attribute__ ((aligned(32)));
static double two_mass_simd2[8] __attribute__ ((aligned(32)));


/*
 * Temporarily revert experiment while 
 * I figure out twisted
template <class Float>
double bfmbase<Float>::Mprec(Fermion_t psi, 
			       Fermion_t chi,
			       Fermion_t tmp, 
			       int dag,int do_nrm)
{
    int dperp_yes=1;
    int axpy_yes=1;
    int axpy_no=0;

    double nrm;
    double     ident = (5.-this->M5);
    double       hop = -0.25/(ident);
    double ident_hop = ident/hop;
    // tmp= d psi
    this->dslash_generic(psi,tmp,tmp,Even,dag,dperp_yes,axpy_no ,0.0);
    // chi= d d psi + ident * psi
    nrm = this->dslash_generic(tmp,psi,chi,Odd ,dag,dperp_yes,axpy_yes,ident_hop,hop);
    return nrm;
}
*/

template <class Float>
void bfmbase<Float>::G5R(Fermion_t input, Fermion_t &result)
{

  int x[4];
  int Nspinco = 12;
  int cb0, cb1, site, bidx, rb_idx;
  double sgn;


#if 0
  // PAB
  // Very inefficient
  for(int s=0;s<Ls;s++){
    for ( x[3]=0; x[3]<node_latt[3];x[3]++ ) { 
      for ( x[2]=0; x[2]<node_latt[2];x[2]++ ) { 
        for ( x[1]=0; x[1]<node_latt[1];x[1]++ ) { 
          for ( x[0]=0; x[0]<node_latt[0];x[0]++ ) { 

            site = x[0]+x[1]+x[2]+x[3];   
            
	    cb0 = ((site)&0x1);
	    cb1 = ((site)&0x1);

	    for ( int co=0;co<Nspinco;co++ ) { 
	      for ( int reim=0;reim<2;reim++ ) { 

		bidx = bagel_idx5d(x,Ls-1-s,reim,co,Nspinco,1);
		rb_idx = bagel_idx5d(x,s,reim,co,Nspinco,1);  
		sgn = 1.0;
		if(co > 5) sgn = -1.0;

		Float * forward = (Float *)input;
		Float * backward = (Float *)result;
		backward[rb_idx] = sgn * forward[bidx]; 
                  
	      }}
          }}}}
  }
#else 
  for(int s=0;s<Ls;s++){
    int ms = Ls-1-s;
    this->axpbg5y_ssp(result,0.0,input,1.0,input,s,ms);
  }
#endif


}

template<class Float>
void bfmbase<Float>::DslashDeriv(Fermion_t X, 
				 Fermion_t Y, 
				 Matrix_t  Force,
				 int cb,
				 int dag)
{
  /*-----------------------------------------------------
   *Divide up work between threads
   *-----------------------------------------------------
   */
  int me, thrvol,throff;
  
  thread_work(this->simd_cbvol,me,thrvol,throff);

  /*Turn global checkerboard into the checkerboard within this node*/
  int lcb = (cb + this->base_parity) & 1;
  cb = lcb;


  /*Comms & D5 */
  this->comm_start(lcb,Y,dag);
  this->comm_complete(lcb,Y);

  integer args[10];
  integer arg_p = (integer)args; 
  args[0] = ((integer)X) + throff
                         * this->cbLs
                         * this->nsimd
                         * 24*sizeof(Float);
  args[1] = (integer)Y;
  Float *Ftmp = (Float *)Force;
  args[2] = (integer)(Ftmp + throff*18*4*this->simd()); 
  args[3] = (integer)thrvol;
  args[4] = (integer)this->cbLs;
  args[5] = (integer) & this->shift_table[cb][throff*16];
  args[6] = (integer)complex_i_simd;
  args[7] = (integer) & this->lebesgue_sites[0];
  args[8] = (integer)recvbufs[0];

  int single =  1;
  if ( sizeof(Float) == sizeof(double) ) single = 0;

  this->thread_barrier();

  if ( dag ) {
    if (single)   vmx_deriv_dag_s(arg_p);
    else          vmx_deriv_dag  (arg_p);
  } else {
    if ( single ) vmx_deriv_s(arg_p);
    else          vmx_deriv  (arg_p);
  }

  this->thread_barrier();

}


template <class Float>
void bfmbase<Float>::dslash(Fermion_t psi, // Input
			    Fermion_t chi, // Result
			    int cb,
			    int dag)

{
  double z=0.0;
  double one=1.0;
  int dperp;

  if ( this->solver == DWF ) dperp = 1;
  else dperp=0;

  if ( this->precision_test ) {
    if ( sizeof(Float) != 4 ) { 
      printf("Bad assumption in sloppy comms word size\n"); exit(-1);
    }
    if ( dperp ) { 
      printf("Bad assumption in sloppy comms dperp\n"); exit(-1);
    }
    dslash_generic_sloppy(psi,chi,chi,cb,dag,0,z,one);
  } else { 
    //dslash_generic(psi,chi,chi,cb,dag,dperp,0,z,one);
    dslash_generic_serial(psi,chi,chi,cb,dag,dperp,0,z,one);
  }
}


template <class Float>
double bfmbase<Float>::dslash_generic(Fermion_t psi, // Input
              	    Fermion_t chi_in,  // Axpy
              	    Fermion_t chi_out, // Result
                    int cb,
                    int dag,
		    int dperp,
		    // if true, chi =b ( a chi + dslash psi )
    	            int addscale ,
  	            double a,
		    double b)
{
  if ( addscale ) exit(-1);
  /*-----------------------------------------------------
   *Divide up work between threads
   *-----------------------------------------------------
   */
  int me, thrvol,throff;
  
  thread_work(this->simd_cbvol,me,thrvol,throff);

  /*Turn global checkerboard into the checkerboard within this node*/
  int lcb = (cb + this->base_parity) & 1;
  cb = lcb;

  uint64_t t1,t2;
  uint64_t t_comm_complete=0;
  uint64_t t_Dslashext    =0;
  uint64_t t_Dslash       =0;

  /*Comms & D5 */
  t1 = GetTimeBase();
  this->comm_start(lcb,psi,dag);
  t2 = GetTimeBase();
  uint64_t t_comm_start=t2-t1;

  //////////////////////////////////////////////////////////
  // Have to put these on stack for dot product even though it is ugly
  // CPS style DWF implementation
  // possible B(AX+Y) operation
  /////////////////////////////////////////////////////////
  int ns = this->nsimd;

  //  double stack_i_simd[ns*8*2] __attribute__ ((aligned(32)));
  double *stack_i_simd = complex_i_simd;

  for(int n=0;n<ns*2*4;n++) {
    stack_i_simd[n] = complex_i_simd[n];
  }
  for(int n=0;n<ns*2;n++) {
    stack_i_simd[ns*2*4+n] = -2*this->mass;
    stack_i_simd[ns*2*5+n] = a;
    stack_i_simd[ns*2*6+n] = a;
  }
  
  /*
   * Setup gauge fields for 5d/4d DWF or Wilson
   */
  Float *gauge_par;
  Float *utmp = (Float *)this->u;
  if ( cb == 1 ) { 
    gauge_par = utmp+ (2*GAUGE_SIZE*this->simd_cbvol*this->simd()); 
  } else { 
    gauge_par = utmp;
  }
  if ( this->precon_5d ) { 
    gauge_par = utmp;
  } 

  bgq_l1p_optimisation(1);

  integer args[10];
  integer arg_p = (integer)args; 
  args[0] = ((integer)chi_in) + throff
                         * this->cbLs
                         * this->nsimd
                         * 24*sizeof(Float);
  args[1] = (integer)psi;
  args[2] = (integer)(gauge_par + throff*18*8*this->simd());
  args[3] = (integer)thrvol;
  args[4] = (integer)this->cbLs;
  args[5] = (integer) & this->shift_table[cb][throff*16];
  args[6] = (integer)stack_i_simd;
  //args[6] = (integer)complex_i_simd;
  args[7] = (integer) & this->lebesgue_sites[0];
   // ((integer)this->two_spinor) + me 
   //                                      * sizeof(Float)
   //                                      * this->cbLs 
   //                                      * HALF_SPINOR_SIZE * this->nsimd * 10 ;

  args[8] = (integer)recvbufs[0];
  args[9] = ((integer)chi_out) + throff
                         * this->cbLs
                         * this->nsimd
                         * 24*sizeof(Float);

  int single =  1;
  if ( sizeof(Float) == sizeof(double) ) single = 0;


  printf("Calling interior mult\n");fflush(stdout);
  this->thread_barrier();

  t1 = GetTimeBase();
  if ( thrvol ) {
    if ( dperp ) {
      if ( addscale ) {
	if ( dag ) {
	  if (single)   vmx_dwf_dp_a_int_dag_s(arg_p);
	  else          vmx_dwf_dp_a_int_dag(arg_p);
	} else {
	  if ( single ) vmx_dwf_dp_a_int_s(arg_p);
	  else          vmx_dwf_dp_a_int(arg_p);
	}
      } else { 
	if ( dag ) {
	  if (single)   vmx_dwf_dp_int_dag_s(arg_p);
	  else          vmx_dwf_dp_int_dag(arg_p);
	} else {
	  if ( single ) vmx_dwf_dp_int_s(arg_p);
	  else          vmx_dwf_dp_int(arg_p);
	}
      }
    } else {
      if ( addscale ) {
	if ( dag ) {
	  if (single)   vmx_wil_a_int_dag_s(arg_p);
	  else          vmx_wil_a_int_dag(arg_p);
	} else {
	  if ( single ) vmx_wil_a_int_s(arg_p);
	  else          vmx_wil_a_int(arg_p);
	}
      } else { 
	if ( dag ) {
	  if (single)   vmx_wil_int_dag_s(arg_p);
	  else          vmx_wil_int_dag(arg_p);
	} else {
	  if ( single ) vmx_wil_int_s(arg_p);
	  else          vmx_wil_int(arg_p);
	}
      }
    }
  }
  this->thread_barrier();
  t2 = GetTimeBase();
  t_Dslash = t2-t1;

  t1 = GetTimeBase();
  this->comm_complete(lcb,psi);
  t2 = GetTimeBase();
  t_comm_complete=t2-t1;
  
  args[0] = args[9];
  t1 = GetTimeBase();

  printf("Calling exterior mult\n");fflush(stdout);
  if ( addscale ) {
    if ( dag ) {
      if (single)   vmx_dwf_ext_scale_dag_s(arg_p);
      else          vmx_dwf_ext_scale_dag(arg_p);
    } else {
      if ( single ) vmx_dwf_ext_scale_s(arg_p);
      else          vmx_dwf_ext_scale(arg_p);
    }
  } else { 
    if ( dag ) {
      if (single)   vmx_dwf_ext_dag_s(arg_p);
      else          vmx_dwf_ext_dag(arg_p);
    } else {
      if ( single ) vmx_dwf_ext_s(arg_p);
      else          vmx_dwf_ext(arg_p);
    }
  }

  t2 = GetTimeBase();
  t_Dslashext = t2-t1;

  double dot=0;
  for(int i=0;i<ns*2;i++){
    dot = dot+stack_i_simd[7*2*ns+i];
  }
  this->thread_sum(dot,me);
  this->comm_gsum(dot);

  // This is painful... all the performance analysis doubles the code length
  //   if ( this->iter == this->time_report_iter ) {
  static int printed;
  if ( !printed  ) {
    if ( this->isBoss() && !me ) { 
    
      printf("comm_start    : %ld cycles\n",t_comm_start);
      printf("comm_complete : %ld cycles\n",t_comm_complete);
      printf("dslash_int    : %ld cycles\n",t_Dslash);
      printf("dslash_ext    : %ld cycles\n",t_Dslashext);
    }
    if ( !me )  printed =1;

  }
     
  if ( this->iter == this->time_report_iter ) {
    if ( this->isBoss() && !me ) { 
	 int allbound = 0;
	 int maxbound = 0;
	 
	 for (int mu=0; mu<4;mu++)
	   if ( !this->local_comm[mu] ) allbound += this->simd_nbound[mu];
	 for (int mu=0; mu<4;mu++) 
	   if ( (!this->local_comm[mu]) && this->simd_nbound[mu]>maxbound ) 
	     maxbound = this->simd_nbound[mu];
	 
	 double XBytes = 2*allbound * this->simd() * 12* sizeof(Float)*this->cbLs;
	 int simd_gm  = this->simd_cbvol*8;
	 int simd_egm = 2*allbound;
	 int simd_igm = simd_gm-simd_egm;

	 this->ThreadBossPerformance("simd_gauge_muls     = %d\n",simd_gm);
	 this->ThreadBossPerformance("simd_int_gauge_muls = %d\n",simd_igm);
	 this->ThreadBossPerformance("simd_ext_gauge_muls = %d\n",simd_egm);
	 double int_flops = (12+132)*simd_igm+7*24*this->simd_cbvol;
	 double ext_flops = (12+132)*simd_egm;
	 double flops = 1344.*this->simd_cbvol*this->simd()*this->cbLs;
	 
	 this->ThreadBossPerformance("int_flops     = %le\n",int_flops);
	 this->ThreadBossPerformance("ext_flops     = %le\n",ext_flops);
	 this->ThreadBossPerformance("flops         = %le\n",flops);

	 this->ThreadBossPerformance("t_comm_start    %lld\n", t_comm_start);
	 this->ThreadBossPerformance("t_comm_complete %lld\n", t_comm_complete);
	 this->ThreadBossPerformance("comm bytes %le completed in %le usecs\n", XBytes*2.0 ,1600./(t_Dslash+t_comm_start+t_comm_complete+t_Dslashext));

	 this->ThreadBossPerformance("t_Dslash        %lld\n", t_Dslash);
	 this->ThreadBossPerformance("t_Dslashext     %lld\n", t_Dslashext);
	
	 this->ThreadBossPerformance("Pass-1+2   Mflop/s %16.3le\n", flops*1600./(t_Dslash+t_Dslashext));
	 this->ThreadBossPerformance("Pass-1     Mflop/s %16.3le\n", int_flops*1600./(t_Dslash));
	 this->ThreadBossPerformance("Pass-2     Mflop/s %16.3le\n", ext_flops*1600./(t_Dslashext));

	 this->ThreadBossPerformance("Aggregate        Mflop/s %16.3le\n", flops*1600./(t_Dslash+t_comm_start+t_comm_complete+t_Dslashext));

	 fflush(stdout);
    }
  }
   this->thread_barrier();
   bgq_l1p_optimisation(0);

  return dot;

}
 
template <class Float>
double bfmbase<Float>::dslash_generic_serial(Fermion_t psi, // Input
              	    Fermion_t chi_in,  // Axpy
              	    Fermion_t chi_out, // Result
                    int cb,
                    int dag,
		    int dperp,
		    // if true, chi =b ( a chi + dslash psi )
    	            int addscale ,
  	            double a,
		    double b)
{

  /*-----------------------------------------------------
   *Divide up work between threads
   *-----------------------------------------------------
   */
  int me, thrvol,throff;

  //  if ( addscale ) exit(-1);
  //  this->ThreadBossLog("Serial dslash addscale %d",addscale);
  
  thread_work(this->simd_cbvol,me,thrvol,throff);

  /*Turn global checkerboard into the checkerboard within this node*/
  int lcb = (cb + this->base_parity) & 1;
  cb = lcb;

  /*Comms & D5 */
  uint64_t t0 = GetTimeBase();
  if ( ! addscale) this->comm_start(lcb,psi,dag);
  uint64_t t_comm_start = GetTimeBase()-t0;

  t0 = GetTimeBase();
  if ( !addscale) this->comm_complete(lcb,psi);
  uint64_t t_comm_complete = GetTimeBase()-t0;
  addscale = 0;


  //////////////////////////////////////////////////////////
  // Have to put these on stack for dot product even though it is ugly
  //
  // CPS style DWF implementation
  // possible B(AX+Y) operation
  /////////////////////////////////////////////////////////


  int ns = this->nsimd;

  double *stack_i_simd = complex_i_simd;

  //  double stack_i_simd[ns*8*2] __attribute__ ((aligned(32)));
  for(int n=0;n<ns*2*4;n++) {
    stack_i_simd[n] = complex_i_simd[n];
  }
  for(int n=0;n<ns*2;n++) {
    stack_i_simd[ns*2*4+n] = -2*this->mass;
    stack_i_simd[ns*2*5+n] = a;
    stack_i_simd[ns*2*6+n] = a;
  }

  
  /*
   * Setup gauge fields for 5d/4d DWF or Wilson
   */
  Float *gauge_par;
  Float *utmp = (Float *)this->u;
  if ( cb == 1 ) { 
    gauge_par = utmp+ (2*GAUGE_SIZE*this->simd_cbvol*this->simd()); 
  } else { 
    gauge_par = utmp;
  }
  if ( this->precon_5d ) { 
    gauge_par = utmp;
  } 

  bgq_l1p_optimisation(1);

  integer args[10];
  integer arg_p = (integer)args; 
  args[0] = (integer)chi_in;
  args[1] = (integer)psi;
  args[2] = (integer)gauge_par;
  args[3] = (integer)thrvol;
  args[4] = (integer)this->cbLs;
  args[5] = (integer)this->shift_table[cb];
  args[6] = (integer)stack_i_simd;
  args[7] = (integer)&this->lebesgue_sites[throff];
  args[8] = (integer)recvbufs[0];
  args[9] = (integer)chi_out;

  int single =  1;
  if ( sizeof(Float) == sizeof(double) ) single = 0;

  this->thread_barrier();
  t0 = GetTimeBase();

  if ( thrvol ) {
    if ( dperp ) {
      if ( addscale ) {
	if ( dag ) {
	  if (single)   vmx_dwf_dp_as_dag_s(arg_p);
	  else          vmx_dwf_dp_as_dag(arg_p);
	} else {
	  if ( single ) vmx_dwf_dp_as_s(arg_p);
	  else          vmx_dwf_dp_as(arg_p);
	}
      } else { 
	if ( dag ) {
	  if (single)   vmx_dwf_dp_dag_s(arg_p);
	  else          vmx_dwf_dp_dag(arg_p);
	} else {
	  if ( single ) vmx_dwf_dp_s(arg_p);
	  else          vmx_dwf_dp(arg_p);
	}
      }
    } else {
      if ( addscale ) {
	if ( dag ) {
	  if (single)   vmx_wil_as_dag_s(arg_p);
	  else          vmx_wil_as_dag(arg_p);
	} else {
	  if ( single ) vmx_wil_as_s(arg_p);
	  else          vmx_wil_as(arg_p);
	}
      } else { 
	if ( dag ) {
	  if (single)   vmx_wil_dag_s(arg_p);
	  else          vmx_wil_dag(arg_p);
	} else {
	  if ( single ) vmx_wil_s(arg_p);
	  else          vmx_wil(arg_p);
	}
      }
    }
  }
  this->thread_barrier();

  uint64_t t_Dslash = GetTimeBase()-t0;

  double dot=0;
  for(int i=0;i<ns*2;i++){
    dot = dot+stack_i_simd[7*2*ns+i];
  }
  this->thread_sum(dot,me);
  this->comm_gsum(dot);
  static int printed;
  if ( !printed  ) {
    if ( this->isBoss() && !me ) { 
      printf("comm_start    : %ld cycles\n",t_comm_start);
      printf("comm_complete : %ld cycles\n",t_comm_complete);
      printf("dslash_int    : %ld cycles\n",t_Dslash);
    }
    if ( !me )  printed =1;

  }


  bgq_l1p_optimisation(0);

  return dot;
}

extern "C" { 
void vmx_compress2b(int N,float  *in, float *out);
void vmx_decompress2b(int N,float  *in, float *out);
}
template <class Float>
double bfmbase<Float>::dslash_generic_sloppy(Fermion_t psi, // Input
					     Fermion_t chi_in,  // Axpy
					     Fermion_t chi_out, // Result
					     int cb,
					     int dag,
					     // if true, chi =b ( a chi + dslash psi )
					     int addscale ,
					     double a,
					     double b)
{
#ifdef BGQ
  /*-----------------------------------------------------
   *Divide up work between threads
   *-----------------------------------------------------
   */
  static int printed=0;

  if ( addscale ) exit(-1);

  int me, thrvol,throff;
  uint64_t t0,t1,t2,t3,t4;
  thread_work(this->simd_cbvol,me,thrvol,throff);

  /*Turn global checkerboard into the checkerboard within this node*/
  int lcb = (cb + this->base_parity) & 1;
  cb = lcb;

  /*Comms & D5 */
  t0 = GetTimeBase();
  this->comm_start(lcb,psi,dag);
  t1 = GetTimeBase();
  this->comm_complete(lcb,psi);
  t2 = GetTimeBase();

  //////////////////////////////////////////////////////////
  // Have to put these on stack for dot product even though it is ugly
  //
  // CPS style DWF implementation
  // possible B(AX+Y) operation
  /////////////////////////////////////////////////////////

  int ns=this->nsimd;
  //  double stack_i_simd[ns*8*2] __attribute__ ((aligned(32)));
  double *stack_i_simd = complex_i_simd;

  for(int n=0;n<ns*2*4;n++) {
    stack_i_simd[n] = complex_i_simd[n];
  }
  for(int n=0;n<ns*2;n++) {
    stack_i_simd[ns*2*4+n] = -2*this->mass;
    stack_i_simd[ns*2*5+n] = a;
    stack_i_simd[ns*2*6+n] = a;
  }
  
  /*
   * Setup gauge fields for 5d/4d DWF or Wilson
   */
  Float *gauge_par;
  Float *utmp = (Float *)this->u;
  if ( cb == 1 ) { 
    gauge_par = utmp+ (2*GAUGE_SIZE*this->simd_cbvol*this->simd()); 
  } else { 
    gauge_par = utmp;
  }
  if ( this->precon_5d ) { 
    gauge_par = utmp;
  } 

  bgq_l1p_optimisation(1);

  integer args[10];
  integer arg_p = (integer)args; 
  args[0] = ((integer)chi_in) + throff
                         * this->cbLs
                         * this->nsimd
                         * 24*sizeof(Float);
  args[1] = (integer)psi;
  args[2] = (integer)(gauge_par + throff*18*8*this->simd());
  args[3] = (integer)thrvol;
  args[4] = (integer)this->cbLs;
  args[5] = (integer) & this->shift_table_half[cb][throff*16];
  args[6] = (integer)stack_i_simd;
  args[7] = (integer) & this->lebesgue_sites[0];

  args[8] = (integer)recvbufs[0];
  args[9] = ((integer)chi_out) + throff
                         * this->cbLs
                         * this->nsimd
                         * 24*sizeof(Float);

  int single =  1;
  if ( sizeof(Float) == sizeof(double) ) single = 0;

  this->thread_barrier();

  if ( thrvol ) {
    if ( addscale ) {
      if ( dag ) {
	vmx_wil_as_dag_hs(arg_p);
      } else {
	vmx_wil_as_hs(arg_p);
      }
    } else { 
      if ( dag ) {
	vmx_wil_dag_hs(arg_p);
      } else {
	vmx_wil_hs(arg_p);
      }
    }
  }
  this->thread_barrier();
  t3 = GetTimeBase();
   
  double dot=0;
  for(int i=0;i<ns*2;i++){
    dot = dot+stack_i_simd[7*2*ns+i];
  }
  this->thread_sum(dot,me);
  this->comm_gsum(dot);

  if ( printed < 10 ) {
    if ( this->isBoss() && !me ) { 
      if ( printed==5 ){
	printf("comm_start    : %ld cycles\n",t1-t0);
	printf("comm_complete : %ld cycles\n",t2-t1);
	printf("dslash        : %ld cycles\n",t3-t2);
      }
      printed++
    }
  }

  bgq_l1p_optimisation(0);

  return dot;
#else 
  return 0.0;
#endif
}




#if 0
// This routine works out Dperp in DWF
template <class Float>
void bfmbase<Float>::Dperp(Fermion_t psi, 
		   Fermion_t chi, 
		   int cb,
		   int dag)
{
  this->Error("Deprecated\n");
  exit(-1);
  //Temporary debug hack 
  two_mass_simd1[0] = two_mass_simd1[1] = 2.0;
  two_mass_simd1[2] = two_mass_simd1[3] =-2.0*this->mass ;

  two_mass_simd2[0] = two_mass_simd2[1] = 2.0;
  two_mass_simd2[2] = two_mass_simd2[3] = 2.0;
  two_mass_simd2[4] = two_mass_simd2[5] =-2.0*this->mass ;
  two_mass_simd2[6] = two_mass_simd2[7] =-2.0*this->mass ;

  int me,thrvol,throff;
  thread_work(this->simd_cbvol,me,thrvol,throff);

  integer args[10];
  integer arg_p = (integer)args; 
  int fermoff =  24*throff
                   * this->nsimd
                   * this->cbLs
                   * sizeof(Float);
  args[0] = (integer) psi + fermoff;
  args[1] = (integer) chi + fermoff;
  args[2] = (integer) thrvol;
  args[3] = (integer) this->cbLs;
  args[4] = (integer) this->cb_table[cb] + throff*sizeof(integer);
  args[5] = (integer) complex_i_simd;

  if(this->simd() == 1 ) args[6] = (integer) two_mass_simd1;
  if(this->simd() == 2 ) args[6] = (integer) two_mass_simd2;

  int single = 1;
  if (sizeof(Float)==sizeof(double) ) single = 0; 

  if ( dag ) {
    if (single ) vmx_dperp_dag_s(arg_p);
    else         vmx_dperp_dag(arg_p);
  } else { 
    if (single ) vmx_dperp_s(arg_p);
    else         vmx_dperp(arg_p);
  }

}
#endif

uint32_t alignup(uint32_t n);
std::vector<short_integer> cache_reorder(const std::vector<uint32_t> &dims);
std::vector<short_integer> cache_interleave(std::vector<short_integer> &reord);

std::vector<short_integer> cache_interleave(std::vector<short_integer> &reord,int nthread)
{
  std::vector<short_integer> ret(0);
  for(int t=0;t<nthread;t++){
    for(int s=0;s<reord.size();s++){
      if ((s%nthread)==(nthread-1-t)) ret.push_back(reord[s]);
    }
  }
  return ret;
}

uint32_t alignup(uint32_t n)
{
  n--;           // 1000 0011 --> 1000 0010
  n |= n >> 1;   // 1000 0010 | 0100 0001 = 1100 0011
  n |= n >> 2;   // 1100 0011 | 0011 0000 = 1111 0011
  n |= n >> 4;   // 1111 0011 | 0000 1111 = 1111 1111
  n |= n >> 8;   // ... (At this point all bits are 1, so further bitwise-or
  n |= n >> 16;  //      operations produce no effect.)
  n++;           // 1111 1111 --> 1 0000 0000
  return n;
}
std::vector<short_integer> cache_reorder(const std::vector<int> &dims)
{
  std::vector<short_integer> ret(0);

  // Align up dimensions to power of two.
  uint32_t ND = dims.size();
  std::vector<uint32_t> adims(ND);
  std::vector<std::vector<uint32_t> > bitlist(ND);
  const uint32_t one=1;

  for(uint32_t mu=0;mu<ND;mu++){
    assert ( dims[mu] != 0 );
    adims[mu] = alignup(dims[mu]);
  }

  // List which bits of padded volume coordinate contribute; this strategy 
  // i) avoids recursion 
  // ii) has loop lengths at most the width of a 32 bit word.
  int sitebit=0;
  for(int bit=0;bit<32;bit++){
    uint32_t mask = one<<bit;
    for(int mu=0;mu<ND;mu++){
      if ( mask&(adims[mu]-1) ){
	bitlist[mu].push_back(sitebit);
	sitebit++;
      }
    }
  }

  // Work out padded and unpadded volumes
  uint32_t avol = 1;
  for(int mu=0;mu<ND;mu++) avol = avol * adims[mu];

  uint32_t vol = 1;
  for(int mu=0;mu<ND;mu++) vol = vol * dims[mu];
  
  // Loop over padded volume, following Lebesgue curve
  // We interleave the bits from sequential "mu".
  std::vector<uint32_t> ax(ND);

  for(uint32_t asite=0;asite<avol;asite++){

    // Start with zero and collect bits
    for(int mu=0;mu<ND;mu++) ax[mu] = 0;

    int contained = 1;
    for(int mu=0;mu<ND;mu++){

      // Build the coordinate on the aligned volume
      for(int bit=0;bit<bitlist[mu].size();bit++){
	int sbit=bitlist[mu][bit];

	if(asite&(one<<sbit)){
	  ax[mu]|=one<<bit;
	}
      }

      // Is it contained in original box
      if ( ax[mu]>dims[mu]-1 ) contained = 0;

    }

    if ( contained ) {
      int site = ax[0]
	+        dims[0]*ax[1]
        +dims[0]*dims[1]*ax[2]
        +dims[0]*dims[1]*dims[2]*ax[3];
      ret.push_back(site);
      //ret.push_back(0);
    }
  }

  assert( ret.size() == vol );
  return ret;
}

template <class Float>
void bfmbase<Float>::lebesgue_order(void)
{
  std::vector<int> dims(4);
  for(int mu=0;mu<4;mu++){
    dims[mu] = this->simd_latt[mu];
    if ( mu==0 ) dims[mu] = dims[mu]>>1; // drop out parity
  }

  lebesgue_sites = cache_reorder(dims);
  //  for(int i=0;i<lebesgue_sites.size();i++){
  //   lebesgue_sites[i]=i;
  //  }
  lebesgue_sites = cache_interleave(lebesgue_sites,nthread);
}

template <class Float>
void bfmbase<Float>::pointers_end(void)
{
  for ( int cb = 0;cb<Ncb;cb++ ) {
    bfm_free(this->shift_table[cb]);
    bfm_free(this->shift_table_half[cb]);
    bfm_free(this->cb_table[cb]);
  }
  for(int cb=0;cb<this->ncb;cb++) {
    for(int pm=0;pm<2;pm++) {
      for(int mu=0; mu<4;mu++) {
	bfm_free(this->face_table[cb][pm][mu]);
      }
    }
  }
}


/*Layout is:                                     */
/* PSI   [site/Nsimd][Ls][spin][co][reim][Nsimd] */
/* Gauge [site/Nsimd][mu][co][co][reim][Nsimd]   */
/* PSI contains attached face buffers at the back of array */

template <class Float>
void bfmbase<Float>::pointers_init(void)
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

  /////////////////////////////////
  //Next nearest wilson checks
  /////////////////////////////////
  if ( this->solver == WilsonNN ) {
    if ( this ->halo_depth != 2 ) {
      this->Error("Halo depth !=2 and WilsonNN\n"); 
      exit(-1);
    }
  } else if ( this ->halo_depth != 1 ) {
    this->Error("Halo depth >1 and not WilsonNN");
  }
  //////////////////////

  for ( cb = 0 ; cb<Ncb ; cb++) {
    table_size = NMinusPlus*(this->simd_cbvol+1)*2*Nmu;
    this->shift_table[cb]=(short_integer *)bfm_alloc( table_size*sizeof(short_integer));

    this->shift_table_half[cb]=(short_integer *)bfm_alloc( table_size*sizeof(short_integer));
    for(int i=0;i<table_size;i++) this->shift_table_half[cb][i]=this->shift_table[cb][i]=0;
  }
  for ( cb = 0 ; cb<Ncb ; cb++) {
    table_size = (this->simd_cbvol+8);
    this->cb_table[cb]=(integer *)bfm_alloc( table_size*sizeof(integer));
  }

  for(cb=0;cb<this->ncb;cb++){
    for(pm=0;pm<2;pm++){
      for(mu=0; mu<4;mu++){

	this->face_table[cb][pm][mu]= (integer *)
	     bfm_alloc(this->simd_allbound*sizeof(integer)*this->halo_depth);
	if ( this->ncb==1 ) this->face_table[1][pm][mu]=NULL;

      }
    }
  }
  
  for( mu = 0 ; mu < Nmu ; mu++ ) {
    bound_index[0][0][mu] = 0;
    bound_index[0][1][mu] = 0;
    bound_index[1][0][mu] = 0;
    bound_index[1][1][mu] = 0;
  }

  /*
   * This bit points to the first element after the end of the
   * body for comm_offset[Minus][0] 
   */

  this->comm_offset[Minus][0]=0;
  /* Now we do the rest of the comms offsets packing densely.
     nbound[mu-1] sites after the previous */
  this->comm_offset[Minus][1]  = this->comm_offset[Minus][0] + this->simd_nbound[0] ;
  this->comm_offset[Minus][2]  = this->comm_offset[Minus][1] + this->simd_nbound[1] ;
  this->comm_offset[Minus][3]  = this->comm_offset[Minus][2] + this->simd_nbound[2] ;

  /* Now the PLUS directions. They are the same as the minus directions 
     pushed by allbound */
  this->comm_offset[Plus][ 0 ] = this->comm_offset[Minus][0] + this->simd_allbound ;
  this->comm_offset[Plus][ 1 ] = this->comm_offset[Minus][1] + this->simd_allbound ;
  this->comm_offset[Plus][ 2 ] = this->comm_offset[Minus][2] + this->simd_allbound ;
  this->comm_offset[Plus][ 3 ] = this->comm_offset[Minus][3] + this->simd_allbound ;


  /* Now work out the shifts */
  for ( pm = 0; pm<2; pm++ ) {

    /* Shift direction */
    if ( pm == Plus ) { /*forwards or backwards*/
      shift = 1;
    } else { 
      shift = -1;
    }
    shift = shift * this->halo_depth; // Two site hopping

      /*Loop naively over local lattice*/
    for ( t=0 ; t<lt ; t++ ) {
     local_addr[3] = t;
     for ( z=0 ; z<lz ; z++ ) {
      local_addr[2] = z;
      for ( y=0 ; y<ly ; y++ ) {
       local_addr[1] = y;
       for ( x=0 ; x<lx ; x++ ) {
        local_addr[0] = x;

	/* local_addr contains the coordinates of the point */
	/* we get the parity of the point                   */
        local_p = this->parity(local_addr);
	/* get the cbsite of the point for within parity */
        cbsite   = this->psite(local_addr,this->simd_latt);


	// Used for DWF Dperp
	this->cb_table[0][cbsite] = (x+y+z+t)&0x1;
	this->cb_table[1][cbsite] = 1-this->cb_table[0][cbsite];

	// Store info in table for the dslash as want to integrate the 5d hopping term
	//
	// Extra info: 0,5 <- Link mask for connects to halo ..  queue_iload_short(mask,TAB_IMM[10],tab);
	// Extra info: 1,4 <- Base pointer    queue_iload_short(Chi_mu,TAB_IMM[9],tab);  
	// Extra info: 1,5 <- checkerboard for 5d prec Dperp inside dslash    queue_iload_short(dperp_XX,TAB_IMM[11],tab); 

	int info_idx;
	info_idx = this->interleave_site(1,5,cbsite);
	this->shift_table_half[0][info_idx] = this->shift_table[0][info_idx] =   this->cb_table[0][cbsite];
	this->shift_table_half[1][info_idx] = this->shift_table[1][info_idx] = 1-this->cb_table[0][cbsite];

	info_idx = this->interleave_site(1,4,cbsite);
	this->shift_table_half[local_p][info_idx] = this->shift_table[local_p][info_idx] = cbsite*24*this->nsimd*sizeof(Float)*this->cbLs;
	

        for ( mu = 0 ; mu < 4; mu ++) {

	  /* Get the coordinates of the point shifted in direction 
	     +/- mu. (+/- is encoded in shift, depending on where 
	     we are in the +/- loop) */

          shift_addr[0] = x;
          shift_addr[1] = y;
          shift_addr[2] = z;
          shift_addr[3] = t;

          shift_addr[mu] += shift;

	  /* Get the interleaved site address of the local point
	     (to which we are shifting data for dslash stencil*/
	  tab = this->interleave_site(pm,mu,cbsite);

          /*we're on the boundary*/
	  int on_boundary = 0;
	  if ( (shift_addr[mu] < 0) || (shift_addr[mu]>=this->simd_latt[mu] ) ) {
	    on_boundary	= 1;
          }

          /*
	   * Get the offset to the neighbour site in the spinor with periodicity
	   */
          shift_addr[mu]= (shift_addr[mu] + this->simd_latt[mu])%this->simd_latt[mu];
          shift_cbsite   = psite(shift_addr,this->simd_latt);

	  info_idx = this->interleave_site(0,5,cbsite);
          if ( on_boundary == 0 ) { 

	    //	    this->shift_table[local_p][tab] = shift_cbsite*2;
	    this->shift_table_half[local_p][tab] = this->shift_table[local_p][tab] = shift_cbsite*24*this->nsimd*sizeof(Float)*this->cbLs;
	    

	  } else {

	    /*
	     * The face_table contains a pointer to the
	     * neighbour site, to be gathered into the comms buffers.
             * local_p is the RESULT checkerboard/parity
	     * This should be sent in the OPPOSITE (1-pm) direction.
	     * Corresponding received data (same if loopback) 
	     * is stored the comm_offset[pm]
	     */
	    int bdx = bound_index[local_p][pm][mu];
	    this->face_table[local_p][pm][mu][bdx] = shift_cbsite;

	    /*
	     * And we use face_table[local_p][pm][mu] in conjunction
	     * with the send buffer[pm][mu] and source parity local_p
	     */
            if ( (!this->local_comm[mu]) || (this->simd_grid[mu]>1) ) {
	      this->shift_table[local_p][tab] = (this->comm_offset[pm][mu]*24 + bdx*12) * this->nsimd*sizeof(Float)*this->cbLs;
	      this->shift_table[local_p][info_idx] |=(1<<(pm+2*mu));

	      this->shift_table_half[local_p][tab] = 
   		this->comm_offset[pm][mu]*24*this->nsimd*sizeof(Float)*this->cbLs
	   	                  + bdx * 12*this->nsimd*sizeof(uint16_t)*this->cbLs; // Half size comms words, halfspinors
	      this->shift_table_half[local_p][info_idx] |=(1<<(pm+2*mu));
            } else { 
	      this->shift_table_half[local_p][tab] = this->shift_table[local_p][tab] = shift_cbsite*24*this->nsimd*sizeof(Float)*this->cbLs;
            }
	    /*
	     *   X X X A
	     * 
	     *   face_table[plus] = 4
	     *   shift_table[4,plus] = this->comm_offset[plus]
	     *
	     *   gather (face_table) -> buffer
             *   This should be copied to -> this->comm_offset[minus]
	     *   Reading the corresponding minus face from shift_table pulls A in the hopping term
	     */
	    bound_index[local_p][pm][mu] ++;
	  }


        } /*mu*/
       }/*x*/
      }/*y*/
     }/*z*/
    }/*t*/
  }/*pm*/


  if ( this->precon_5d ) { // Checkerboarding for 5d is different
    for(int s=0;s<this->simd_cbvol;s++){
      int idx=Nmu*NMinusPlus*2*s;
      for(int i=0;i<Nmu*NMinusPlus*2;i++){
	this->shift_table[1][idx+i]=this->shift_table[0][idx+i];
	this->shift_table_half[1][idx+i]=this->shift_table_half[0][idx+i];
      }
      this->shift_table_half[1][idx+11] = this->shift_table[1][idx+11] =1-this->shift_table[0][idx+11]; // Flip the parity
      if ( s==this->simd_cbvol-1 ) {
	this->shift_table_half[1][idx+8] = this->shift_table[1][idx+8] =0;
	
      } else {
	this->shift_table_half[1][idx+8] = this->shift_table[1][idx+8] = this->shift_table[1][idx+16];
      }
    }
  }

#if 0
      /*Loop naively over local lattice*/
  for ( t=0 ; t<lt ; t++ ) {
  for ( z=0 ; z<lz ; z++ ) {
  for ( y=0 ; y<ly ; y++ ) {
  for ( x=0 ; x<lx ; x++ ) {

    local_addr[3] = t;
    local_addr[2] = z;
    local_addr[1] = y;
    local_addr[0] = x;

    /* local_addr contains the coordinates of the point */
    /* we get the parity of the point                   */
    local_p = this->parity(local_addr);
    /* get the cbsite of the point for within parity */
    cbsite   = this->psite(local_addr,this->simd_latt);
    this->Error("%d %d %d %d Site %d\n",x,y,z,t,cbsite);
    for(int i=0;i<16;i++){
      this->Error("%d\t",this->shift_table[local_p][cbsite*16+i]);
    }
    this->Error("\n");

  }}}}
#endif


  //  int bound_index[Ncb][NMinusPlus][ND];
  for( cb = 0 ; cb < 2 ; cb++ ) {
  for( pm = 0 ; pm < 2 ; pm++ ) {
  for( mu = 0 ; mu < Nmu ; mu++ ) {
    if ( (!this->precon_5d) || cb == 0 ) {
      if ( bound_index[cb][pm][mu] != this->simd_nbound[mu] ) {
	this->Error("bfm<Float>::pointers Boundary size mismatch[%d][%d][mu=%d] : %d != %d \n",
	       cb,pm,mu,bound_index[cb][pm][mu],this->simd_nbound[mu]);
	exit(-1);
      }
    }
  }
  }
  }

  return;
}

template <class Float>
int bfmbase<Float>::interleave_site(int mp,int mu, int site)
{
  return(mp+mu*NMinusPlus+site*NMinusPlus*ND*2);
}

//
// To switch to 5d preconditioning, I think I just need to
// make parity always return zero, and make psite not divide by two.
// introduce a flat precon_5d in the bfmarg class
//

template <class Float>
int bfmbase<Float>::parity(int x[4])
{
  if ( this->precon_5d ) return 0;
  else return((x[0]+x[1]+x[2]+x[3]) % 2);
}

template <class Float>
int bfmbase<Float>::psite(int addr[4],integer latt[4])
{
  int site;
  site =  addr[0] + latt[0]*(addr[1]
                             + latt[1]*(addr[2]
                                        + latt[2]*addr[3]));
  if ( this->precon_5d ) return site;
  else return(site/2);
}

template <class Float>
void bfmbase<Float>::bgq_l1p_optimisation(int mode)
{
#ifdef BGQ
#define L1P_OPT
#endif
#ifdef L1P_OPT
#undef L1P_CFG_PF_USR
#define L1P_CFG_PF_USR  (0x3fde8000108ll)   /*  (64 bit reg, 23 bits wide, user/unpriv) */

  uint64_t cfg_pf_usr;
  if ( mode ) { 
    cfg_pf_usr =
        L1P_CFG_PF_USR_ifetch_depth(0)       
      | L1P_CFG_PF_USR_ifetch_max_footprint(1)   
      | L1P_CFG_PF_USR_pf_stream_est_on_dcbt 
      | L1P_CFG_PF_USR_pf_stream_establish_enable
      | L1P_CFG_PF_USR_pf_stream_optimistic
      | L1P_CFG_PF_USR_pf_adaptive_throttle(0xF) ;
    if ( sizeof(Float) == sizeof(double) ) {
      cfg_pf_usr |=  L1P_CFG_PF_USR_dfetch_depth(2)| L1P_CFG_PF_USR_dfetch_max_footprint(3)   ;
    } else {
      cfg_pf_usr |=  L1P_CFG_PF_USR_dfetch_depth(1)| L1P_CFG_PF_USR_dfetch_max_footprint(2)   ;
    }
  } else { 
    cfg_pf_usr = L1P_CFG_PF_USR_dfetch_depth(1)
      | L1P_CFG_PF_USR_dfetch_max_footprint(2)   
      | L1P_CFG_PF_USR_ifetch_depth(0)       
      | L1P_CFG_PF_USR_ifetch_max_footprint(1)   
      | L1P_CFG_PF_USR_pf_stream_est_on_dcbt 
      | L1P_CFG_PF_USR_pf_stream_establish_enable
      | L1P_CFG_PF_USR_pf_stream_optimistic
      | L1P_CFG_PF_USR_pf_stream_prefetch_enable;
  }

  *((uint64_t *)L1P_CFG_PF_USR) = cfg_pf_usr;

#endif

}
/*
template <class Float>
void bfmbase<Float>::gather_proj(integer *sites,
			     Float     *psi,
			     integer comm_off, 
			     Float  *comm_buf,
			     integer nface,
			     int permute, 
			     int L,
			     int pm,
			     int mu)
{

  int me,thrvol,throff;

  this->thread_work(nface,me,thrvol,throff);

  integer args[8];
  args[0] = (integer) &sites[throff];
  args[1] = (integer) psi;
  args[2] = (integer) comm_off+throff;
  args[3] = (integer) comm_buf; // Fix me -- need offset into this
  args[4] = (integer) thrvol;
  args[5] = (integer) L;
  args[6] = pm+2*mu;
  args[7] = (integer)complex_i_simd;

  integer arg_p = (integer) args;

  if ( permute ) { 
    if ( sizeof(Float) == sizeof(double ) ) vmx_gather_proj_perm(arg_p);
    else                                    vmx_gather_proj_perm_s(arg_p);
  } else { 
    if ( sizeof(Float) == sizeof(double ) ) vmx_gather_proj(arg_p);
    else                                    vmx_gather_proj_s(arg_p);
  }  
  this->thread_barrier();
}
*/
template <class Float>
void bfmbase<Float>::applyCloverCpp(Fermion_t psi,Fermion_t chi,int cb,int inv)
{
#ifdef BFM_QDP_BINDING
    int Nelem=42; /*6x6 at each site ; number of complex*/
    int Nspinco=12;
    int Nc=3;
    Float * spchi = (Float *)chi;
    Float * sppsi = (Float *)psi;

    int vol = this->node_latt[0]*this->node_latt[1]*this->node_latt[2]*this->node_latt[3];

#pragma omp parallel
  {
#pragma omp for
    for (int site=0; site<vol/2; site++ )
    {
      for(int cbsite = site;cbsite<vol;cbsite+=vol/2)
      {
        int x[4];
        int s=cbsite;
            x[0]=s%this->node_latt[0];
            s=s/this->node_latt[0];
            x[1]=s%this->node_latt[1];
            s=s/this->node_latt[1];
            x[2]=s%this->node_latt[2];
            s=s/this->node_latt[2];
            x[3]=s%this->node_latt[3];	        

        int ccb   = ((x[0]+x[1]+x[2]+x[3])&0x1); /*Work out local checkerboard of site*/  
   
        if(ccb!=cb)
           continue;
        
        Float* ppsi = &(sppsi[bagel_idx(x,0,0,Nspinco,cb)]);
        Float* pA;

        if(inv==0)
        {
           pA   = &(this->A[bagel_idx(x,0,0,Nelem,0)]);
        }
        else if(inv==1)
        {
           pA   = &(this->Ainv[bagel_idx(x,0,0,Nelem,0)]);
        }

        complex<Float> soln[12];

	soln[0] =  complex<Float>(pA[0],pA[1])         * complex<Float>(ppsi[0],ppsi[1])
               +   std::conj(complex<Float>(pA[24],pA[25])) * complex<Float>(ppsi[4],ppsi[5])
               +   std::conj(complex<Float>(pA[28],pA[29])) * complex<Float>(ppsi[8],ppsi[9])
               +   std::conj(complex<Float>(pA[36],pA[37])) * complex<Float>(ppsi[12],ppsi[13])
               +   std::conj(complex<Float>(pA[48],pA[49])) * complex<Float>(ppsi[16],ppsi[17])
               +   std::conj(complex<Float>(pA[64],pA[65])) * complex<Float>(ppsi[20],ppsi[21]);      

       soln[1] =   complex<Float>(pA[24],pA[25])       * complex<Float>(ppsi[0],ppsi[1])
               +   complex<Float>(pA[4] ,pA[5] )       * complex<Float>(ppsi[4],ppsi[5])
               +   std::conj(complex<Float>(pA[32],pA[33])) * complex<Float>(ppsi[8],ppsi[9])
               +   std::conj(complex<Float>(pA[40],pA[41])) * complex<Float>(ppsi[12],ppsi[13])
               +   std::conj(complex<Float>(pA[52],pA[53])) * complex<Float>(ppsi[16],ppsi[17])
               +   std::conj(complex<Float>(pA[68],pA[69])) * complex<Float>(ppsi[20],ppsi[21]);

       soln[2] =   complex<Float>(pA[28],pA[29])       * complex<Float>(ppsi[0],ppsi[1])
               +   complex<Float>(pA[32],pA[33])       * complex<Float>(ppsi[4],ppsi[5])
               +   complex<Float>(pA[8] ,pA[9] )       * complex<Float>(ppsi[8],ppsi[9])
               +   std::conj(complex<Float>(pA[44],pA[45])) * complex<Float>(ppsi[12],ppsi[13])
               +   std::conj(complex<Float>(pA[56],pA[57])) * complex<Float>(ppsi[16],ppsi[17])
               +   std::conj(complex<Float>(pA[72],pA[73])) * complex<Float>(ppsi[20],ppsi[21]);

       soln[3] =   complex<Float>(pA[36],pA[37])       * complex<Float>(ppsi[0],ppsi[1])
               +   complex<Float>(pA[40],pA[41])       * complex<Float>(ppsi[4],ppsi[5])
               +   complex<Float>(pA[44],pA[45])       * complex<Float>(ppsi[8],ppsi[9])
               +   complex<Float>(pA[12],pA[13])       * complex<Float>(ppsi[12],ppsi[13])
               +   std::conj(complex<Float>(pA[60],pA[61])) * complex<Float>(ppsi[16],ppsi[17])
               +   std::conj(complex<Float>(pA[76],pA[77])) * complex<Float>(ppsi[20],ppsi[21]);

       soln[4] =   complex<Float>(pA[48],pA[49])       * complex<Float>(ppsi[0],ppsi[1])
               +   complex<Float>(pA[52],pA[53])       * complex<Float>(ppsi[4],ppsi[5])
               +   complex<Float>(pA[56],pA[57])       * complex<Float>(ppsi[8],ppsi[9])
               +   complex<Float>(pA[60],pA[61])       * complex<Float>(ppsi[12],ppsi[13])
               +   complex<Float>(pA[16],pA[17])       * complex<Float>(ppsi[16],ppsi[17])
               +   std::conj(complex<Float>(pA[80],pA[81])) * complex<Float>(ppsi[20],ppsi[21]);

       soln[5] =   complex<Float>(pA[64],pA[65])       * complex<Float>(ppsi[0],ppsi[1])
               +   complex<Float>(pA[68],pA[69])       * complex<Float>(ppsi[4],ppsi[5])
               +   complex<Float>(pA[72],pA[73])       * complex<Float>(ppsi[8],ppsi[9])
               +   complex<Float>(pA[76],pA[77])       * complex<Float>(ppsi[12],ppsi[13])
               +   complex<Float>(pA[80],pA[81])       * complex<Float>(ppsi[16],ppsi[17])
               +   complex<Float>(pA[20],pA[21])       * complex<Float>(ppsi[20],ppsi[21]);

        ppsi = &(sppsi[bagel_idx(x,0,Nspinco/2,Nspinco,cb)]);
        
        if(inv==0)
        {
           pA   = &(this->A[bagel_idx(x,0,Nelem/2,Nelem,0)]);
        }
        else if(inv==1)
        {
           pA   = &(this->Ainv[bagel_idx(x,0,Nelem/2,Nelem,0)]);
        }


	soln[6] =  complex<Float>(pA[0],pA[1])         * complex<Float>(ppsi[0],ppsi[1])
               +   std::conj(complex<Float>(pA[24],pA[25])) * complex<Float>(ppsi[4],ppsi[5])
               +   std::conj(complex<Float>(pA[28],pA[29])) * complex<Float>(ppsi[8],ppsi[9])
               +   std::conj(complex<Float>(pA[36],pA[37])) * complex<Float>(ppsi[12],ppsi[13])
               +   std::conj(complex<Float>(pA[48],pA[49])) * complex<Float>(ppsi[16],ppsi[17])
               +   std::conj(complex<Float>(pA[64],pA[65])) * complex<Float>(ppsi[20],ppsi[21]);

	soln[7] =  complex<Float>(pA[24],pA[25])       * complex<Float>(ppsi[0],ppsi[1])
               +   complex<Float>(pA[4] ,pA[5] )       * complex<Float>(ppsi[4],ppsi[5])
               +   std::conj(complex<Float>(pA[32],pA[33])) * complex<Float>(ppsi[8],ppsi[9])
               +   std::conj(complex<Float>(pA[40],pA[41])) * complex<Float>(ppsi[12],ppsi[13])
               +   std::conj(complex<Float>(pA[52],pA[53])) * complex<Float>(ppsi[16],ppsi[17])
               +   std::conj(complex<Float>(pA[68],pA[69])) * complex<Float>(ppsi[20],ppsi[21]);

	soln[8] =  complex<Float>(pA[28],pA[29])       * complex<Float>(ppsi[0],ppsi[1])
               +   complex<Float>(pA[32],pA[33])       * complex<Float>(ppsi[4],ppsi[5])
               +   complex<Float>(pA[8] ,pA[9] )       * complex<Float>(ppsi[8],ppsi[9])
               +   std::conj(complex<Float>(pA[44],pA[45])) * complex<Float>(ppsi[12],ppsi[13])
               +   std::conj(complex<Float>(pA[56],pA[57])) * complex<Float>(ppsi[16],ppsi[17])
               +   std::conj(complex<Float>(pA[72],pA[73])) * complex<Float>(ppsi[20],ppsi[21]);

        soln[9] =  complex<Float>(pA[36],pA[37])       * complex<Float>(ppsi[0],ppsi[1])
               +   complex<Float>(pA[40],pA[41])       * complex<Float>(ppsi[4],ppsi[5])
               +   complex<Float>(pA[44],pA[45])       * complex<Float>(ppsi[8],ppsi[9])
               +   complex<Float>(pA[12],pA[13])       * complex<Float>(ppsi[12],ppsi[13])
               +   std::conj(complex<Float>(pA[60],pA[61])) * complex<Float>(ppsi[16],ppsi[17])
               +   std::conj(complex<Float>(pA[76],pA[77])) * complex<Float>(ppsi[20],ppsi[21]);

	soln[10] =  complex<Float>(pA[48],pA[49])       * complex<Float>(ppsi[0],ppsi[1])
               +   complex<Float>(pA[52],pA[53])       * complex<Float>(ppsi[4],ppsi[5])
               +   complex<Float>(pA[56],pA[57])       * complex<Float>(ppsi[8],ppsi[9])
               +   complex<Float>(pA[60],pA[61])       * complex<Float>(ppsi[12],ppsi[13])
               +   complex<Float>(pA[16],pA[17])       * complex<Float>(ppsi[16],ppsi[17])
               +   std::conj(complex<Float>(pA[80],pA[81])) * complex<Float>(ppsi[20],ppsi[21]);
        
       soln[11] =   complex<Float>(pA[64],pA[65])       * complex<Float>(ppsi[0],ppsi[1])
               +   complex<Float>(pA[68],pA[69])       * complex<Float>(ppsi[4],ppsi[5])
               +   complex<Float>(pA[72],pA[73])       * complex<Float>(ppsi[8],ppsi[9])
               +   complex<Float>(pA[76],pA[77])       * complex<Float>(ppsi[12],ppsi[13])
               +   complex<Float>(pA[80],pA[81])       * complex<Float>(ppsi[16],ppsi[17])
               +   complex<Float>(pA[20],pA[21])       * complex<Float>(ppsi[20],ppsi[21]);
 
     	Float* pchi = &(spchi[bagel_idx(x,0,0,Nspinco,cb)]);
    	for(int k=0;k<12;k++)
        {	   
	   pchi[k*4] = real(soln[k]);
	   pchi[k*4+1] = imag(soln[k]);
        }
      }

    }

  }
  
           
#endif
}


template <class Float>
void bfmbase<Float>::applyClover(Fermion_t psi,Fermion_t chi,int cb,int inv)
{

     /*-----------------------------------------------------
     *Divide up work between threads
     *-----------------------------------------------------
     */
     int me, thrvol,throff;
     thread_work(this->simd_cbvol,me,thrvol,throff);     

     /*Turn global checkerboard into the checkerboard within this node*/
     int lcb = (cb + this->base_parity) & 1;
     cb = lcb;
     
     /*
     * Setup Clover fields 
     */
    Float *clover_par;
    Float *clovtmp;
    if(inv==0)
       clovtmp = (Float *)this->A;
    else
       clovtmp = (Float *)this->Ainv;
    if ( cb == 1 )
    {
        clover_par = clovtmp+ (2*CLOVER_MAT_SIZE*this->simd_cbvol*this->simd());
    }
    else
    {
        clover_par = clovtmp;
    }

    thrvol =  2*thrvol;
    throff =  2*throff;

    //bgq_l1p_optimisation(1);

    integer args[6];
    integer arg_p = (integer)args;
    args[0] = ((integer)chi) + throff * this->cbLs  * this->nsimd * 12 * sizeof(Float);
    args[1] = ((integer)psi) + throff * this->cbLs  * this->nsimd * 12 * sizeof(Float);
    args[2] = ((integer)clover_par) + throff * 42  * sizeof(Float) * this->simd(); 
    args[3] = (integer)thrvol;
    args[4] = (integer)cb;                             //do I need this

    int single =  1;
    if ( sizeof(Float) == sizeof(double) ) single = 0;

    uint64_t t1 = 0,t2 = 0, tclov = 0;

    this->thread_barrier();

    t1 = GetTimeBase();
    if (single)   
	vmx_clov_apply_s(arg_p);
    else
        vmx_clov_apply(arg_p);

    this->thread_barrier();
    t2 = GetTimeBase();

    tclov = t2-t1;
      
    //bgq_l1p_optimisation(0);
    
}
template <class Float> 
void bfmbase<Float>::gather(int result_cb, Fermion_t psi,int dag)
{
  for(int mu=0;mu<4;mu++) gather(mu,result_cb,psi,dag);
}
template <class Float> 
void bfmbase<Float>::gather(int mu,int result_cb, Fermion_t psi,int dag)
{
  int nspinco = 12;
  int cb = 1-result_cb;

  int me = this->thread_barrier();

  for(int pm=0;pm<2;pm++){

      int permute=0;
      int extract=0;
      int ftcb =0;
      int dir=pm+2*mu;
      int sgn =pm;
      Float *comm_buf = this->sendbufs[dir];

      if ( dag ) sgn = 1-pm;

      if ( this->precon_5d ) ftcb = 0;
      else ftcb = result_cb;

      if ( (this->simd_grid[mu]>1) && (this->local_comm[mu]) ){

	/******************************************************
	 * Direction with SIMD interleave is special
	 * Gather into recvbuf with permute set if local.
	 ******************************************************/
        permute = 1;
	comm_buf = this->recvbufs[dir];

      }

      if ( (this->simd_grid[mu]>1) && (!this->local_comm[mu]) ){
	/******************************************************
	 * Direction with SIMD interleave is special
	 * Gather into sendbuf extracting subset of SIMD vec 
	 ******************************************************/
	// Here we must take array : (AE BF CG DH) (IM JN KO LP)
	// Pull AE in, A-> sendbuf
	//
	// After receive of "I" from next node into the receive buf.
	// Should receive into high part to enable in place merge.
	// This may come later, for now can do it poorly prefetched descending order
        extract = 1;
      }

      int do_comm = 0;
      if ( !this->local_comm[mu]   ) do_comm=1;
      if ( this->simd_grid[mu] > 1 ) do_comm=1;
      if ( this->simd_grid[mu] > 2 ) {
	this->Error("simd_grid[%d] > 2 not supported\n",mu);
	exit(0); // Not supported yet
      }
      if ( this->simd()        > 4 ) {
	this->Error("simd > 4 not supported\n");
	exit(0); // Not supported yet
      }


      if (do_comm) {

	int thrvol,throff;
	this->thread_work_nobarrier(this->simd_nbound[mu],me,thrvol,throff);

	integer args[9];
	args[0] = (integer) &this->face_table[ftcb][pm][mu][throff];
	args[1] = (integer) psi;
	args[2] = (integer) throff;
	args[3] = (integer) comm_buf; 
	args[4] = (integer) thrvol;
	args[5] = (integer) this->cbLs;
	args[6] = sgn+2*mu;
	args[7] = (integer) this->complex_i_simd;

	integer arg_p = (integer) args;

	if ( permute ) { 

	  if ( sizeof(Float) == sizeof(double ) ) vmx_gather_proj_perm(arg_p);
#ifdef BGQ
	  else if (this->precision_test)          vmx_gather_proj_perm_hs(arg_p);
#endif
	  else                                    vmx_gather_proj_perm_s(arg_p);

	} else if ( extract ) { 
	  args[3] = (integer) this->sendbufs[dir]; // Off node
	  args[8] = (integer) this->sendbufs[dir]  // On node
	    + this->simd_nbound[mu]*this->cbLs*12*sizeof(Float); // NB no nsimd
	  if ( pm == 0 ) {
	    integer tmp =args[3];
	    args[3]=args[8];
	    args[8]=tmp;
	  }

	  if ( sizeof(Float) == sizeof(double ) ) vmx_gather_proj_extr_l(arg_p);
#ifdef BGQ
	  else if ( this->precision_test ) vmx_gather_proj_extr_l_hs(arg_p);
#endif
	  else  vmx_gather_proj_extr_l_s(arg_p);
	  
	} else { 
	  if ( sizeof(Float) == sizeof(double ) ) vmx_gather_proj(arg_p);
#ifdef BGQ
	  else if (this->precision_test)          vmx_gather_proj_hs(arg_p);
#endif
	  else                                    vmx_gather_proj_s(arg_p);
	}  
      }
  }

  this->thread_barrier();
}

		   
template class bfmbase<float>;
template class bfmbase<double>;
