
#include "bfm.h"
#include "bfm_vmx.h"
#include "bfm_mprec.h"
#include "bfm_linalg.h"
#include "bfm_cg.h"
#include "bfm_mcr.h"
#include "bfm_qdp.h"
#include "bfm_import_export.h"


#include <sys/time.h>



template <class Float>
void bfmbase<Float>::comm_gsum(double *val,int N) { return ; };


template <class Float>
void bfmbase<Float>::comm_gsum(double &val) { return ; };

template <class Float>
bool bfmbase<Float>::isBoss() {return 1;};

// Default values
double bfmActionParams::mobius_scale = 1.7;
int bfmarg::threads = 1;
int bfmarg::reproduce = 0;
int bfmarg::reproduce_checksum = 0;
int bfmarg::reproduce_master_check= 0;
int bfmarg::CGdiagonalMee= 0;
int bfmarg::part_frac_chroma_convention = 1;
char *bfmarg::watchfile=NULL;
int bfmarg::onepluskappanorm;
int bfmarg::default_verbose= BfmMessage|BfmError;

template <class Float>
void bfmbase<Float>::check_is_same(double d,char *s)
{
  int me = this->thread_barrier();
  double checker =0;
  if (this->isBoss() && !me ) { 
    checker=d;
  }
  this->thread_sum(checker,me);
  this->comm_gsum(checker);
  if ( checker != d ) {
    Error("check_is_same %s found discrepancy %d me %le node 0 %le\n",s,me, d,checker);
  }
  return;
}

template <class Float>
void bfmbase<Float>::ThreadBossDebug(const char *fmt,...)
{ 
  va_list args;
  va_start(args,fmt);
  int me = this->thread_barrier();
  if ( this->isBoss() && (verbose&BfmDebug) && !me ) {
    printf("Bfm:   Debug: ");
    vprintf(fmt,args);
    fflush(stdout);
  }
  va_end(args);
}
template <class Float>
void bfmbase<Float>::ThreadBossMessage(const char *fmt,...)
{
  va_list args;
  va_start(args,fmt);
  int me = this->thread_barrier();
  if ( this->isBoss() && (verbose&BfmMessage) && !me) {
    printf("Bfm: Message: ");
    vprintf(fmt,args);
    fflush(stdout);
  }
  va_end(args);
}
template <class Float>
void bfmbase<Float>::ThreadBossPerformance(const char *fmt,...)
{
  va_list args;
  va_start(args,fmt);
  int me = this->thread_barrier();
  if ( this->isBoss() && (verbose&BfmPerformance) && !me) {
    printf("Bfm:    Perf: ");
    vprintf(fmt,args);
    fflush(stdout);
  }
  va_end(args);
}
template <class Float>
void bfmbase<Float>::ThreadBossLog(const char *fmt,...)
{
  va_list args;
  va_start(args,fmt);
  int me = this->thread_barrier();
  if ( this->isBoss() && !me) {
    printf("Bfm:     Log: ");
    vprintf(fmt,args);
    fflush(stdout);
  }
  va_end(args);
}
template <class Float>
void bfmbase<Float>::BossDebug(const char *fmt,...)
{ 
  va_list args;
  va_start(args,fmt);
  if ( this->isBoss() && (verbose&BfmDebug) ) {
    vprintf(fmt,args);
    fflush(stdout);
  }
  va_end(args);
}
template <class Float>
void bfmbase<Float>::BossLog(const char *fmt,...)
{ 
  va_list args;
  va_start(args,fmt);
  if ( this->isBoss() ) {
    printf("Bfm:     Log: ");
    vprintf(fmt,args);
    fflush(stdout);
  }
  va_end(args);
}
template <class Float>
void bfmbase<Float>::BossMessage(const char *fmt,...)
{
  va_list args;
  va_start(args,fmt);
  if ( this->isBoss() && (verbose&BfmMessage) ) {
    printf("Bfm: Message: ");
    vprintf(fmt,args);
    fflush(stdout);
  }
  va_end(args);
}
template <class Float>
void bfmbase<Float>::NodeDebug(const char *fmt,...)
{
  va_list args;
  va_start(args,fmt);
  if ( (verbose&BfmDebug) ) {
    printf("Bfm:   Debug: ");
    vprintf(fmt,args);
    fflush(stdout);
  }
  va_end(args);
}
template <class Float>
void bfmbase<Float>::NodeMessage(const char *fmt,...)
{
  va_list args;
  va_start(args,fmt);
  if ( (verbose&BfmMessage) ) {
    printf("Bfm: Message: ");
    vprintf(fmt,args);
    fflush(stdout);
  }
  va_end(args);
}
template <class Float>
void bfmbase<Float>::Error(const char *fmt,...)
{
  va_list args;
  va_start(args,fmt);
  if ( (verbose&BfmError) ) {
    printf("Bfm:   Error: ");
    vprintf(fmt,args);
    fflush(stdout);
  }
  va_end(args);
}

template <class Float>
int bfmbase<Float>::CheckStopWatchFile(void)
{
  if ( watchfile ) {
    FILE * fp = fopen(watchfile,"r");
    if ( fp ) {
      fclose(fp);
      Error("Found watchfile %s\n",watchfile);
      return 1;
    }
  }
  return 0;
}

#define MAX_ENTRIES (50000)

#ifdef BFM_REPRO
template <class Float>
void bfmbase<Float>::repro_control(int state) { 

  int me = thread_barrier();
  repro_state = state; 

  thread_barrier();

  if ( repro_state == ReproIdle ) return;

  repro_entry = 0;
  csum_entry = 0;

  if ( (repro_state == ReproRecord) && (me == 0) ){
    repro_log.resize(nthread+1);
    for(int tid=0;tid<nthread+1;tid++){
      repro_log[tid].resize(MAX_ENTRIES);
    }
    csum_log.resize(MAX_ENTRIES);
  }
  thread_barrier();
};

template <class Float>
void bfmbase<Float>::repro(uint64_t csum) {

  if ( csum_entry >= MAX_ENTRIES-1 )  return;

  if ( repro_state == ReproRecord ) csum_log[csum_entry] = csum;

  if ( repro_state == ReproCheck  ) {
    uint64_t pop = csum_log[csum_entry];
    if( csum != pop ) { 
    
      Error("Repro mismatch %lx e=%d csum : %16.16lx %16.16lx\n",
	      (uint64_t)this,csum_entry,csum,pop);
      
      
      fflush(stderr);
      abort();
    }
  }

  csum_entry++;

}
template <class Float>
void bfmbase<Float>::repro(double d,int who) {

  if ( repro_entry >= MAX_ENTRIES-1 ) return;

  if ( (who < 0) || (who > nthread) ) { 
    Error("repro: bad who %d\n",who);
    exit(-1);
  }

  if ( repro_state == ReproRecord ) repro_log[who][repro_entry] = d;
  if ( repro_state == ReproCheck  ) {
    double pop = repro_log[who][repro_entry];
    if( d != pop ) { 
      Error("Repro mismatch %lx e=%d : %16.8le %16.8le %16.16lx %16.16lx lane %d\n",
	      (uint64_t)this,repro_entry,d,pop,*((uint64_t*)&d),*((uint64_t *)&pop),who);
      fflush(stderr);
      fflush(stdout);
      abort();
    } 
  }
  if ( who ==nthread ) {
    repro_entry++;
  }
};
#endif


template <class Float>
int bfmbase<Float>::doubleStored(void) { 
  return 1;
};

uint64_t NumInversions;
uint64_t FirstInversionUsec;
uint64_t InversionTimeToDateUsec;
int      NestedInversionDepth;

template <class Float>
void bfmbase<Float>::InverterEnter(void)
{
  NestedInversionDepth++;
  if ( NumInversions == 0 ) {
    FirstInversionUsec = Microsecond();
  }
  if ( NestedInversionDepth == 1 ) {
    InversionEnteredUsec=Microsecond();
    NumInversions++;
  }
}
template <class Float>
void bfmbase<Float>::InverterExit (void)
{
  NestedInversionDepth--;

  if ( NestedInversionDepth == 0 ) { 

    uint64_t now = Microsecond();
    uint64_t elapsed = now-InversionEnteredUsec;

    InversionTimeToDateUsec += elapsed;
    BossDebug("Inversion :            %.6f s \n",1.0e-6*elapsed);
    BossDebug("Inversion time to date %.6f s \n",1.0e-6*InversionTimeToDateUsec);

    uint64_t wall_clock =  now-FirstInversionUsec;
    BossDebug("Elapsed   time to date %.6f s \n",1.0e-6*wall_clock);
    BossDebug("Percentage time in inverter %.3f\n",100.0*InversionTimeToDateUsec/wall_clock);
  }
}
template <class Float>
void bfmbase<Float>::InverterLoggingBegin(std::string file)
{
  /*
  std::string InverterLogFile = file;
  if ( this->isBoss() ) {
    InverterLogFP = fopen(InverterLogFile.c_str(),"w");
    fprintf(InverterLogFP,"# <iter> <computed_residual> <true_residual> <true_error>\n");
  }
  LoggingSuppressed=0;
  */
}
template <class Float>
void bfmbase<Float>::InverterLoggingEnd(void)
{
  //
  //  if ( InverterLogFP ) {
  //    fclose(InverterLogFP);
  //  }
}

template <class Float>
void bfmbase<Float>::InverterRegisterExactSolution(Fermion_t solution)
{
  //  LoggingSuppressed=0;
  //  InverterCompareExactSolution=1;
  //  InverterExactSolution = solution;
  //  InverterLogTrueResidual=1;
  //  printf("Logging enaled\n");
}
template <class Float>
void bfmbase<Float>::InverterLogIteration(int iter, double residual, 
					  Fermion_t src,
					  Fermion_t current_solution,
					  Fermion_t tmp1,
					  Fermion_t tmp2,
					  Fermion_t Mtmp)
{
  return;

  if ( LoggingSuppressed ) return;

  int me = this->thread_barrier();


  uint64_t now = Microsecond();
  uint64_t elapsed = now-InversionEnteredUsec;

  double true_residual=0;
  double true_error=0;

  if ( InverterCompareExactSolution ) {
    this->axpy(tmp1,InverterExactSolution,current_solution,-1.0);
    true_error=this->norm(tmp1);
  }

  if ( InverterLogTrueResidual ) {
    this->Mprec(current_solution,tmp1,Mtmp,0);
    this->Mprec(tmp1,tmp2,Mtmp,1);
    this->axpy(tmp2,tmp2,src,-1.0);
    true_residual = this->norm(tmp2);
  }


  ThreadBossMessage("%d %le %le %le %le\n", iter,1.0e-6*elapsed,
	      residual,
	      true_residual,
	      true_error
	      );

  // Want to get timeslice by timeslice difference
  // Use block inner routine to get the spin trace
  // Implement a gamma multiply

}


template class bfmbase<double>;
template class bfmbase<float>;
template void bfmbase<double>::cps_importGauge<double>(double *importme);
template void bfmbase<float>::cps_importGauge<double>(double *importme);
