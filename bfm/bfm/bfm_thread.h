#ifndef BFM_THREAD_H
#define BFM_THREAD_H


/*
 *Thread model - threads created EXTERNAL to Bagel/BFM.
 *
 *             - Number of threads passed in.
 *
 *             - Threads register themselves with the bagel (get threadid)
 *
 *             - Boss thread 
 *               - does the alloc's, table inits etc...
 *
 *             - All threads call all numerical routines.
 *               - linalg Mprec etc... all to a per thread chunk of work
 *                 & contain barriers()
 *
 *             - CG contains allocs. These should be replaced with AllocBcast.
 *
 *             - Simply need a barrier & gsum across threads & bcast
 *
 *      Thus:
 *           Setup is single thread
 *           Fork and join is done at the ENTIRE CG level.
 *           Registration is tricky. Maps pthread_self to a work unit.
 */

void SpiL2AtomicInit(void);

class ThreadModelSingle {
 public:
  virtual void thread_init      (int nthr) { 
     if (nthr != 1) exit(0); 
     nthread=nthr; 
  };
  virtual void thread_end      (void) { 
  }
  virtual int  thread_barrier   (void) { return 0; };
  virtual void thread_sum       (double &val, int me) { return; };
  virtual void thread_csum   (uint64_t &val,int me) {return;};
  virtual void  *thread_bcast(int me,void * val) { return val; };

  void thread_work_nobarrier(int nwork, int me, int & mywork, int & myoff){
    int threads = nthread;
    int basework = nwork/threads;
    int backfill = threads-(nwork%threads);
    if ( me >= threads ) { 
      mywork = myoff = 0;
    } else { 
      mywork = (nwork+me)/threads;
      myoff  = basework * me;
      if ( me > backfill ) 
	myoff+= (me-backfill);
    }
    return;
  }
  void thread_work(int nwork, int &me, int & mywork, int & myoff){
    me     = thread_barrier();
    thread_work_nobarrier(nwork,me,mywork,myoff);
  }
  int nthread;
};

#ifdef THREAD_MODEL_SINGLE

typedef ThreadModelSingle ThreadModel;

#endif

#ifdef THREAD_MODEL_PTHREAD
#error PTHREAD not supported
#endif

#ifdef THREAD_MODEL_OPENMP
#include <omp.h>
class ThreadModelOpenMP : public ThreadModelSingle {
 public:
  virtual void thread_init      (int nthr) {
    nthread=nthr;
    omp_set_num_threads(nthread);
    sum_array = (double *)malloc(nthr*sizeof(double)); 
  };
  virtual void thread_end       (void){
    free(sum_array);
  };
  virtual int  thread_barrier(void) {
#pragma omp barrier
    return omp_get_thread_num();
  };
  virtual void thread_sum    (double &val,int me){
    sum_array[me] = val;
    val=0;
    thread_barrier();
    for(int i=0;i<nthread;i++) val+= sum_array[i];
    thread_barrier();
  }
  virtual void thread_csum   (uint64_t &val,int me){
    uint64_t *csum_array = (uint64_t *)sum_array;
    csum_array[me] = val;
    val=0;
    thread_barrier();
    for(int i=0;i<nthread;i++) val^= csum_array[i];
    thread_barrier();
  }
  virtual void  *thread_bcast(int me,void * val)
  {
    if (me == 0) bcast_ptr = val;
    thread_barrier();
    val = bcast_ptr;
    thread_barrier();
    return val; 
  }
 public:
  double * sum_array;
  void   * bcast_ptr;
};
typedef ThreadModelOpenMP ThreadModel;

#endif

#ifdef THREAD_MODEL_SPI

#include <spi/include/l2/atomic.h>
//#include <spi/include/l2/barrier.h>

#define L1D_CACHE_LINE_SIZE 64

typedef struct {
  volatile __attribute__((aligned(L1D_CACHE_LINE_SIZE)))
    uint64_t start;  /*!< Thread count at start of current round. */
  volatile __attribute__((aligned(L1D_CACHE_LINE_SIZE)))
    uint64_t count;  /*!< Current thread count. */
} L2_Barrier_t;


#define MAX_HW_THREADS (512)
class ThreadModelSPI : public ThreadModelSingle {
 public:
  virtual void thread_init      (int nthread);
  virtual void thread_end       (void);
  virtual int  thread_barrier   (void) ;
  virtual void thread_sum       (double &val, int me) ;
  virtual void thread_csum   (uint64_t &val,int me) ;
  virtual void  *thread_bcast(int me,void * val);
  void thread_init_noalloc      (int nthread);

  //SPI specific ; prevent work migration between threads
  int  bgq_hwid(void);
  void bgq_hwid_register(int hwid);

  int ID_map[MAX_HW_THREADS];
  int ID_map_initialized;

 protected:
 public:
  double * sum_array;
  void   * bcast_ptr;
};
void SPI_PhysicalThreadMapping(int yesno);
typedef ThreadModelSPI ThreadModel;
#endif


#endif
