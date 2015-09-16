/*
 * Can I come up with a lock-free implementation that is the
 * same for pthread and BGQ? Need N+1
 */

#include <bfm.h>
#include <stdio.h>
#include <spi/include/kernel/location.h>
#include <pthread.h>

L2_Barrier_t init_b = {0,0};

/*Locate in atomic op space*/
L2_Barrier_t  b ;

int  PhysicalThreadMapping=0;
void SPI_PhysicalThreadMapping(int yesno)
{
  PhysicalThreadMapping = yesno;
}

int L2_BarrierWithTicket(L2_Barrier_t *b, int numthreads);

int TLBmapped=0;

void SpiL2AtomicInit(void);

void SpiL2AtomicInit(void)
{
    Kernel_L2AtomicsAllocate(&b,sizeof(L2_Barrier_t));
}

int  ThreadModelSPI::bgq_hwid(void)
{
  int hwid = Kernel_PhysicalProcessorID() | (Kernel_PhysicalHWThreadID() <<4);
  if ( hwid >= nthread ) exit(-1);
  return hwid;
}


int L2_BarrierWithTicket(L2_Barrier_t *b, int numthreads)
{
  if( !TLBmapped)  SpiL2AtomicInit();
  TLBmapped=1;

  uint64_t start = b->start;
  ppc_msync();  // make sure we pick up the correct start for this round
  uint64_t count = L2_AtomicLoadIncrement(&b->count);

  uint64_t target = start + numthreads;
  uint64_t current = count + 1;

  if (current == target) {
    b->start = current;  // advance to next round
  } else {
    int count = 0;
    while (b->start < current) {
      count++;
      if ( count > 10*1024 ) {
	count =0; 
	pthread_yield();
      }
    }  // wait for advance to next round
  }

  return (int) (count - start);
}


void ThreadModelSPI::thread_init_noalloc  (int nthr) 
{
  nthread = nthr;
  b = init_b;
}
void ThreadModelSPI::thread_end  (void)  
{
  free(sum_array);
}

void ThreadModelSPI::thread_init  (int nthr)  
{
  ID_map_initialized = 0;
  thread_init_noalloc  ( nthr) ;
  sum_array = (double *)malloc(nthr*sizeof(double)); // Missing free
}
 
int  ThreadModelSPI::thread_barrier(void) 
{
  int me ;
  
  if ( PhysicalThreadMapping ) {
    L2_BarrierWithTicket(&b,nthread);
    me = bgq_hwid();
  } else { 
    me = L2_BarrierWithTicket(&b,nthread);
  }

  //  Delay(me*200);

  return me;
}

#undef MASTER_CHECK
/*This could be shared*/
void ThreadModelSPI::thread_sum    (double &val,int me) 
{
  sum_array[me] = val;
  val=0;
  thread_barrier();
  for(int i=0;i<nthread;i++) val+= sum_array[i];
  thread_barrier();

#ifdef MASTER_CHECK
    double * master = thread_bcast(me,&val);
    double mval = *master;
    if ( val != mval ) {
      printf("thread_sum : tid %d not all threads agree on sum %16.8le %16.8le\n",me,val,mval);
      for(int i=0;i<nthread;i++){
	printf(" sum_array[%d]=%le\n",i,sum_array[i]);
      }
      exit(-1);
    }
  thread_barrier();
#endif


}

void ThreadModelSPI::thread_csum    (uint64_t &val,int me) 
{
  uint64_t *csum_array = (uint64_t *)sum_array;
  csum_array[me] = val;
  val=0;

  ppc_msync();
  thread_barrier();
  for(int i=0;i<nthread;i++) val^= csum_array[i];
  thread_barrier();
  
#ifdef MASTER_CHECK
    uint64_t * master = thread_bcast(me,&val);
    if ( val != *master ) {
      printf("thread_csum : tid %d not all threads agree on sum %16.16lx %16.16lx\n",me,val,*master);
      for(int i=0;i<nthread;i++){
	printf(" csum_array[%d]=%lx\n",i,csum_array[i]);
      }
      exit(-1);
    }
    thread_barrier();
#endif
}



void  *ThreadModelSPI::thread_bcast(int me,void * val) 
{ 
  if (me == 0) bcast_ptr = val;
  thread_barrier();
  val = bcast_ptr;
  thread_barrier();
  return val; 
};

