#include <bfm.h>

#define BFM_ALIGN_ARG (128)
#define BFM_LEAK_CHECK

extern "C" { 
   void* memalign(size_t, size_t);
};

#ifdef BFM_LEAK_CHECK
struct CallTree { 
  uint64_t caller[4];
};
#include <map>
std::map<uint64_t, CallTree> Pointers;
#endif


void *bfm_alloc(size_t size,int pooh)
{
  void * ptr;

  ptr = memalign(BFM_ALIGN_ARG,size+512);

  if ( (uint64_t)ptr &0x1F ) {
    printf("bad align\n"); fflush(stdout);
    exit(0);
  }
  if( ptr==NULL) {
    printf("bad alloc\n"); fflush(stdout);
    exit(0);
  }

  // Write 128 byte fenceposts
  char *cp = (char *)ptr;
  for(int i=0;i<128;i++){
    cp[i]=0x5A;
  }
  for(int i=size+128;i<size+512;i++){
    cp[i]=0xA5;
  }  
  ptr = (void *) & cp[128];


#ifdef BFM_LEAK_CHECK
  CallTree trace;
  trace.caller[0] = (uint64_t)__builtin_return_address(0);
  trace.caller[1] = (uint64_t)__builtin_return_address(1);
  trace.caller[2] = (uint64_t)__builtin_return_address(2);
  trace.caller[3] = (uint64_t)0;
  Pointers[(uint64_t)ptr] = trace;
#ifdef BFM_ALLOC_VERBOSE
  printf("BFM heap: %lx-%lx allocated by %lx %lx %lx %lx\n",
	 ptr,ptr+size,
	 trace.caller[0],trace.caller[1],trace.caller[2],trace.caller[3]);
#endif
#endif

  return ptr;
}

#ifdef BFM_LEAK_CHECK
void bfm_alloc_memcheck(void)
{
  std::map<uint64_t,CallTree>::iterator it;
  for(it = Pointers.begin();it!=Pointers.end();it++){
    CallTree trace = it->second;
    printf("BFM heap: unfreed pointer %lx allocated by %lx %lx %lx %lx\n",it->first,trace.caller[0],trace.caller[1],trace.caller[2],trace.caller[3]);
  }
}
#endif

void  bfm_ptr_check (void *ptr,int size)
{
  void *optr=ptr;

  integer iptr = (integer)ptr;
  iptr-=128;
  int bad = 0;

  ptr      = (void *) iptr;
  unsigned char *cp = (unsigned char *)ptr;
  unsigned char hi = 0xA5;
  unsigned char lo = 0x5A;
  for(int i=0;i<128;i++){
    if ( cp[i]!=lo ) {
      printf("Low Fence post overwritten pointer %lx \n",optr);
      bad = 1 ;
      break;
    };
  }

  if ( size ) { 
    for(int i=size+128;i<size+512;i++){
      if ( cp[i]!=hi ) {
	printf("High Fence post @ %d overwritten (%2.2x %2.2x) pointer %lx \n",i-size-128,cp[i],hi,optr);
	bad = 1 ;
	break;
      }
    }  
  }

  if ( bad ) {
    printf("Free call trace %lx %lx %lx %lx\n", 
	   (uint64_t)__builtin_return_address(0),
	   (uint64_t)__builtin_return_address(1),
	   (uint64_t)__builtin_return_address(2),
	   (uint64_t)__builtin_return_address(3));

#ifdef BFM_LEAK_CHECK
    CallTree trace = Pointers[(uint64_t)optr];
    printf("Allocator call trace %lx %lx %lx %lx\n", 
	   trace.caller[0],
	   trace.caller[1],
	   trace.caller[2],
	   trace.caller[3]);
      
#endif    
    exit(-1);
  }

}

void  bfm_free (void *ptr,int size)
{
  bfm_ptr_check(ptr,size);

#ifdef BFM_LEAK_CHECK
  Pointers.erase((uint64_t)ptr);
#endif

  integer iptr = (integer)ptr;
  iptr-=128;
  ptr      = (void *) iptr;
  free(ptr);
}


