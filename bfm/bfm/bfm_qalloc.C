#include <bfm.h>
#include <qalloc.h>

#define BFM_ALIGN_ARG (128)

void *bfm_alloc(size_t size,int mem_type)
{
  int mask;
  if ( mem_type == mem_fast    ) mask = QCOMMS|QFAST;
  if ( mem_type == mem_slow    ) mask = QCOMMS;
  if ( mem_type == mem_sendbuf ) mask = QFAST|QNONCACHE;

  void *ptr = qalloc(mask,size);
  if (ptr == NULL ) { 
    if ( mem_type == mem_fast ) { 
      return bfm_alloc(size,mem_slow);
    }
    printf("allocate failed 0x%x bytes of type %d\n",size,mem_type);
    fflush(stdout);
  }
  return ptr;
}

void  bfm_free (void *ptr)
{
  qfree(ptr);
}


