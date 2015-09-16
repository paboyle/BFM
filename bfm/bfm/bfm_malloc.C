#include <bfm.h>

void *bfm_alloc(size_t size,int pooh)
{
  void *ptr = malloc(size);
  if ( ptr == NULL ) {
    printf("bad bfm_alloc\n");
    exit(-1);
  }
  return ptr;
}

void  bfm_free (void *ptr)
{
  free(ptr);
}


