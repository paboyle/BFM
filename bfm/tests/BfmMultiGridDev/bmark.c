#include <stdio.h>
#include <malloc.h>
#include <hwi/include/bqc/A2_inlines.h>

int main (int argc, char **argv)
{

  int bsize = 64*1024*1024;
  void * buffer = memalign(bsize,bsize);
  uint64_t base_addr = (uint64_t)buffer;
  uint64_t t1,t1;

  for( int mb=1;mb<64;mb++){
    double d=0;
    for(int i=0;i<3;i++){
      t1=GetTimeBase();
#pragma omp parallel 
      {
#pragma omp for 
	for(int o=0;o<mb*1024*1024;o+=32){
	  d+=*( (double *) (base_addr+o));
	}
      }
      t2=GetTimeBase();
      printf("%d MB  %ld cycles %ld MB/s d=%le\n",mb,t2-t1,mb*1600.*1024.*1024./(t2-t1),d);
    }
    printf("**** \n");
  }
}
