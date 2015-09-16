#ifndef _BFM_COMM_SPI_BAGEL2_H_
#define _BFM_COMM_SPI_BAGEL2_H_


#include <utility>
#include <vector>
#include <algorithm>


#define NBRS (8)
#define INJ_MEMORY_FIFO_SIZE (2*1024)
#define NBALL (3*3*3*3)

#warning BFM USING SPI COMMUNICATIONS

// Must be called prior to MPI. Ideally make it a static object initialiser
int comm_spi_reserve_hardware(int log2mb) ;

// NB no free. Allocate hardware related stuff in reserve_hardware & save heap state
// Restore heap state each time we allocate volume dependent arrays in comm_init.
void * comm_spi_malloc(size_t bytes); 
void comm_spi_heap_save(void) ;
void comm_spi_heap_restore(void) ;

template <class Float>
class bfmcommspi : public bfmcommQMP<Float> {  // Inherit from  commQMP?
public:

  int allBytes;
  std::vector<std::pair<int,int> > haloes ;

  void comm_ldop_init(void) ;
  virtual void comm_init(void) ;
  virtual void comm_end(void) ;

  virtual void comm_start(int result_cb, Fermion_t psi,int dag) ;
  virtual void comm_complete(int result_cb, Fermion_t psi) ;

  void comm_start_precision_test(void);

  virtual void comm_spi_start(int mu) ;
  virtual void comm_spi_complete(void) ;
  virtual void comm_spi_barrier(void) ;
  virtual void comm_gsum(double &val) ;// override the gsum
  virtual void comm_gsum(double *val,int N) ;// override the gsum
  virtual int  SPIcomms(void) { return 1; };

};
#endif
