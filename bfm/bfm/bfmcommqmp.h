#ifndef _BFM_COMM_QMP_BAGEL2_H_
#define _BFM_COMM_QMP_BAGEL2_H_
#include <qmp.h>
template <class Float>
class bfmcommQMP : public bfmbase<Float> { 
public:

  // QMP thingy-me-bob's
  QMP_msghandle_t multi_handle[2];

  QMP_msgmem_t send_ops_msgmem_t[8];
  QMP_msgmem_t recv_ops_msgmem_t[8];

  QMP_msghandle_t send_multi_handle;
  QMP_msghandle_t send_handles[8];
 

  QMP_msghandle_t recv_handles[8];
  QMP_msghandle_t all_handles;

  Float *simd_rbuf[8];
  Float *receive_area;

  int num_op;

  virtual bool isBoss() ;


  virtual void recv_init(void);
  virtual void recv_end(void) ;

  virtual void comm_init(void) ;
  virtual void comm_end(void) ;
  virtual void comm(int result_cb, Fermion_t psi,int dag) ;
  virtual void comm_start(int result_cb, Fermion_t psi,int dag) ;
  virtual void comm_qmp_start(Fermion_t psi) ;
  virtual void comm_qmp_complete(void) ;
  virtual void comm_merge(void);
  virtual void comm_complete(int result_cb, Fermion_t psi) ;
  virtual void comm_gsum(double &val) ;
  virtual void comm_gsum(double *val,int N) ;

};
#endif
