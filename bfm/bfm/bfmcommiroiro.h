#ifndef _BFM_COMM_IROIRO_BAGEL2_H_
#define _BFM_COMM_IROIRO_BAGEL2_H_


template <class Float>
class bfmcommIroIro : public bfmbase<Float> { 
public:

  Float *simd_rbuf[8];

  virtual bool isBoss() ;

  virtual void recv_init(Fermion_t psi);
  virtual void recv_end(void) ;
  virtual void comm_init(void) ;
  virtual void comm_end(void) ;
  virtual void comm(int result_cb, Fermion_t psi,int dag) ;
  virtual void comm_start(int result_cb, Fermion_t psi,int dag) ;
  virtual void comm_merge(void);
  virtual void comm_complete(int result_cb, Fermion_t psi) ;
  virtual void comm_gsum(double &val) ;

};
#endif
