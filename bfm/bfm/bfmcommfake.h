#ifndef BFM_COMM_FAKE_BAGEL2_H
#define BFM_COMM_FAKE_BAGEL2_H

template <class Float>
class bfmcommfake : public bfmbase<Float> { 
public:

  Float *simd_rbuf[8];
  Float *receive_area;

  virtual void comm_init(void) ;
  virtual void comm_end(void)  ;
  virtual void comm(int result_cb, Fermion_t psi,int dag) ; 
  virtual void comm_start(int result_cb, Fermion_t psi, int dag) ; 
  virtual void comm_complete(int result_cb, Fermion_t psi) ; 

};

#endif
