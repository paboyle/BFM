#include <bfm.h>
#include <stdio.h>
#include <stdlib.h>

#include <bfmcommiroiro.h>
#include <Communicator/communicator.hpp>

extern int time_report;

/*
 * BGNET/IroIro comms implementation
 */

template <class Float> 
void bfmcommIroIro<Float>::recv_end(void)
{
  for(int mu=0; mu < ND; mu++) { 
    if (!this->local_comm[mu]) { 
      if ( this->simd_grid[mu]>1 ){
	bfm_free(simd_rbuf[2*mu]);
	bfm_free(simd_rbuf[2*mu+1]);
      }
    }
  }
}
template <class Float> 
void bfmcommIroIro<Float>::recv_init(Fermion_t psi_t)
{
  Float *psi = (Float *) psi_t;  
  for(int mu=0;mu<4;mu++){
    for(int pm=0;pm<2;pm++ ) {

      int idx=2*mu+pm; 
      this->recvbufs[idx] = &psi[this->comm_offset[pm][mu]*SPINOR_SIZE*this->simd()*this->cbLs];

      if( !this->local_comm[mu] ) {
	
	// receive buffers
	// For a SIMD direction should reduce SIMD() factor....
	int simd_factor = this->simd();
	simd_factor = simd_factor/this->simd_grid[mu];

	int words  = this->simd_nbound[mu] * HALF_SPINOR_SIZE * simd_factor * this->cbLs;
	int delta;

	if ( pm == Plus ) delta = 1;
	else delta = -1;

	Float *rbuf = this->recvbufs[idx];
	if ( this->simd_grid[mu] > 1 ) {
	  rbuf       = (Float *)bfm_alloc(words*sizeof(Float));
	} 
	simd_rbuf[idx] = rbuf;

      } // if ( !local_comm )
    }      
    } // for(mu=... )

  return;

}
template <class Float> 
void bfmcommIroIro<Float>::comm_init(void)
{
                                                                            
  for(int mu=0;mu<4;mu++){
    for(int pm=0;pm<2;pm++ ) {

      int words, idx, alloc_words;
      Float * ptr;

      // Allocate send buffers
      idx=2*mu+pm; 

      int simd_factor = this->simd();
      simd_factor = simd_factor/this->simd_grid[mu];

      words = this->simd_nbound[mu] * HALF_SPINOR_SIZE * simd_factor * this->cbLs;
      alloc_words = this->simd_nbound[mu] * HALF_SPINOR_SIZE * this->simd() * this->cbLs;
      ptr= this->sendbufs[idx] = (Float *)bfm_alloc(alloc_words*sizeof(Float)); 

    }      
  } // for(mu=... )

  return;
}

template <class Float> 
void bfmcommIroIro<Float>::comm_end(void)
{
  /*---------------------------------------------------------------------*
   * I need to free the associated msgmem_t-s here.                      *
   *---------------------------------------------------------------------*/
  for(int mu=0; mu < ND; mu++) { 
    bfm_free(this->sendbufs[2*mu]);
    bfm_free(this->sendbufs[2*mu+1]);
  }
  return;
}

template <class Float> 
void bfmcommIroIro<Float>::comm  (int result_cb, Fermion_t psi,int dag)  
{
  comm_start(result_cb,psi,dag);
  comm_complete(result_cb,psi);
}


template <class Float> 
void bfmcommIroIro<Float>::comm_start  (int result_cb, Fermion_t psi,int dag)  
{
  // gather the faces. Routines here are threaded.
  // All threads cooperate in gathering.

  int me = this->thread_barrier();
  if ( me == 0 ) {
    recv_init((Fermion_t)psi);
  }
  this->thread_barrier();

  gather (result_cb,psi,dag);

  this->thread_barrier();
  if ( me == 0 ) {

    // Slow initial implementation 
    // Ideally insert BGNET calls and exit here after the BGNET_Put starts 
    // running asynchronously.
    for(int mu=0;mu<4;mu++){
      int simd_factor = this->simd();
      simd_factor = simd_factor/this->simd_grid[mu];
      int words  = this->simd_nbound[mu] * HALF_SPINOR_SIZE * simd_factor * this->cbLs;
      int idx=2*mu;
      int idxp=2*mu+1;
      Communicator::instance()->transfer_bk((double *)this->simd_rbuf[idx],(double *)this->sendbufs[idx],words,mu);
      Communicator::instance()->transfer_fw((double *)this->simd_rbuf[idxp],(double *)this->sendbufs[idxp],words,mu);
    }
  }
  this->thread_barrier();

  return;
}

template <class Float> 
void bfmcommIroIro<Float>::comm_merge  (void)  
{
  int mynode = Communicator::instance()->nodeid();

  for (int mu=0;mu<4;mu++) { 
    for(int pm=0;pm<2;pm++){
      int dir=pm+2*mu;
      if ( !this->local_comm[mu] && this->simd_grid[mu]>1 ) {
	int length = this->simd_nbound[mu]*this->cbLs;
	int skip   = length*12*sizeof(Float); // NB no nsimd
	Float *off_node = simd_rbuf[dir];
	Float * on_node = (Float *)((integer)this->sendbufs[dir]   + skip);

	if( pm==0) {
	  this->merge(this->recvbufs[dir],off_node,on_node,length);
	} else { 
	  this->merge(this->recvbufs[dir],on_node,off_node,length);
	}
      }
    }
  }
}

template <class Float> 
void bfmcommIroIro<Float>::comm_complete  (int result_cb, Fermion_t psi)  
{
  int me=this->thread_barrier();
  this->comm_merge();
  if ( me == 0) {
    this->recv_end();
  }
  return;
}

template <class Float> 
bool bfmcommIroIro<Float>::isBoss(void)
{
  return Communicator::instance()->primaryNode();
}

template <class Float> 
void bfmcommIroIro<Float>::comm_gsum(double &val) 
{
  int me = this->thread_barrier();
  double myval=0;
  if ( me == 0 ) { 
    myval=Communicator::instance()->reduce_sum(val);
  }
  this->thread_sum(myval,me);
  val = myval;
}

template class bfmcommIroIro <float>;
template class bfmcommIroIro <double>;

