#include <bfm.h>
#include <stdio.h>
#include <stdlib.h>

#include <bfmcommqmp.h>
#include <qmp.h>

extern int time_report;


/*
 * QMP comms implementation
 * If ( local_comm[mu]== 1 ) these routines do not get called
 * If ( local_comm[mu]== 0 ) these routines get called but loop back
 *
 * For now a brain dead implementation with no spin projection.
 */

#define recv_multi_handle multi_handle[0]
#define send_multi_handle multi_handle[1]

template <class Float> 
void bfmcommQMP<Float>::recv_end(void)
{
  if ( num_op>0 ) {
    QMP_free_msghandle(recv_multi_handle);
  }
  for(int mu=0; mu < ND; mu++) { 
    if (!this->local_comm[mu]) { 
      QMP_free_msgmem(recv_ops_msgmem_t[2*mu]);
      QMP_free_msgmem(recv_ops_msgmem_t[2*mu+1]);
      if ( this->simd_grid[mu]>1 ){
	bfm_free(simd_rbuf[2*mu]);
	bfm_free(simd_rbuf[2*mu+1]);
      }
    }
  }
  bfm_free(receive_area);
}
template <class Float> 
void bfmcommQMP<Float>::recv_init(void)
{
  int words = this->simd_allbound*SPINOR_SIZE*this->simd()*this->cbLs*2;
  receive_area = bfm_alloc(words * sizeof(Float));
  num_op = 0;

  for(int mu=0;mu<4;mu++){
    for(int pm=0;pm<2;pm++ ) {

      int idx=2*mu+pm; 
      
      this->recvbufs[idx] = &receive_area[this->comm_offset[pm][mu]*SPINOR_SIZE*this->simd()*this->cbLs];

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
	  rbuf       = simd_rbuf[idx] = (Float *)bfm_alloc(words*sizeof(Float));
	}
	recv_ops_msgmem_t[idx]  = QMP_declare_msgmem((void *)rbuf,words*sizeof(Float));
	recv_handles[num_op++]  = QMP_declare_receive_relative(recv_ops_msgmem_t[idx], mu, delta, 0);

      } // if ( !local_comm )
    }      
    } // for(mu=... )
    
    if( num_op > 0 ) { 
      recv_multi_handle = QMP_declare_multiple(recv_handles, num_op);
    }

  return;

}
template <class Float> 
void bfmcommQMP<Float>::comm_init(void)
{
                                                                            
  num_op = 0;
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
      if( !this->local_comm[mu] ) {

	// Declare send & msgmem
	int delta;
	if ( pm == Plus ) delta = 1;
	else delta = -1;

	// We pull the "Plus" face_table and send in the Minus direction
	// This reflects an oddity of the A
	send_ops_msgmem_t[idx] = QMP_declare_msgmem((void *)ptr,words*sizeof(Float));
	send_handles[num_op++] = QMP_declare_send_relative(send_ops_msgmem_t[idx], mu, -delta, 0);

      } 
    }      
    } // for(mu=... )
    
    if( num_op > 0 ) { 
      send_multi_handle = QMP_declare_multiple(send_handles, num_op);
    }

    recv_init();

  return;
}

template <class Float> 
void bfmcommQMP<Float>::comm_end(void)
{
  /*-----------------------------------------------------------------------* 
   * I need to free the multiple msghandles                                *
   * The individual handles (tmp_handles) used should have been freed      *
   * when the multiple was created and according to the spec I don't need  *
   * to free them                                                          *
   *---------------------------------------------------------------------- */
  if ( num_op>0 ) QMP_free_msghandle(send_multi_handle);
  /*---------------------------------------------------------------------*
   * I need to free the associated msgmem_t-s here.                      *
   *---------------------------------------------------------------------*/
  for(int mu=0; mu < ND; mu++) { 
    bfm_free(this->sendbufs[2*mu]);
    bfm_free(this->sendbufs[2*mu+1]);
    if (!this->local_comm[mu]) { 
      QMP_free_msgmem(send_ops_msgmem_t[2*mu]);
      QMP_free_msgmem(send_ops_msgmem_t[2*mu+1]);
    }
  }

  return;
}

template <class Float> 
void bfmcommQMP<Float>::comm  (int result_cb, Fermion_t psi,int dag)  
{
  comm_start(result_cb,psi,dag);
  comm_complete(result_cb,psi);
}


template <class Float> 
void bfmcommQMP<Float>::comm_start  (int result_cb, Fermion_t psi,int dag)  
{
  // gather the faces. Routines here are threaded.
  // All threads cooperate in gathering.
  

  //  if ( this->isBoss() && (me == 0)) {
  //    this->Error("comm_start checking heap\n"); fflush(stdout);
  //    mcheck_check_all();
  //    this->Error("comm_start mcheck_all");fflush(stdout);
  //  }

  int me = this->thread_barrier();
  
  this->thread_barrier();

  gather (result_cb,psi,dag);

  this->thread_barrier();
  if ( me == 0 ) {
    comm_qmp_start(psi);
  }
  this->thread_barrier();

  return;
}
template <class Float> 
void bfmcommQMP<Float>::comm_qmp_start(Fermion_t psi) 
{
    if ( num_op > 0 ) {
      QMP_start(recv_multi_handle);
      QMP_start(send_multi_handle);
    }
}
template <class Float> 
void bfmcommQMP<Float>::comm_qmp_complete(void) 
{
  if ( num_op > 0 ) { 
    
    if( QMP_wait(recv_multi_handle) != QMP_SUCCESS ) { 
      if ( isBoss() ) this->Error("QMP_wait returned unsuccessful value: \n"); fflush(stdout);
      exit(-1);
    }
    
    if( QMP_wait(send_multi_handle) != QMP_SUCCESS ) { 
      if ( isBoss() ) this->Error("QMP_wait returned unsuccessful value: \n");fflush(stdout);
      exit(-1);
    }

  }
}

template <class Float> 
void bfmcommQMP<Float>::comm_merge  (void)  
{
  for (int mu=0;mu<4;mu++) { 
    for(int pm=0;pm<2;pm++){
      int dir=pm+2*mu;
      if ( !this->local_comm[mu] && this->simd_grid[mu]>1 ) {
	int length = this->simd_nbound[mu]*this->cbLs;
	int skip   = length*12*sizeof(Float); // NB no nsimd
	//	if (this->precision_test) skip = length*12*2;// 2byte truncation ??
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
void bfmcommQMP<Float>::comm_complete  (int result_cb, Fermion_t psi)  
{
  uint64_t t1 = GetTimeBase();
  int me=this->thread_barrier();
  if ( me == 0) { 
    comm_qmp_complete() ;
  }
  this->thread_barrier();
  uint64_t t2 = GetTimeBase();

  this->comm_merge();
  uint64_t t3 = GetTimeBase();

  this->thread_barrier();

  uint64_t t4 = GetTimeBase();

  if ( this->iter == this->time_report_iter+1 ) {
    this->ThreadBossMessage("qmp_complete %ld cyc\n",t2-t1);
    this->ThreadBossMessage("merge        %ld cyc\n",t3-t2);
  }
  this->thread_barrier();
  
  return;
}

template <class Float> 
bool bfmcommQMP<Float>::isBoss(void)
{
  return QMP_is_primary_node() ;
}

template <class Float> 
void bfmcommQMP<Float>::comm_gsum(double *val,int N) 
{
  for(int i=0;i<N;i++){
    comm_gsum(val[i]);
  }
}
template <class Float> 
void bfmcommQMP<Float>::comm_gsum(double &val) 
{
  int me = this->thread_barrier();
  double myval=0;
  if ( me == 0 ) { 
    QMP_sum_double_array(&val,1);
    myval = val;
  }
  this->thread_sum(myval,me);
  val = myval;
}

template class bfmcommQMP <float>;
template class bfmcommQMP <double>;

