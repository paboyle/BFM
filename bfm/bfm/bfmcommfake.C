/****************************************************************************/
/* PAB Jul 2007                                                             */
/****************************************************************************/

#include <stdlib.h>
#include <stdio.h>

#include <bfm.h>
#include <bfmcommfake.h>

template<class Float>
void bfmcommfake<Float>::comm_init (void)
{
  int words = this->simd_allbound*SPINOR_SIZE*this->simd()*this->cbLs*2;
  receive_area = bfm_alloc(words * sizeof(Float));

  for(int mu=0;mu<4;mu++){
    for(int pm=0;pm<2;pm++ ) {

      int idx=2*mu+pm; 
      // receive buffers
      // For a SIMD direction should reduce SIMD() factor....
      int simd_factor = this->simd();
      simd_factor = simd_factor/this->simd_grid[mu];
      
      this->recvbufs[idx] = &receive_area[this->comm_offset[pm][mu]*SPINOR_SIZE*this->simd()*this->cbLs];
#if 0
      printf("recvbufs[%d] = %lx - %lx\n",idx,this->recvbufs[idx],&this->recvbufs[idx][ this->simd_nbound[mu] * HALF_SPINOR_SIZE * simd_factor * this->cbLs]);
#endif
      if( !this->local_comm[mu] ) {
	

	int words  = this->simd_nbound[mu] * HALF_SPINOR_SIZE * simd_factor * this->cbLs;
	int delta;

	if ( pm == Plus ) delta = 1;
	else delta = -1;

	Float *rbuf = this->recvbufs[idx];
	if ( this->simd_grid[mu] > 1 ) {
	  rbuf       = simd_rbuf[idx] = (Float *)bfm_alloc(words*sizeof(Float));
	}
      } // if ( !local_comm )
    }      
    } // for(mu=... )

  return;
}
template<class Float>
void bfmcommfake<Float>::comm_end (void)
{
  for(int mu=0; mu < ND; mu++) { 
    if (!this->local_comm[mu]) { 
      if ( this->simd_grid[mu]>1 ){
	bfm_free(simd_rbuf[2*mu]);
	bfm_free(simd_rbuf[2*mu+1]);
      }
    }
  }
  bfm_free(receive_area);
}

// cb is cb of result, 1-cb is cb of input field, and it is this
// we communicate
template<class Float>
void bfmcommfake<Float>::comm_start(int result_cb, Fermion_t psi,int dag)
{
  comm(result_cb,psi,dag);
}
template<class Float>
void bfmcommfake<Float>::comm_complete(int result_cb, Fermion_t psi )
{
}
template<class Float>
void bfmcommfake<Float>::comm(int result_cb, Fermion_t psi, int dag)
{

  int nspinco = 12;
  int cb = 1-result_cb;

  int me = this->thread_barrier();
  for(int mu=0;mu<ND;mu++){
    for(int pm=0;pm<2;pm++ ) {

      int ftcb =0;
      int dir=pm+2*mu;
      Float *comm_buf = this->sendbufs[dir];

      int sgn =pm;
      if ( dag ) sgn = 1-pm;

      if ( this->precon_5d ) ftcb = 0;
      else ftcb = result_cb;
      int idx=2*mu+pm; 
      
      Float *psi_t = (Float *) psi;  

      if ( (!this->local_comm[mu])  ) {
	Error("comm fake is broken; should use gather_proj\n");
	exit(-1);
      } 

      if (this->simd_grid[mu]>1) {

	int thrvol,throff;
	this->thread_work_nobarrier(this->simd_nbound[mu],me,thrvol,throff);

	integer args[9];
	args[0] = (integer) &this->face_table[ftcb][pm][mu][throff];
	args[1] = (integer) psi;
	args[2] = (integer) throff;
	args[3] = (integer) this->recvbufs[dir];
	args[4] = (integer) thrvol;
	args[5] = (integer) this->cbLs;
	args[6] = sgn+2*mu;
	args[7] = (integer) this->complex_i_simd;

	integer arg_p = (integer) args;

	if ( sizeof(Float) == sizeof(double ) ) vmx_gather_proj_perm(arg_p);
	else                                    vmx_gather_proj_perm_s(arg_p);

#if 0
	Float *ptr = this->recvbufs[dir];
	for(int b=0;b<this->simd_nbound[mu];b++){
	  for(int sc=0;sc<6;sc++){
	    printf("b=%d sc=%d site %d",b,sc,this->face_table[ftcb][pm][mu][b]);
	    for(int rr=0;rr<2*this->simd();rr++){
	      printf(" %f ",*(ptr++));
	    }
	    printf("\n");
	  }
	}
#endif
      }
    }
  }

}

template class bfmcommfake<float>;
template class bfmcommfake<double>;
