#ifndef BFM_LINALG_H
#define BFM_LINALG_H

/*  2007
 *  Calls PAB's assembler routines to give an implementation
 *  of the dwf dslash. Aim is to give very high performance
 *  in a reasonably portable/retargettable library.
 */
#include "bfm.h"
#include <stdio.h>
#include <stdlib.h>
#ifdef BFM_G5SIGN
#define G5SIGN(A) (-A)
#else
#define G5SIGN(A) (A)
#endif

extern "C" { 
void   qpx_inner4(double *dot,int N, double *A,  double *B);
void   qpx_inner4_f(double *dot,int N,float *A, float *B);
void qpx_site_inner_dd (int N,double *inner, double *X,double *Y);
void qpx_site_inner_df (int N,double *inner, float  *X  , float * Y);
}


#ifndef BGQ
uint64_t GetTimeBase(void) {
  struct timeval tv;
  gettimeofday(&tv,NULL);
  return tv.tv_usec+1000*1000*(tv.tv_sec%1000);
}
#endif

template <class Float>
double bfmbase<Float>::axpby_norm(Fermion_t r, Fermion_t x, Fermion_t y,double a,double b)
{
  double nrm = 0;
  double *aa  = &a;
  double *bb  = &b;

  int me,thrlen,throff;
  thread_work(2*cbLs*simd_cbvol,me,thrlen,throff);// units of 6 Floats

  Float *xx = ((Float * )x) + throff*12*nsimd;
  Float *yy = ((Float * )y) + throff*12*nsimd;
  Float *rr = ((Float * )r) + throff*12*nsimd;

  if ( sizeof(Float) == sizeof(double) ) {
    vmx_vaxpby_norm((double *)rr,aa,bb,(double *)xx,(double *)yy,thrlen,&nrm);
  } else { 
    vmx_vaxpby_norm_s((float *)rr,aa,bb,(float *)xx,(float *)yy,thrlen,&nrm);
  }
  double dot = nrm;
  double mydot = dot;

  thread_sum(dot,me);

  if ( reproduce ) {
    repro(mydot,me);
    if ( reproduce_checksum ){
      uint64_t sum=checksum(r);
      if ( me == 0 ) repro(sum);
    }
    thread_barrier();
    if ( me == 0 ) repro(dot,nthread);
  }

  comm_gsum(dot);

  return dot;
}
template <class Float>
double bfmbase<Float>::axpy_norm(Fermion_t r, Fermion_t x, Fermion_t y,double a)
{
  double nrm=0;

  int me,thrlen,throff;
  thread_work(2*cbLs*simd_cbvol,me,thrlen,throff);// units of 6 Floats


  Float *xx = ((Float * )x) + throff*12*nsimd;
  Float *yy = ((Float * )y) + throff*12*nsimd;
  Float *rr = ((Float * )r) + throff*12*nsimd;

  if ( sizeof(Float) == sizeof(double) ) {
    vmx_vaxpy_norm((double *)rr,&a,(double *)xx,(double *)yy,thrlen,&nrm);
  } else { 
    vmx_vaxpy_norm_s((float *)rr,&a,(float *)xx,(float *)yy,thrlen,&nrm);
  }
  double dot = nrm;
  double mydot = dot;
  thread_sum(dot,me);

  if ( reproduce ) { 
    repro(mydot,me);
    if ( reproduce_checksum  ){
      uint64_t sum=checksum(r);
      if ( me == 0 ) repro(sum);
    }
    thread_barrier();
    if ( me == 0 ) repro(dot,nthread);
  }
  comm_gsum(dot);

  return dot;
}

template <class Float>
void bfmbase<Float>::axpibg5x(Fermion_t r, Fermion_t x, double a,double b)
{
  //r = a x + i b g5 x
  double aa,bb;

  int me,thrlen,throff;
  thread_work(cbLs*simd_cbvol,me,thrlen,throff);// units of 24 Floats
  Float *xx = ((Float * )x) + throff*24*nsimd;
  Float *rr = ((Float * )r) + throff*24*nsimd;

  aa=a;
  bb=G5SIGN(b);
  
  if ( sizeof(Float) == sizeof(double) ) {
    vmx_tmass((double *)rr,
               (double *)&aa,
               (double *)&bb,
               (double *)xx,
               thrlen);
  } else { 
    vmx_tmass_s((float *)rr,
		(double *)&aa,
		(double *)&bb,
		(float *)xx,
                thrlen);
  }
  thread_barrier();
}

template <class Float>
void bfmbase<Float>::merge(Float *r, Float *lo, Float *hi, int nsite)
{
  // nsite should be simd_nbound[mu]*cbLs
  int me,thrlen,throff;
  thread_work(nsite,me,thrlen,throff);// units of  Floats


  if ( sizeof(Float) == sizeof(double) ) {
    Float *llo = ((Float * )lo) + throff*6*nsimd;
    Float *hhi = ((Float * )hi) + throff*6*nsimd;
    Float *rr  = ((Float * )r)  + throff*12*nsimd;
    vmx_merge ((double *)rr,
               (double *)llo,
               (double *)hhi,
               thrlen);
#ifdef BGQ
  } else if ( this->precision_test )  {
    uint16_t *llo = ((uint16_t * )lo) + throff*6*nsimd;
    uint16_t *hhi = ((uint16_t * )hi) + throff*6*nsimd;
    uint16_t *rr  = ((uint16_t * )r)  + throff*12*nsimd;
    vmx_merge_hs((float *)rr,
		(float *)llo,
		(float *)hhi,
                thrlen);
#endif
  } else{ 
    Float *llo = ((Float * )lo) + throff*6*nsimd;
    Float *hhi = ((Float * )hi) + throff*6*nsimd;
    Float *rr  = ((Float * )r)  + throff*12*nsimd;
    vmx_merge_s((float *)rr,
		(float *)llo,
		(float *)hhi,
                thrlen);
  }

  thread_barrier();
}

template <class Float>
void bfmbase<Float>::axpby(Fermion_t r, Fermion_t x, Fermion_t y,double a,double b)
{
  // write_vaxpy 2 way unrolled operating on 3 colors * complex * nsimd
  double aa = a;
  double bb = b;

  int me,thrlen,throff;
  thread_work(2*cbLs*simd_cbvol,me,thrlen,throff);// units of 6 Floats

  Float *xx = ((Float * )x) + throff*12*nsimd;
  Float *yy = ((Float * )y) + throff*12*nsimd;
  Float *rr = ((Float * )r) + throff*12*nsimd;


  if ( sizeof(Float) == sizeof(double) ) {
    vmx_vaxpby((double *)rr,
               (double *)&aa,
               (double *)&bb,
               (double *)xx,
               (double *)yy,thrlen);
  } else { 
    vmx_vaxpby_s((float *)rr,
                      (double *)&aa,
                      (double *)&bb,
                      (float *)xx,
                      (float *)yy,thrlen);
  }

  thread_barrier();
}




template <class Float>
double bfmbase<Float>::cg_psi_p(Fermion_t psi, Fermion_t p, Fermion_t r,Fermion_t mmp,double a,double b)
{
  // 2 way unrolled operating on 3 colors * complex * nsimd
  double nrm = 0;

  int me,thrlen,throff;
  thread_work(2*cbLs*simd_cbvol,me,thrlen,throff);// units of 6 Floats

  Float *ppsi = ((Float * )psi) + throff*12*nsimd;
  Float *pp   = ((Float * )p)   + throff*12*nsimd;
  Float *rr   = ((Float * )r)   + throff*12*nsimd;
  Float *mmmp = ((Float * )mmp) + throff*12*nsimd;

  integer iargs[10];
  iargs[0] = (integer)&a;
  iargs[1] = (integer)&b;
  iargs[2] = (integer)ppsi;
  iargs[3] = (integer)pp;
  iargs[4] = (integer)rr;
  iargs[5] = (integer)mmmp;
  iargs[6] = (integer)thrlen;
  iargs[7] = (integer)&nrm;

  if ( sizeof(Float) == sizeof(double) ) {
    vmx_cg(iargs);
  } else { 
    vmx_cg_s(iargs);
  }
  double dot = nrm;
  thread_sum(dot,me);
  comm_gsum (dot);
  return(dot);
}

template <class Float>
void bfmbase<Float>::axpy(Fermion_t r, Fermion_t x, Fermion_t y,double a)
{
  int me,thrlen,throff;

  thread_work(2*cbLs*simd_cbvol,me,thrlen,throff);// units of 6 complex

  Float *xx = ((Float * )x) + throff*12*nsimd;
  Float *yy = ((Float * )y) + throff*12*nsimd;
  Float *rr = ((Float * )r) + throff*12*nsimd;

  if ( sizeof(Float) == sizeof(double) ) {
    vmx_vaxpy((double *)rr,(double *)&a,(double *)xx,(double *)yy,thrlen);
  } else { 
    vmx_vaxpy_s((float *)rr,(double *)&a,(float *)xx,(float *)yy,thrlen);
  } 
  thread_barrier();
  return;
}
template <class Float>
void bfmbase<Float>::caxpy(Fermion_t r, Fermion_t x, Fermion_t y,double a, double b)
{
  int me,thrlen,throff;

  thread_work(2*cbLs*simd_cbvol,me,thrlen,throff);// units of 6 complex

  Float *xx = ((Float * )x) + throff*12*nsimd;
  Float *yy = ((Float * )y) + throff*12*nsimd;
  Float *rr = ((Float * )r) + throff*12*nsimd;

  std::complex<double> aa(a,b);

  if ( sizeof(Float) == sizeof(double) ) {
    vmx_caxpy((double *)rr ,(double *)&aa,(double *)xx,(double *)yy,thrlen);
  } else { 
    vmx_caxpy_s((float *)rr,(double *)&aa,(float *)xx,(float *)yy,thrlen);
  } 
  thread_barrier();
  return;
}
template <class Float>
void bfmbase<Float>::caxpby(Fermion_t r, Fermion_t x, Fermion_t y,double a,double b,double c, double d)
{
  int me,thrlen,throff;
  thread_work(2*cbLs*simd_cbvol,me,thrlen,throff);// units of 6 Floats

  Float *xx = ((Float * )x) + throff*12*nsimd;
  Float *yy = ((Float * )y) + throff*12*nsimd;
  Float *rr = ((Float * )r) + throff*12*nsimd;
  
  std::complex<double> aa(a,b);
  std::complex<double> bb(c,d);

  if ( sizeof(Float) == sizeof(double) ) {
    vmx_caxpby((double *)rr,
               (double *)&aa,
               (double *)&bb,
               (double *)xx,
               (double *)yy,thrlen);
  } else { 
    vmx_caxpby_s((float *)rr,
                      (double *)&aa,
                      (double *)&bb,
                      (float *)xx,
                      (float *)yy,thrlen);
  }

  thread_barrier();
}
template <class Float>
std::complex<double> bfmbase<Float>::dot(Fermion_t x, Fermion_t y)
{
  std::complex<double> inner;

  int me,thrlen,throff;
  thread_work(2*cbLs*simd_cbvol,me,thrlen,throff);// units of 6 Floats


  Float *xx = ((Float * )x) + throff*12*nsimd;
  Float *yy = ((Float * )y) + throff*12*nsimd;

  if ( sizeof(Float) == sizeof(double) ) {
    vmx_inner((double *)xx,(double *)yy,thrlen,(double *)&inner);
  } else { 
    vmx_inner_s((float *)xx,(float *)yy,thrlen,(double *)&inner);
  }

  double rr = real(inner);
  thread_sum(rr,me);
  comm_gsum(rr);

  double ii = imag(inner);
  thread_sum(ii,me);
  comm_gsum(ii);

  std::complex<double> mdot(rr,ii);
  return mdot;
}
/* 
 * Performance non-critical
 */
template <class Float>
void bfmbase<Float>::scale(Fermion_t x_t, double a)
{
  int nspinco = 12;
  Float *x = (Float *) x_t;
  Float *xx;
  Float nrm = 0.0;

  int me,thrlen,throff;
  thread_work(cbLs*simd_cbvol*nsimd,me,thrlen,throff);

  for(int s=0;s<thrlen;s++){

    xx = &x[(s+throff)*nspinco*2 ];

    for(int spinco=0;spinco<nspinco;spinco++){
    for(int reim=0;reim<2;reim++){
      int idx = reim + spinco*2;
      xx[idx] = a *xx[idx];
    }}

  }  
  
  thread_barrier();
}


template <class Float>
void bfmbase<Float>::scale(Fermion_t x_t, double ar, double ai)
{
  int nspinco = 12;
  Float *x = (Float *) x_t;
  Float *xx;
  Float nrm = 0.0;

  int me,thrlen,throff;
  thread_work(cbLs*simd_cbvol*nsimd,me,thrlen,throff);

  for(int s=0;s<thrlen;s++){

    xx = &x[(s+throff)*nspinco*2 ];

    for(int spinco=0;spinco<nspinco;spinco++){
      int idx = spinco*2;
      double xr = xx[idx];
      double xi = xx[idx+1];
      xx[idx]   = ar*xr - ai*xi;
      xx[idx+1] = ai*xr + ar*xi;
    }

  }  
  thread_barrier();
}

template <class Float>
void bfmbase<Float>::copy(Fermion_t x_t, Fermion_t copyme)
{
  axpy(x_t,copyme,copyme,0.0);
}

 
template <class Float>
void bfmbase<Float>::set_zero(Fermion_t x_t)
{
  this->fill(x_t,0.0);
}

template <class Float>
complex<double> bfmbase<Float>::inner(Fermion_t x_t,Fermion_t y_t)
{
  return dot(x_t,y_t);
}

template <class Float>
double bfmbase<Float>::inner_real(Fermion_t x_t,Fermion_t y_t)
{
  std::complex<double> in = dot(x_t,y_t);
  return real(in);
}

template <class Float>
void bfmbase<Float>::fill(Fermion_t x_t, double a)
{
  int nspinco = 12;
  Float *x = (Float *) x_t;
  Float *xx;
  Float nrm = 0.0;

  int me,thrlen,throff;
  thread_work(cbLs*simd_cbvol*nsimd,me,thrlen,throff);

  for(int s=0;s<thrlen;s++){
    xx = &x[(s+throff)*nspinco*2 ];

    for(int spinco=0;spinco<nspinco;spinco++){
    for(int reim=0;reim<2;reim++){
      int idx = reim + spinco*2;
      xx[idx] = a ;
    }}
  }  
  
  thread_barrier();
}

template <class Float>
void bfmbase<Float>::master_fill(Fermion_t x_t, double a)
{
  int nspinco = 12;
  Float *x = (Float *) x_t;
  Float *xx;
  Float nrm = 0.0;

  int throff=0;
  int thrlen=cbLs*simd_cbvol*nsimd;
  int me=0;

  for(int s=0;s<thrlen;s++){
    xx = &x[(s+throff)*nspinco*2 ];

    for(int spinco=0;spinco<nspinco;spinco++){
    for(int reim=0;reim<2;reim++){
      int idx = reim + spinco*2;
      xx[idx] = a ;
    }}

  }  
  
}

template <class Float>
double bfmbase<Float>::local_norm(Fermion_t x_t)
{
  int nspinco = 12;

  Float *x = (Float *) x_t;
  Float *xx;
  double nrm = 0.0;

  int me,thrlen,throff;
  thread_work(cbLs*simd_cbvol*nsimd,me,thrlen,throff);
  for(int s=throff;s<throff+thrlen;s++){

    xx = &x[s*nspinco*2];

    for(int spinco=0;spinco<nspinco;spinco++){
    for(int reim=0;reim<2;reim++){
      int idx = reim + spinco*2;
      nrm = nrm+ xx[idx]*xx[idx];
    }}
  }
  double mydot = nrm;
  thread_sum(nrm,me);

  return nrm;
}

template <class Float>
double bfmbase<Float>::norm(Fermion_t x_t)
{
  int nspinco = 12;

  Float *x = (Float *) x_t;
  Float *xx;
  double nrm = 0.0;


  int me,thrlen,throff;
  thread_work(cbLs*simd_cbvol*nsimd,me,thrlen,throff);

  for(int s=throff;s<throff+thrlen;s++){

    xx = &x[s*nspinco*2];

    for(int spinco=0;spinco<nspinco;spinco++){
    for(int reim=0;reim<2;reim++){
      int idx = reim + spinco*2;
      nrm = nrm+ xx[idx]*xx[idx];
    }}
  }

  double mydot = nrm;
  thread_sum(nrm,me);

  if ( reproduce ) { 
    repro(mydot,me);
    if ( reproduce_checksum ) { 
      uint64_t sum=checksum(x_t);
      if ( me == 0 ) repro(sum);
    }
    thread_barrier();
    if ( me == 0 ) repro(nrm,nthread);
  }

  comm_gsum(nrm);

  return nrm;
}

template <class Float>
void bfmbase<Float>::conj(Fermion_t x_t)
{
  int nspinco = 12;
  Float *x = (Float *) x_t;
  Float *xx;
  Float nrm = 0.0;


  int me,thrlen,throff;
  thread_work(cbLs*simd_cbvol*nsimd,me,thrlen,throff);

  for(int s=0;s<thrlen;s++){

    xx = &x[(s+throff)*nspinco*2 ];

    for(int spinco=0;spinco<nspinco;spinco++){

      int idx_re = 0 + spinco*2;
      int idx_im = 1 + spinco*2;

      xx[idx_re] = xx[idx_re];
      xx[idx_im] = -xx[idx_im];
    
   }
  }  
  
  thread_barrier();
}

template <class Float>
uint64_t bfmbase<Float>::checksum(Fermion_t x_t)
{
  int nspinco = 12;

  Float *x = (Float *) x_t;
  Float *xx;
  union { double d ; uint64_t u ; } conv;
        
  uint64_t csum = 0;

  int me,thrlen,throff;
  thread_work(cbLs*simd_cbvol*nsimd,me,thrlen,throff);

  for(int s=throff;s<throff+thrlen;s++){

    xx = (Float *)&x[s*nspinco*2];

    for(int spinco=0;spinco<nspinco;spinco++){
    for(int reim=0;reim<2;reim++){
      int idx = reim + spinco*2;
      conv.d = xx[idx];
      csum = csum^conv.u;
    }}
  }

  thread_csum(csum,me);

  return csum;
}


template <class Float>
uint64_t bfmbase<Float>::checksumGauge(void)
{
  int ncoco = 9;
  // 8 links per site      => vol * 8
  // 2 sites on inner most => /2
  // dag+undag             => /2
  int vol = node_latt[0]*node_latt[1]*node_latt[2]*node_latt[3]*2;

  Float *xx;
  union { double d ; uint64_t u ; } conv;
        
  uint64_t csum = 0;
  uint64_t csumdag = 0;

  int me,thrlen,throff;
  thread_work(vol,me,thrlen,throff);

  for(int s=throff;s<throff+thrlen;s++){

    xx = (Float *)&u[s*ncoco*2*2*2];
    uint64_t gsum = 0;
    uint64_t gsumdag= 0;
    for(int coco=0;coco<ncoco;coco++){
    for(int reim=0;reim<4;reim++){
      int idx = reim + coco*4;
      conv.d = xx[idx];
      gsum = gsum^conv.u;
      conv.d = xx[idx+36];
      gsumdag = gsumdag^conv.u;
    }}
    csum = csum ^ gsum;
    csumdag = csumdag ^ gsumdag;
  }
  

  thread_csum(csum,me);
  thread_csum(csumdag,me);

  if ( csumdag != csum ){
    Error("Bad dag checksum comparison %16.16lx %16.16lx %lx\n",csum,csumdag, this);
  }

  return csum;
}

template <class Float>
void bfmbase<Float>::copy_slice(Fermion_t in_t ,int si,int Lsi,
				Fermion_t out_t,int so,int Lso)
{
  int nspinco = 12;
  Float *in = (Float *) in_t;
  Float *out= (Float *)out_t;
  Float *ii;
  Float *oo;

  if ( this->precon_5d ) { 
    Error("copy_slice does not work for 5d preconditioning\n");
    exit(0);
  }
  int me,thrlen,throff;
  thread_work(simd_cbvol,me,thrlen,throff);
  for(int site=0;site<thrlen;site++){

    // Nsimd x 24 x Ls x 4dSite ordering
    ii = & in[((site+throff)*Lsi +si)*nspinco*2*this->nsimd];
    oo = &out[((site+throff)*Lso +so)*nspinco*2*this->nsimd];

    for(int idx=0;idx<2*nspinco*this->nsimd;idx++){
      oo[idx]=ii[idx];
    }

  }  
  thread_barrier();
}
template <class Float>
double bfmbase<Float>::norm(Fermion_t x_t,int s)
{
  int nspinco = 12;
  Float *x = (Float *) x_t;
  Float *xx;
  int me,thrlen,throff;
  thread_work(simd_cbvol,me,thrlen,throff);
  double nrm=0;
  for(int site=0;site<thrlen;site++){

    // Nsimd x 24 x Ls x 4dSite ordering
    xx = &x[((site+throff)*this->cbLs +s)*nspinco*2*this->nsimd];

    for(int idx=0;idx<2*nspinco*this->nsimd;idx++){
      nrm+=xx[idx]*xx[idx];
    }

  }  
  thread_barrier();
  thread_sum(nrm,me);
  comm_gsum(nrm);
  thread_barrier();
  return nrm;

}

template <class Float>
void bfmbase<Float>::UniformRandom(Fermion_t x_t)
{
  int nspinco = 12;
  Float *x = (Float *) x_t;
  Float *xx;

  for(int site=0;site<simd_cbvol*cbLs;site++){

    // Nsimd x 24 x Ls x 4dSite ordering
    xx = &x[site*nspinco*2*this->nsimd];

    for(int idx=0;idx<2*nspinco*this->nsimd;idx++){
      xx[idx] = drand48()-0.5;
    }
  }  
}

template <class Float>
std::complex<double> bfmbase<Float>::dot(Fermion_t x_t,Fermion_t y_t,int sx,int sy)
{
  int nspinco = 12;
  Float *x = (Float *) x_t;
  Float *xx;
  Float *y = (Float *) y_t;
  Float *yy;
  int me,thrlen,throff;
  thread_work(simd_cbvol,me,thrlen,throff);
  double rr = 0.0;
  double ii = 0.0;
  for(int site=throff;site<throff+thrlen;site++){


    xx = &x[(site*this->cbLs +sx)*nspinco*2*this->nsimd];
    yy = &y[(site*this->cbLs +sy)*nspinco*2*this->nsimd];

    for(int idx=0;idx<2*nspinco*this->nsimd;idx+=2){ // ReIm x simd x 3 x 4
      rr = rr + xx[idx] *yy[idx]    + xx[idx+1] *yy[idx+1] ;
      ii = ii + xx[idx] *yy[idx+1]  - xx[idx]   *yy[idx+1] ;
    }
  }

  thread_sum(rr,me);
  comm_gsum(rr);

  thread_sum(ii,me);
  comm_gsum(ii);

  std::complex<double> mdot(rr,ii);
  return mdot;
}


template <class Float>
double bfmbase<Float>::normMat(Matrix_t x_t)
{
  int nmucoco = 9*4;
  Float *x = (Float *) x_t;
  Float *xx;
  int me,thrlen,throff;
  thread_work(simd_cbvol,me,thrlen,throff);
  double nrm=0;
  for(int site=0;site<thrlen;site++){

    // Nsimd x 36 x Ls x 4dSite ordering
    xx = &x[((site+throff))*nmucoco*2*this->nsimd];

    for(int idx=0;idx<2*nmucoco*this->nsimd;idx++){
      nrm+=xx[idx]*xx[idx];
    }

  }  
  thread_barrier();
  thread_sum(nrm,me);
  comm_gsum(nrm);
  thread_barrier();
  return nrm;

}

template <class Float>
void bfmbase<Float>::axpby_ssp_internal(Fermion_t out,
					double al,double ah,Fermion_t x_t,
					double bl,double bh,Fermion_t y_t,
					int sxo,int sy)
{

#if 1
  int nspinco = 12;
  Float *x = (Float *) x_t;
  Float *y = (Float *) y_t;
  Float *o = (Float *) out;

  // write_vaxpy 2 way unrolled operating on 3 colors * complex * nsimd

  int me,thrlen,throff;
  thread_work(simd_cbvol,me,thrlen,throff);

  int ssxo= sxo*24*nsimd;
  int ssy = sy *24*nsimd;

  Float *xx = ((Float * )x) + throff*24*nsimd*cbLs +ssxo;
  Float *yy = ((Float * )y) + throff*24*nsimd*cbLs +ssy;
  Float *oo = ((Float * )o) + throff*24*nsimd*cbLs +ssxo;

  double aa[4];
  aa[0] = al;
  aa[1] = ah;
  aa[2] = bl;
  aa[3] = bh;

  integer args[6];
  args[0]=(integer)oo;
  args[1]=(integer)xx;
  args[2]=(integer)yy;
  args[3]=thrlen;
  args[4]=cbLs*nsimd*24*sizeof(Float);
  args[5]=(integer)aa;

  if ( sizeof(Float) == sizeof(double) ) {
    vmx_vaxpby_ssp((integer *)args);
  } else { 
    vmx_vaxpby_ssp_s((integer *)args);
  }

  thread_barrier();

#else
  int nspinco = 12;
  Float *x = (Float *) x_t;
  Float *y = (Float *) y_t;
  Float *o = (Float *) out;
  Float *xx;
  Float *yy;
  Float *oo;

  int me,thrlen,throff;
  thread_work(simd_cbvol,me,thrlen,throff);

  int skip = 12*this->nsimd;

  for(int site=0;site<thrlen;site++){

    // Nsimd x 24 x Ls x 4dSite ordering
    xx = &x[((site+throff)*this->cbLs +sxo)*nspinco*2*this->nsimd];
    oo = &o[((site+throff)*this->cbLs +sxo)*nspinco*2*this->nsimd];
    yy = &y[((site+throff)*this->cbLs +sy )*nspinco*2*this->nsimd];

    for(int idx=0;idx<2*6*this->nsimd;idx++){
      oo[idx]      = al*xx[idx]+bl*yy[idx];
      oo[idx+skip] = ah*xx[idx+skip]+bh*yy[idx+skip];
    }

  }  

  thread_barrier();
#endif
}
template <class Float>
void bfmbase<Float>::axpby_ssp_proj   (Fermion_t out,double a,Fermion_t x,double b,Fermion_t y,int sxo,int sy,int psign)
{
  psign = G5SIGN(psign);
  if ( psign==1 )          axpby_ssp_internal(out,a,a,x,b,0,y,sxo,sy);
  else if ( psign==-1 )    axpby_ssp_internal(out,a,a,x,0,b,y,sxo,sy);
  else { 
    Error("axpby_ssp_proj bad psign\n");
    exit(-1);
  }

}

template <class Float>
void bfmbase<Float>::ag5xpby_ssp_proj   (Fermion_t out,double a,Fermion_t x,double b,Fermion_t y,int sxo,int sy,int psign)
{
  psign = G5SIGN(psign);
  if ( psign==1 )          axpby_ssp_internal(out,G5SIGN(a),G5SIGN(-a),x,b,0,y,sxo,sy);
  else if ( psign==-1 )    axpby_ssp_internal(out,G5SIGN(a),G5SIGN(-a),x,0,b,y,sxo,sy);
  else { 
    Error("ag5xpby_ssp_proj bad psign\n");
    exit(-1);
  }
}

template <class Float>
void bfmbase<Float>::axpby_ssp   (Fermion_t out,double a,Fermion_t x,double b,Fermion_t y,int sxo,int sy)
{
  axpby_ssp_internal(out,a,a,x,b,b,y,sxo,sy);
}
template <class Float>
void bfmbase<Float>::ag5xpby_ssp   (Fermion_t out,double a,Fermion_t x,double b,Fermion_t y,int sxo,int sy)
{
  axpby_ssp_internal(out,G5SIGN(a),G5SIGN(-a),x,b,b,y,sxo,sy);
}
template <class Float>
void bfmbase<Float>::axpbg5y_ssp   (Fermion_t out,double a,Fermion_t x,double b,Fermion_t y,int sxo,int sy)
{
  axpby_ssp_internal(out,a,a,x,G5SIGN(b),G5SIGN(-b),y,sxo,sy);
}

template <class Float>
void bfmbase<Float>::ag5xpbg5y_ssp (Fermion_t out,double a,Fermion_t x,double b,Fermion_t y,int sxo,int sy)
{
  axpby_ssp_internal(out,G5SIGN(a),G5SIGN(-a),x,G5SIGN(b),G5SIGN(-b),y,sxo,sy);

}

template <class Float>
void bfmbase<Float>::four_to_five(Fermion_t input[2], Fermion_t result5d[2]){
  Error("four_to_five not implmementd\n");
  exit(0);
}
template <class Float>
void bfmbase<Float>::five_to_four(Fermion_t input5d[2], Fermion_t result[2], int sgn){
  Error("five_to_four not implmementd\n");
  exit(0);
}


#endif
