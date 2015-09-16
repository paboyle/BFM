
/*Peter Boyle Dec 2012*/
#include "bfm.h"
#include <omp.h>

#include <stdio.h>
#include <stdlib.h>

#define ALIGNIT(A) (double *)( (((uint64_t)A) + 31)& (~0x1FUL) );

extern "C" {
  double qpx_site_norm  (int N, double *  A);

  void   qpx_site_zaxpy (int N,double *  R,  double *  X,  double * Y,double *c);
  void   qpx_site_zscale(int N,double *  R,  double *  X, double *c);

  void   qpx_site_zaxpy_f (int N,float *  R,  float  *  X, float * Y,double *c);
  void   qpx_site_zscale_f(int N,float *  R,  float  *  X, double *c);

  void qpx_site_inner_dd (int N,double *inner, double *X,double *Y);
  void qpx_site_inner_df (int N,double *inner, float *X,float *Y);
  void qpx_site_inner_fd (int N,float *inner, double *X,double *Y);
  void qpx_site_inner_ff (int N,float *inner, float *X,float *Y);

}

template<class Float>
void bfm_site_zaxpy(int N, Float *  R,  Float *  X,  Float * Y,double *c)  __attribute__ ((optimize("-O6")));
template<class Float, class cFloat>
void bfm_site_inner (int N,cFloat *inner, Float *X,Float *Y)   __attribute__ ((optimize("-O6")));
template<class Float>
void bfm_site_zaxpy(int N, Float *  R,  Float *  X,  Float * Y,double *c)
{
  if ( sizeof(Float) == sizeof(double) ) {
    qpx_site_zaxpy(N,(double *)R,(double *)X,(double *)Y,c);
  } else { 
    qpx_site_zaxpy_f(N,(float *)R,(float *)X,(float *)Y,c);
  }
}


template<class Float, class cFloat>
void bfm_site_inner (int N,cFloat *inner, Float *X,Float *Y)
{
  if(sizeof(Float)==sizeof(double)){
    if ( sizeof(cFloat)==sizeof(double) ) qpx_site_inner_dd (N,(double *)inner,(double *)X,(double *)Y);
    else                                  qpx_site_inner_fd (N,(float *) inner,(double *)X,(double *)Y);
  } else {
    if ( sizeof(cFloat)==sizeof(double) ) qpx_site_inner_df (N,(double *)inner,(float *)X,(float *)Y);
    else                                  qpx_site_inner_ff (N,(float *) inner,(float *)X,(float *)Y);
  }
}
template<class Float>
void   bfm_site_zscale(int N,Float *  R,  Float  *  X, double *c)
{
  if(sizeof(Float)==sizeof(double)){
    qpx_site_zscale(N,(double *)R,(double *)X,c);
  } else {
    qpx_site_zscale_f(N,(float *)R,(float *)X,c);
  }
}


template<class Float>
void bfmbase<Float>::block_pick (Fermion_t out,Fermion_t match,Fermion_t nomatch,int pickblock,int *blockids)
{
  int Nspinco=12;

  Float *matchp   = (Float *)match;
  Float *nomatchp = (Float *)nomatch;
  Float *outp     = (Float *)out;

  int me,thrlen,throff;
  thread_work(this->node_cbvol,me,thrlen,throff);
  for (int cbsite=throff;cbsite<throff+thrlen;cbsite++ ) { 
    for (int s=0;s<this->cbLs;s++ ) { 

	int cidx = bagel_idx5d_tmp(cbsite,s,0,0,1); cidx=cidx>>1;
	int block  = blockids[cidx];

	for ( int co=0;co<Nspinco;co++ ) { 
	  int bidx = bagel_idx5d_tmp(cbsite,s,0,co,Nspinco);
	  if ( block == pickblock ) { 
	    outp[bidx]   = matchp[bidx];
	    outp[bidx+1] = matchp[bidx+1];
	  } else { 
	    outp[bidx]   = nomatchp[bidx];
	    outp[bidx+1] = nomatchp[bidx+1];
	  }
	}
    }
  }
  thread_barrier();
}

template<class Float>
void bfmbase<Float>::block_inner(Fermion_t in1,Fermion_t in2,int nblock,double *result,double *reduce,int *blockids)
{
  block_inner(in1,in2,nblock,result,reduce,blockids,0,this->cbLs-1);
}
template<class Float>
template<class cFloat>
void bfmbase<Float>::block_project(Fermion_t *basis,int nbasis,Fermion_t in,
				   int nblock,cFloat *result,
				   double *reduce,int *blockids,int smin,int smax)
{
  int Nspinco=12;
  Float *in2p = (Float *)in;
  int me,thrlen,throff;

  //////////////////
  // Lots of thread contract for each 4d site
  //////////////////
  if ( ((uint64_t)reduce) & 0x1F ) {
    Error("Bfm::block_inner -- Unaligned pointer\n");
    exit(0);
  }

  int fifth = smax-smin+1;
  if ( smin > smax ) {
    fifth = smax+1+ this->cbLs-smin;
  }
  uint64_t bytes = nbasis*this->node_cbvol*24*sizeof(Float)*fifth;

  int unit=this->simd();
  thread_work(this->node_cbvol/unit,me,thrlen,throff);
  thrlen*=unit;
  throff*=unit;

  //  static int print;
  uint64_t t0 = GetTimeBase();
  for (int cbsite=throff;cbsite<thrlen+throff;cbsite+=this->simd() ) { 
    
    for(int vec=0;vec<nbasis;vec++){

      int vo = vec*this->node_cbvol*2;
      Float *in1p = (Float *)basis[vec];

      for(int v=0;v<2*this->simd();v++) reduce[vo+cbsite*2+v]=0;
    
      if ( smin > smax ) { 
	int s,bidx,sslices;
	sslices=this->cbLs-smin;
	bidx = bagel_idx5d_tmp(cbsite,smin,0,0,Nspinco);
	bfm_site_inner<Float>(sslices*Nspinco,(double *)&reduce[vo+cbsite*2],&in1p[bidx],&in2p[bidx]);
	s=0; sslices=smax+1;
	bidx = bagel_idx5d_tmp(cbsite,s,0,0,Nspinco);
	bfm_site_inner<Float>(sslices*Nspinco,(double *)&reduce[vo+cbsite*2],&in1p[bidx],&in2p[bidx]);
	//	print=1;
      } else { 
	int bidx = bagel_idx5d_tmp(cbsite,smin,0,0,Nspinco);
	int sslices=smax-smin+1;
	__builtin_prefetch(&in1p[bidx],1);
	__builtin_prefetch(&in2p[bidx],1);
	__builtin_prefetch(&reduce[vo+cbsite*2],1);
	bfm_site_inner<Float>(sslices*Nspinco,
			      (double *)&reduce[vo+cbsite*2],&in1p[bidx],&in2p[bidx]);
      }
    }
  }

  thread_barrier();
  uint64_t t1 = GetTimeBase();

  //////////////////////////
  // Divide blockid's across threads
  //////////////////////////
  thread_work(nblock,me,thrlen,throff);

  for(int vec=0;vec<nbasis;vec++){
    for(int b=throff;b<throff+thrlen;b++){
      result[vec*2+b*2*nbasis] = result[vec*2+b*2*nbasis+1] = 0.0;
    }
  }
  for (int cbsite=0;cbsite<this->node_cbvol;cbsite++ ) { 
    int cidx = bagel_idx5d_tmp(cbsite,0,0,0,1); cidx=cidx>>1;
    int b    = blockids[cidx];
    if ( (b>=throff)&&(b<throff+thrlen) ) {
      for(int vec=0;vec<nbasis;vec++){
	int vo = vec*this->node_cbvol*2;
	result[vec*2+b*2*nbasis]   += reduce[vo+cbsite*2];
	result[vec*2+b*2*nbasis+1] += reduce[vo+cbsite*2+1];
      }
    }  
  }
  thread_barrier();

  uint64_t t2 = GetTimeBase();
  /*
  if ( print < 20 && !me ) { 
    uint64_t micsec = (t2-t0)/MHz();
    double MBps =1.0*bytes/micsec;
    BossLog("Project<%d> %d/%d/%d microsecs smin,smax=%d,%d %f MB/s\n",sizeof(Float),(t1-t0)/MHz(),(t2-t1)/MHz(),(t2-t0)/MHz(),smin,smax,MBps);
    print++;
  }
  */
}

template<class Float>
void bfmbase<Float>::block_inner(Fermion_t in1,Fermion_t in2,int nblock,double *result,double *reduce,int *blockids,int smin,int smax)
{
  int Nspinco=12;
  Float *in1p = (Float *)in1;
  Float *in2p = (Float *)in2;

  int me,thrlen,throff;

  //////////////////
  // Lots of thread contract for each 4d site
  //////////////////
  if ( ((uint64_t)reduce) & 0x1F ) {
    Error("Bfm::block_inner -- Unaligned pointer\n");
    exit(0);
  }

  int unit = this->simd();
  thread_work(this->node_cbvol/unit,me,thrlen,throff);
  thrlen*=unit;
  throff*=unit;
  
  uint64_t t1 = GetTimeBase();
  for (int cbsite=throff;cbsite<thrlen+throff;cbsite+=unit ) { 

    for(int v=0;v<2*this->simd();v++) reduce[cbsite*2+v]=0;

    if ( smin > smax ) { 

      int s,bidx,sslices;
      sslices=this->cbLs-smin;
      bidx = bagel_idx5d_tmp(cbsite,smin,0,0,Nspinco);
      bfm_site_inner<Float>(sslices*Nspinco,(double *)&reduce[cbsite*2],&in1p[bidx],&in2p[bidx]);
 
      s=0; sslices=smax+1;
      bidx = bagel_idx5d_tmp(cbsite,s,0,0,Nspinco);
      bfm_site_inner<Float>(sslices*Nspinco,(double *)&reduce[cbsite*2],&in1p[bidx],&in2p[bidx]);

    } else { 
      int bidx = bagel_idx5d_tmp(cbsite,smin,0,0,Nspinco);
      int sslices=smax-smin+1;
      __builtin_prefetch(&in2p[bidx],1);
      __builtin_prefetch(&in1p[bidx],1);
      __builtin_prefetch(&reduce[cbsite*2],1);
      bfm_site_inner<Float>(sslices*Nspinco,(double *)&reduce[cbsite*2],&in1p[bidx],&in2p[bidx]);
    }
    
  }
  thread_barrier();
  uint64_t t2 = GetTimeBase();

  //////////////////////////
  // Divide blockid's across threads
  //////////////////////////
  thread_work(nblock,me,thrlen,throff);
  for(int b=throff;b<throff+thrlen;b++){
    result[b*2] = result[b*2+1] = 0.0;
  }
  for (int cbsite=0;cbsite<this->node_cbvol;cbsite++ ) { 
    int cidx = bagel_idx5d_tmp(cbsite,0,0,0,1); cidx=cidx>>1;
    int b    = blockids[cidx];
    if ( (b>=throff)&&(b<throff+thrlen) ) {
      result[b*2]   += reduce[cbsite*2];
      result[b*2+1] += reduce[cbsite*2+1];
    }
  }  
  thread_barrier();
  uint64_t t3 = GetTimeBase();
}

template<class Float>
template<class cFloat>
void   bfmbase<Float>::block_zaxpy(Fermion_t result,Fermion_t x,Fermion_t y, double a, cFloat *zap,int *blockids)
{
  block_zaxpy(result,x,y,a,zap,blockids,0,this->cbLs-1);
}

template<class Float>
template<class cFloat>
void bfmbase<Float>::block_promote(Fermion_t result,Fermion_t *basis, int nbasis,cFloat *zap,int *blockids,int smin,int smax)
{
  int Nspinco=12;
  Float *rp = (Float *)result;

  int me,thrlen,throff;
  int unit=this->simd();
  thread_work(this->node_cbvol/unit,me,thrlen,throff);
  thrlen*=unit;
  throff*=unit;
  
  int fifth = smax-smin+1;
  if ( smin > smax ) {
    fifth = smax+1+ this->cbLs-smin;
  }

  uint64_t bytes = nbasis*this->node_cbvol*fifth*24*sizeof(Float);

  std::vector<double> zz(2*this->simd());

  if( (((uint64_t)&zz[0])&0x1F) 
  ||( ((uint64_t)rp)&0x1F ) 
  ||( thrlen & 0x1 )        ) { 
    Error("bfm::block_promote -- unaligned pointer %lx %lx\n",&zz[0],rp);
    exit(0);
  }

  //  static int print;
  uint64_t t0 = GetTimeBase();
  for (int cbsite=throff;cbsite<throff+thrlen;cbsite+=this->simd() ) { 

    int cidx = bagel_idx5d_tmp(cbsite,smin,0,0,1); cidx=cidx>>1;

    for(int vec=0;vec<nbasis;vec++) {

      Float *xp = (Float *)basis[vec];

      for(int v=0;v<this->simd();v++){
	zz[2*v+0] = zap[2*blockids[cidx+v]*nbasis+vec*2+0];
	zz[2*v+1] = zap[2*blockids[cidx+v]*nbasis+vec*2+1];
      }

      if ( smin > smax ) {
	int bidx = bagel_idx5d_tmp(cbsite,smin,0,0,Nspinco);
	int sslices=this->cbLs-smin;
	//	print=1;
	if ( vec == 0 ) {
	  bfm_site_zscale<Float>(sslices*Nspinco,&rp[bidx],&xp[bidx],(double *)&zz[0]);
	} else {
	  bfm_site_zaxpy<Float>(sslices*Nspinco,&rp[bidx],&xp[bidx],&rp[bidx],(double *)&zz[0]);
	}
	int s=0;
	bidx = bagel_idx5d_tmp(cbsite,s,0,0,Nspinco);
	sslices=smax+1;
	bfm_site_zaxpy(sslices*Nspinco,&rp[bidx],&xp[bidx],&rp[bidx],(double *)&zz[0]);
      } else {
	int bidx = bagel_idx5d_tmp(cbsite,smin,0,0,Nspinco);
	int sslices=smax-smin+1;
	__builtin_prefetch(&xp[bidx],1);
	__builtin_prefetch(&rp[bidx],1);
	if ( vec == 0 ) {
	  bfm_site_zscale<Float>(sslices*Nspinco,&rp[bidx],&xp[bidx],(double *)&zz[0]);
	} else { 
	  bfm_site_zaxpy<Float>(sslices*Nspinco,&rp[bidx],&xp[bidx],&rp[bidx],(double *)&zz[0]);
	}
      }
    }
  }
  uint64_t t1 = GetTimeBase();
  /*
  if ( print < 20 && !me ) { 
    uint64_t micsec = (t1-t0)/MHz();
    double MBps =1.0*bytes/micsec;
    BossLog("Promote<%d> %d microsecs smin,smax=%d,%d %f MB/s\n",
	    sizeof(Float),(t1-t0)/MHz(),smin,smax,MBps);
    print++;
  }
  */

  thread_barrier();


}

template<class Float>
template<class cFloat>
void   bfmbase<Float>::block_zaxpy(Fermion_t result,Fermion_t x,Fermion_t y, double a, cFloat *zap,
				   int *blockids,int smin,int smax)
{
  int Nspinco=12;
  Float *xp = (Float *)x;
  Float *yp = (Float *)y;
  Float *rp = (Float *)result;


  int me,thrlen,throff;
  int unit=this->simd();
  thread_work(this->node_cbvol/unit,me,thrlen,throff);
  thrlen*=unit;
  throff*=unit;

  int block[4];
  double zz_b[4*4] = {0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0};
  double *zz = ALIGNIT(zz_b);

  if ( ((uint64_t)&zz[0])&0x1F ) { 
    Error("bfm::block_zaxpy -- unaligned pointer\n");
    exit(0);
  }
  if ( ((uint64_t)xp)&0x1F ) { 
    Error("bfm::block_zaxpy -- unaligned pointer\n");
    exit(0);
  }
  if ( ((uint64_t)yp)&0x1F ) { 
    Error("bfm::block_zaxpy -- unaligned pointer\n");
    exit(0);
  }
  if ( ((uint64_t)rp)&0x1F ) { 
    Error("bfm::block_zaxpy -- unaligned pointer\n");
    exit(0);
  }
  if ( thrlen & 0x1 ) { 
    Error("odd thread work unit \n");
    exit(0);
  }
  for (int cbsite=throff;cbsite<throff+thrlen;cbsite+=this->simd() ) { 

    int cidx = bagel_idx5d_tmp(cbsite,smin,0,0,1); cidx=cidx>>1;

    for(int v=0;v<this->simd();v++){
      block[v]= blockids[cidx+v];
      zz[2*v+0] = a*zap[2*block[v]+0];
      zz[2*v+1] = a*zap[2*block[v]+1];
    }

    if ( smin > smax ) {
      int bidx = bagel_idx5d_tmp(cbsite,smin,0,0,Nspinco);
      int sslices=this->cbLs-smin;
      bfm_site_zaxpy<Float>(sslices*Nspinco,&rp[bidx],&xp[bidx],&yp[bidx],(double *)&zz[0]);

      int s=0;
      bidx = bagel_idx5d_tmp(cbsite,s,0,0,Nspinco);
      sslices=smax+1;
      bfm_site_zaxpy<Float>(sslices*Nspinco,&rp[bidx],&xp[bidx],&yp[bidx],(double *)&zz[0]);
    } else {
      int bidx = bagel_idx5d_tmp(cbsite,smin,0,0,Nspinco);
      int sslices=smax-smin+1;
      bfm_site_zaxpy<Float>(sslices*Nspinco,&rp[bidx],&xp[bidx],&yp[bidx],(double *)&zz[0]);
    }
  }
  thread_barrier();
}
/*
template<class Float>
void bfmbase<Float>::block_norm(Fermion_t in,int nblock,double *result,int *blockids)
{
  int Nspinco=12;
 
  for(int i=0;i<nblock;i++){
    result[i]=0;
  }

  Float *inp = (Float *)in;

  int me,thrlen,throff;
  thread_work(this->node_cbvol,me,thrlen,throff);
  for (int cbsite=throff;cbsite<throff+thrlen;cbsite++ ) { 
    for (int s=0;s<this->cbLs;s++ ) { 

	int cidx = bagel_idx5d_tmp(cbsite,s,0,0,1); cidx=cidx>>1;
	int block  = blockids[cidx];

	Float fr=0;
	for ( int co=0;co<Nspinco;co++ ) { 
	  int bidx = bagel_idx5d_tmp(cbsite,s,0,co,Nspinco);
	  fr = fr + inp[bidx]*inp[bidx]   + inp[bidx+1]*inp[bidx+1];   //rr+ii
	}

#pragma omp critical
	{
	  result[block]  +=fr;
	}
    } 
  }
  thread_barrier();
}
*/
template<class Float>
void bfmbase<Float>::block_norm(Fermion_t in,int nblock,double *result,int *blockids)
{
  int Nspinco=12;

  int me,thrlen,throff;
  thread_work(nblock,me,thrlen,throff);
 
  for(int i=throff;i<thrlen+throff;i++){
    result[i]=0;
  }

  Float *inp = (Float *)in;
  for (int cbsite=0;cbsite<node_cbvol;cbsite++ ) { 
    int cidx = bagel_idx5d_tmp(cbsite,0,0,0,1); cidx=cidx>>1;
    int block  = blockids[cidx];
    if ( (block >= throff)  && (block < thrlen+throff) ) {
      Float fr=0;
      for (int s=0;s<this->cbLs;s++ ) { 
	for ( int co=0;co<Nspinco;co++ ) { 
	  int bidx = bagel_idx5d_tmp(cbsite,s,0,co,Nspinco);
	  fr = fr + inp[bidx]*inp[bidx]   + inp[bidx+1]*inp[bidx+1];   //rr+ii
	}
      }
      result[block]  +=fr;
    }
  }
  thread_barrier();
}


template<class Float>
void bfmbase<Float>::block_normalise(Fermion_t in,int nblock,double *scale,int *blockids)
{
  int Nspinco=12;
  Float *inp = (Float *)in;

  int me,thrlen,throff;
  thread_work(this->node_cbvol,me,thrlen,throff);
  for (int cbsite=throff;cbsite<throff+thrlen;cbsite++ ) { 
    for (int s=0;s<this->cbLs;s++ ) { 

	int cidx = bagel_idx5d_tmp(cbsite,s,0,0,1); cidx=cidx>>1;
	int block  = blockids[cidx];
	double cc   = 1.0/sqrt(scale[block]);

	for ( int co=0;co<Nspinco;co++ ) { 
	  int bidx = bagel_idx5d_tmp(cbsite,s,0,co,Nspinco);
	  inp[bidx]= cc*inp[bidx];
	  inp[bidx+1]= cc*inp[bidx+1];
	}
    } 
  }
  thread_barrier();
}


template class bfmbase<double>;
template class bfmbase<float>;
template void bfmbase<double>::block_promote<double>(Fermion_t result,Fermion_t *basis, int nbasis,double *zap,int *blockids,int smin,int smax);
template void bfmbase<double>::block_promote<float>(Fermion_t result,Fermion_t *basis, int nbasis,float *zap,int *blockids,int smin,int smax);
template void bfmbase<float>::block_promote<double>(Fermion_t result,Fermion_t *basis, int nbasis,double *zap,int *blockids,int smin,int smax);
template void bfmbase<float>::block_promote<float>(Fermion_t result,Fermion_t *basis, int nbasis,float *zap,int *blockids,int smin,int smax);

template void bfmbase<float>::block_project<float>(Fermion_t *basis,int nbasis,Fermion_t in,
						   int nblock,float *result,
						   double *reduce,int *blockids,int smin,int smax);
template void bfmbase<float>::block_project<double>(Fermion_t *basis,int nbasis,Fermion_t in,
						    int nblock,double *result,
						   double *reduce,int *blockids,int smin,int smax);

template void bfmbase<double>::block_project<float>(Fermion_t *basis,int nbasis,Fermion_t in,
						   int nblock,float *result,
						   double *reduce,int *blockids,int smin,int smax);
template void bfmbase<double>::block_project<double>(Fermion_t *basis,int nbasis,Fermion_t in,
						    int nblock,double *result,
						   double *reduce,int *blockids,int smin,int smax);

template void bfmbase<float>::block_zaxpy<float>(Fermion_t result,Fermion_t x,Fermion_t y, double a, float  *zap, int *blockids,int smin,int smax);
template void bfmbase<float>::block_zaxpy<double>(Fermion_t result,Fermion_t x,Fermion_t y, double a, double *zap, int *blockids,int smin,int smax);
template void bfmbase<double>::block_zaxpy<float>(Fermion_t result,Fermion_t x,Fermion_t y, double a, float  *zap, int *blockids,int smin,int smax);
template void bfmbase<double>::block_zaxpy<double>(Fermion_t result,Fermion_t x,Fermion_t y, double a, double *zap, int *blockids,int smin,int smax);

template void bfmbase<float>::block_zaxpy<float>(Fermion_t result,Fermion_t x,Fermion_t y, double a, float  *zap, int *blockids);
template void bfmbase<float>::block_zaxpy<double>(Fermion_t result,Fermion_t x,Fermion_t y, double a, double *zap, int *blockids);
template void bfmbase<double>::block_zaxpy<float>(Fermion_t result,Fermion_t x,Fermion_t y, double a, float  *zap, int *blockids);
template void bfmbase<double>::block_zaxpy<double>(Fermion_t result,Fermion_t x,Fermion_t y, double a, double *zap, int *blockids);
