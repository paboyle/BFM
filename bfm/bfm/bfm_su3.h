#ifndef BFM_SU3_H
#define BFM_SU3_H

#include "bfm.h"
#include <stdio.h>
#include <stdlib.h>

////////////////////////////////////////////////////////////
// Complex vector arithmetic to accelerate C++
// Would love to find a way to generate intrinsics for this
// that is portable
////////////////////////////////////////////////////////////
template<class Float, int NSIMD> class bfmComplex
{
  typedef std::complex<Float> Fcomplex;
 public:
  Fcomplex v[NSIMD] ;
  Fcomplex operator [] (int n) { return v[n]; };
};

template<class Float, int NSIMD> Float norm2 (const bfmComplex<Float,NSIMD> &l){
  Float f = 0;
  bfmComplex<Float,NSIMD> ll = l * conj(l) ;
  for(int i=0;i<NSIMD;i++){
    f += ll[i].real() ;
  }
  return f;
};
template<class Float, int NSIMD> bfmComplex<Float,NSIMD> conj(const bfmComplex<Float,NSIMD> &l)
{
  bfmComplex<Float,NSIMD>  rr;
  for(int i=0;i<NSIMD;i++)
    rr.v[i] = conj(l.v[i]);
  return rr;
};

template<class Float, int NSIMD> bfmComplex<Float,NSIMD> operator * (const bfmComplex<Float,NSIMD> &l,const bfmComplex<Float,NSIMD> &r){
  bfmComplex<Float,NSIMD>  rr;
  for(int i=0;i<NSIMD;i++)
    rr.v[i] = l.v[i]*r.v[i];
  return rr;
};

typedef bfmComplex<float,BAGEL_NSIMD> bfmFcomplex;
typedef bfmComplex<double,BAGEL_NSIMD> bfmDcomplex;

///////////////////////////////////////////////////////////////////
// SU3 covariant deriv support
///////////////////////////////////////////////////////////////////
template <class Float>
class bfmsu3 : public bfmbase<Float>{

  typedef bfmComplex<Float,BAGEL_NSIMD> ComplexVec;

  void impexMatrix  (QDPdouble *mat, Matrix_t handle, int import) ;
  void      importMatrix  (QDPdouble *mat, Matrix_t handle);
  void      exportMatrix  (QDPdouble *mat, Matrix_t handle);

  Matrix_t allocMatrix   (void); 
  void      freeMatrix    (Matrix_t handle);

  // Always treat matrix as both parities as no need for checkerboarding and this is simpler
  void CovariantShift(Matrix_t out, Matrix_t Umu, Matrix_t in,int mu,int dir);
  void          Shift(Matrix_t out, Matrix_t in,int mu,int dir);

  void CovariantShift(Fermion_t out[2], Matrix_t Umu, Fermion_t in[2],int mu,int dir);
  void          Shift(Fermion_t out[2], Fermion_t in[2],int mu,int dir);

  double Plaquette(int mu,int nu, Matrix_t U[4]);
  double Rectangle(int mu,int nu, Matrix_t U[4]);
  void   Staple(Matrix_t stap, Matrix_t U[4],int mu);

  void      axpby(Matrix_t result, Matrix_t x, Matrix_t y,double a,double b);
  void      axpy(Matrix_t result, Matrix_t x, Matrix_t y,double a,double b);
  void      conj(Matrix_t src);
  void      mult(Matrix_t result, Matrix left, Matrix_t right,int num);
  void      mult(Fermion_t result, Matrix left, Fermion_t right);
}

template<class Float> void bfmsu3<Float>::impexMatrix  (QDPdouble *mat, Matrix_t handle, int import)
{
  int Ncoco=9;

  omp_set_num_threads(this->nthread);

#pragma omp parallel
  {
#pragma omp for 
  for (int site=0;site<node_latt[0]*node_latt[1]*node_latt[2]*node_latt[3];site++ ) { 
    int x[4] ;
    int si=site;
    x[0]=si%node_latt[0];    si=si/node_latt[0];
    x[1]=si%node_latt[1];    si=si/node_latt[1];
    x[2]=si%node_latt[2];    si=si/node_latt[2];
    x[3]=si%node_latt[3];

    for(int cb=0;cb<2;cb++){
    if ( ((x[0]+x[1]+x[2]+x[3])&0x1) == cb ) {      
      for ( int co=0;co<Ncoco;co++ ) { 
      for ( int reim=0;reim<2;reim++ ) { 
        Float * bagel = (Float *)handle;
        int cidx = chroma_idx(x,reim,co,Ncoco);
        int bidx = bagel_idx(x,reim,co,Ncoco,0); // 0 here means both parities stored in single vector

	if ( doimport ) bagel[bidx] = mat[cidx]; 
	else mat[cidx] = bagel[bidx] ;
      }}
    }      
    }
  }
  }
}
template<class Float> void bfmsu3<Float>::importMatrix  (QDPdouble *mat, Matrix_t handle)
{
  this->impexMatrix(mat,handle,1);
}
template<class Float> void bfmsu3<Float>::exportMatrix  (QDPdouble *mat, Matrix_t handle)
{
  this->impexMatrix(mat,handle,0);
}
template<class Float> Matrix_t bfmsu3<Float>::allocMatrix   (void)
{
  int words =    this->node_latt[0]*this->node_latt[1]
                *this->node_latt[2]*this->node_latt[3];
  words = words * 18*sizeof(Float);
  Float *mat = (Float *)bfm_alloc(words*cbLs*sizeof(Float));
  if(mat == 0){
    if ( isBoss() ) printf("bfmbase::allocMatrix\n");
    fflush(stdout);
    exit(-1);
  }
  return (Matrix_t)mat;
}
template<class Float> void     bfmsu3<Float>::freeMatrix    (Matrix_t handle)
{
  Float * mat = ( Float * ) handle;
  bfm_free(mat);
}

template<class Float>   void bfmsu3<Float>::CovariantShift(Matrix_t out, Matrix_t Umu,Matrix_t in,int mu,int dir)
{
  ComplexVec *out_p = (ComplexVec *)     out;
  ComplexVec * in_p = (ComplexVec *)      in;
  ComplexVec *  u_p = (ComplexVec *)     Umu;

  // Need to gather face, communicate, fill halo
  for(int cb=0;cb<this->ncb;cb++){
    for(int cbsite=0;cbsite<this->simd_cbvol;cbsite++){

      int tab_idx  = this->interleave_site(pm,mu,cbsite);
      int info_idx = this->interleave_site(0,5,cbsite);

      int nbr = this->shift_table[cb][tab_idx];
      int msk = this->shift_table[cb][info_idx];

      ComplexVec * nbr_p = in_p + nbr*9;
      // nbr should be the _other_ checkerboard

      for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
	out_p[i*3+j] = u_p[i*3+0]*nbr_p[0*3+j]
                     + u_p[i*3+1]*nbr_p[1*3+j]
  	             + u_p[i*3+2]*nbr_p[2*3+j];
      }}
      u_p = u_p + 9;
    }
  }
}

template<class Float>   double Plaquette(int mu,int nu, Matrix_t U[4])
{
}
template<class Float>   double Rectangle(int mu,int nu, Matrix_t U[4])
{
}
template<class Float>   void   Staple(Matrix_t stap, Matrix_t U[4],int mu)
{
}
template<class Float>   void bfmsu3<Float>::Shift         (Matrix_t out, Matrix_t Umu,Matrix_t in,int mu,int dir)
{
  // Easier to implement once finish covariant shift & simplify
}
template<class Float>   void bfmsu3<Float>::CovariantShift(Fermion_t out[2], Matrix_t Umu,Fermion_t in[2],int mu,int dir)
{
  // Easier to implement once finish covariant shift & simplify
}
template<class Float>   void bfmsu3<Float>::Shift         (Fermion_t out[2], Fermion_t in[2],int mu,int dir)
{
  // Easier to implement once finish covariant shift & simplify
}
template<class Float>   void bfmsu3<Float>::axpby(Matrix_t result, Matrix_t x, Matrix_t y,double a,double b)
{
  // Use existing axpy etc
}
template<class Float>   void bfmsu3<Float>::axpy (Matrix_t result, Matrix_t x, Matrix_t y,double a,double b)
{
  // Use existing axpy etc
}
template<class Float>   void bfmsu3<Float>::conj (Matrix_t src)
{
  // Easy to implement
}
template<class Float>   void bfmsu3<Float>::mult (Matrix_t result, Matrix left, Matrix_t right)
{
  // Easy to implement
}
