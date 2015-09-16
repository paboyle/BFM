#undef BFM_USE_ESSL
#define BFM_USE_GSL

#ifdef BFM_USE_ESSL
#include <essl.h>
#endif
#ifdef BFM_USE_GSL
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif
#ifdef BFM_USE_ESSL
void LapackHermitianDiagonalise(int Dim,double *mat,double *evecs,double *evals)
{
  int info;					
  double ddum=1.0;
  int    idum=1;
  double abstol=1.0e-7;
  _ESVCOM work[Dim*2];
  double  rwork[Dim*7];
  int     iwork[Dim*5];
  int     ifail[Dim];
  /*
  void   esvzheevx(const char *, const char *, const char *, _ESVINT, void *,
		   _ESVINT, double, double, _ESVINT, _ESVINT, double, _ESVINT *,
		   double *, void *, _ESVINT, _ESVCOM *, _ESVINT, double *,
		   _ESVINT *, _ESVINT *, _ESVI);
  */
  zheevx("V","A","U",Dim,mat,
	 Dim,ddum,ddum,idum,idum,abstol,&idum,
	 evals,evecs,Dim,work,Dim*2,rwork,
	 iwork,ifail,&info);
}

void LapackHermitianInvert(int Dim,double *mat,int isboss)
{
  int info;
  int pivots[Dim];
  _ESVCOM work[1];

  zgetrf(Dim,Dim,mat,Dim,pivots,&info);
  if ( info ) 
    printf("Blas factorisation fail %d\n",info);
  zgetri(Dim,mat,Dim,pivots,work,0,&info);
  if ( info ) 
    printf("Blas zgetri fail %d\n",info);
  
  return;
}


void LapackHermitianDiagonalise_f(int Dim,float *mat,float *evecs,float *evals)
{
  int info;					
  float ddum=1.0;
  int    idum=1;
  float abstol=1.0e-7;
  _ESVCOM work[Dim*2];
  float   rwork[Dim*7];
  int     iwork[Dim*5];
  int     ifail[Dim];
  /*
  void   esvzheevx(const char *, const char *, const char *, _ESVINT, void *,
		   _ESVINT, double, double, _ESVINT, _ESVINT, double, _ESVINT *,
		   double *, void *, _ESVINT, _ESVCOM *, _ESVINT, double *,
		   _ESVINT *, _ESVINT *, _ESVI);
  */
  cheevx("V","A","U",Dim,mat,
	 Dim,ddum,ddum,idum,idum,abstol,&idum,
	 evals,evecs,Dim,work,Dim*2,rwork,
	 iwork,ifail,&info);
}

void LapackHermitianInvert_f(int Dim,float *mat,int isboss)
{
  int info;
  int pivots[Dim];
  _ESVCOM work[1];

  cgetrf(Dim,Dim,mat,Dim,pivots,&info);

  cgetri(Dim,mat,Dim,pivots,work,0,&info);
  
  return;
}
#endif


#ifdef BFM_USE_GSL
void LapackHermitianDiagonalise(int Dim,double *mat,double *evecs,double *evals)
{
  if ( Dim <= 0 ) {
    printf("Oops diagonalise_d dimension %d %x %x\n",Dim,__builtin_return_address(1),__builtin_return_address(2));
    exit(0);
  }
  gsl_matrix_complex *   m = gsl_matrix_complex_alloc(Dim,Dim);
  gsl_matrix_complex *evec = gsl_matrix_complex_alloc(Dim,Dim);
  gsl_vector *eval = gsl_vector_alloc(Dim);
  gsl_eigen_hermv_workspace *WORK= gsl_eigen_hermv_alloc (Dim);

  double re,im;
  int i,j;
  gsl_complex z;
  for(i=0;i<Dim;i++){
    for(j=0;j<Dim;j++){
     re = mat[2*i*Dim+2*j];
     im = mat[2*i*Dim+2*j+1];
     GSL_SET_COMPLEX(&z,re,im);
     gsl_matrix_complex_set(m,i,j,z);
   }
  }

  gsl_eigen_hermv(m,eval,evec,WORK);
  gsl_eigen_hermv_sort(eval,evec,GSL_EIGEN_SORT_ABS_ASC);

  for(i=0;i<Dim;i++){
    evals[i] =  gsl_vector_get(eval,i);
    for(j=0;j<Dim;j++){
      z = gsl_matrix_complex_get(evec,i,j);
      evecs[2*j*Dim+2*i+0] = GSL_REAL(z); 
      evecs[2*j*Dim+2*i+1] =-GSL_IMAG(z); 
    }
  }  

  gsl_matrix_complex_free(m);
  gsl_matrix_complex_free(evec);
  gsl_vector_free(eval);
  gsl_eigen_hermv_free (WORK);
  return;
}

void LapackHermitianInvert(int Dim,double *mat,int isboss)
{
  if ( Dim <= 0 ) {
    printf("Oops dimension %d %x %x\n",Dim,__builtin_return_address(1),__builtin_return_address(2));
    exit(0);
  }
  gsl_matrix_complex *   m = gsl_matrix_complex_alloc(Dim,Dim);
  gsl_matrix_complex * minv= gsl_matrix_complex_alloc(Dim,Dim);
  gsl_permutation *      p = gsl_permutation_alloc(Dim);
  int signum;

  gsl_complex z;
  double re,im;
  int i,j;
  for(i=0;i<Dim;i++){
    for(j=0;j<Dim;j++){
     re = mat[2*i*Dim+2*j];
     im = mat[2*i*Dim+2*j+1];
     GSL_SET_COMPLEX(&z,re,im);
     gsl_matrix_complex_set(m,i,j,z);
   }
  }

  gsl_linalg_complex_LU_decomp(m,p,&signum);
  gsl_linalg_complex_LU_invert(m,p,minv);

  for(i=0;i<Dim;i++){
    for(j=0;j<Dim;j++){
      z = gsl_matrix_complex_get(minv,i,j);
      mat[2*i*Dim+2*j+0] = GSL_REAL(z); 
      mat[2*i*Dim+2*j+1] = GSL_IMAG(z); 
    }
  }  

  gsl_matrix_complex_free(m);
  gsl_matrix_complex_free(minv);
  gsl_permutation_free(p);
  return;
}


void LapackHermitianDiagonalise_f(int Dim,float *mat,float *evecs,float *evals)
{
  if ( Dim <= 0 ) {
    printf("Oops diagonalis_f dimension %d %x %x\n",Dim,__builtin_return_address(1),__builtin_return_address(2));
    exit(0);
  }
  gsl_matrix_complex *   m = gsl_matrix_complex_alloc(Dim,Dim);
  gsl_matrix_complex *evec = gsl_matrix_complex_alloc(Dim,Dim);
  gsl_vector  *eval = gsl_vector_alloc(Dim);
  gsl_eigen_hermv_workspace *WORK= gsl_eigen_hermv_alloc (Dim);

  gsl_complex z;
  double re,im;
  int i,j;
  for(i=0;i<Dim;i++){
    for(j=0;j<Dim;j++){
     re = mat[2*i*Dim+2*j];
     im = mat[2*i*Dim+2*j+1];
     GSL_SET_COMPLEX(&z,re,im);
     gsl_matrix_complex_set(m,i,j,z);
    }
  }

  gsl_eigen_hermv(m,eval,evec,WORK);
  gsl_eigen_hermv_sort(eval,evec,GSL_EIGEN_SORT_ABS_ASC);

  for(i=0;i<Dim;i++){
    
    evals[i] =gsl_vector_get(eval,i);
    for(j=0;j<Dim;j++){
      z = gsl_matrix_complex_get(evec,i,j);
      evecs[2*i*Dim+2*j+0] = GSL_REAL(z); 
      evecs[2*i*Dim+2*j+1] = GSL_IMAG(z); 
    }
  }  

  gsl_matrix_complex_free(m);
  gsl_matrix_complex_free(evec);
  gsl_vector_free(eval);
  gsl_eigen_hermv_free (WORK);
  return;
}

void LapackHermitianInvert_f(int Dim,float *mat,int isboss)
{
  if ( Dim <= 0 ) {
    printf("Oops dimension_f %d %x %x\n",Dim,__builtin_return_address(1),__builtin_return_address(2));
    exit(0);
  }

  gsl_matrix_complex *   m = gsl_matrix_complex_alloc(Dim,Dim);
  gsl_matrix_complex * minv= gsl_matrix_complex_alloc(Dim,Dim);
  gsl_permutation *      p = gsl_permutation_alloc(Dim);
  int signum;

  gsl_complex z;
  double re,im;
  int i,j;
  for(i=0;i<Dim;i++){
   for(j=0;j<Dim;j++){
     re = mat[2*i*Dim+2*j];
     im = mat[2*i*Dim+2*j+1];
     GSL_SET_COMPLEX(&z,re,im);
     gsl_matrix_complex_set(m,i,j,z);
   }
  }

  gsl_linalg_complex_LU_decomp(m,p,&signum);
  gsl_linalg_complex_LU_invert(m,p,minv);

  for(i=0;i<Dim;i++){
    for(j=0;j<Dim;j++){
      z = gsl_matrix_complex_get(minv,i,j);
      mat[2*i*Dim+2*j+0] = GSL_REAL(z); 
      mat[2*i*Dim+2*j+1] = GSL_IMAG(z); 
    }
  }  

  gsl_matrix_complex_free(m);
  gsl_matrix_complex_free(minv);
  gsl_permutation_free(p);
  return;
}
#endif
#ifdef __cplusplus
}
#endif
