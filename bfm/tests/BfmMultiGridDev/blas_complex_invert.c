#include <essl.h>

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

  zgetri(Dim,mat,Dim,pivots,work,0,&info);
  
  return;
}
