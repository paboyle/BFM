#include <essl.h>

void LapackHermitianInvert(int Dim,double *mat)
{
  int info;
  int pivots[Dim];
  _ESVCOM work[1];

  zgetrf(Dim,Dim,mat,Dim,pivots,&info);

  zgetri(Dim,mat,Dim,pivots,work,0,&info);
  
  return;
}
