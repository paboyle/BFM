void qpx_zgemv(int N,_Complex double * __restrict A, _Complex double * __restrict in, double c, _Complex double * __restrict out);
void qpx_zgemv(int N,_Complex double * __restrict A, _Complex double * __restrict in, double c, _Complex double * __restrict out)
{
  vector4double tmp;
  vector4double AA;
  vector4double II;
  vector4double cc = {c,c,c,c};

  for(int i=0;i<N;i++){
    tmp = vec_ld(i,out);
    tmp = vec_mul(tmp,cc);
    for(int j=0;j<N;j+=2){

      AA =vec_ld(N*i+j,A);
      II =vec_ld(j,in);
      tmp = vec_xmadd(AA,II,tmp);
      tmp = vec_xxnpmadd(II,AA,tmp);

    }
    vec_st(tmp,i,out);
  }
}
