
//////////////////////////////////////////////////////////////////////////////
// Routines used by MultGrid to operate on coarse grid vectors
//////////////////////////////////////////////////////////////////////////////

void   qpx_zgemv(int N,_Complex double *  A, _Complex double *  in, double c, _Complex double *  out);
void   qpx_cgemv(int N,_Complex float   *  A, _Complex double *  in, double c, _Complex double *  out);
void   qpx_axpy(int N, double *  R,  double *  X,  double * Y,double c);
void   qpx_zaxpy(int N, double *  R,  double *  X,  double * Y,double *c);
void   qpx_axpby(int N, double *  R,  double *  X,  double * Y,double a,double b);

double qpx_inner_real(int N, double *A,  double *B);
void   qpx_inner(double *dot,int N, double *A,  double *B);
void   qpx_inner_f(double *dot,int N,float *A, float *B);
void   qpx_inner4(double *dot,int N, double *A,  double *B);
void   qpx_inner4_f(double *dot,int N,float *A, float *B);
double qpx_inner_real_f(int N,float *A, float *B);


double qpx_norm     (int N, double *  A);


void   qpx_axpy_f(int N,float *  R, float *  X, float * Y,double c);
void   qpx_zaxpy_f(int N, float *  R,  float *  X,  float * Y,double *c);
void   qpx_axpby_f(int N, float *  R,  float *  X,  float* Y,double a,double b);
double qpx_norm_f      (int N,float *A);
void   qpx_cgemv_f(int N,_Complex float   *  A, _Complex float *  in, double c, _Complex float *  out);

void qpx_gather(int Nvec,
		int *sbuf_gather,
		int sbuf_gather_size,
		double *sendbuf,
		double *vec);
void qpx_gather_f(int Nvec,
		int *sbuf_gather,
		int sbuf_gather_size,
		float *sendbuf,
		float *vec);

//////////////////////////////////////////////////////////////////////////////
// Routines used by MultGrid to operate on Bagel order fine grid vectors
// These contain a stride to stay on single site
//////////////////////////////////////////////////////////////////////////////
double qpx_site_norm  (int N, double *  A);
void   qpx_site_zaxpy (int N, double *  R,  double *  X,  double * Y,double *c);
void   qpx_site_zscale(int N, double *  R,  double *  X, double *c);

void   qpx_site_zaxpy_f (int N, float *  R, float *X, float * Y,double *c);

void   qpx_single_to_double (int N, float *  in, double *out);
void   qpx_double_to_single (int N, double *  in, float *out);

void qpx_site_inner_dd (int N,double *inner, double *X,double *Y);
void qpx_site_inner_df (int N,double *inner, float  *X  , float * Y);
void qpx_site_inner_fd (int N,float *inner, double *X,double *Y);
void qpx_site_inner_ff (int N,float *inner, float  *X  , float * Y);

void   qpx_site_zscale_f(int N, float * R, float *  X, double *c);

#define COMPLEX_SIMD_WIDTH (2)

void   qpx_single_to_double (int N, float *  in, double *out)
{
  __vector4double XX0;
  __vector4double XX1;
  __vector4double XX2;
  __vector4double XX3;
  __vector4double XX4;
  __vector4double XX5;
  for(int i=0;i<N;i+=6) { 
    XX0=vec_lda( COMPLEX_SIMD_WIDTH*(i+0)*2*sizeof(float),in); 
    XX4=vec_lda( COMPLEX_SIMD_WIDTH*(i+4)*2*sizeof(float),in); 
    XX1=vec_lda( COMPLEX_SIMD_WIDTH*(i+1)*2*sizeof(float),in); 
    XX2=vec_lda( COMPLEX_SIMD_WIDTH*(i+2)*2*sizeof(float),in); 
    XX3=vec_lda( COMPLEX_SIMD_WIDTH*(i+3)*2*sizeof(float),in); 
    XX5=vec_lda( COMPLEX_SIMD_WIDTH*(i+5)*2*sizeof(float),in); 
    vec_sta(XX0,COMPLEX_SIMD_WIDTH*(i+0)*2*sizeof(double),out);
    vec_sta(XX1,COMPLEX_SIMD_WIDTH*(i+1)*2*sizeof(double),out);
    vec_sta(XX2,COMPLEX_SIMD_WIDTH*(i+2)*2*sizeof(double),out);
    vec_sta(XX3,COMPLEX_SIMD_WIDTH*(i+3)*2*sizeof(double),out);
    vec_sta(XX4,COMPLEX_SIMD_WIDTH*(i+4)*2*sizeof(double),out);
    vec_sta(XX5,COMPLEX_SIMD_WIDTH*(i+5)*2*sizeof(double),out);
  }
}
void   qpx_double_to_single (int N, double *  in, float *out)
{
  __vector4double XX0;
  __vector4double XX1;
  __vector4double XX2;
  __vector4double XX3;
  __vector4double XX4;
  __vector4double XX5;
  for(int i=0;i<N;i+=6) { 
    XX0=vec_lda( COMPLEX_SIMD_WIDTH*(i+0)*2*sizeof(double),in); 
    XX2=vec_lda( COMPLEX_SIMD_WIDTH*(i+2)*2*sizeof(double),in); 
    XX4=vec_lda( COMPLEX_SIMD_WIDTH*(i+4)*2*sizeof(double),in); 
    XX1=vec_lda( COMPLEX_SIMD_WIDTH*(i+1)*2*sizeof(double),in); 
    XX3=vec_lda( COMPLEX_SIMD_WIDTH*(i+3)*2*sizeof(double),in); 
    XX5=vec_lda( COMPLEX_SIMD_WIDTH*(i+5)*2*sizeof(double),in); 
    vec_sta(XX0,COMPLEX_SIMD_WIDTH*(i+0)*2*sizeof(float),out);
    vec_sta(XX1,COMPLEX_SIMD_WIDTH*(i+1)*2*sizeof(float),out);
    vec_sta(XX2,COMPLEX_SIMD_WIDTH*(i+2)*2*sizeof(float),out);
    vec_sta(XX3,COMPLEX_SIMD_WIDTH*(i+3)*2*sizeof(float),out);
    vec_sta(XX4,COMPLEX_SIMD_WIDTH*(i+4)*2*sizeof(float),out);
    vec_sta(XX5,COMPLEX_SIMD_WIDTH*(i+5)*2*sizeof(float),out);
  }
}


void qpx_axpy(int N,double *  R,  double *  X,  double * Y,double c)
{
  __vector4double XX;
  __vector4double YY;
  __vector4double RR;
  __vector4double CC ={c,c,c,c};
  for(int i=0;i<N;i++) { 
      XX= vec_ld2a( i *2*sizeof(double),X);
      YY= vec_ld2a( i *2*sizeof(double),Y);
      RR = vec_madd(CC,XX,YY);
      vec_st2a(RR,    i *2*sizeof(double),R);
  }
}


void qpx_axpby(int N,double *  R,  double *  X,  double * Y,double a,double b)
{
  __vector4double XX;
  __vector4double YY;
  __vector4double RR;
  __vector4double AA ={a,a,a,a};
  __vector4double BB ={b,b,b,b};
  for(int i=0;i<N;i++) { 
      XX= vec_ld2a( i *2*sizeof(double),X);
      YY= vec_ld2a( i *2*sizeof(double),Y);
      RR = vec_mul (AA,XX);
      RR = vec_madd(BB,YY,RR);
      vec_st2a(RR,    i *2*sizeof(double),R);
  }
}

void   qpx_zaxpy(int N, double *  R,  double *  X,  double * Y,double *c)
{
  __vector4double XX;
  __vector4double YY;
  __vector4double RR;

  __vector4double CC = {c[0],c[1],c[0],c[1]};

  for(int i=0;i<N;i++) { 


    XX= vec_ld2a( COMPLEX_SIMD_WIDTH*(i)*sizeof(double),X);
    YY= vec_ld2a( COMPLEX_SIMD_WIDTH*(i)*sizeof(double),Y);

    RR = vec_xmadd(CC,XX,YY);
    RR = vec_xxnpmadd(XX,CC,RR);

    vec_st2a(RR, COMPLEX_SIMD_WIDTH*(i)*sizeof(double),R);
  }
}

void qpx_inner4(double *dot,int N, double *X,  double *Y)
{
  __vector4double XX0;
  __vector4double XX1;
  __vector4double XX2;
  __vector4double XX3;
  __vector4double YY0;
  __vector4double YY1;
  __vector4double YY2;
  __vector4double YY3;
  __vector4double tmp0 = { 0,0,0,0};
  __vector4double tmp1 = { 0,0,0,0};
  __vector4double tmp2 = { 0,0,0,0};
  __vector4double tmp3 = { 0,0,0,0};

  for(int i=0;i<N;i+=4) { 
    XX0= vec_ld2a( (i+0)*2*sizeof(double),X);
    XX1= vec_ld2a( (i+1)*2*sizeof(double),X);
    XX2= vec_ld2a( (i+2)*2*sizeof(double),X);
    XX3= vec_ld2a( (i+3)*2*sizeof(double),X);

    YY0= vec_ld2a( (i+0)*2*sizeof(double),Y);
    YY1= vec_ld2a( (i+1)*2*sizeof(double),Y);
    YY2= vec_ld2a( (i+2)*2*sizeof(double),Y);
    YY3= vec_ld2a( (i+3)*2*sizeof(double),Y);

    tmp0 = vec_xmadd(XX0,YY0,tmp0);
    tmp1 = vec_xmadd(XX1,YY1,tmp1);
    tmp2 = vec_xmadd(XX2,YY2,tmp2);
    tmp3 = vec_xmadd(XX3,YY3,tmp3);

    tmp0 = vec_xxcpnmadd(YY0,XX0,tmp0);
    tmp1 = vec_xxcpnmadd(YY1,XX1,tmp1);
    tmp2 = vec_xxcpnmadd(YY2,XX2,tmp2);
    tmp3 = vec_xxcpnmadd(YY3,XX3,tmp3);
  }
  dot[0] = vec_extract(tmp0,0)+vec_extract(tmp2,0)+
           vec_extract(tmp1,0)+vec_extract(tmp3,0);
  dot[1] = vec_extract(tmp0,1)+vec_extract(tmp2,1)+
           vec_extract(tmp1,1)+vec_extract(tmp3,1);
}

double qpx_norm(int N, double *  X)
{
  __vector4double XX;
  __vector4double tmp = { 0,0,0,0};

  for(int i=0;i<N;i++) { 
      XX= vec_ld2a( i *2*sizeof(double),X);
      tmp = vec_madd(XX,XX,tmp);
  }
  double d=vec_extract(tmp,0)+vec_extract(tmp,1);
  return d;
}


void qpx_gather(int Nvec,
		int *sbuf_gather,
		int sbuf_gather_size,
		double *sendbuf,
		double *vec)
{
  __vector4double XX;
  for(int j=0;j<sbuf_gather_size;j++){
    int m=sbuf_gather[j]*Nvec;
    int o=j*Nvec;
    for(int i=0;i<Nvec;i++){
      // sendbuf[j*Nvec+i] = vec[i+sbuf_gather[j]*Nvec];
      XX= vec_ld2a((i+m)*2*sizeof(double),vec);
      vec_st2a(XX,(i+o)*2*sizeof(double),sendbuf);
    }
  }
}
void qpx_gather_f(int Nvec,
		int *sbuf_gather,
		int sbuf_gather_size,
		float *sendbuf,
		float *vec)
{
  __vector4double XX;
  for(int j=0;j<sbuf_gather_size;j++){
    int m=sbuf_gather[j]*Nvec;
    int o=j*Nvec;
    for(int i=0;i<Nvec;i++){
      // sendbuf[j*Nvec+i] = vec[i+sbuf_gather[j]*Nvec];
      XX= vec_ld2a((i+m)*2*sizeof(float),vec);
      vec_st2a(XX,(i+o)*2*sizeof(float),sendbuf);
    }
  }
}


void qpx_zgemv(int N,_Complex double *  A, _Complex double *  in, double c, _Complex double *  out)
{
  __vector4double AA1;
  __vector4double AA2;
  __vector4double AA3;
  __vector4double AA4;
  __vector4double II;

  if ( N&0x3 ) exit(0);

  for(int i=0;i<N;i+=4){

    __vector4double tmp1 = { 0,0,0,0};
    __vector4double tmp2 = { 0,0,0,0};
    __vector4double tmp3 = { 0,0,0,0};
    __vector4double tmp4 = { 0,0,0,0};


    for(int j=0;j<N;j+=2){

      AA1 =vec_lda((N*i    +j)*2*sizeof(double),A);
      AA3 =vec_lda((N*(i+2)+j)*2*sizeof(double),A);

      II  =vec_lda(         j *2*sizeof(double),in);

      AA2 =vec_lda((N*(i+1)+j)*2*sizeof(double),A);
      AA4 =vec_lda((N*(i+3)+j)*2*sizeof(double),A);

      tmp1 = vec_xmadd(AA1,II,tmp1);
      tmp2 = vec_xmadd(AA2,II,tmp2);
      tmp3 = vec_xmadd(AA3,II,tmp3);
      tmp4 = vec_xmadd(AA4,II,tmp4);

      tmp1 = vec_xxnpmadd(II,AA1,tmp1);
      tmp2 = vec_xxnpmadd(II,AA2,tmp2);
      tmp3 = vec_xxnpmadd(II,AA3,tmp3);
      tmp4 = vec_xxnpmadd(II,AA4,tmp4);

    }
    
    double *dp = (double *)out;
    if ( c==0.0 ){
      dp[i    *2  ]=vec_extract(tmp1,0)+vec_extract(tmp1,2);
      dp[i    *2+1]=vec_extract(tmp1,1)+vec_extract(tmp1,3);
      dp[(i+1)*2  ]=vec_extract(tmp2,0)+vec_extract(tmp2,2);
      dp[(i+1)*2+1]=vec_extract(tmp2,1)+vec_extract(tmp2,3);
      dp[(i+2)*2  ]=vec_extract(tmp3,0)+vec_extract(tmp3,2);
      dp[(i+2)*2+1]=vec_extract(tmp3,1)+vec_extract(tmp3,3);
      dp[(i+3)*2  ]=vec_extract(tmp4,0)+vec_extract(tmp4,2);
      dp[(i+3)*2+1]=vec_extract(tmp4,1)+vec_extract(tmp4,3);
    } else { 
      dp[i    *2  ]=c*dp[i*2       ] +vec_extract(tmp1,0)+vec_extract(tmp1,2);
      dp[i    *2+1]=c*dp[i*2+1     ] +vec_extract(tmp1,1)+vec_extract(tmp1,3);
      dp[(i+1)*2  ]=c*dp[(i+1)*2   ] +vec_extract(tmp2,0)+vec_extract(tmp2,2);
      dp[(i+1)*2+1]=c*dp[(i+1)*2 +1] +vec_extract(tmp2,1)+vec_extract(tmp2,3);
      dp[(i+2)*2  ]=c*dp[(i+2)*2   ] +vec_extract(tmp3,0)+vec_extract(tmp3,2);
      dp[(i+2)*2+1]=c*dp[(i+2)*2+1 ] +vec_extract(tmp3,1)+vec_extract(tmp3,3);
      dp[(i+3)*2  ]=c*dp[(i+3)*2   ] +vec_extract(tmp4,0)+vec_extract(tmp4,2);
      dp[(i+3)*2+1]=c*dp[(i+3)*2+1 ] +vec_extract(tmp4,1)+vec_extract(tmp4,3);
    }
  }
}


void qpx_cgemv(int N,_Complex float *  A, _Complex double *  in, double c, _Complex double *  out)
{
  __vector4double AA1;
  __vector4double AA2;
  __vector4double AA3;
  __vector4double AA4;
  __vector4double II;

  if ( N&0x3 ) exit(0);

  for(int i=0;i<N;i+=4){

    __vector4double tmp1 = { 0,0,0,0};
    __vector4double tmp2 = { 0,0,0,0};
    __vector4double tmp3 = { 0,0,0,0};
    __vector4double tmp4 = { 0,0,0,0};


    for(int j=0;j<N;j+=2){

      AA1 =vec_lda((N*i    +j)*2*sizeof(float),A);
      AA3 =vec_lda((N*(i+2)+j)*2*sizeof(float),A);

      II  =vec_lda(         j *2*sizeof(double),in);

      AA2 =vec_lda((N*(i+1)+j)*2*sizeof(float),A);
      AA4 =vec_lda((N*(i+3)+j)*2*sizeof(float),A);

      tmp1 = vec_xmadd(AA1,II,tmp1);
      tmp2 = vec_xmadd(AA2,II,tmp2);
      tmp3 = vec_xmadd(AA3,II,tmp3);
      tmp4 = vec_xmadd(AA4,II,tmp4);

      tmp1 = vec_xxnpmadd(II,AA1,tmp1);
      tmp2 = vec_xxnpmadd(II,AA2,tmp2);
      tmp3 = vec_xxnpmadd(II,AA3,tmp3);
      tmp4 = vec_xxnpmadd(II,AA4,tmp4);

    }
    
    double *dp = (double *)out;
    if ( c==0.0 ){
      dp[i    *2  ]=vec_extract(tmp1,0)+vec_extract(tmp1,2);
      dp[i    *2+1]=vec_extract(tmp1,1)+vec_extract(tmp1,3);
      dp[(i+1)*2  ]=vec_extract(tmp2,0)+vec_extract(tmp2,2);
      dp[(i+1)*2+1]=vec_extract(tmp2,1)+vec_extract(tmp2,3);
      dp[(i+2)*2  ]=vec_extract(tmp3,0)+vec_extract(tmp3,2);
      dp[(i+2)*2+1]=vec_extract(tmp3,1)+vec_extract(tmp3,3);
      dp[(i+3)*2  ]=vec_extract(tmp4,0)+vec_extract(tmp4,2);
      dp[(i+3)*2+1]=vec_extract(tmp4,1)+vec_extract(tmp4,3);
    } else { 
      dp[i    *2  ]=c*dp[i*2       ] +vec_extract(tmp1,0)+vec_extract(tmp1,2);
      dp[i    *2+1]=c*dp[i*2+1     ] +vec_extract(tmp1,1)+vec_extract(tmp1,3);
      dp[(i+1)*2  ]=c*dp[(i+1)*2   ] +vec_extract(tmp2,0)+vec_extract(tmp2,2);
      dp[(i+1)*2+1]=c*dp[(i+1)*2 +1] +vec_extract(tmp2,1)+vec_extract(tmp2,3);
      dp[(i+2)*2  ]=c*dp[(i+2)*2   ] +vec_extract(tmp3,0)+vec_extract(tmp3,2);
      dp[(i+2)*2+1]=c*dp[(i+2)*2+1 ] +vec_extract(tmp3,1)+vec_extract(tmp3,3);
      dp[(i+3)*2  ]=c*dp[(i+3)*2   ] +vec_extract(tmp4,0)+vec_extract(tmp4,2);
      dp[(i+3)*2+1]=c*dp[(i+3)*2+1 ] +vec_extract(tmp4,1)+vec_extract(tmp4,3);
    }
  }
}

void qpx_cgemv_f(int N,_Complex float *  A, _Complex float *  in, double c, _Complex float *  out)
{
  float *Ap = (float *)A;
  float *Ip = (float *)in;
  float *Op = (float *)out;

  __vector4double AA1;
  __vector4double AA2;
  __vector4double AA3;
  __vector4double AA4;
  __vector4double II;

  if ( N&0x3 ) exit(0);

  for(int i=0;i<N;i+=4){

    __vector4double tmp1 = { 0,0,0,0};
    __vector4double tmp2 = { 0,0,0,0};
    __vector4double tmp3 = { 0,0,0,0};
    __vector4double tmp4 = { 0,0,0,0};


    for(int j=0;j<N;j+=2){

      AA1 =vec_lda((N*i    +j)*2*sizeof(float),Ap);
      AA3 =vec_lda((N*(i+2)+j)*2*sizeof(float),Ap);

      II  =vec_lda(         j *2*sizeof(float),Ip);

      AA2 =vec_lda((N*(i+1)+j)*2*sizeof(float),A);
      AA4 =vec_lda((N*(i+3)+j)*2*sizeof(float),A);

      tmp1 = vec_xmadd(AA1,II,tmp1);
      tmp2 = vec_xmadd(AA2,II,tmp2);
      tmp3 = vec_xmadd(AA3,II,tmp3);
      tmp4 = vec_xmadd(AA4,II,tmp4);

      tmp1 = vec_xxnpmadd(II,AA1,tmp1);
      tmp2 = vec_xxnpmadd(II,AA2,tmp2);
      tmp3 = vec_xxnpmadd(II,AA3,tmp3);
      tmp4 = vec_xxnpmadd(II,AA4,tmp4);

    }
    
    float *dp = (float *)out;
    if ( c==0.0 ){
      dp[i    *2  ]=vec_extract(tmp1,0)+vec_extract(tmp1,2);
      dp[i    *2+1]=vec_extract(tmp1,1)+vec_extract(tmp1,3);
      dp[(i+1)*2  ]=vec_extract(tmp2,0)+vec_extract(tmp2,2);
      dp[(i+1)*2+1]=vec_extract(tmp2,1)+vec_extract(tmp2,3);
      dp[(i+2)*2  ]=vec_extract(tmp3,0)+vec_extract(tmp3,2);
      dp[(i+2)*2+1]=vec_extract(tmp3,1)+vec_extract(tmp3,3);
      dp[(i+3)*2  ]=vec_extract(tmp4,0)+vec_extract(tmp4,2);
      dp[(i+3)*2+1]=vec_extract(tmp4,1)+vec_extract(tmp4,3);
    } else { 
      dp[i    *2  ]=c*dp[i*2       ] +vec_extract(tmp1,0)+vec_extract(tmp1,2);
      dp[i    *2+1]=c*dp[i*2+1     ] +vec_extract(tmp1,1)+vec_extract(tmp1,3);
      dp[(i+1)*2  ]=c*dp[(i+1)*2   ] +vec_extract(tmp2,0)+vec_extract(tmp2,2);
      dp[(i+1)*2+1]=c*dp[(i+1)*2 +1] +vec_extract(tmp2,1)+vec_extract(tmp2,3);
      dp[(i+2)*2  ]=c*dp[(i+2)*2   ] +vec_extract(tmp3,0)+vec_extract(tmp3,2);
      dp[(i+2)*2+1]=c*dp[(i+2)*2+1 ] +vec_extract(tmp3,1)+vec_extract(tmp3,3);
      dp[(i+3)*2  ]=c*dp[(i+3)*2   ] +vec_extract(tmp4,0)+vec_extract(tmp4,2);
      dp[(i+3)*2+1]=c*dp[(i+3)*2+1 ] +vec_extract(tmp4,1)+vec_extract(tmp4,3);
    }
  }
}


/////////////////////////////////////////////////////////////////////
// Remaining routines operate on a single site for bfm_block support
/////////////////////////////////////////////////////////////////////


void qpx_site_zaxpy(int N,double *  R,  double *  X,  double * Y,double *c)
{
  __vector4double XX0;
  __vector4double XX1;
  __vector4double XX2;
  __vector4double XX3;
  __vector4double XX4;
  __vector4double XX5;

  __vector4double YY0;
  __vector4double YY1;
  __vector4double YY2;
  __vector4double YY3;
  __vector4double YY4;
  __vector4double YY5;

  __vector4double RR0;
  __vector4double RR1;
  __vector4double RR2;
  __vector4double RR3;
  __vector4double RR4;
  __vector4double RR5;

  __vector4double CC = vec_lda(0,c);

  for(int i=0;i<N;i+=6) { 

    XX0= vec_lda( COMPLEX_SIMD_WIDTH*(i+0)*2*sizeof(double),X);
    YY0= vec_lda( COMPLEX_SIMD_WIDTH*(i+0)*2*sizeof(double),Y);

    XX1= vec_lda( COMPLEX_SIMD_WIDTH*(i+1)*2*sizeof(double),X);
    YY1= vec_lda( COMPLEX_SIMD_WIDTH*(i+1)*2*sizeof(double),Y);

    XX2= vec_lda( COMPLEX_SIMD_WIDTH*(i+2)*2*sizeof(double),X);
    YY2= vec_lda( COMPLEX_SIMD_WIDTH*(i+2)*2*sizeof(double),Y);

    XX3= vec_lda( COMPLEX_SIMD_WIDTH*(i+3)*2*sizeof(double),X);
    YY3= vec_lda( COMPLEX_SIMD_WIDTH*(i+3)*2*sizeof(double),Y);

    XX4= vec_lda( COMPLEX_SIMD_WIDTH*(i+4)*2*sizeof(double),X);
    YY4= vec_lda( COMPLEX_SIMD_WIDTH*(i+4)*2*sizeof(double),Y);

    XX5= vec_lda( COMPLEX_SIMD_WIDTH*(i+5)*2*sizeof(double),X);
    YY5= vec_lda( COMPLEX_SIMD_WIDTH*(i+5)*2*sizeof(double),Y);

    RR0 = vec_xmadd(CC,XX0,YY0);
    RR1 = vec_xmadd(CC,XX1,YY1);
    RR2 = vec_xmadd(CC,XX2,YY2);
    RR3 = vec_xmadd(CC,XX3,YY3);
    RR4 = vec_xmadd(CC,XX4,YY4);
    RR5 = vec_xmadd(CC,XX5,YY5);

    RR0 = vec_xxnpmadd(XX0,CC,RR0);
    RR1 = vec_xxnpmadd(XX1,CC,RR1);
    RR2 = vec_xxnpmadd(XX2,CC,RR2);
    RR3 = vec_xxnpmadd(XX3,CC,RR3);
    RR4 = vec_xxnpmadd(XX4,CC,RR4);
    RR5 = vec_xxnpmadd(XX5,CC,RR5);

    vec_sta(RR0, COMPLEX_SIMD_WIDTH*(i+0)*2*sizeof(double),R);
    vec_sta(RR1, COMPLEX_SIMD_WIDTH*(i+1)*2*sizeof(double),R);
    vec_sta(RR2, COMPLEX_SIMD_WIDTH*(i+2)*2*sizeof(double),R);
    vec_sta(RR3, COMPLEX_SIMD_WIDTH*(i+3)*2*sizeof(double),R);
    vec_sta(RR4, COMPLEX_SIMD_WIDTH*(i+4)*2*sizeof(double),R);
    vec_sta(RR5, COMPLEX_SIMD_WIDTH*(i+5)*2*sizeof(double),R);
  }
}


void qpx_site_zscale(int N,double *  R,  double *  X, double *c)
{
  __vector4double XX0;
  __vector4double XX1;
  __vector4double XX2;
  __vector4double XX3;
  __vector4double XX4;
  __vector4double XX5;

  __vector4double YY0 = {0,0,0,0};

  __vector4double RR0;
  __vector4double RR1;
  __vector4double RR2;
  __vector4double RR3;
  __vector4double RR4;
  __vector4double RR5;

  __vector4double CC = vec_lda(0,c);

  for(int i=0;i<N;i+=6) { 

    XX0= vec_lda( COMPLEX_SIMD_WIDTH*(i+0)*2*sizeof(double),X);
    XX1= vec_lda( COMPLEX_SIMD_WIDTH*(i+1)*2*sizeof(double),X);
    XX2= vec_lda( COMPLEX_SIMD_WIDTH*(i+2)*2*sizeof(double),X);
    XX3= vec_lda( COMPLEX_SIMD_WIDTH*(i+3)*2*sizeof(double),X);
    XX4= vec_lda( COMPLEX_SIMD_WIDTH*(i+4)*2*sizeof(double),X);
    XX5= vec_lda( COMPLEX_SIMD_WIDTH*(i+5)*2*sizeof(double),X);

    RR0 = vec_xmadd(CC,XX0,YY0);
    RR1 = vec_xmadd(CC,XX1,YY0);
    RR2 = vec_xmadd(CC,XX2,YY0);
    RR3 = vec_xmadd(CC,XX3,YY0);
    RR4 = vec_xmadd(CC,XX4,YY0);
    RR5 = vec_xmadd(CC,XX5,YY0);

    RR0 = vec_xxnpmadd(XX0,CC,RR0);
    RR1 = vec_xxnpmadd(XX1,CC,RR1);
    RR2 = vec_xxnpmadd(XX2,CC,RR2);
    RR3 = vec_xxnpmadd(XX3,CC,RR3);
    RR4 = vec_xxnpmadd(XX4,CC,RR4);
    RR5 = vec_xxnpmadd(XX5,CC,RR5);

    vec_sta(RR0, COMPLEX_SIMD_WIDTH*(i+0)*2*sizeof(double),R);
    vec_sta(RR1, COMPLEX_SIMD_WIDTH*(i+1)*2*sizeof(double),R);
    vec_sta(RR2, COMPLEX_SIMD_WIDTH*(i+2)*2*sizeof(double),R);
    vec_sta(RR3, COMPLEX_SIMD_WIDTH*(i+3)*2*sizeof(double),R);
    vec_sta(RR4, COMPLEX_SIMD_WIDTH*(i+4)*2*sizeof(double),R);
    vec_sta(RR5, COMPLEX_SIMD_WIDTH*(i+5)*2*sizeof(double),R);
  }
}


double qpx_site_norm(int N, double *  X)
{
  __vector4double XX;
  __vector4double tmp = { 0,0,0,0};

  for(int i=0;i<N;i++) { 
      XX= vec_ld2a( COMPLEX_SIMD_WIDTH*i *2*sizeof(double),X);
      tmp = vec_madd(XX,XX,tmp);
  }
  double d=vec_extract(tmp,0)+vec_extract(tmp,1);
  return d;
}


void qpx_site_zaxpy_f(int N,float *  R,  float *  X,  float * Y,double *c)
{
  __vector4double XX0;
  __vector4double XX1;
  __vector4double XX2;
  __vector4double XX3;
  __vector4double XX4;
  __vector4double XX5;

  __vector4double YY0;
  __vector4double YY1;
  __vector4double YY2;
  __vector4double YY3;
  __vector4double YY4;
  __vector4double YY5;

  __vector4double RR0;
  __vector4double RR1;
  __vector4double RR2;
  __vector4double RR3;
  __vector4double RR4;
  __vector4double RR5;

  __vector4double CC = vec_lda(0,c);

  for(int i=0;i<N;i+=6) { 

    XX0= vec_lda( COMPLEX_SIMD_WIDTH*(i+0)*2*sizeof(float),X);
    YY0= vec_lda( COMPLEX_SIMD_WIDTH*(i+0)*2*sizeof(float),Y);

    XX1= vec_lda( COMPLEX_SIMD_WIDTH*(i+1)*2*sizeof(float),X);
    YY1= vec_lda( COMPLEX_SIMD_WIDTH*(i+1)*2*sizeof(float),Y);

    XX2= vec_lda( COMPLEX_SIMD_WIDTH*(i+2)*2*sizeof(float),X);
    YY2= vec_lda( COMPLEX_SIMD_WIDTH*(i+2)*2*sizeof(float),Y);

    XX3= vec_lda( COMPLEX_SIMD_WIDTH*(i+3)*2*sizeof(float),X);
    YY3= vec_lda( COMPLEX_SIMD_WIDTH*(i+3)*2*sizeof(float),Y);

    XX4= vec_lda( COMPLEX_SIMD_WIDTH*(i+4)*2*sizeof(float),X);
    YY4= vec_lda( COMPLEX_SIMD_WIDTH*(i+4)*2*sizeof(float),Y);

    XX5= vec_lda( COMPLEX_SIMD_WIDTH*(i+5)*2*sizeof(float),X);
    YY5= vec_lda( COMPLEX_SIMD_WIDTH*(i+5)*2*sizeof(float),Y);

    RR0 = vec_xmadd(CC,XX0,YY0);
    RR1 = vec_xmadd(CC,XX1,YY1);
    RR2 = vec_xmadd(CC,XX2,YY2);
    RR3 = vec_xmadd(CC,XX3,YY3);
    RR4 = vec_xmadd(CC,XX4,YY4);
    RR5 = vec_xmadd(CC,XX5,YY5);

    RR0 = vec_xxnpmadd(XX0,CC,RR0);
    RR1 = vec_xxnpmadd(XX1,CC,RR1);
    RR2 = vec_xxnpmadd(XX2,CC,RR2);
    RR3 = vec_xxnpmadd(XX3,CC,RR3);
    RR4 = vec_xxnpmadd(XX4,CC,RR4);
    RR5 = vec_xxnpmadd(XX5,CC,RR5);

    vec_sta(RR0, COMPLEX_SIMD_WIDTH*(i+0)*2*sizeof(float),R);
    vec_sta(RR1, COMPLEX_SIMD_WIDTH*(i+1)*2*sizeof(float),R);
    vec_sta(RR2, COMPLEX_SIMD_WIDTH*(i+2)*2*sizeof(float),R);
    vec_sta(RR3, COMPLEX_SIMD_WIDTH*(i+3)*2*sizeof(float),R);
    vec_sta(RR4, COMPLEX_SIMD_WIDTH*(i+4)*2*sizeof(float),R);
    vec_sta(RR5, COMPLEX_SIMD_WIDTH*(i+5)*2*sizeof(float),R);
  }
}

void qpx_site_zscale_f(int N,float *  R,  float *  X, double *c)
{
  __vector4double XX0;
  __vector4double XX1;
  __vector4double XX2;
  __vector4double XX3;
  __vector4double XX4;
  __vector4double XX5;

  __vector4double YY0 = {0,0,0,0};

  __vector4double RR0;
  __vector4double RR1;
  __vector4double RR2;
  __vector4double RR3;
  __vector4double RR4;
  __vector4double RR5;

  __vector4double CC = vec_lda(0,c);

  for(int i=0;i<N;i+=6) { 

    XX0= vec_lda( COMPLEX_SIMD_WIDTH*(i+0)*2*sizeof(float),X);
    XX1= vec_lda( COMPLEX_SIMD_WIDTH*(i+1)*2*sizeof(float),X);
    XX2= vec_lda( COMPLEX_SIMD_WIDTH*(i+2)*2*sizeof(float),X);
    XX3= vec_lda( COMPLEX_SIMD_WIDTH*(i+3)*2*sizeof(float),X);
    XX4= vec_lda( COMPLEX_SIMD_WIDTH*(i+4)*2*sizeof(float),X);
    XX5= vec_lda( COMPLEX_SIMD_WIDTH*(i+5)*2*sizeof(float),X);

    RR0 = vec_xmadd(CC,XX0,YY0);
    RR1 = vec_xmadd(CC,XX1,YY0);
    RR2 = vec_xmadd(CC,XX2,YY0);
    RR3 = vec_xmadd(CC,XX3,YY0);
    RR4 = vec_xmadd(CC,XX4,YY0);
    RR5 = vec_xmadd(CC,XX5,YY0);

    RR0 = vec_xxnpmadd(XX0,CC,RR0);
    RR1 = vec_xxnpmadd(XX1,CC,RR1);
    RR2 = vec_xxnpmadd(XX2,CC,RR2);
    RR3 = vec_xxnpmadd(XX3,CC,RR3);
    RR4 = vec_xxnpmadd(XX4,CC,RR4);
    RR5 = vec_xxnpmadd(XX5,CC,RR5);

    vec_sta(RR0, COMPLEX_SIMD_WIDTH*(i+0)*2*sizeof(float),R);
    vec_sta(RR1, COMPLEX_SIMD_WIDTH*(i+1)*2*sizeof(float),R);
    vec_sta(RR2, COMPLEX_SIMD_WIDTH*(i+2)*2*sizeof(float),R);
    vec_sta(RR3, COMPLEX_SIMD_WIDTH*(i+3)*2*sizeof(float),R);
    vec_sta(RR4, COMPLEX_SIMD_WIDTH*(i+4)*2*sizeof(float),R);
    vec_sta(RR5, COMPLEX_SIMD_WIDTH*(i+5)*2*sizeof(float),R);
  }
}


void   qpx_zaxpy_f(int N, float *  R,  float *  X,  float * Y,double *c)
{
  __vector4double XX;
  __vector4double YY;
  __vector4double RR;

  __vector4double CC = {c[0],c[1],c[0],c[1]};

  for(int i=0;i<N;i++) { 

    // 2x possible?
    // Paired Complex??
    // Prefetch & unroll ??

    XX= vec_ld2a( COMPLEX_SIMD_WIDTH*(i)*sizeof(float),X);
    YY= vec_ld2a( COMPLEX_SIMD_WIDTH*(i)*sizeof(float),Y);

    RR = vec_xmadd(CC,XX,YY);
    RR = vec_xxnpmadd(XX,CC,RR);

    vec_st2a(RR, COMPLEX_SIMD_WIDTH*(i)*sizeof(float),R);
  }
}
void qpx_inner(double *dot,int N, double *X,  double *Y)
{
  __vector4double XX;
  __vector4double YY;
  __vector4double tmp = { 0,0,0,0};

  for(int i=0;i<N;i++) { 
      XX= vec_ld2a( i *2*sizeof(double),X);
      YY= vec_ld2a( i *2*sizeof(double),Y);
      tmp = vec_xmadd(XX,YY,tmp);
      tmp = vec_xxcpnmadd(YY,XX,tmp);
  }
  dot[0] = vec_extract(tmp,0);
  dot[1] = vec_extract(tmp,1);
}

double qpx_inner_real(int N, double *X,  double *Y)
{
  __vector4double XX;
  __vector4double YY;
  __vector4double tmp = { 0,0,0,0};

  for(int i=0;i<N;i++) { 
      XX= vec_ld2a( i *2*sizeof(double),X);
      YY= vec_ld2a( i *2*sizeof(double),Y);
      tmp = vec_xmadd(XX,YY,tmp);
      tmp = vec_xxcpnmadd(YY,XX,tmp);
  }
  return vec_extract(tmp,0);
}
void   qpx_inner_f(double *dot,int N,float *X, float *Y)
{
  __vector4double XX;
  __vector4double YY;
  __vector4double tmp = { 0,0,0,0};

  for(int i=0;i<N;i++) { 
      XX= vec_ld2a( i *2*sizeof(float),X);
      YY= vec_ld2a( i *2*sizeof(float),Y);
      tmp = vec_xmadd(XX,YY,tmp);
      tmp = vec_xxcpnmadd(YY,XX,tmp);
  }
  dot[0] = vec_extract(tmp,0);
  dot[1] = vec_extract(tmp,1);
}
double qpx_inner_real_f(int N,float *X, float *Y)
{
  __vector4double XX;
  __vector4double YY;
  __vector4double tmp = { 0,0,0,0};

  for(int i=0;i<N;i++) { 
      XX= vec_ld2a( i *2*sizeof(float),X);
      YY= vec_ld2a( i *2*sizeof(float),Y);
      tmp = vec_xmadd(XX,YY,tmp);
      tmp = vec_xxcpnmadd(YY,XX,tmp);
  }
  return vec_extract(tmp,0);
}

void   qpx_inner4_f(double *dot,int N,float *X, float *Y)
{
  __vector4double XX0;
  __vector4double XX1;
  __vector4double XX2;
  __vector4double XX3;
  __vector4double YY0;
  __vector4double YY1;
  __vector4double YY2;
  __vector4double YY3;
  __vector4double tmp0 = { 0,0,0,0};
  __vector4double tmp1 = { 0,0,0,0};
  __vector4double tmp2 = { 0,0,0,0};
  __vector4double tmp3 = { 0,0,0,0};

  for(int i=0;i<N;i+=4) { 
    XX0= vec_ld2a( (i+0)*2*sizeof(float),X);
    XX1= vec_ld2a( (i+1)*2*sizeof(float),X);
    XX2= vec_ld2a( (i+2)*2*sizeof(float),X);
    XX3= vec_ld2a( (i+3)*2*sizeof(float),X);

    YY0= vec_ld2a( (i+0)*2*sizeof(float),Y);
    YY1= vec_ld2a( (i+1)*2*sizeof(float),Y);
    YY2= vec_ld2a( (i+2)*2*sizeof(float),Y);
    YY3= vec_ld2a( (i+3)*2*sizeof(float),Y);

    tmp0 = vec_xmadd(XX0,YY0,tmp0);
    tmp1 = vec_xmadd(XX1,YY1,tmp1);
    tmp2 = vec_xmadd(XX2,YY2,tmp2);
    tmp3 = vec_xmadd(XX3,YY3,tmp3);

    tmp0 = vec_xxcpnmadd(YY0,XX0,tmp0);
    tmp1 = vec_xxcpnmadd(YY1,XX1,tmp1);
    tmp2 = vec_xxcpnmadd(YY2,XX2,tmp2);
    tmp3 = vec_xxcpnmadd(YY3,XX3,tmp3);
  }
  dot[0] = vec_extract(tmp0,0)+vec_extract(tmp1,0)+
           vec_extract(tmp2,0)+vec_extract(tmp3,0);
  dot[1] = vec_extract(tmp0,1)+vec_extract(tmp2,1)+
           vec_extract(tmp1,1)+vec_extract(tmp3,1);
}
double qpx_norm_f      (int N,float *  X)
{
  __vector4double XX;
  __vector4double tmp = { 0,0,0,0};

  for(int i=0;i<N;i++) { 
      XX= vec_ld2a( i *2*sizeof(float),X);
      tmp = vec_madd(XX,XX,tmp);
  }
  double d=vec_extract(tmp,0)+vec_extract(tmp,1);
  return d;
}
void qpx_axpy_f(int N,float *  R,  float *  X,  float * Y,double c)
{
  __vector4double XX;
  __vector4double YY;
  __vector4double RR;
  __vector4double CC ={c,c,c,c};
  for(int i=0;i<N;i++) { 
      XX= vec_ld2a( i *2*sizeof(float),X);
      YY= vec_ld2a( i *2*sizeof(float),Y);
      RR = vec_madd(CC,XX,YY);
      vec_st2a(RR,    i *2*sizeof(float),R);
  }
}
void qpx_axpby_f(int N,float *  R,  float *  X,  float * Y,double a,double b)
{
  __vector4double XX;
  __vector4double YY;
  __vector4double RR;
  __vector4double AA ={a,a,a,a};
  __vector4double BB ={b,b,b,b};
  for(int i=0;i<N;i++) { 
      XX= vec_ld2a( i *2*sizeof(float),X);
      YY= vec_ld2a( i *2*sizeof(float),Y);
      RR = vec_mul (AA,XX);
      RR = vec_madd(BB,YY,RR);
      vec_st2a(RR,    i *2*sizeof(float),R);
  }
}


void qpx_site_inner_dd (int N,double *inner, double *X,double *Y)
{
  __vector4double XX0;
  __vector4double XX1;
  __vector4double XX2;
  __vector4double XX3;
  __vector4double XX4;
  __vector4double XX5;

  __vector4double YY0;
  __vector4double YY1;
  __vector4double YY2;
  __vector4double YY3;
  __vector4double YY4;
  __vector4double YY5;


  __vector4double tmp0 = vec_lda(0,inner);
  __vector4double tmp1 = { 0,0,0,0};
  __vector4double tmp2 = { 0,0,0,0};
  __vector4double tmp3 = { 0,0,0,0};
  __vector4double tmp4 = { 0,0,0,0};
  __vector4double tmp5 = { 0,0,0,0};

  for(int i=0;i<N;i+=6) { 
    XX0= vec_lda( COMPLEX_SIMD_WIDTH*i*2*sizeof(double),X);   
    YY0= vec_lda( COMPLEX_SIMD_WIDTH*i*2*sizeof(double),Y); 

    XX2= vec_lda( COMPLEX_SIMD_WIDTH*(i+2)*2*sizeof(double),X);
    YY2= vec_lda( COMPLEX_SIMD_WIDTH*(i+2)*2*sizeof(double),Y); 

    XX4= vec_lda( COMPLEX_SIMD_WIDTH*(i+4)*2*sizeof(double),X);
    YY4= vec_lda( COMPLEX_SIMD_WIDTH*(i+4)*2*sizeof(double),Y); 

    XX1= vec_lda( COMPLEX_SIMD_WIDTH*(i+1)*2*sizeof(double),X);
    YY1= vec_lda( COMPLEX_SIMD_WIDTH*(i+1)*2*sizeof(double),Y); 

    XX3= vec_lda( COMPLEX_SIMD_WIDTH*(i+3)*2*sizeof(double),X);
    YY3= vec_lda( COMPLEX_SIMD_WIDTH*(i+3)*2*sizeof(double),Y); 

    XX5= vec_lda( COMPLEX_SIMD_WIDTH*(i+5)*2*sizeof(double),X);
    YY5= vec_lda( COMPLEX_SIMD_WIDTH*(i+5)*2*sizeof(double),Y); 

    tmp0 = vec_xmadd(XX0,YY0,tmp0);
    tmp1 = vec_xmadd(XX1,YY1,tmp1);
    tmp2 = vec_xmadd(XX2,YY2,tmp2);
    tmp3 = vec_xmadd(XX3,YY3,tmp3);
    tmp4 = vec_xmadd(XX4,YY4,tmp4);
    tmp5 = vec_xmadd(XX5,YY5,tmp5);

    __dcbt(&X[COMPLEX_SIMD_WIDTH*(i+6)*2]);
    __dcbt(&Y[COMPLEX_SIMD_WIDTH*(i+6)*2]);
    __dcbt(&X[COMPLEX_SIMD_WIDTH*(i+8)*2]);
    __dcbt(&Y[COMPLEX_SIMD_WIDTH*(i+8)*2]);
    __dcbt(&X[COMPLEX_SIMD_WIDTH*(i+10)*2]);
    __dcbt(&Y[COMPLEX_SIMD_WIDTH*(i+10)*2]);

    tmp0 = vec_xxcpnmadd(YY0,XX0,tmp0);
    tmp1 = vec_xxcpnmadd(YY1,XX1,tmp1);
    tmp2 = vec_xxcpnmadd(YY2,XX2,tmp2);
    tmp3 = vec_xxcpnmadd(YY3,XX3,tmp3);
    tmp4 = vec_xxcpnmadd(YY4,XX4,tmp4);
    tmp5 = vec_xxcpnmadd(YY5,XX5,tmp5);

  }
  tmp0 = vec_add(tmp0,tmp1);
  tmp0 = vec_add(tmp0,tmp2);
  tmp0 = vec_add(tmp0,tmp3);
  tmp0 = vec_add(tmp0,tmp4);
  tmp0 = vec_add(tmp0,tmp5);
  vec_sta(tmp0,0,inner);
}
void qpx_site_inner_df (int N,double *inner, float *X,float *Y)
{
  __vector4double XX0;
  __vector4double XX1;
  __vector4double XX2;
  __vector4double XX3;
  __vector4double XX4;
  __vector4double XX5;

  __vector4double YY0;
  __vector4double YY1;
  __vector4double YY2;
  __vector4double YY3;
  __vector4double YY4;
  __vector4double YY5;

  __vector4double tmp0 = vec_lda(0,inner);
  __vector4double tmp1 = { 0,0,0,0};
  __vector4double tmp2 = { 0,0,0,0};
  __vector4double tmp3 = { 0,0,0,0};
  __vector4double tmp4 = { 0,0,0,0};
  __vector4double tmp5 = { 0,0,0,0};

  for(int i=0;i<N;i+=6) { 
    XX0= vec_lda( COMPLEX_SIMD_WIDTH*i*2*sizeof(float),X);   
    YY0= vec_lda( COMPLEX_SIMD_WIDTH*i*2*sizeof(float),Y); 

    XX4= vec_lda( COMPLEX_SIMD_WIDTH*(i+4)*2*sizeof(float),X);
    YY4= vec_lda( COMPLEX_SIMD_WIDTH*(i+4)*2*sizeof(float),Y); 

    XX1= vec_lda( COMPLEX_SIMD_WIDTH*(i+1)*2*sizeof(float),X);
    YY1= vec_lda( COMPLEX_SIMD_WIDTH*(i+1)*2*sizeof(float),Y); 

    XX2= vec_lda( COMPLEX_SIMD_WIDTH*(i+2)*2*sizeof(float),X);
    YY2= vec_lda( COMPLEX_SIMD_WIDTH*(i+2)*2*sizeof(float),Y); 

    XX3= vec_lda( COMPLEX_SIMD_WIDTH*(i+3)*2*sizeof(float),X);
    YY3= vec_lda( COMPLEX_SIMD_WIDTH*(i+3)*2*sizeof(float),Y); 

    XX5= vec_lda( COMPLEX_SIMD_WIDTH*(i+5)*2*sizeof(float),X);
    YY5= vec_lda( COMPLEX_SIMD_WIDTH*(i+5)*2*sizeof(float),Y); 


    tmp0 = vec_xmadd(XX0,YY0,tmp0);
    tmp1 = vec_xmadd(XX1,YY1,tmp1);
    tmp2 = vec_xmadd(XX2,YY2,tmp2);
    tmp3 = vec_xmadd(XX3,YY3,tmp3);
    tmp4 = vec_xmadd(XX4,YY4,tmp4);
    tmp5 = vec_xmadd(XX5,YY5,tmp5);

    __dcbt(&X[COMPLEX_SIMD_WIDTH*(i+6)*2]);
    __dcbt(&Y[COMPLEX_SIMD_WIDTH*(i+6)*2]);
    __dcbt(&X[COMPLEX_SIMD_WIDTH*(i+10)*2]);
    __dcbt(&Y[COMPLEX_SIMD_WIDTH*(i+10)*2]);

    tmp0 = vec_xxcpnmadd(YY0,XX0,tmp0);
    tmp1 = vec_xxcpnmadd(YY1,XX1,tmp1);
    tmp2 = vec_xxcpnmadd(YY2,XX2,tmp2);
    tmp3 = vec_xxcpnmadd(YY3,XX3,tmp3);
    tmp4 = vec_xxcpnmadd(YY4,XX4,tmp4);
    tmp5 = vec_xxcpnmadd(YY5,XX5,tmp5);

  }
  tmp0 = vec_add(tmp0,tmp1);
  tmp0 = vec_add(tmp0,tmp2);
  tmp0 = vec_add(tmp0,tmp3);
  tmp0 = vec_add(tmp0,tmp4);
  tmp0 = vec_add(tmp0,tmp5);
  vec_sta(tmp0,0,inner);
}








void qpx_site_inner_fd (int N,float *inner, double *X,double *Y)
{
  __vector4double XX0;
  __vector4double XX1;
  __vector4double XX2;
  __vector4double XX3;
  __vector4double XX4;
  __vector4double XX5;

  __vector4double YY0;
  __vector4double YY1;
  __vector4double YY2;
  __vector4double YY3;
  __vector4double YY4;
  __vector4double YY5;

  __vector4double tmp0 = vec_lda(0,inner);
  __vector4double tmp1 = { 0,0,0,0};
  __vector4double tmp2 = { 0,0,0,0};
  __vector4double tmp3 = { 0,0,0,0};
  __vector4double tmp4 = { 0,0,0,0};
  __vector4double tmp5 = { 0,0,0,0};

  for(int i=0;i<N;i+=6) { 
    XX0= vec_lda( COMPLEX_SIMD_WIDTH*i*2*sizeof(double),X);   
    YY0= vec_lda( COMPLEX_SIMD_WIDTH*i*2*sizeof(double),Y); 

    XX2= vec_lda( COMPLEX_SIMD_WIDTH*(i+2)*2*sizeof(double),X);
    YY2= vec_lda( COMPLEX_SIMD_WIDTH*(i+2)*2*sizeof(double),Y); 

    XX4= vec_lda( COMPLEX_SIMD_WIDTH*(i+4)*2*sizeof(double),X);
    YY4= vec_lda( COMPLEX_SIMD_WIDTH*(i+4)*2*sizeof(double),Y); 

    XX1= vec_lda( COMPLEX_SIMD_WIDTH*(i+1)*2*sizeof(double),X);
    YY1= vec_lda( COMPLEX_SIMD_WIDTH*(i+1)*2*sizeof(double),Y); 

    XX3= vec_lda( COMPLEX_SIMD_WIDTH*(i+3)*2*sizeof(double),X);
    YY3= vec_lda( COMPLEX_SIMD_WIDTH*(i+3)*2*sizeof(double),Y); 

    XX5= vec_lda( COMPLEX_SIMD_WIDTH*(i+5)*2*sizeof(double),X);
    YY5= vec_lda( COMPLEX_SIMD_WIDTH*(i+5)*2*sizeof(double),Y); 

    tmp0 = vec_xmadd(XX0,YY0,tmp0);
    tmp1 = vec_xmadd(XX1,YY1,tmp1);
    tmp2 = vec_xmadd(XX2,YY2,tmp2);
    tmp3 = vec_xmadd(XX3,YY3,tmp3);
    tmp4 = vec_xmadd(XX4,YY4,tmp4);
    tmp5 = vec_xmadd(XX5,YY5,tmp5);

    __dcbt(&X[COMPLEX_SIMD_WIDTH*(i+6)*2]);
    __dcbt(&Y[COMPLEX_SIMD_WIDTH*(i+6)*2]);
    __dcbt(&X[COMPLEX_SIMD_WIDTH*(i+8)*2]);
    __dcbt(&Y[COMPLEX_SIMD_WIDTH*(i+8)*2]);
    __dcbt(&X[COMPLEX_SIMD_WIDTH*(i+10)*2]);
    __dcbt(&Y[COMPLEX_SIMD_WIDTH*(i+10)*2]);

    tmp0 = vec_xxcpnmadd(YY0,XX0,tmp0);
    tmp1 = vec_xxcpnmadd(YY1,XX1,tmp1);
    tmp2 = vec_xxcpnmadd(YY2,XX2,tmp2);
    tmp3 = vec_xxcpnmadd(YY3,XX3,tmp3);
    tmp4 = vec_xxcpnmadd(YY4,XX4,tmp4);
    tmp5 = vec_xxcpnmadd(YY5,XX5,tmp5);

  }
  tmp0 = vec_add(tmp0,tmp1);
  tmp0 = vec_add(tmp0,tmp2);
  tmp0 = vec_add(tmp0,tmp3);
  tmp0 = vec_add(tmp0,tmp4);
  tmp0 = vec_add(tmp0,tmp5);
  vec_sta(tmp0,0,inner);
}
void qpx_site_inner_ff (int N,float *inner, float *X,float *Y)
{
  __vector4double XX0;
  __vector4double XX1;
  __vector4double XX2;
  __vector4double XX3;
  __vector4double XX4;
  __vector4double XX5;

  __vector4double YY0;
  __vector4double YY1;
  __vector4double YY2;
  __vector4double YY3;
  __vector4double YY4;
  __vector4double YY5;

  __vector4double tmp0 = vec_lda(0,inner);
  __vector4double tmp1 = { 0,0,0,0};
  __vector4double tmp2 = { 0,0,0,0};
  __vector4double tmp3 = { 0,0,0,0};
  __vector4double tmp4 = { 0,0,0,0};
  __vector4double tmp5 = { 0,0,0,0};

  for(int i=0;i<N;i+=6) { 
    XX0= vec_lda( COMPLEX_SIMD_WIDTH*i*2*sizeof(float),X);   
    YY0= vec_lda( COMPLEX_SIMD_WIDTH*i*2*sizeof(float),Y); 

    XX4= vec_lda( COMPLEX_SIMD_WIDTH*(i+4)*2*sizeof(float),X);
    YY4= vec_lda( COMPLEX_SIMD_WIDTH*(i+4)*2*sizeof(float),Y); 

    XX1= vec_lda( COMPLEX_SIMD_WIDTH*(i+1)*2*sizeof(float),X);
    YY1= vec_lda( COMPLEX_SIMD_WIDTH*(i+1)*2*sizeof(float),Y); 

    XX2= vec_lda( COMPLEX_SIMD_WIDTH*(i+2)*2*sizeof(float),X);
    YY2= vec_lda( COMPLEX_SIMD_WIDTH*(i+2)*2*sizeof(float),Y); 

    XX3= vec_lda( COMPLEX_SIMD_WIDTH*(i+3)*2*sizeof(float),X);
    YY3= vec_lda( COMPLEX_SIMD_WIDTH*(i+3)*2*sizeof(float),Y); 

    XX5= vec_lda( COMPLEX_SIMD_WIDTH*(i+5)*2*sizeof(float),X);
    YY5= vec_lda( COMPLEX_SIMD_WIDTH*(i+5)*2*sizeof(float),Y); 

    tmp0 = vec_xmadd(XX0,YY0,tmp0);
    tmp1 = vec_xmadd(XX1,YY1,tmp1);
    tmp2 = vec_xmadd(XX2,YY2,tmp2);
    tmp3 = vec_xmadd(XX3,YY3,tmp3);
    tmp4 = vec_xmadd(XX4,YY4,tmp4);
    tmp5 = vec_xmadd(XX5,YY5,tmp5);

    __dcbt(&X[COMPLEX_SIMD_WIDTH*(i+6)*2]);
    __dcbt(&Y[COMPLEX_SIMD_WIDTH*(i+6)*2]);
    __dcbt(&X[COMPLEX_SIMD_WIDTH*(i+10)*2]);
    __dcbt(&Y[COMPLEX_SIMD_WIDTH*(i+10)*2]);

    tmp0 = vec_xxcpnmadd(YY0,XX0,tmp0);
    tmp1 = vec_xxcpnmadd(YY1,XX1,tmp1);
    tmp2 = vec_xxcpnmadd(YY2,XX2,tmp2);
    tmp3 = vec_xxcpnmadd(YY3,XX3,tmp3);
    tmp4 = vec_xxcpnmadd(YY4,XX4,tmp4);
    tmp5 = vec_xxcpnmadd(YY5,XX5,tmp5);

  }
  tmp0 = vec_add(tmp0,tmp1);
  tmp0 = vec_add(tmp0,tmp2);
  tmp0 = vec_add(tmp0,tmp3);
  tmp0 = vec_add(tmp0,tmp4);
  tmp0 = vec_add(tmp0,tmp5);
  vec_sta(tmp0,0,inner);
}
