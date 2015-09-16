#include "littleDiracOp.h"
#include <iostream>
#include <iomanip>


void LittleDiracOperator::OrthogonaliseSubspace(void)
{
  // Remove earlier (normalised) vectors
  for(int v=0;v<Nvec;v++){

    for(int sb=0;sb<NblockS;sb++){

      for(int u=0;u<v;u++){

	multi1d<Complex> bip(Nblock4);
	for(int i=0;i<Nblock4;i++) bip[i]= zero;

	//Inner product over whole sblock
	for(int ssb=0;ssb<BlockSize[4];ssb++){
	  int s = sb*BlockSize[4]+ssb;
	  LatticeComplex lip   = localInnerProduct(subspace[u][s],subspace[v][s]);
	  bip= bip+sumMulti(lip,Blocks);
	}

	LatticeComplex linearcomb;
	//	QDPIO::cout << "Inner product <"<<u<< "|" <<v<<">";
	for(int b=0;b<bip.size();b++){
	  linearcomb[Blocks[b]]=bip[b];
	  //	  QDPIO::cout << bip[b]<<" ";
	}
	//	QDPIO::cout << endl;
	for(int ssb=0;ssb<BlockSize[4];ssb++){
	  int s = sb*BlockSize[4]+ssb;
	  subspace[v][s] = subspace[v][s]-linearcomb* subspace[u][s];
	}
      }
    }
    
    // Normalise this vector
    for(int sb=0;sb<NblockS;sb++){
      multi1d<Real> bn(Nblock4); for(int i=0;i<bn.size();i++) bn[i] = zero;
      for(int ssb=0;ssb<BlockSize[4];ssb++){
	int s = sb*BlockSize[4]+ssb;
	LatticeReal ln   = localNorm2(subspace[v][s]);
	bn = bn+sumMulti(ln,Blocks);
      }
      LatticeReal scale;
      for(int b=0;b<bn.size();b++){
	scale[Blocks[b]] = 1.0/sqrt(bn[b]);
      }
      for(int ssb=0;ssb<BlockSize[4];ssb++){
	int s = ssb+sb*BlockSize[4];
	subspace[v][s] = subspace[v][s]*scale;
      }
    }
  }

  for(int v=0;v<Nvec;v++){
  for(int u=0;u<Nvec;u++){
    for(int sb=0;sb<NblockS;sb++){

      multi1d<Complex> bip(Nblock4);
      for(int i=0;i<Nblock4;i++) bip[i]=zero;
      for(int ssb=0;ssb<BlockSize[4];ssb++){
	int s = sb*BlockSize[4]+ssb;
	LatticeComplex lip   = localInnerProduct(subspace[u][s],subspace[v][s]);
	bip=bip+sumMulti(lip,Blocks);
      }

      Complex c=zero;
      if ( u== v ) { 
	c=1;
      }
      for(int i=0;i<bip.size();i++) { 
	if ( toDouble(norm2(bip[i]-c )) > 1.0e-9 ) {  
	  QDPIO::cout << u<< " " << v << " " << bip[i]<<endl;
	  QDPIO::cerr<<"Doesnt look orthogonal and should be"<<endl;
	  QMP_barrier();
	  exit(0);
	}
      }

    }
  }}
  for(int v=0;v<Nvec;v++){
    QDPIO::cout << norm2(subspace[v]) << endl;
  }
  QDPIO::cout<<"Orthogonalised and check orthogonal now"<<endl;
}
void LittleDiracOperator::PickBlock (const LatticeFermion &was, LatticeFermion &is,int b)
{
  is=zero;
  is[Blocks[b]] = was;
}
void LittleDiracOperator::PickBlock (const multi1d<LatticeFermion> &was, multi1d<LatticeFermion> &is,int b,int sb)
{
  for(int i=0;i<was.size();i++) is[i]=zero;
  for(int ss=0;ss<BlockSize[4];ss++){
    int s = sb*BlockSize[4]+ss;
    PickBlock(was[s],is[s],b);
  }
}
void LittleDiracOperator::ComputeLittleMatrix (
					       bfm_qdp<double> &dwf)
{
  multi1d<LatticeFermion> phi_i(N5);
  multi1d<LatticeFermion> Dphi_i(N5);
  multi1d<LatticeFermion> tmp(N5);
  std::vector<std::vector<std::complex<double> > > Acopy;
    for(int s=0;s<N5;s++){
      Dphi_i[s]=zero;
    }
  
  Fermion_t phi_t = dwf.allocFermion();
  Fermion_t tmp_t   = dwf.allocFermion();
  Fermion_t mmp   = dwf.allocFermion();
  Fermion_t mp    = dwf.allocFermion();
  for(int v_i=0;v_i<Nvec;v_i++){
    //    QDPIO::cout << "subspace " <<norm2(subspace[v_i],rb[0]) << " " << norm2(subspace[v_i],rb[1]) << endl;
    for(int s_i=0;s_i<NblockS;s_i++){
      for(int b_i=0;b_i<Nblock4;b_i++){

	//	QDPIO::cout<<"Compute for little matrix for block  v"<< v_i <<" s"<< s_i << "  b" << b_i<<endl;
	int i=SubspaceIndex(b_i,s_i,v_i); // subvector number

	// Inefficient first cut. Should use coloring to apply D to 
	// multiple blocks in parallel. Modest factor easily obtain
	// By doing this in the s-dimension only.
	// Presently applies D^dagD Nblock times
	// Big perf jump if simply throw this off to Bagel.
	PickBlock(subspace[v_i],phi_i,b_i,s_i);

#if 0
	//MdagM
	(*linop)(tmp,phi_i,PLUS);
	(*linop)(Dphi_i,tmp,MINUS);
#else
	dwf.importFermion(phi_i,phi_t,1);
#pragma omp parallel 
  {
    omp_set_num_threads(dwf.nthread);
#pragma omp for 
    for(int i=0;i<dwf.nthread;i++) {
	dwf.Mprec(phi_t,mp,tmp_t,0);
	dwf.Mprec(mp,mmp,tmp_t,1); 
    }
  }
	dwf.exportFermion(Dphi_i,mmp,1);
#endif
	//	QDPIO::cout << "phi " <<norm2(phi_i) << " " << norm2(Dphi_i) << endl;

	std::vector<std::complex<double> > proj(Nsubspace);
	ProjectToSubspace(Dphi_i,proj);
	for(int j=0;j<Nsubspace;j++){
	  A[j][i] = proj[j];
	}
      }
    }
  }
  QDPIO::cout << "A matrix"<<endl;
  for(int i=0;i<0;i++){
  for(int j=0;j<Nsubspace;j++){
    QDPIO::cout <<"A[][] "<< i<<" "<<j<<" "<<real(A[i][j]) << " "<< imag(A[i][j])<<endl;
  }
  //  QDPIO::cout << endl;
  } 
  QDPIO::cout << "A matrix (im)"<<endl;
  for(int i=0;i<Nsubspace;i++){
  for(int j=0;j<Nsubspace;j++){
    //    QDPIO::cout << imag(A[i][j]) << " ";
  }
  //  QDPIO::cout << endl;
  }
  Acopy=A;

  QDPIO::cout << "Asparse matrix "<<endl;
  for(int b=0;b<Nblock;b++){   // One per 5d-block
    for(int n=0;n<Nball;n++){
      int nn=nbr[b][n];
      for(int i=0;i<Nvec;i++){
	int si=b*Nvec+i;
	for(int j=0;j<Nvec;j++){
	  int sj=nn*Nvec+j;
	  Asparse[b][n][Nvec*i+j] = A[si][sj];
	  Acopy[si][sj]=0.0;
	}
      }
    }
    if ( b == 0 ) {
      for(int n=0;n<Nball;n++){
      for(int i=0;i<Nvec;i++){
	for(int j=0;j<Nvec;j++){
	  QDPIO::cout << "A block 0, n " << n<< " ij "<<i<<","<<j<< real(Asparse[b][n][j+Nvec*i]) <<" "<<imag(Asparse[b][n][j+Nvec*i])<<endl;
	}
      }
      }
    }
  }    

  if ( QMP_is_primary_node() ) { 
    FILE * fp = fopen("A.dat","w");
    for(int b=0;b<Nblock;b++){   // One per 5d-block
      for(int n=0;n<Nball;n++){
	std::complex<double> ip = Asparse[b][n][0];
	fprintf(fp,"A %d %d %le %le\n",b,n,real(ip),imag(ip));
      }
    }
    fclose(fp);
  }

  for(int i=0;i<Nsubspace;i++){
    for(int j=0;j<Nsubspace;j++){
      if ( abs(Acopy[i][j])>0 ) {
	QDPIO::cout << "Element of A not taken ij " << i<< " " << j<<" "<<real(Acopy[i][j]) << " " << imag(Acopy[i][j])<< endl;
	int gb[5];
	int nb[5];
	int idx=j;
	for(int mu=0;mu<4;mu++){
	  nb[mu]=Layout::lattSize()[mu]/BlockSize[mu];
	}
	nb[4]=N5/BlockSize[4];

	gb[0] = idx%nb[0]; idx = idx/nb[0];
	gb[1] = idx%nb[1]; idx = idx/nb[1];
	gb[2] = idx%nb[2]; idx = idx/nb[2];
	gb[3] = idx%nb[3]; idx = idx/nb[3];
	gb[4] = idx%nb[4];
	
	QDPIO::cout << "Other block is ";
	for(int mu=0;mu<5;mu++){
	  QDPIO::cout << gb[mu] << " ";
	}
	QDPIO::cout << endl;
	
      }
    }
  }

  dwf.freeFermion( phi_t );
  dwf.freeFermion( tmp_t );
  dwf.freeFermion( mmp   );
  dwf.freeFermion( mp    );

}
void LittleDiracOperator::ProjectToSubspace(multi1d<LatticeFermion> &vec, std::vector<std::complex<double> > &proj)
{
  proj.resize(Nsubspace);
  for(int i=0;i<Nsubspace;i++) proj[i]=0;
  for(int v_j=0;v_j<Nvec;v_j++){
    for(int s_j=0;s_j<NblockS;s_j++){ 

      // Local inner product route is reasonably efficient
      // Block ops performed on complex
      multi1d<Complex> bip(Nblock4) ;
      for(int i=0;i<Nblock4;i++) bip[i]=zero;

      for(int ss=0;ss<BlockSize[4];ss++){
	int s = ss+s_j*BlockSize[4];
	LatticeComplex lip   = localInnerProduct(subspace[v_j][s],vec[s]);
	bip = bip + sumMulti(lip,Blocks);
      }

      for(int b_j=0;b_j<Nblock4;b_j++){
	int j = SubspaceIndex(b_j,s_j,v_j);
	std::complex<double> cc (toDouble(real(bip[b_j])),toDouble(imag(bip[b_j])));
	proj[j] = cc;
      }
    }
  }
}
void LittleDiracOperator::PromoteFromSubspace(std::vector<std::complex<double> > &v, multi1d<LatticeFermion> &prom)
{
  LatticeFermion phi_i;
  prom=zero;
  for(int v_i=0;v_i<Nvec;v_i++){
    for(int s_i=0;s_i<NblockS;s_i++){
      LatticeComplex coeff;
      for(int b_i=0;b_i<Nblock4;b_i++){
	int i=SubspaceIndex(b_i,s_i,v_i);
	Complex cc = cmplx(Real(std::real(v[i])),Real(std::imag(v[i])));
	coeff[Blocks[b_i]] = cc;
      }
      for(int ss=0;ss<BlockSize[4];ss++){
	int s = ss+s_i*BlockSize[4];
	prom[s]+=coeff*subspace[v_i][s]; // Reasonably efficient way
      }
    }
  }
}

void myzgemv(int N,double * __restrict A, double * __restrict in,double c, double * __restrict out);
void myzgemv(int N,double * __restrict A, double * __restrict in,double c, double * __restrict out)
{
  register double tmp_r;
  register double tmp_i;
  register double A_r;
  register double A_i;
  register double in_r;
  register double in_i;
  register double *ap;
  register double *op;
  register double *ip;

  op = out;
  for(int i=0;i<N;i++){
    tmp_r = c*op[0];
    tmp_i = c*op[1];
    ap = &A[2*N*i];
    ip = &in[0];
    for(int j=0;j<N;j++){
      
      A_r = ap[0];
      A_i = ap[1];

      in_r= ip[0];
      in_i= ip[1];
      
      tmp_r = tmp_r + A_r * in_r - A_i * in_i;
      tmp_i = tmp_i + A_r * in_i + A_i * in_r;

      ip =ip+2;
      ap = ap+2;
    }
    op[0]=tmp_r;
    op[1]=tmp_i;
    op=op+2;
  }
}

void LittleDiracOperator::ApplySparse(std::vector<std::complex<double> > &in,std::vector<std::complex<double> > &out )
{

#pragma omp parallel 
  {
#pragma omp for 
  for(int b=0;b<Nblock;b++){
    for(int mu=0;mu<Nball;mu++){
      int nn = nbr[b][mu];
      double coeff=1.0;
      if ( mu==0) {
	coeff=0.0;
      }
      myzgemv(Nvec,(double *)&Asparse[b][mu][0],(double *)&in[nn*Nvec],coeff,(double *)&out[b*Nvec]);
    }
  }
  }
  /*
  for(int b=0;b<Nblock;b++){
    for(int i=0;i<Nvec;i++){
      out[i+b*Nvec]=0.0;
    }
    for(int mu=0;mu<Nball;mu++){
      int nn = nbr[b][mu];
      for(int i=0;i<Nvec;i++){
	for(int j=0;j<Nvec;j++){
	  if (b==0){
	    QDPIO::cout << " Asr " << real(Asparse[b][mu][j+Nvec*i]) << " phi " << real(in[j+nn*Nvec]) <<endl;
	    QDPIO::cout << " Asi " << imag(Asparse[b][mu][j+Nvec*i]) << " phi " << imag(in[j+nn*Nvec]) <<endl;
	  }

	  out[i+b*Nvec]+= Asparse[b][mu][j+Nvec*i]*in[j+nn*Nvec];
	}
      }
    }
  }
  */
}




void LittleDiracOperator::Apply(std::vector<std::complex<double> > &in,std::vector<std::complex<double> > &out)
{
  ApplySparse(in,out);
  return;
    
  if(0){

    for(int i=0;i<Nsubspace;i++){
      out[i]=0;
      for(int j=0;j<Nsubspace;j++){
	if ((i==0) && (real(A[i][j]) != 0.0) ) { 
	  QDPIO::cout << " Ar " << real(A[i][j]) << " phi " << real(in[j]) << "ij "<< i<< " " << j<<endl;
	  QDPIO::cout << " Ai " << imag(A[i][j]) << " phi " << imag(in[j]) <<endl;
	}
	out[i] += A[i][j] * in[j];
      }
    }
  
    std::vector<std::complex<double> > outtmp(Nsubspace);
    ApplySparse(in,outtmp);
    for(int i=0;i<Nsubspace;i++){
      if(abs(outtmp[i]-out[i])>1.0e-10 ){
	QDPIO::cout << "ApplySparse check : "<<i<<" "<<real(outtmp[i])<< " " << real(out[i])<<"     " <<imag(outtmp[i])<< " " << imag(out[i])<<endl;
      }
    }
  }
}
void LittleDiracOperator::ApplyInverse(std::vector<std::complex<double> > &v,std::vector<std::complex<double> > &vinv)
{
  // Conjugate gradient on A. A is hermitian posdef.
  std::vector<std::complex<double> >  p(Nsubspace);
  std::vector<std::complex<double> > Ap(Nsubspace);
  std::vector<std::complex<double> >  r(Nsubspace);

  for(int i=0;i<Nsubspace;i++){
    vinv[i]=0;
    r[i]=v[i];
    p[i]=r[i];
  }
  double a;
  double b;
  double c;
  double d;
  double cp;
  double ssq;
  double rsq;

  ssq=norm_vec(v);
  a  =norm_vec(p);
  cp =norm_vec(r);
  rsq = 1.0e-10*ssq;
 
  QDPIO::cout << "Source^2 = " << ssq<<endl;
  QDPIO::cout << "Target^2 = " << rsq<<endl;

  for(int k=1;k<10000;k++){
    c=cp;
    Apply(p,Ap);
    d=innerProductReal(p,Ap);
    a=c/d;
    for(int i=0;i<Nsubspace;i++){
      r[i]=r[i]-a*Ap[i];
    }
    cp=norm_vec(r);
    b=cp/c;
    for(int i=0;i<Nsubspace;i++){
      vinv[i]=vinv[i]+a*p[i];
      p[i]=b*p[i]+r[i];
    }
    QDPIO::cout << "LittleDiracInversion ["<<k<<"] "<<cp<<"/"<<rsq<<endl;
    if(cp<rsq) {
      QDPIO::cout << "LittleDiracInversion converged "<<endl;
      Apply(vinv,Ap);
      for(int i=0;i<Nsubspace;i++){
	Ap[i]=Ap[i]-v[i];
      }
      QDPIO::cout << "LittleDiracInversion true residual "<<norm_vec(Ap) <<endl;
      return;
    }
  }
  QDP_error_exit("Little Dirac Matrix inversion failed\n");
}
