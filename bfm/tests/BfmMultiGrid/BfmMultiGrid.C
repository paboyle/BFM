#include "BfmMultiGrid.h"
#include <iostream>
#include <iomanip>

int MGreport;
double SubspaceRationalLo=0.001;
double SubspaceRationalResidual=1.0e-6;

double SubspaceChebyLo = 0.001;
double SubspaceChebyHi = 100.0;
int    SubspaceChebyOrder = 800; 
double SubspaceChebyPower = -4.0;

int SubspaceSurfaceDepth = 256;

double CoarseSubspaceChebyLo = 0.002;
double CoarseSubspaceChebyHi = 200.0;
int    CoarseSubspaceChebyOrder = 800; 
double CoarseSubspaceChebyPower = -2.0;

int  SubspaceRelaxChebyshev=0;


int SolverControl=SolverDeflatedMdagMpc;
double LittleDopSolverResidual=1.0e-11;
int    LittleDopSolver=LittleDopSolverDeflCG;

double InnerKrylovResidual = 1.0e-6;
int    InnerKrylovIterMax  = 10;
double InnerKrylovShift= 0.2;


double PreconChebyLo = 0.3;
double PreconChebyHi = 100.0;
int    PreconChebyOrder = 10; 
double PreconChebyPower = 0.0;

//double PreconChebyPower = 0.0;
void BfmMultiGrid::ReadAlgorithmParms(void)
{
  FILE * fp = fopen("MultiGridParms.txt","r");
  fscanf(fp,"%le %le %le %d",
	  &InnerKrylovResidual,
	  &InnerKrylovShift,
	  &SubspaceRationalLo,
	  &InnerKrylovIterMax);
  fclose(fp);
  if ( linop->isBoss() ) { 
    printf("**** Agorithm parms ******\n");
    printf("**** InnerKrylovResidual       %le\n",InnerKrylovResidual);
    printf("**** InnerKrylovIterMax        %d \n",InnerKrylovIterMax);
    printf("**** InnerKrylovShift          %le\n",InnerKrylovShift);
    printf("**** SubspaceRationalLo        %le\n",SubspaceRationalLo);
    printf("**** SubspaceRationalResidual  %le\n",SubspaceRationalResidual);
  }
}

// Pass a linop and solve in single
template<class Float>
void BfmMultiGrid::RationalSubspace(bfm_qdp<Float> * rop,Fermion_t src,Fermion_t sol)
{
  const int nshift=4;
  Fermion_t sol_guess[nshift];

  double Lo = SubspaceRationalLo;
  //
  // [ 1/6(x+Lo)  - 1/2(x+2Lo) + 1/2(x+3Lo)  -1/6(x+4Lo) = Lo^3 /[ (x+1Lo)(x+2Lo)(x+3Lo)(x+4Lo) ]
  //
  // 1/(x+Lo)  - 1/(x+2 Lo)
  //
  double alpha[4]     = {1.0,1.0,1.0,1.0};
  double shifts[4]    = {Lo,2.0*Lo,3.0*Lo,4.0*Lo};
  double mresidual[4] = {3.0*SubspaceRationalResidual,SubspaceRationalResidual,SubspaceRationalResidual,SubspaceRationalResidual};

  for(int i=0;i<nshift;i++){
    sol_guess[i] = rop->threadedAllocFermion();
  }
  int power=1;
  for(int i=0;i<power;i++){
    int single = 0;
    rop->CGNE_prec_MdagM_multi_shift(sol_guess,
				       src,
				       shifts,
				       alpha,
				       nshift,
				       mresidual,
				       single);

    rop->axpby(sol,sol_guess[0],sol_guess[1],1.0/6.0,-1.0/2.0);
    rop->axpy(sol,sol_guess[2],sol,1.0/2.0);
    rop->axpy(sol,sol_guess[3],sol,-1.0/6.0);

    rop->axpy(src,sol,sol,0.0);
  }

  for(int i=0;i<nshift;i++){
    rop->threadedFreeFermion(sol_guess[i]);
  }

}

void BfmMultiGrid::AugmentSubspace(void)
{
  int Ls=N5;
  multi1d<LatticeFermion> gauss(Ls);
  int vec_all =Nvec;
  int vec_skip=1;

  if ( linop->SPIcomms() ) linop->comm_init();

  subspace = (Fermion_t *) malloc(vec_all*sizeof(Fermion_t));
  
  int threads   = linop->threads;
  Fermion_t src = linop->allocFermion();
  Fermion_t sol = linop->allocFermion();
  Fermion_t res = linop->allocFermion();
  double restore_residual = linop->residual;

  for(int v=0;v<Nvec;v++){
    subspace[v] = linop->allocFermion();
    for(int s=0;s<Ls;s++) gaussian(gauss[s]);
    linop->importFermion(gauss,subspace[v],1); // Fill with noise to begin with
  }
  Fermion_t tmp[vec_skip];
  for(int v=0;v<vec_skip;v++) tmp[v] = linop->allocFermion();


  for(int nv=0;nv<vec_all;nv+=vec_skip){

    for(int v=0;v<vec_skip;v++){

      int vi = v+nv;

      if (nv==0) {

	linop->residual = 1.0e-6;

	for(int s=0;s<Ls;s++) gauss[s]=zero;
	linop->importFermion(gauss,sol,1);
	
	gaussian(gauss[0]);
	gaussian(gauss[Ls-1]);
	linop->importFermion(gauss,src,1);

#pragma omp parallel 
	{
#pragma omp for 
	  for(int i=0;i<threads;i++) {
	    RationalSubspace(linop,src,sol);
	    linop->axpy(tmp[v],sol,sol,0.0);
	  }
	}

      } else { 

	linop->residual = 5.0e-5;

	for(int s=0;s<Ls;s++) gauss[s]=zero;
	gaussian(gauss[0]);
	gaussian(gauss[Ls-1]);
	linop->importFermion(gauss,src,1);
	
	for(int s=0;s<Ls;s++) gaussian(gauss[s]);
	linop->importFermion(gauss,sol,1);

#pragma omp parallel 
	{
#pragma omp for 
	  for(int i=0;i<threads;i++) {
	    MCR_PolyPrec(sol,src,res);
	    linop->axpy(tmp[v],res,res,0.0);
	  }
	}

      }
    }

#pragma omp parallel 
	{
#pragma omp for 
	  for(int i=0;i<threads;i++) {
	    for(int v=0;v<vec_skip;v++){
	      linop->axpy(subspace[v+nv],tmp[v],tmp[v],0.0);
	    }
	  }
	}

    OrthogonaliseSubspace();

    if ( linop->SPIcomms() ) linop->comm_init();

#pragma omp parallel 
    {
#pragma omp for 
      for(int i=0;i<threads;i++) {
	ComputeLittleMatrixColored();
      }
    }

  }
  linop->residual = restore_residual;
  for(int v=0;v<vec_skip;v++) linop->freeFermion(tmp[v]);
  linop->freeFermion(res);
  linop->freeFermion(sol);
  linop->freeFermion(src);
}

#if 0
template<class Float>
void BfmMultiGrid::RelaxSubspace(bfm_qdp<Float> *rop)
{
  int Ls=N5;
  multi1d<LatticeFermion> gauss(Ls);

  Fermion_t sol = linop->allocFermion();
  Fermion_t src = linop->allocFermion();
  subspace = (Fermion_t *) malloc(2*Nvec*sizeof(Fermion_t));

  if ( linop->SPIcomms() ) linop->comm_init();

  QDPIO::cout << "Chebyshev: Building subspace using rational  filter 1/(x+lambda)^2 for lambda= "
		<<SubspaceRationalLo << " tol " << SubspaceRationalResidual<<endl;

  QDPIO::cout << "Subdividing blocks into "<<Nquadrant<<" quadrants"<<endl;

  for(int v=0;v<Nvec;v+=Nquadrant){

    for(int s=0;s<Ls;s++) gauss[s]=zero;
    gaussian(gauss[0]);
    gaussian(gauss[Ls-1]);
    linop->importFermion(gauss,src,1);

    double nn;
#pragma omp parallel 
    {
#pragma omp for 
      for(int i=0;i<linop->threads;i++) {
	RationalSubspace(linop,src,sol);
	nn=linop->norm(sol);
      }
    }
    if ( linop->isBoss() ) printf("Relax rational norm %le\n",nn);

    for(int i=0;i<Nquadrant;i++){
      int vi=v+i;
      subspace[vi] = linop->allocFermion();
      PickQuadrant(sol,subspace[vi],i);
      s_min[vi]=0; s_max[vi]=Ls-1;
    }
  }

  linop->freeFermion(src);
  linop->freeFermion(sol);

  if ( linop->isBoss() ) {
    printf("Relax Subspace Orthogonalising\n"); fflush(stdout);
  }

  OrthogonaliseSubspace();
}
#else
template<class Float>
void BfmMultiGrid::RelaxSubspace(bfm_qdp<Float> *rop)
{
  int Ls=N5;
  multi1d<LatticeFermion> gauss(Ls);

  Fermion_t sol = rop->allocFermion();
  Fermion_t src = rop->allocFermion();
  subspace = (Fermion_t *) malloc(2*Nvec*sizeof(Fermion_t));

  if ( rop->SPIcomms() ) rop->comm_init();

  if ( SubspaceRelaxChebyshev ) {
    QDPIO::cout << "Chebyshev: Building subspace using chebyshev filter order "<<SubspaceChebyOrder << " x^"
		<< SubspaceChebyPower<<"  ["<<SubspaceChebyLo<<","<<SubspaceChebyHi<<"]" <<endl;
  } else {
    QDPIO::cout << "Chebyshev: Building subspace using rational  filter 1/(x+lambda)^2 for lambda= "
		<<SubspaceRationalLo << " tol " << SubspaceRationalResidual<<endl;
  }

#define VSKIP (1)

  for(int v=0;v<Nvec;v+=VSKIP){

    for(int s=0;s<Ls;s++) gauss[s]=zero;
    gaussian(gauss[0]);
    gaussian(gauss[Ls-1]);

    rop->importFermion(gauss,src,1);

#pragma omp parallel 
    {
#pragma omp for 
      for(int i=0;i<rop->threads;i++) {
	if ( SubspaceRelaxChebyshev ) { 

	  ChebyshevInit(SubspaceChebyLo,
			SubspaceChebyHi,
			SubspaceChebyPower,
			SubspaceChebyOrder);

	  PolyMdagMprec(src,sol);
	} else { 
	  RationalSubspace(rop,src,sol);
	}

	double nn=rop->norm(sol);
	if ( rop->isBoss() && (i==0)) {
	  printf("Relax: vector %d Polynomial norm %le\n",v,nn);
	}
      }
    }

    for(int i=0;i<VSKIP;i++){

      for(int s=0;s<Ls;s++) gauss[s]=zero;

      int vi=v+i;
      subspace[vi] = linop->allocFermion();
      rop->exportFermion(gauss,sol,1);    // + rest


      if ( VSKIP == 4 ) {

	switch(i) {
	case 0:
	  s_min[vi]=0; s_max[vi]=0;
	  break;
	case 1:
	  s_min[vi]=Ls-1; s_max[vi]=Ls-1;
	  break;
	case 2:
	  s_min[vi]=1;    s_max[vi]=2;
	  break;
	case 3:
	  s_min[vi]=Ls-3; s_max[vi]=Ls-2;
	  break;
	}
	for(int s=0;s<s_min[vi];s++)    gauss[s]=zero;
	for(int s=s_max[vi]+1;s<Ls;s++) gauss[s]=zero;

      } else if ( VSKIP == 2 ) { 

	if (i==0) {
	  s_min[vi]=0;
	  s_max[vi]=Ls/2-1;
	  for(int s=s_max[vi+1];s<Ls-1;s++) gauss[s]=zero;
	} else { 	    
	  s_min[vi]=Ls/2;
	  s_max[vi]=Ls-1;
	  for(int s=0;s<s_min[vi];s++) gauss[s]=zero;
	}


      } else { 
	if ( SubspaceSurfaceDepth > Ls/2 ) { 
	  s_min[vi]=0; s_max[vi]=Ls-1;
	} else { 
	  s_min[vi]=Ls-SubspaceSurfaceDepth; s_max[vi]=SubspaceSurfaceDepth-1;
	  for(int s=s_max[vi]+1;s<s_min[vi];s++) gauss[s]=zero;
	}
      }
      for(int s=0;s<Ls;s++) { 
	QDPIO::cout << "Vec["<<vi<<"][s=" <<s<<"] == "<<norm2(gauss[s])<<endl;
      }
      linop->importFermion(gauss,subspace[vi],1);
    }
  }

  rop->freeFermion(src);
  rop->freeFermion(sol);

  if ( rop->isBoss() ) {
    printf("Relax Subspace Orthogonalising\n"); fflush(stdout);
  }

  OrthogonaliseSubspace();

  if ( linop->SPIcomms() ) linop->comm_init();
  
}
#endif

void BfmMultiGrid::AddToSubspace(Fermion_t vec)
{	

  int Ls=N5;
  int vi=Nvec-1;// Replace last vec

#pragma omp parallel 
    {
#pragma omp for 
      for(int i=0;i<linop->threads;i++) {
	linop->axpy(subspace[vi],vec,vec,0.0);
      }
    }

    s_min[vi]=0; s_max[vi]=Ls-1;
    OrthogonaliseSubspace();

#pragma omp parallel 
    {
#pragma omp for 
      for(int i=0;i<linop->threads;i++) {
	ComputeLittleMatrixColored();  // Expensive recompute matrix
      }
    }
}

template<class Float> 
void BfmMultiGrid::PolyMdagMprec(bfm_qdp<Float> *lop,Fermion_t in,Fermion_t out)
{

    Fermion_t y=    lop->threadedAllocFermion();
    Fermion_t Mtmp= lop->threadedAllocFermion();
    Fermion_t Tnm = lop->threadedAllocFermion();
    Fermion_t Tn  = lop->threadedAllocFermion();
    Fermion_t Tnp = lop->threadedAllocFermion();
    Fermion_t tmp=Tnp;

    int me = lop->thread_barrier();

    double xscale = 2.0/(hi-lo);
    double mscale = -(hi+lo)/(hi-lo);
    
    lop->axpy(Tnm,in,in,0.0);              // Tnm=T0=in

    lop->Mprec(in ,tmp,Mtmp,DaggerNo);  
    lop->Mprec(tmp,y  ,Mtmp,DaggerYes);    
    lop->axpby(Tn ,y,in,xscale,mscale);    // TN=T1 = (xscale MdagM + mscale)in 

    lop->axpby(out,Tnm,Tn,coeffs[0]/2.0,coeffs[1]); // sum = .5c[0]T0 + c[1]T1
    
    for(int i=2;i<order;i++){
      
      lop->Mprec(Tn,tmp,Mtmp,DaggerNo);  
      lop->Mprec(tmp, y,Mtmp,DaggerYes);    
      lop->axpby(y,y,Tn,xscale,mscale);    // y Tn = [xscale MdagM+mscale] Tn
      lop->axpby(Tnp,y,Tnm,2.0,-1.0);      // Tnp=2yTn - Tnm

      lop->axpy(Tnm,Tn,Tn,0.0); 
      lop->axpy(Tn,Tnp,Tnp,0.0);
      lop->axpy(out,Tn,out,coeffs[i]);//Accumulate

    }
    lop->threadedFreeFermion(y);
    lop->threadedFreeFermion(Mtmp);
    lop->threadedFreeFermion(Tnm);
    lop->threadedFreeFermion(Tn);
    lop->threadedFreeFermion(Tnp);
};

template void BfmMultiGrid::PolyMdagMprec<float>(bfm_qdp<float> *lop,Fermion_t in,Fermion_t out);
template void BfmMultiGrid::RelaxSubspace<float>(bfm_qdp<float> *rop);
template void BfmMultiGrid::RationalSubspace<float>(bfm_qdp<float> * rop,Fermion_t src,Fermion_t sol);

template void BfmMultiGrid::PolyMdagMprec<double>(bfm_qdp<double> *lop,Fermion_t in,Fermion_t out);
template void BfmMultiGrid::RelaxSubspace<double>(bfm_qdp<double> *rop);
template void BfmMultiGrid::RationalSubspace<double>(bfm_qdp<double> * rop,Fermion_t src,Fermion_t sol);

void BfmMultiGrid::SolverMatrix(Fermion_t in,Fermion_t out)
{


 switch ( SolverControl ) { 
 case SolverDeflatedMdagMpc: 
    DeflatedSolverMdagMpc(in,out);
    break;
 case SolverUndeflatedMdagMpc:  
   UndeflatedSolverMdagMpc(in,out);
   break;
 case SolverUndeflatedMpc:  
    UndeflatedSolverMpc(in,out);
   break;
 case SolverDeflatedMpc:
    DeflatedSolverMpc(in,out);
   break;
 default: 
   exit(-1);
   break;
 }

}



BfmMultiGrid::BfmMultiGrid(int _N5,
			   int _Nvec,
			   int _BlockSize[5],
			   int _QuadrantSize[4],
			   bfm_t * _linop) 
{
  Nvec=_Nvec;
  N5=_N5;
  s_min.resize(Nvec);
  s_max.resize(Nvec);
  for(int i=0;i<Nvec;i++){
    s_min[i]=0;
    s_max[i]=N5-1;
  }
  
  for(int mu=0;mu<5;mu++){
    BlockSize[mu]=_BlockSize[mu];
  }
  for(int mu=0;mu<4;mu++){
    QuadrantSize[mu]=_QuadrantSize[mu];
  }
  linop=_linop;
  preclinop=_linop; // By default use same linop

  QDPIO::cout << "Building subspace from " << Nvec<<" vectors"<<endl;


  // Block layout
  for(int mu=0;mu<4;mu++){
    BlockSize[mu] = _BlockSize[mu];
    glatt[mu]     = Layout::lattSize()[mu];
    global_nb[mu] = glatt[mu]/BlockSize[mu];
    if ( global_nb[mu]*BlockSize[mu] !=  Layout::lattSize()[mu] ) {
      QDPIO::cout<<"global BlockSize err"<<endl;
      exit(0);
    }
    local_nb[mu] = QDP::Layout::subgridLattSize()[mu]/BlockSize[mu];
    if ( local_nb[mu]*BlockSize[mu] !=  QDP::Layout::subgridLattSize()[mu] ) {
      QDPIO::cout<<"local BlockSize err "<< mu << " " << local_nb[mu] << " "<< BlockSize[mu] << " "  <<endl;
      exit(0);
    }
  }
  QDPIO::cout << "4d block size"<<endl;

  // Fifth dimension
  NblockS= N5/BlockSize[4];  
  global_nb[4]=NblockS;
  local_nb[4] =NblockS;
  if ( NblockS*BlockSize[4]!=N5 ) {
    QDP_error_exit("BlockSize err");
  }

  // Dimensions of subspace blocks (globally)
  Nblock4= global_nb[0]*global_nb[1]*global_nb[2]*global_nb[3];
  Nblock =NblockS*Nblock4;

  Nball=1;
  for(int d=0;d<5;d++){
    if ( global_nb[d] == 1 ) {
      stencil_lo[d] = 0;
      stencil_hi[d] = 0;
      stencil_size[d]= 1;
    } else if ( global_nb[d] == 2 ) {
      stencil_lo[d] = -1;
      stencil_hi[d] = 0;
      stencil_size[d]= 2;
    } else if ( global_nb[d] > 2 ) {
      stencil_lo[d] = -1;
      stencil_hi[d] =  1;
      stencil_size[d]= 3;
    }
    Nball=Nball*stencil_size[d];  // Keep matrix representation local
  }
  StencilNonZero.resize(Nball);


  // Dimensions of subspace blocks (locally)
  LocalNblock4= local_nb[0]*local_nb[1]*local_nb[2]*local_nb[3];
  LocalNblock = NblockS*LocalNblock4;
  Nsubspace=Nblock*Nvec;
  LocalNsubspace=LocalNblock*Nvec;

  ProjectInnerProduct.resize(LocalNblock);
  PromoteCoeff.resize(LocalNblock);
  ComputeProj.resize(Nball);
  for(int b=0;b<Nball;b++){
    ComputeProj[b].resize(LocalNsubspace);
  }
  phases.resize(LocalNblock);

  slatt[0]=QDP::Layout::subgridLattSize()[0];
  slatt[1]=QDP::Layout::subgridLattSize()[1];
  slatt[2]=QDP::Layout::subgridLattSize()[2];
  slatt[3]=QDP::Layout::subgridLattSize()[3];
  slatt[4]=N5;

  ncoor[0] = QDP::Layout::nodeCoord()[0];
  ncoor[1] = QDP::Layout::nodeCoord()[1];
  ncoor[2] = QDP::Layout::nodeCoord()[2];
  ncoor[3] = QDP::Layout::nodeCoord()[3];
  ncoor[4] = 0;

  int vol   = slatt[0]*slatt[1]*slatt[2]*slatt[3]*slatt[4];
  int cbvol = vol/2;

  inner_reduce.resize(cbvol);
  QDPIO::cout << "BlockSize "
	      << BlockSize[0] << "x"
	      << BlockSize[1] << "x"
	      << BlockSize[2] << "x"
	      << BlockSize[3] << "x"
	      << BlockSize[4] << endl;

  QDPIO::cout << "Making integer tables"<<endl;

  Asparse.resize(LocalNblock);
  for(int b=0;b<LocalNblock;b++){   // One per 5d-block
    Asparse[b].resize(Nball);       // One per neigbour
    for(int n=0;n<Nball;n++){
      Asparse[b][n].resize(Nvec*Nvec);
    }
  }    

  // Preallocate vectors needed in threaded routines
  Krylov_p.resize(LocalNsubspace) ;
  Krylov_Ap.resize(LocalNsubspace) ;
  Krylov_Ar.resize(LocalNsubspace) ;
  Krylov_r.resize(LocalNsubspace) ;
  Krylov_mu.resize(LocalNsubspace) ;
  Krylov_Amu.resize(LocalNsubspace) ;
  Krylov_Atmp.resize(LocalNsubspace*Nball) ;
  PleftProj.resize(LocalNsubspace) ;
  PleftMss_proj.resize(LocalNsubspace) ;

  block_id.resize(2);
  quadrant_id.resize(2);
  local_block_id.resize(2);
  for(int cb=0;cb<2;cb++){
    block_id[cb].resize(cbvol);
    local_block_id[cb].resize(cbvol);
    quadrant_id[cb].resize(cbvol);
  }

  std::vector<std::vector< std::vector<int> > > coordinate(5);
  for(int mu=0;mu<5;mu++){
    coordinate[mu].resize(2);
    coordinate[mu][0].resize(cbvol); 
    coordinate[mu][1].resize(cbvol);
  }
    

  int nq[4];
  int qs[4];
  for(int mu=0;mu<4;mu++){
    nq[mu] = BlockSize[mu]/QuadrantSize[mu];
    qs[mu] = QuadrantSize[mu];
    if( nq[mu]*QuadrantSize[mu] != BlockSize[mu] ) { 
      QDPIO::cout << "Quadrant does not divide block"<<endl;
      nq[mu]++;
    }
  }
  Nquadrant = nq[0]*nq[1]*nq[2]*nq[3];

  QDPIO::cout << "Mapping blockids"<<endl;
  for(int site=0;site<vol;site++){
      
    int ss=site;
    int x[5];
    x[0]=ss%slatt[0];    ss=ss/slatt[0];
    x[1]=ss%slatt[1];    ss=ss/slatt[1];
    x[2]=ss%slatt[2];    ss=ss/slatt[2];
    x[3]=ss%slatt[3];    ss=ss/slatt[3];
    x[4]=ss%slatt[4];
	
    int sp;
    if ( linop->precon_5d ) sp=x[4];
    else sp=0;
      
    int cb = ((x[0]+x[1]+x[2]+x[3]+sp)&0x1) ;
    int cbsite = linop->bagel_idx5d(x,x[4],0,0,1,1); cbsite=cbsite>>1; //wipes out reim factor of 2

    int b[5];
    int lb[5];
    int xb[5];

    for(int mu=0;mu<5;mu++){
      coordinate[mu][cb][cbsite] = x[mu]+ncoor[mu]*slatt[mu]; //Global coordinate
      b[mu]=coordinate[mu][cb][cbsite]/BlockSize[mu];         //Global block
      lb[mu]=x[mu]/BlockSize[mu];                         //Local block
      xb[mu]=x[mu]%BlockSize[mu];
    }

    // Needed for constructing little dirac op by picking subvectors
    block_id[cb][cbsite] = b[0]+global_nb[0]*(b[1]+global_nb[1]*(b[2]+global_nb[2]*(b[3]+global_nb[3]*b[4])));

    quadrant_id[cb][cbsite] =                   (xb[0]/QuadrantSize[0]) ;
    quadrant_id[cb][cbsite]+=             nq[0]*(xb[1]/QuadrantSize[1]) ;
    quadrant_id[cb][cbsite]+=       nq[1]*nq[0]*(xb[2]/QuadrantSize[2]) ;
    quadrant_id[cb][cbsite]+= nq[2]*nq[1]*nq[0]*(xb[3]/QuadrantSize[3]) ;

    //      ball_id_cb[cb][cbsite]= ball_id[cb][cbsite] + ball_cb[cb][cbsite]*32;

    local_block_id[cb][cbsite]=lb[0]+local_nb[0]*(lb[1]+local_nb[1]*(lb[2]+local_nb[2]*(lb[3]+local_nb[3]*lb[4])));
    
  }

  QDPIO::cout << "Reading algorithm parameters"<<endl;
  ReadAlgorithmParms();

  QDPIO::cout << "Halo Init"<<endl;
  HaloInit();

  QDPIO::cout << "Initialised Little Dirac Operator"<<endl;
}

int BfmMultiGrid::GetLocalBlockIdx(int d[5])
{
  return d[0] + local_nb[0]*(d[1]+local_nb[1]*(d[2]+local_nb[2]*(d[3]+local_nb[3]*d[4])));
}


void BfmMultiGrid::GlobalToLocalBlock(int gb[5],int &proc,int lb[5])
{
  multi1d<int> gcoor(4);
  for(int mu=0;mu<5;mu++)   lb[mu]=gb[mu]%local_nb[mu];
  for(int mu=0;mu<4;mu++)   gcoor[mu] = gb[mu]*BlockSize[mu];
  proc = QDP::Layout::nodeNumber(gcoor);  
}
int BfmMultiGrid::ReverseDirection(int mu)
{
  int d[5];
  int rd[5];
  d[0] = mu % stencil_size[0]; mu = mu/stencil_size[0];
  d[1] = mu % stencil_size[1]; mu = mu/stencil_size[1];
  d[2] = mu % stencil_size[2]; mu = mu/stencil_size[2];
  d[3] = mu % stencil_size[3]; mu = mu/stencil_size[3];
  d[4] = mu % stencil_size[4]; 
  for(int dir=0;dir<5;dir++){
    if ( stencil_size[dir]==2 )
      rd[dir] = stencil_size[dir]-d[dir];
    else 
      rd[dir] = d[dir];
  }
  mu = rd[4]+stencil_size[4]*(rd[3]+stencil_size[3]*(rd[2]+stencil_size[2]*(rd[1]+stencil_size[1]*rd[0])));
  return mu;
}

void BfmMultiGrid::GetDelta(int delta[5],int _ballidx,int rev)
{
  int d[5];
  int ballidx=0;

  for(d[0]=stencil_lo[0];d[0]<=stencil_hi[0];d[0]++){
    for(d[1]=stencil_lo[1];d[1]<=stencil_hi[1];d[1]++){
      for(d[2]=stencil_lo[2];d[2]<=stencil_hi[2];d[2]++){
	for(d[3]=stencil_lo[3];d[3]<=stencil_hi[3];d[3]++){
	  for(d[4]=stencil_lo[4];d[4]<=stencil_hi[4];d[4]++){
	    if ( _ballidx == ballidx ) {
	      if( rev ) 
		for(int mu=0;mu<5;mu++) {
		  if ( stencil_size[mu]==2 ) delta[mu] = -d[mu];
		  else delta[mu] = d[mu];
		}
	      else for(int mu=0;mu<5;mu++) delta[mu] = d[mu];
	      return;
	    }
	    ballidx++;
	  }
	}
      }
    }
  }
  exit(0);
}

int BfmMultiGrid::NeighbourBlockId(int myblock,int _mu, int &me,
				   int &from,int &xblock,int &to,int rev)
{
  int idx = myblock;
  
  int b[5];// Local block coord of this site
  b[0] = idx%local_nb[0]; idx = idx/local_nb[0];
  b[1] = idx%local_nb[1]; idx = idx/local_nb[1];
  b[2] = idx%local_nb[2]; idx = idx/local_nb[2];
  b[3] = idx%local_nb[3]; idx = idx/local_nb[3];
  b[4] = idx%local_nb[4];

  int gb[5]; // Global block coordinate
  for(int mu=0;mu<5;mu++)    gb[mu] = b[mu]+ncoor[mu]*local_nb[mu];

  int lwb[5]; 
  int seeker[5]; 

  int delta[5]; 
  GetDelta(delta,_mu,rev); 

  // Block plus delta
  int gbpd[5];
  for(int mu=0;mu<5;mu++){
    // Periodic wrap at global lattice boundaries
    gbpd[mu] = gb[mu]+delta[mu];
    if ( gbpd[mu] >= global_nb[mu] ) gbpd[mu]-=global_nb[mu];
    if ( gbpd[mu] <  0             ) gbpd[mu]+=global_nb[mu];

    // Locally wrap internal to node. Someone must be also seeking this point
    lwb[mu]= b[mu]+delta[mu];
    if ( lwb[mu] >= local_nb[mu] ) lwb[mu]-=local_nb[mu];
    if ( lwb[mu] <  0            ) lwb[mu]+=local_nb[mu];
    
    // Who would be looking for that
    seeker[mu] = lwb[mu]+ncoor[mu]*local_nb[mu] - delta[mu];
    if ( seeker[mu]>= global_nb[mu] ) seeker[mu]-=global_nb[mu];
    if ( seeker[mu]<0 )               seeker[mu]+=global_nb[mu];

  }

  int fb[5];
  int tb[5];
  int mb[5];
  GlobalToLocalBlock(gbpd,from,fb);
  GlobalToLocalBlock(seeker,to,tb);

  GlobalToLocalBlock(gb,me,mb);
  int mecheck= QDP::Layout::nodeNumber();
  if( me != mecheck) QDP_error_exit("oops");

  xblock   =GetLocalBlockIdx(fb);
  int lwblock =GetLocalBlockIdx(lwb);
  if ( xblock != lwblock ) QDP_error_exit("block oops");

  int meblock= GetLocalBlockIdx(mb);
  if ( meblock != myblock ) {
    cout << "logic error meblock = "<<meblock<<" by myblock="<< myblock<<endl;
    exit(0);
  }

}
void BfmMultiGrid::HaloEnd(void)
{
  int Nmsg = tproc.size();
  for(int m=0;m<Nmsg;m++) {
    QMP_free_msghandle(xmit_msghandle[m]);
    QMP_free_msghandle(recv_msghandle[m]);
    QMP_free_msgmem(xmit_msgmem[m]);
    QMP_free_msgmem(recv_msgmem[m]);
#ifndef BFM_BGQ_SPI_DSLASH
    bfm_free(sendbufs[m]);
    bfm_free(recvbufs[m]);
#endif
  }
}
void BfmMultiGrid::HaloInit(void)
{
  tproc.resize(0);
  fproc.resize(0);

  sendbufs.resize(Nball);
  recvbufs.resize(Nball);
  
  sendbuf_bytes.resize(Nball);
  recvbuf_bytes.resize(Nball);
  for(int mu=0;mu<Nball;mu++){
    sendbuf_bytes[mu]=0;
    recvbuf_bytes[mu]=0;
  }
 
  sbuf_id.resize (LocalNblock*Nball); 
  sbuf_idx.resize(LocalNblock*Nball); 
  rbuf_id.resize (LocalNblock*Nball); 
  rbuf_idx.resize(LocalNblock*Nball); 
 
  //Build and Fill halo buffers
  QDPIO::cout << "Initialising haloes "<<endl;
  for(int idx=0;idx<LocalNblock;idx++){
    for(int mu=0;mu<Nball;mu++){
 
      int me,from,fromblock,to;
      int rev=0;
 
      NeighbourBlockId(idx,mu,me,from, fromblock,to,rev);
 
      rbuf_id [idx*Nball+mu] = -1; // Indicates local to node
      rbuf_idx[idx*Nball+mu] = fromblock*Nvec;
 
      int fidx=-1;
      if ( me != from ) {
	for(int p=0;p<fproc.size();p++){
	  if (fproc[p]==from) fidx=p;
	}
	if(fidx==-1) { 
	  fproc.push_back(from);
	  fidx = fproc.size()-1;
	}
 	rbuf_id [idx*Nball+mu] = fidx;
 	rbuf_idx[idx*Nball+mu] = recvbuf_bytes[fidx]/(2*sizeof(double));
 	for(int i=0;i<Nvec;i++){
	  recvbuf_bytes[fidx]+=2*sizeof(double);
	}
      }
      
      int tidx=-1;
      if ( me != to ) {
	for(int p=0;p<tproc.size();p++){
 	  if (tproc[p]==to) tidx=p;
	}
	if(tidx==-1) {
	  std::vector<int> empty(0); // gather table for filling buffer
	  sbuf_gather.push_back(empty);
	  tproc.push_back(to);
	  tidx = tproc.size()-1;
	}
 	sbuf_id [idx*Nball+mu] = tidx;
 	sbuf_idx[idx*Nball+mu] = sendbuf_bytes[tidx]/(2*sizeof(double));
	sbuf_gather[tidx].push_back(fromblock);
	for(int i=0;i<Nvec;i++){
	  sendbuf_bytes[tidx]+=2*sizeof(double);
	}      
      }
    } //mu
  } // idx


  if ( tproc.size() != fproc.size() ) {
    QDP_error_exit("Halo exchange logic bomb");
  }

  HaloInitBuffers(this);
 
  // QMP init
  int Nmsg = tproc.size();
  xmit_msgmem.resize(Nmsg);
  recv_msgmem.resize(Nmsg);
  xmit_msghandle.resize(Nmsg);
  recv_msghandle.resize(Nmsg);

  for(int m=0;m<Nmsg;m++){
    xmit_msgmem[m] = QMP_declare_msgmem((void *) &sendbufs[m][0],sendbuf_bytes[m]) ;  // send buf
    recv_msgmem[m] = QMP_declare_msgmem((void *) &recvbufs[m][0],recvbuf_bytes[m]) ; // receive buf
    
    xmit_msghandle[m] = QMP_declare_send_to(xmit_msgmem[m],tproc[m],0);
    recv_msghandle[m] = QMP_declare_receive_from(recv_msgmem[m],fproc[m],0);
  }
  QDPIO::cout << "Initialised haloes "<<endl;
}
std::complex<double> * BfmMultiGrid::HaloGetStencilData(int _idx,int _ballidx,
							std::vector<std::complex<double> > &my_data)
{
  // Locate the neighbour
  int me, from, to;
  int fromblock;
  int mu;

  std::complex<double> *nbr_data;
  mu =_ballidx;

  int bid = rbuf_id [_idx*Nball+mu] ;
  int bidx= rbuf_idx [_idx*Nball+mu] ;

  if (bid==-1){ 
    nbr_data = &my_data[bidx];
  } else { 
    nbr_data = &recvbufs[bid][bidx];
  }
  return nbr_data;
}
void BfmMultiGrid::HaloExchange(std::vector<std::complex<double> >&my_data)
{
  int Nmsg = tproc.size();

  int nwork = sbuf_gather.size();
  int me, thrlen, throff;
  static int pooh;
  static unsigned long totbytes=0;

  linop->thread_work(nwork,me,thrlen,throff);
  if ( linop->isBoss() && (!me)  && (!pooh)) {
    totbytes=0;
    printf("Dividing %d gather tasks between %d threads\n",nwork,linop->threads);
    for(int t=0;t<nwork;t++){
      int sz =sbuf_gather[t].size();
      int bytes = sz*Nvec*2*sizeof(double);
      totbytes+=bytes;
      printf("Unit for %d is size() %d, %d bytes\n",t,bytes,totbytes);
    }
    pooh=1;
  }

  // Better to run this loop and make a single flat table, single flat sendbuf
  // Could then load balance across threads => ~4x saving in gather time
  uint64_t t0 = GetTimeBase();
  linop->thread_barrier();
  for(int tidx=throff;tidx<throff+thrlen;tidx++){
#if 0
    for(int j=0;j<sbuf_gather[tidx].size();j++){
      for(int i=0;i<Nvec;i++){
	sendbufs[tidx][j*Nvec+i] = my_data[i+sbuf_gather[tidx][j]*Nvec];
      }
    }
#else
    qpx_gather(Nvec,
	       &sbuf_gather[tidx][0],
	        sbuf_gather[tidx].size(),
	       (double *)&sendbufs[tidx][0],
	       (double *)&my_data[0]);
#endif
  }
  linop->thread_barrier();
  uint64_t t1 = GetTimeBase();
  HaloExchangeCommStart(this);
  uint64_t tw = GetTimeBase();
  HaloExchangeCommComplete(this);
  uint64_t t2 = GetTimeBase();
  if ( MGreport &&(!me) ) {
    printf("Halo: Gather %ld cyc\n",t1-t0);
    printf("Halo: Start  %ld cyc\n",tw-t1);
    printf("Halo: Wait   %ld cyc\n",t2-tw);
    printf("Halo: Gather %f Mbyte/s\n",1600.*totbytes/(t1-t0) );
    printf("Halo: Start  %f Mbyte/s\n",1600.*totbytes/(tw-t1) );
    printf("Halo: Wait   %f Mbyte/s\n",1600.*totbytes/(t2-tw) );
  }
  linop->thread_barrier();
}
#ifndef BFM_BGQ_SPI_DSLASH
  ////////////////////////////////////////////////
  // QMP implementation
  ////////////////////////////////////////////////
void BfmMultiGrid::HaloExchangeCommComplete(BfmMultiGrid * BMG)
{
  int Nmsg=BMG->tproc.size();
  int me = BMG->linop->thread_barrier();
  if ( me == 0){
    for(int m=0;m<Nmsg;m++) QMP_wait (BMG->xmit_msghandle[m]);
    for(int m=0;m<Nmsg;m++) QMP_wait (BMG->recv_msghandle[m]);
  }
  linop->thread_barrier();
}
void HaloExchangeCommStart(BfmMultiGrid * BMG)
{
  int Nmsg=BMG->tproc.size();
  int me = BMG->linop->thread_barrier();
  if ( me == 0){
    for(int m=0;m<Nmsg;m++) QMP_start(BMG->xmit_msghandle[m]);
    for(int m=0;m<Nmsg;m++) QMP_start(BMG->recv_msghandle[m]);
  }
  BMG->linop->thread_barrier();
}
void BfmMultiGrid::HaloInitBuffers(BfmMultiGrid * BMG)
{
  int Nmsg = tproc.size();

  for(int m=0;m<Nmsg;m++){
    BMG->sendbufs[m] = (std::complex<double> *)bfm_alloc(BMG->sendbuf_bytes[m]);
    BMG->recvbufs[m] = (std::complex<double> *)bfm_alloc(BMG->recvbuf_bytes[m]);
  }

  // QMP init
  BMG->xmit_msgmem.resize(Nmsg);
  BMG->recv_msgmem.resize(Nmsg);
  BMG->xmit_msghandle.resize(Nmsg);
  BMG->recv_msghandle.resize(Nmsg);

  for(int m=0;m<Nmsg;m++){
    BMG->xmit_msgmem[m] = QMP_declare_msgmem((void *) &BMG->sendbufs[m][0],2*sizeof(double)*BMG->sendbufs[m].size()) ;  // send buf
    BMG->recv_msgmem[m] = QMP_declare_msgmem((void *) &BMG->recvbufs[m][0],2*sizeof(double)*BMG->recvbufs[m].size()) ; // receive buf
    
    BMG->xmit_msghandle[m] = QMP_declare_send_to(BMG->xmit_msgmem[m],BMG->tproc[m],0);
    BMG->recv_msghandle[m] = QMP_declare_receive_from(BMG->recv_msgmem[m],BMG->fproc[m],0);
  }
  QDPIO::cout << "Initialised haloes "<<endl;
}
#endif

void BfmMultiGrid::HaloTest(std::vector<std::complex<double> >&my_data)
{
  QDPIO::cout << "Halo Init"<<endl;
  HaloExchange(my_data);
  QDPIO::cout << "Halo Stencil cross check"<<endl;

  std::complex<double> * halo_data;
  std::vector<std::complex<double> > ref_data(Nvec);
  for(int idx=0;idx<LocalNblock;idx++){
    for(int mu=0;mu<Nball;mu++){
      int rev=0;
      GetStencilData(idx,mu,ref_data,my_data,rev);
      halo_data= HaloGetStencilData(idx,mu,my_data);
      for(int i=0;i<Nvec;i++){
	if ( halo_data[i]!=ref_data[i] ) {
	  QDPIO::cout << "Stencil Mismatch " << idx<< " mu " << mu << " g "<< real(ref_data[i]) << " b "<< real(halo_data[i]) <<endl;
	}
      }
    }
  }
  QDPIO::cout << "Check complete\n" <<endl;
}

void BfmMultiGrid::GetStencilData(int _idx,int _ballidx,
				  std::vector<std::complex<double> >& nbr_data,
				  std::vector<std::complex<double> >&my_data,int rev)
{
  // Locate the neighbour
  int me, from, to;
  int fromblock;

  NeighbourBlockId(_idx,_ballidx,me,from, fromblock,to,rev);
  
  if ( me == from ) {
    for(int i=0;i<Nvec;i++){
      nbr_data[i] = my_data[i+fromblock*Nvec];
    }
  } else { 

    QMP_msgmem_t msgmem[2];
    msgmem[0] = QMP_declare_msgmem((void *) &my_data[fromblock*Nvec],2*Nvec*sizeof(double));  // send buf
    msgmem[1] = QMP_declare_msgmem((void *) &nbr_data[0],2*Nvec*sizeof(double)); // receive buf
    
    QMP_msghandle_t msghandle[2];
    msghandle[0] = QMP_declare_send_to(msgmem[0],to,0);
    msghandle[1] = QMP_declare_receive_from(msgmem[1],from,0);
    
    for(int m=0;m<2;m++)  QMP_start(msghandle[m]);
    for(int m=0;m<2;m++)  QMP_wait (msghandle[m]);

    for(int m=0;m<2;m++)  QMP_free_msghandle(msghandle[m]);
    for(int m=0;m<2;m++)  QMP_free_msgmem(msgmem[m]);

  }

}

void BfmMultiGrid::OrthogonaliseSubspace(void)
{
  // Remove earlier (normalised) vectors
  std::vector<std::complex<double> > bip(LocalNblock);
  std::vector<double> bn(LocalNblock);

#pragma omp parallel 
  {
#pragma omp for 
    for(int t=0;t<linop->nthread;t++) {

      for(int v=0;v<Nvec;v++){
	for(int u=0;u<v;u++){
	  //Inner product & remove component
	  linop->block_inner(subspace[u],subspace[v],LocalNblock,(double *)&bip[0],(double *)&inner_reduce[0],&local_block_id[Odd][0]);
	  linop->block_zaxpy(subspace[v],subspace[u],subspace[v],-1,(double *)&bip[0],&local_block_id[Odd][0]);
	}
	// Normalise this vector
	linop->block_norm(subspace[v],LocalNblock,(double *)&bn[0],&local_block_id[Odd][0]);
	linop->block_normalise(subspace[v],LocalNblock,(double *)&bn[0],&local_block_id[Odd][0]);
      }
    }
  }

#if 1
  for(int v=0;v<Nvec;v++){
  for(int u=0;u<Nvec;u++){
    
#pragma omp parallel 
  {
#pragma omp for 
    for(int t=0;t<linop->nthread;t++) {
      linop->block_inner(subspace[u],subspace[v],LocalNblock,(double *)&bip[0],(double *)&inner_reduce[0],&local_block_id[Odd][0]);
    }
  }
      double c = (u==v);

      for(int i=0;i<bip.size();i++) { 
	if ( abs(bip[i]-c) > 1.0e-6 ) {  
	  cout << "Block "<< i<<" vec "<< u<< "," << v << " inner= " << bip[i] << " : " 
	       <<"Doesnt look orthogonal and should be "<<c<<endl;
	}
      }
  }}
#endif
}
void BfmMultiGrid::PickBlock (Fermion_t was,Fermion_t is,int b)
{
#pragma omp parallel 
    {
#pragma omp for 
      for(int i=0;i<linop->threads;i++) {
	linop->fill(is,0.0);
	linop->block_pick (is,was,is,b,&block_id[1][0]);
      }
    }
    exit(0);
}
void BfmMultiGrid::PickQuadrant(Fermion_t was, Fermion_t is,int q)
{
  
  double nn;
#pragma omp parallel 
    {
#pragma omp for 
      for(int i=0;i<linop->threads;i++) {
	linop->fill(is,0.0);
	linop->block_pick (is,was,is,q,&quadrant_id[1][0]);
	nn=linop->norm(is);
      }
    }
    if(linop->isBoss()) { 
      printf("Quadrant %d picked norm is %le\n",q,nn);
    }
    exit(0);
}

void BfmMultiGrid::ComputeLittleMatrixColored (void)
{
  Fermion_t zero_t= linop->threadedAllocFermion();
  Fermion_t phi_t = linop->threadedAllocFermion();
  Fermion_t tmp_t = linop->threadedAllocFermion();
  Fermion_t mmp   = linop->threadedAllocFermion();
  Fermion_t mp    = linop->threadedAllocFermion();

  int me, thrlen,throff;
  me = linop->thread_barrier();
  linop->thread_barrier();
  linop->fill(zero_t,0.0);

  for(int i=0;i<Nvec;i++){

    // Problematic for threading
    if ( !me ) { 
      FlightLog("ComputeColored: Vector "<<i);
    }
    for(int b=0;b<Nball;b++){  // Loop over momenta (Nball)

      if ( (i==0)&&(b==1)&&!me ) { 
	FlightLog("BallId "<<b);
      }
      /////////////////////////////////////////////////////
      // Stick a different phase on every block
      /////////////////////////////////////////////////////
      int lb[5];// Local Block
      int gb[5];// Global Block
      int    imom[5];
      double dmom[5];

      GetDelta(imom,b,0);  
      for(int mu=0;mu<5;mu++){
 	dmom[mu] = imom[mu]*2*M_PI/global_nb[mu];
      }
      
      linop->thread_barrier();
      for(lb[0]=0;lb[0]<local_nb[0];lb[0]++){ // Redundant computation alert
      for(lb[1]=0;lb[1]<local_nb[1];lb[1]++){
      for(lb[2]=0;lb[2]<local_nb[2];lb[2]++){
      for(lb[3]=0;lb[3]<local_nb[3];lb[3]++){
      for(lb[4]=0;lb[4]<local_nb[4];lb[4]++){

	double pdotx = 0.0;
	for(int mu=0;mu<5;mu++) gb[mu] = ncoor[mu]*local_nb[mu]+lb[mu];
	for(int mu=0;mu<5;mu++) pdotx+= gb[mu]*dmom[mu];

	std::complex<double> pha(cos(pdotx),sin(pdotx));

	int lbidx=GetLocalBlockIdx(lb);
	phases[lbidx] = pha;

      }}}}}
      linop->thread_barrier();

      if ( (i==0)&&(b==1)&&!me ) { 
	FlightLog("Calculated phases");
      }

      ///////////////////////////////////////////////////////
      // We apply the matrix Nball times using these phases
      ///////////////////////////////////////////////////////
      linop->block_zaxpy(phi_t,subspace[i],zero_t,1.0,(double *)&phases[0],&local_block_id[Odd][0]);

      if ( (i==0)&&(b==1)&&!me ) { 
	FlightLog("zaxpy done");
      }

      linop->Mprec(phi_t,mp,tmp_t,0);
      linop->Mprec(mp,mmp,tmp_t,1); 
      double mmp_norm= linop->norm(mmp);
      double phi_norm= linop->norm(phi_t);
      double ss_norm = linop->norm(subspace[i]);
      if ( (i==0)&&(b==1)&&!me ) { 
	FlightLog("matmul done "<< mmp_norm<< " phi "<<phi_norm <<" subspace "<<ss_norm);
      }

      ///////////////////////////////////////////////////////
      // Project this into subspace
      ///////////////////////////////////////////////////////
      ProjectToSubspace(mmp,ComputeProj[b]);

      if ( (i==0) && (b==1) && !me ) { 
	FlightLog("Projection done");
      }

      linop->thread_work(LocalNblock,me,thrlen,throff);
      for(int lbidx=throff;lbidx<throff+thrlen;lbidx++) { 
	for(int j=0;j<Nvec;j++){
	  ComputeProj[b][lbidx*Nvec+j] = ComputeProj[b][lbidx*Nvec+j]*conj(phases[lbidx]);
	}
      }
      linop->thread_barrier();
    }

    if ( i==0 && !me ) { 
      FlightLog("APPLIED ALL PROJECTIONS");
    }


    //////////////////////////////////////////////////////////
    // Solve system of equations to get Aij
    //////////////////////////////////////////////////////////
    /*
     *     Here, k,l index which possible shift within the 3^Nd "ball" connected by MdagM.
     *
     *     conj(phases[block]) proj[k][ block*Nvec+j ] =  \sum_ball  e^{i q_k . delta} < phi_{block,j} | MdagM | phi_{(block+delta),i} > 
     *                                                 =  \sum_ball e^{iqk.delta} A_ji
     *
     *     Must invert matrix M_k,l = e^[i q_k . delta_l]
     *
     *     Where q_k = delta_k . (2*M_PI/global_nb[mu])
     */
    {
      
      linop->thread_work(LocalNblock,me,thrlen,throff);
      for(int lbidx=throff;lbidx<throff+thrlen;lbidx++){
	for(int j=0;j<Nvec;j++){

	for(int mu=0;mu<5;mu++){

	  double pmu=2*M_PI/global_nb[mu];      
	  std::complex<double> pha(cos(pmu),-sin(pmu));
	  std::vector<std::complex<double> > FT(Nball);
	  std::complex<double> a = pha+1.0;
	  std::complex<double> b = (pha-1.0)*(pha-1.0);
	  std::complex<double> ab=a*b;

	  ///////////////////////////////////////////////////////////////////////////////
	  // Stencil == 3
	  //  pha* 1 pha                 pha/a -pha        pha^2/a
	  //  1    1  1     -> inv ->    -pha  (1+pha^2)   -pha           / [ 1-pha]^2     where a=1+pha
 	  //  pha  1 pha*               pha^2/a -pha       pha/a
	  ///////////////////////////////////////////////////////////////////////////////
	  ///////////////////////////////////////////////////////////////////////////////
	  // Stencil == 2
	  // -1 1                  1/2   -1 1 
	  //  1 1      -> inv ->          1 1
 	  //                  
	  ///////////////////////////////////////////////////////////////////////////////

	  std::vector< std::vector< std::complex<double> > > imat; // Inverse matrix for this mu.
	  imat.resize(stencil_size[mu]);
	  for(int i=0;i<stencil_size[mu];i++) imat[i].resize(stencil_size[mu]);

	  if ( stencil_size[mu]==3 ) {
	    imat[0][0] = pha/ab; 	   imat[0][1] =-pha/b;	                  imat[0][2] = pha*pha/ab;
	    imat[1][0] =-pha/b;	           imat[1][1] = (1.0+pha*pha)/b;	  imat[1][2] =-pha/b;
	    imat[2][0] = pha*pha/ab;       imat[2][1] =-pha/b;         	          imat[2][2] = pha/ab;
	  } else if (stencil_size[mu]==2) { 
	    imat[0][0] = -0.5; 	   imat[0][1] = 0.5;
	    imat[1][0] =  0.5;	   imat[1][1] = 0.5;
	  } else if (stencil_size[mu]==1) { 
	    imat[0][0]=1.0;
	  }

	  int d[5];
	  int nd[5];

	  for(d[0]=0;d[0]<stencil_size[0];d[0]++){
	  for(d[1]=0;d[1]<stencil_size[1];d[1]++){
	  for(d[2]=0;d[2]<stencil_size[2];d[2]++){
	  for(d[3]=0;d[3]<stencil_size[3];d[3]++){
	  for(d[4]=0;d[4]<stencil_size[4];d[4]++){

	    int  d_mu;
	    d_mu =   d[4]+stencil_size[4]*(d[3]+stencil_size[3]*(d[2]+stencil_size[2]*(d[1]+stencil_size[1]*d[0])));
	    for(int i=0;i<5;i++) nd[i]=d[i];
	    FT[d_mu] = 0;
	    
	    for(nd[mu]=0;nd[mu]<stencil_size[mu];nd[mu]++){
	      int nd_mu;
	      nd_mu = nd[4]+stencil_size[4]*(nd[3]+stencil_size[3]*(nd[2]+stencil_size[2]*(nd[1]+stencil_size[1]*nd[0])));
	      
	      FT[d_mu]+= imat[d[mu]][nd[mu]] * ComputeProj[nd_mu][lbidx*Nvec+j];
	      
	    }
	  }}}}}

	  for(int b=0;b<Nball;b++) ComputeProj[b][lbidx*Nvec+j] = FT[b];

	}//5-Directions FT
	////////////////////////
	// Copy back solution of system of eqns
	////////////////////////
	for(int b=0;b<Nball;b++)  Asparse[lbidx][b][Nvec*j+i] = ComputeProj[b][lbidx*Nvec+j];
	}//j
      }// Block
    }// Scope
    if ( i==0 && !me ) { 
      FlightLog("APPLIED INVERSE FT");
    }

  }//Nvec

  // Find non-zero elements of stencil
  if( !me ) {
    for(int b=0;b<Nball;b++){
      double nrm =0.0;
      for(int j=0;j<Nvec*Nvec;j++){
	std::complex<double> e = Asparse[0][b][j];
	nrm=nrm+real(e)*real(e)+imag(e)*imag(e);
      }
      if(nrm < 1.0e-10 ) {
	if ( linop->isBoss() ) { 
	  printf("Stencil %d is zero %le\n",b,nrm);
	}
	StencilNonZero[b] = 0;
      } else { 
	if ( linop->isBoss() ) { 
	  printf("Stencil %d is nonzero %le\n",b,nrm);
	}
	StencilNonZero[b] = 1;
      }
    }
  }

  linop->threadedFreeFermion( zero_t ); 
  linop->threadedFreeFermion( phi_t ); 
  linop->threadedFreeFermion( tmp_t ); 
  linop->threadedFreeFermion( mmp   ); 
  linop->threadedFreeFermion( mp    );

}
void BfmMultiGrid::ComputeLittleMatrix (void)
{
  exit(0);
  Fermion_t phi_t = linop->allocFermion();
  Fermion_t tmp_t = linop->allocFermion();
  Fermion_t mmp   = linop->allocFermion();
  Fermion_t mp    = linop->allocFermion();
	
  for(int i=0;i<Nvec;i++){
    for(int b=0;b<Nblock;b++){            // Loop over global blocks

      // Apply DdagD to phi_i[b] -> mmp
#pragma omp parallel 
	{
#pragma omp for 
	  for(int t=0;t<linop->nthread;t++) {
	    PickBlock(subspace[i],phi_t,b);
	    linop->Mprec(phi_t,mp,tmp_t,0);
	    linop->Mprec(mp,mmp,tmp_t,1); 
	  }
	}
	
	// Need to get the local block coord for phi_i[b]
	int gb[5];
	int idx=b;
	int lb[5];
	int lb_idx;
	int proc;
	int me = QDP::Layout::nodeNumber();  

	gb[0] = idx%global_nb[0]; idx = idx/global_nb[0];
	gb[1] = idx%global_nb[1]; idx = idx/global_nb[1];
	gb[2] = idx%global_nb[2]; idx = idx/global_nb[2];
	gb[3] = idx%global_nb[3]; idx = idx/global_nb[3];
	gb[4] = idx%global_nb[4];

	GlobalToLocalBlock(gb,proc,lb);
	lb_idx= GetLocalBlockIdx(lb); // Block index within node "proc" of "b"

	std::vector<std::complex<double> > proj(LocalNsubspace);
	std::vector<std::complex<double> > row (Nvec);


	// A[j,i] = <phi_j| D |phi_i> 
	ProjectToSubspace(mmp,proj); // phi_j^dag M phi_i
	for(int n=0;n<Nball;n++){
	  int me;
	  int from;
	  int to;
	  int xblock;
	  int reverse=1;
	  NeighbourBlockId(lb_idx,n,me,from,xblock,to,reverse);

	  // I am the neighbour of the site containing the non-zero block
	  if( proc == to ) { 
	    for(int j=0;j<Nvec;j++){
	      Asparse[xblock][n][Nvec*j+i] = proj[xblock*Nvec+j];
	    }	    
	  }

	}
    }
  }
  linop->freeFermion( phi_t ); 
  linop->freeFermion( tmp_t ); 
  linop->freeFermion( mmp   ); 
  linop->freeFermion( mp    );
}
// Threaded
void BfmMultiGrid::ProjectToSubspace(Fermion_t vec, std::vector<std::complex<double> > &proj)
{
  int me = linop->thread_barrier();
  if ( proj.size() != LocalNsubspace ){
    printf("ProjectToSubspace size wrong\n");
    exit(0);
  }
  linop->thread_barrier();

  for(int v_j=0;v_j<Nvec;v_j++){

     linop->thread_barrier();
     linop->block_inner(subspace[v_j],vec,
			LocalNblock,
			(double *)&ProjectInnerProduct[0],
			(double *)&inner_reduce[0],
			&local_block_id[Odd][0],
			s_min[v_j],s_max[v_j]);
     linop->thread_barrier();

     for(int b_j=0;b_j<LocalNblock;b_j++){
       int j = SubspaceIndex(b_j,v_j);
       proj[j] = ProjectInnerProduct[b_j];
     }
     linop->thread_barrier();
  }
}
// Threaded
void BfmMultiGrid::PromoteFromSubspace(std::vector<std::complex<double> > &v,Fermion_t prom)
{
  linop->thread_barrier();
  uint64_t tf1=GetTimeBase();
  linop->fill(prom,0.0); // Wasteful
  uint64_t tf2=GetTimeBase();
  uint64_t tc=0;
  uint64_t tz=0;

  int thrlen,throff,me;
  for(int v_i=0;v_i<Nvec;v_i++){
    uint64_t t1=GetTimeBase();
    linop->thread_work(LocalNblock,me,thrlen,throff);
    for(int b_i=throff;b_i<throff+thrlen;b_i++){
      PromoteCoeff[b_i]=v[v_i+b_i*Nvec];
    }
    linop->thread_barrier();
    uint64_t t2=GetTimeBase();
    linop->block_zaxpy(prom,subspace[v_i],prom,1.0,(double *)&PromoteCoeff[0],&local_block_id[Odd][0],s_min[v_i],s_max[v_i]);
    uint64_t t3=GetTimeBase();
    tc+=t2-t1;
    tz+=t3-t2;
  }
  if(linop->isBoss()&&!me) {
    static int print_once;
    if ( print_once == 0 ) { 
        printf("Promote : Fill/coeff/zaxpy \n\t%.16ld\n\t%.16ld\n\t%.16ld\n",tf2-tf1,tc,tz);
	print_once=1;
    }
  }
}
// Threaded
void myzgemv(int N,double * __restrict A, double * __restrict in,double c, double * __restrict out);

void BfmMultiGrid::Apply(std::vector<std::complex<double> > &in,std::vector<std::complex<double> > &out )
{
  // This is tricky to thread because comms is interleaved.
  // Easiest way to thread is to serially
  // fill out a set of vectors for the neighbours, and then
  // apply the matrix in parallel.
  std::vector<std::complex<double> > nbr_data(Nvec);
  std::vector<std::complex<double> > all_nbr_data(LocalNsubspace);
  for(int mu=0;mu<Nball;mu++){
    for(int b=0;b<LocalNblock;b++){
      int nodag=0;
      GetStencilData(b,mu,nbr_data,in,nodag);
      double coeff=1.0;
      if ( mu==0) {
	coeff=0.0;
      }
      
      myzgemv(Nvec,(double *)&Asparse[b][mu][0],(double *)&nbr_data[0],coeff,(double *)&out[b*Nvec]);

    }
  }
}


void BfmMultiGrid::ApplyInverse(std::vector<std::complex<double> > &v,std::vector<std::complex<double> > &vinv)
{
  if ( LittleDopSolver == LittleDopSolverDeflCG ) {
    ApplyInverseDeflCG(v,vinv);
    return;
  }
  if ( LittleDopSolver == LittleDopSolverMCR ) {
    ApplyInverseMCR(v,vinv);
    return;
  }  
  if ( LittleDopSolver == LittleDopSolverCG ) {
    ApplyInverseCG(v,vinv);
    return;
  }
  exit(0);
}
void BfmMultiGrid::ApplyInverseDeflCG(std::vector<std::complex<double> > &src,
				      std::vector<std::complex<double> > &x)
{
  double rtzp,rtz,a,d,b;

  int me = linop->thread_barrier();

  double ssq =  norm_vec(src); // Forces sync of threads & nodes

  // Deflated precond CG (3.6 in Saad umsi-98-97)
  ///////////////////////////////////
  // x_{-1} = src
  ///////////////////////////////////

  axpy(x,src,src,0.0);
  ApplyThread(x,Krylov_Ap);
  axpy(Krylov_r, Krylov_Ap, src,-1.0);        // r_{-1} = src - A x

  ///////////////////////////////////
  // Choose x_0 such that 
  // x_0 = guess +  (A_ss^inv) r_s
  // 
  // W^T (src - A x_0) = src_s - A guess_s - r_s  
  //                   = src_s - (A guess)_s - src_s  + (A guess)_s 
  //                   = 0 
  ///////////////////////////////////

  LdopDeflationProject(Krylov_r,DeflKrylovProj);
  LdopDeflationMatrixInverseMult(DeflKrylovProj,DeflKrylovMss);
  LdopDeflationPromote(DeflKrylovMss,Krylov_Ap);
  axpy(x,x,Krylov_Ap,1.0);


  // Recomputes r=src-x0
  ApplyThread(x,Krylov_Ap);
  axpy (Krylov_r,Krylov_Ap, src,-1.0);  

  rtzp =innerProductReal(Krylov_r,Krylov_r);

  // Check orthogonal
  LdopDeflationProject(Krylov_r,DeflKrylovProj);     
  double nv=norm_vec(DeflKrylovProj);
  
  ///////////////////////////////////////
  // Solve for Mss mu = P A r and set p = r-mu
  ///////////////////////////////////////
  ApplyThread(Krylov_r,Krylov_Ap);
  LdopDeflationProject(Krylov_Ap,DeflKrylovProj);
  LdopDeflationMatrixInverseMult(DeflKrylovProj,DeflKrylovMss);
  LdopDeflationPromote(DeflKrylovMss,Krylov_mu);
  LdopDeflationMatPromote(DeflKrylovMss,Krylov_Amu);

  axpy (Krylov_p,Krylov_mu,Krylov_r,-1.0);
  axpy (Krylov_Ap,Krylov_Amu,Krylov_Ap,-1.0);

  double rsq =  LittleDopSolverResidual*LittleDopSolverResidual*ssq;

  uint64_t t0= GetTimeBase();

  for (int k=1;k<=10000;k++){

    uint64_t mat_0= GetTimeBase();

    rtz=rtzp;
    d  = innerProductReal(Krylov_p,Krylov_Ap);
    a = rtz/d;
    axpy(x,Krylov_p,x,a);
    axpy(Krylov_r,Krylov_Ap,Krylov_r,-a);

    rtzp =norm_vec(Krylov_r);
    b = rtzp/rtz;

    // WAW mu = W A z
    uint64_t mat_2= GetTimeBase();
    ApplyThread(Krylov_r,Krylov_Ar);
    uint64_t mat_3= GetTimeBase();

    uint64_t defl_0= GetTimeBase();
    LdopDeflationProject(Krylov_Ar,DeflKrylovProj);
    uint64_t defl_1= GetTimeBase();
    LdopDeflationMatrixInverseMult(DeflKrylovProj,DeflKrylovMss);
    uint64_t defl_2= GetTimeBase();
    LdopDeflationPromote(DeflKrylovMss,Krylov_mu);
    LdopDeflationMatPromote(DeflKrylovMss,Krylov_Amu);
    uint64_t defl_3= GetTimeBase();

    //=>    Krylov_p = b*Krylov_p + Krylov_r -  Krylov_mu
    axpy(Krylov_p,Krylov_p,Krylov_r,b);     // Overlap with comms if lift
    axpy(Krylov_p,Krylov_mu,Krylov_p,-1.0);

    //=>    Krylov_Ap = b*Krylov_p + Krylov_Ar +  Krylov_Amu
    axpy(Krylov_Ap,Krylov_Ap,Krylov_Ar,b);
    axpy(Krylov_Ap,Krylov_Amu,Krylov_Ap,-1.0);

    

    if (linop->isBoss() && !me ){
      if ( k%250 == 0 ){
	printf("bfmbase::ApplyInverseDelfCG: k= %d residual = %le \n",k,sqrt(rtzp/ssq)); fflush(stdout);
      }
    }

    // Stopping condition
    if ( rtzp <= rsq ) { 
      uint64_t t1= GetTimeBase();
      ApplyThread(x,Krylov_Ap);
      axpy(Krylov_r,src,Krylov_Ap,-1.0);
      double tmpnorm=sqrt(norm_vec(Krylov_r));
      double srcnorm=sqrt(ssq);
      double true_residual = tmpnorm/srcnorm;

      if ( linop->isBoss() && !me ) {
	printf("BfmMultiGrid::ApplyInverseDeflCG: true residual is %le after %d iterations\n",true_residual,k);
	printf("BfmMultiGrid::ApplyInverseDeflCG: %f millisec\n",(t1-t0)/1600.0/1000.0);
	printf("*** Last Iter Stats: *** \n");
	printf("Cycles total %ld\n",t1-mat_0);
	printf("Cycles mat %ld\n",mat_3-mat_2);
	printf("Cycles defl1 %ld\n",defl_1-defl_0);
	printf("Cycles defl2 %ld\n",defl_2-defl_1);
	printf("Cycles defl3 %ld\n",defl_3-defl_2);
	fflush(stdout);
      }
      return;
    }

  }
  if ( linop->isBoss() && !me ) printf("BfmMultiGrid::ApplyInverseDeflCG: CG not converged \n"); fflush(stdout);
  exit(0);
}

void BfmMultiGrid::ApplyInverseMCR(std::vector<std::complex<double> > &v,std::vector<std::complex<double> > &vinv)
{
  int me, thrlen,throff;
  double a;
  double b;
  double c;
  double d;
  double cp;
  double ssq;
  double rsq;

  double rAr;
  double rArp;
  double pAAp;

  // Conjugate gradient on A. A is hermitian posdef.
  // Need these to be shared across threads
  me = linop->thread_barrier();

  linop->thread_work(LocalNsubspace,me,thrlen,throff);
  for(int i=throff;i<throff+thrlen;i++){
    vinv[i]=0;
    Krylov_r[i]=v[i];
    Krylov_p[i]=Krylov_r[i];
  }
  linop->thread_barrier();

  ApplyThread(Krylov_p,Krylov_Ap); //Threaded
  ApplyThread(Krylov_r,Krylov_Ar); //Threaded
  rAr = innerProductReal(Krylov_r,Krylov_Ar);
  pAAp= norm_vec(Krylov_Ap);

  cp =norm_vec(Krylov_r);
  ssq=norm_vec(v);
  rsq = LittleDopSolverResidual*LittleDopSolverResidual*ssq;
  if ( me == 0 ) { 
    QDPIO::cout << "target rsq "<<rsq<< " " << LittleDopSolverResidual << " ssq " << ssq<<endl;
  }

  uint64_t t1 = GetTimeBase();
  for(int k=1;k<10000;k++){


    MGreport=0;
    if(linop->isBoss() && (k==10) ) MGreport=1;
    a   = rAr/pAAp;

    uint64_t axpy_0 = GetTimeBase();
    axpy(vinv,Krylov_p,vinv,a);           // overlap with comms
    axpy(Krylov_r,Krylov_Ap,Krylov_r,-a);
    uint64_t axpy_1 = GetTimeBase();


    uint64_t mat_0 = GetTimeBase();
    ApplyThread(Krylov_r,Krylov_Ar); //          Ar_i+1
    uint64_t mat_1 = GetTimeBase();

    rArp=rAr;
    uint64_t inner_0 = GetTimeBase();
    rAr = innerProductReal(Krylov_r,Krylov_Ar); // Fuse into ApplyThread
    uint64_t inner_1 = GetTimeBase();
    b = rAr/rArp;
 
    uint64_t axpy1_0 = GetTimeBase();
    axpy(Krylov_p,Krylov_p,Krylov_r,b);
    axpy(Krylov_Ap,Krylov_Ap,Krylov_Ar,b); // Axpy_norm
    uint64_t axpy1_1 = GetTimeBase();

    uint64_t norm_0 = GetTimeBase();
    pAAp= norm_vec(Krylov_Ap);// Fuse into axpy norm
    cp  = norm_vec(Krylov_r); // Fuse into earlier ApplyThread
    uint64_t norm_1 = GetTimeBase();

    if ( (me==0) && k == 10 ) { 
      FlightLog("Iteration 10 - linear alg done");
    }

    if(cp<rsq) {
      ApplyThread(vinv,Krylov_Ap);
      linop->thread_work(LocalNsubspace,me,thrlen,throff);
      for(int i=throff;i<throff+thrlen;i++){
	Krylov_Ap[i]=Krylov_Ap[i]-v[i];
      }
      linop->thread_barrier();

      double nv = norm_vec(Krylov_Ap);
      uint64_t t2 = GetTimeBase();
      if ( me == 0 ) { 

	QDPIO::cout<< "LittleDiracInversion via MCR converged : "<<k<<" iterations ,"
		   << " residual "<<sqrt(nv/ssq)
		   << " millisec "<< (t2-t1)/1600.0/1000.0
		   <<endl;

	QDPIO::cout << "***  Last Iter Stats: *** "<< axpy_1-axpy_0 <<endl;
	QDPIO::cout << "Cycles axpy  "<< axpy_1-axpy_0 <<endl;
	QDPIO::cout << "Cycles mat   "<< mat_1-mat_0 <<endl;
	QDPIO::cout << "Cycles inner "<< inner_1-inner_0 <<endl;
	QDPIO::cout << "Cycles axpy1 "<< axpy1_1-axpy1_0 <<endl;
	QDPIO::cout << "Cycles norm  "<< norm_1-norm_0 <<endl;
      }

      return;
    }
  }
  printf("Little Dirac Matrix inversion failed\n");
  fflush(stdout);
  exit(0);
}

// Threaded
void BfmMultiGrid::ApplyInverseCG(std::vector<std::complex<double> > &v,std::vector<std::complex<double> > &vinv)
{
  int me, thrlen,throff;

  // Conjugate gradient on A. A is hermitian posdef.
  // Need these to be shared across threads
  me = linop->thread_barrier();

  linop->thread_work(LocalNsubspace,me,thrlen,throff);

  for(int i=throff;i<throff+thrlen;i++){
    vinv[i]=0;
    Krylov_r[i]=v[i];
    Krylov_p[i]=Krylov_r[i];
  }
  linop->thread_barrier();

  double a;
  double b;
  double c;
  double d;
  double cp;
  double ssq;
  double rsq;


  ssq=norm_vec(v);
  a  =norm_vec(Krylov_p);
  cp =norm_vec(Krylov_r);
  rsq = LittleDopSolverResidual*LittleDopSolverResidual*ssq;
  if ( me == 0 ) { 
    QDPIO::cout << "target rsq "<<rsq<< " " << LittleDopSolverResidual<<endl;
  }

  uint64_t t0 = GetTimeBase();
  for(int k=1;k<10000;k++){
    c=cp;

    MGreport=0;
    if(linop->isBoss() && (k==10) ) MGreport=1;
    
    uint64_t mat_0 = GetTimeBase();
    ApplyThread(Krylov_p,Krylov_Ap); //Threaded
    uint64_t mat_1 = GetTimeBase();

    uint64_t inner_0 = GetTimeBase();
    d=innerProductReal(Krylov_p,Krylov_Ap);
    a=c/d;
    uint64_t inner_1 = GetTimeBase();

    uint64_t axpy_norm_0 = GetTimeBase();
    axpy(Krylov_r,Krylov_Ap,Krylov_r,-a); // Fuse Axpy-Norm
    cp=norm_vec(Krylov_r);
    uint64_t axpy_norm_1 = GetTimeBase();
    b=cp/c;

    uint64_t axpy_0 = GetTimeBase();
    axpy(vinv,Krylov_p,vinv,a);
    axpy(Krylov_p,Krylov_p,Krylov_r,b);
    uint64_t axpy_1 = GetTimeBase();

    if(cp<rsq) {
      uint64_t t1 = GetTimeBase();
      ApplyThread(vinv,Krylov_Ap);
      linop->thread_work(LocalNsubspace,me,thrlen,throff);
      for(int i=throff;i<throff+thrlen;i++){
	Krylov_Ap[i]=Krylov_Ap[i]-v[i];
      }
      linop->thread_barrier();

      double nv = norm_vec(Krylov_Ap);
      if ( me == 0 ) { 
	FlightLog("CG LittleDiracInversion converged : "<<k<<" iterations , residual "<<sqrt(nv/ssq)
		  << " millisec "<< (t1-t0)/1600.0/1000.0)
	QDPIO::cout << "***  Last Iter Stats: *** " <<endl;
	QDPIO::cout << "Cycles mat   "<< mat_1-mat_0 <<endl;
	QDPIO::cout << "Cycles inner "<< inner_1-inner_0 <<endl;
	QDPIO::cout << "Cycles axpy_norm "<< axpy_norm_1-axpy_norm_0 <<endl;
	QDPIO::cout << "Cycles axpy  "<< axpy_1-axpy_0 <<endl;
      }

      return;
    }
  }
  QDP_error_exit("Little Dirac Matrix inversion failed \n");
}
void myzgemv(int N,double * __restrict A, double * __restrict in,double c, double * __restrict out)
{
#ifdef USE_XLC_OPTIMISED_CODE
  if (((uint64_t) A) &0x1f) { printf("unaligned A\n"); fflush(stdout); while(1){};}
  if (((uint64_t)in) &0x1f) { printf("unaligned in\n"); fflush(stdout);while(1){};}
  if (((uint64_t)out)&0x1f) { printf("unaligned out\n"); fflush(stdout);while(1){};}

  qpx_zgemv(N,A,in,c,out);
  
#else

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
#endif
}

void BfmMultiGrid::BenchmarkMultiply(int NumTest)
{
  double *Ap= (double *)&Asparse[0][0][0];
  double *oo= (double *)&Krylov_Atmp[0];
  double *ii= (double *)&Krylov_p[0];

  double flops = Nvec*Nvec*8*NumTest;
  uint64_t t1 = GetTimeBase();
  for (int i=0;i<NumTest;i++){
    myzgemv(Nvec,Ap,ii,0.0,oo);
  }
  uint64_t t2 = GetTimeBase();
  if ( linop->isBoss() ) {
    printf("Bench MatMul [%d]: t2-t1 cycles %ld, %f Mflop/s %f us\n",NumTest,t2-t1,flops*1600.0/(t2-t1),(t2-t1)/1600.0);
  }
}

int gcd(int a, int b)
{
  int c = a % b;
  while(c != 0)
    {
      a = b;
      b = c;
      c = a % b;
    }
  return b;
}

void BfmMultiGrid::ApplyThread(std::vector<std::complex<double> > &in,std::vector<std::complex<double> > &out )
{
  int me, thrlen,throff;

  me = linop->thread_barrier();

  uint64_t t0 = GetTimeBase();
  HaloExchange(in);

  double flops = Nvec*Nvec*Nball*8*LocalNblock;
  int nwork = LocalNblock*Nball;
  uint64_t t1 = GetTimeBase();

  //
  // Need to divide Blocks between threads
  // 
  int threads   = linop->threads;

  int block_threads= gcd(LocalNblock,threads);
  int ball_threads = threads/block_threads;

  ThreadModelSingle WorkAllocator; 
  WorkAllocator.nthread=ball_threads;

  int block_me =me/ball_threads;
  int block_len=LocalNblock/block_threads;
  int block_off=block_me*block_len;
  
  int ball_me  =me%ball_threads;
  int ball_len;
  int ball_off;
  WorkAllocator.thread_work_nobarrier(Nball,ball_me,ball_len,ball_off);

  linop->thread_barrier();
  
  for(int b =block_off;b<block_off+block_len;b++){
  for(int mu= ball_off;mu<ball_off+ball_len;mu++){
    
    int w = ball_me*LocalNblock+b;
    double *nbr_data = (double *)HaloGetStencilData(b,mu,in);
    double *Ap= (double *)&Asparse[b][mu][0];
    double *oo= (double *)&Krylov_Atmp[w*Nvec];
    double coeff=1.0;
    if ( mu==ball_off ) coeff=0.0;
    myzgemv(Nvec,Ap,nbr_data,coeff,oo);

  }}

  uint64_t t3 = GetTimeBase();
  nwork = LocalNblock;
  linop->thread_work(nwork,me,thrlen,throff);
  for(int b=throff;b<throff+thrlen;b++){
    for(int i=0;i<Nvec;i++){
      out[b*Nvec+i]= Krylov_Atmp[b*Nvec+i];
    }
    for(int mu=1;mu<ball_threads;mu++){
      for(int i=0;i<Nvec;i++){
	out[b*Nvec+i]+=	Krylov_Atmp[mu*Nvec*LocalNblock+b*Nvec+i];
      }
    }
  }
  linop->thread_barrier();
  uint64_t t2 = GetTimeBase();
  if ( MGreport && (!me) ) {
    printf("ApplyThread -- MatMul : t2-t1 cycles %ld, %f Mflop/s %f us\n",t2-t1,flops*1600.0/(t2-t1),(t2-t1)/1600.0);
    printf("ApplyThread -- MatMul : reduce cycles %ld, mult cycles %ld \n",t2-t3,t3-t1);
  }
}

/////////////////////////////////////////
// Deflated Solver support
/////////////////////////////////////////

void BfmMultiGrid::ChebyshevInit(double _lo,double _hi,double _alpha,int _order)
{
  lo=_lo;
  hi=_hi;
  order=_order;
  alpha=_alpha;
  
  if( alpha == -0.5 ) {
    RightPrecondition=1;
    LeftPrecondition=1;
  } else if ( alpha == -1 ) {
    RightPrecondition=0;
    LeftPrecondition=1;
  } else {
    RightPrecondition=0;
    LeftPrecondition=0;
  }

  
  if(order < 2) exit(-1);
  coeffs = new double [order];
  for(int j=0;j<order;j++){
    double s=0;
    for(int k=0;k<order;k++){
      double y=cos(M_PI*(k+0.5)/order);
      double x=0.5*(y*(hi-lo)+(hi+lo));
      double f=pow(x,alpha);
      s=s+f*cos( j*M_PI*(k+0.5)/order );
    }
    coeffs[j] = s * 2.0/order;
  }
}

void BfmMultiGrid::Pleft (Fermion_t in,Fermion_t out)
{
  // P_L  = [ 1  -Mbs Mss^-1] 
  //        [ 0   0         ] 
  Fermion_t in_sbar = linop->threadedAllocFermion();
  Fermion_t tmp2= linop->threadedAllocFermion();
  Fermion_t Mtmp= linop->threadedAllocFermion();

  int me = linop->thread_barrier();

  uint64_t t1=GetTimeBase();
  ProjectToSubspace(in,PleftProj);     
  PromoteFromSubspace(PleftProj,out);  
  linop->axpy(in_sbar,out,in,-1.0);      // in_sbar = in - in_s

  uint64_t t2=GetTimeBase();
  ApplyInverse(PleftProj,PleftMss_proj); // Mss^{-1} in_s
  uint64_t t3=GetTimeBase();
  PromoteFromSubspace(PleftMss_proj,out);
  uint64_t t4=GetTimeBase();

  linop->Mprec(out,tmp2,Mtmp,DaggerNo);  // M Mss^{-1} in_s
  linop->Mprec(tmp2,out,Mtmp,DaggerYes);    
  uint64_t t5=GetTimeBase();

  ProjectToSubspace(out,PleftProj);      // Msbar s Mss^{-1}
  PromoteFromSubspace(PleftProj,tmp2);
  linop->axpy(out,tmp2,out,-1.0);

  linop->axpy(out,out,in_sbar,-1.0);     // in_sbar - Msbars Mss^{-1} in_s
  uint64_t t6=GetTimeBase();

  if ( linop->isBoss() && !me ) { 
    printf("Pleft cyc ProjProm/AInv/Prom/MdagM/ProjProm\n\t%10ld\n\t%10ld\n\t%10ld\n\t%10ld\n\t%10ld\n",t2-t1,t3-t2,t4-t3,t5-t4,t6-t5);
    fflush(stdout);
  }

  linop->threadedFreeFermion(in_sbar);
  linop->threadedFreeFermion(Mtmp);
  linop->threadedFreeFermion(tmp2);

}
void BfmMultiGrid::Pright(Fermion_t in,Fermion_t out)
{
  // P_R  = [ 1              0 ] 
  //        [ -Mss^-1 Msb    0 ] 
  Fermion_t in_sbar = linop->threadedAllocFermion();
  Fermion_t tmp = linop->threadedAllocFermion();
  Fermion_t Mtmp= linop->threadedAllocFermion();

  ProjectToSubspace(in,PleftProj);     
  PromoteFromSubspace(PleftProj,out);  
  linop->axpy(in_sbar,out,in,-1.0);       // in_sbar = in - in_s 

  linop->Mprec(in_sbar,tmp,Mtmp,DaggerNo);// M in_sbar
  linop->Mprec(tmp,out,Mtmp,DaggerYes);  
  ProjectToSubspace(out,PleftProj);           // Mssbar in_sbar  (project)

  ApplyInverse     (PleftProj,PleftMss_proj); // Mss^{-1} Mssbar 
  PromoteFromSubspace(PleftMss_proj,out);     // 

  linop->axpy(out,out,in_sbar,-1.0);     // in_sbar - Mss^{-1} Mssbar in_sbar

  linop->threadedFreeFermion(in_sbar);
  linop->threadedFreeFermion(Mtmp);
  linop->threadedFreeFermion(tmp);
}
void BfmMultiGrid::OneMinusPleft (Fermion_t in,Fermion_t out)
{
  Pleft(in,out);
  linop->axpy(out,out,in,-1.0);
}
void BfmMultiGrid::OneMinusPright(Fermion_t in,Fermion_t out)
{
  Pright(in,out);
  linop->axpy(out,out,in,-1.0);
}
void BfmMultiGrid::MPright (Fermion_t in,Fermion_t out)
{
  int me = linop->thread_barrier();
  Fermion_t tmp  = linop->threadedAllocFermion();
  Fermion_t Mtmp  = linop->threadedAllocFermion();

  Pright(in,out);
  linop->Mprec(out,tmp,Mtmp,0);
  linop->Mprec(tmp,out,Mtmp,1); 

  linop->threadedFreeFermion(Mtmp);
  linop->threadedFreeFermion(tmp);
}
void BfmMultiGrid::PleftM (Fermion_t in,Fermion_t out)
{
  int me = linop->thread_barrier();
  Fermion_t tmp  = linop->threadedAllocFermion();
  Fermion_t Mtmp  = linop->threadedAllocFermion();

  uint64_t t1=GetTimeBase();
  linop->Mprec(in,out,Mtmp,0);
  linop->Mprec(out,tmp,Mtmp,1); 
  uint64_t t2=GetTimeBase();
  Pleft(tmp,out);
  uint64_t t3=GetTimeBase();

  if( linop->isBoss() && !me ) {
    printf("MdagM/Pleft cyc = %ld/%ld\n",t2-t1,t3-t2);
  }

  linop->threadedFreeFermion(Mtmp);
  linop->threadedFreeFermion(tmp);
}


int BfmMultiGrid::MCR_PolyPrec(Fermion_t psi, Fermion_t src_uprec,Fermion_t true_res_vec) // Internal
{
  SolverControl=SolverDeflatedMdagMpc;

  KrylovPrecondition=0;

  ChebyshevInit(PreconChebyLo,
		PreconChebyHi,
		PreconChebyPower,
		PreconChebyOrder);
  double f;
  double cp,c,a,d,b,ci;
  complex<double> num;
  int me = linop->thread_barrier();

  if ( linop->isBoss() && (!me) ) { 
    FILE * fp = fopen("mcr.log","w"); // truncate log file
    fclose(fp);
    linop->InverterEnter();
  }

  Fermion_t src     = linop->threadedAllocFermion();  //(Pleft)src
  Fermion_t src_left= linop->threadedAllocFermion();  //(Pleft)src
  Fermion_t p   = linop->threadedAllocFermion(); 
  Fermion_t tmp = linop->threadedAllocFermion(); 
  Fermion_t mp  = linop->threadedAllocFermion(); 
  Fermion_t mmp = linop->threadedAllocFermion(); 
  Fermion_t r   = linop->threadedAllocFermion(); 
  Fermion_t Ap  = linop->threadedAllocFermion(); 
  Fermion_t Ar  = linop->threadedAllocFermion(); 

  //Initial residual computation & set up
  double guess = linop->norm(psi);
  double usrc  = linop->norm(src_uprec);
  if ( linop->isBoss() && !me ) {
      printf("bfmbase::MCR_PolyPrec guess %le \n",guess);
      printf("bfmbase::MCR_PolyPrec usrc  %le \n",usrc);
      fflush(stdout);
  }


  if ( SolverControl == SolverDeflatedMdagMpc ) {
    Pleft(src_uprec,src_left);
    PolyMdagMprecLeft(src_left,src);
  } else { 
    PolyMdagMprecLeft(src_uprec,src);
  }
  //  Mprec(psi,mp,tmp,DaggerNo);
  //  Mprec(mp,Ar,tmp,DaggerYes);

  SolverMatrix(psi,Ar);
  linop->axpy (r, Ar, src,-1.0);
  linop->axpy (p, Ar, src,-1.0);

  a =linop->norm(p);
  cp=linop->norm(r);

  double ssq =  linop->norm(src);
  double residual=linop->residual;
  double rsq =  residual* residual*ssq;

  if ( linop->isBoss() && !me ) {
      printf("bfmbase::MCR_PolyPrec guess %le \n",guess);
      printf("bfmbase::MCR_PolyPrec ssq %le rsq %le\n",ssq,rsq);
      printf("bfmbase::MCR_PolyPrec a %le cp %le\n",a,cp);
      fflush(stdout);
  }

  //Check if guess is really REALLY good :)
  if ( cp <= rsq ) {
    linop->threadedFreeFermion(src);
    linop->threadedFreeFermion(tmp);
    linop->threadedFreeFermion(p);
    linop->threadedFreeFermion(mp);
    linop->threadedFreeFermion(r);
    linop->threadedFreeFermion(Ap);
    linop->threadedFreeFermion(Ar);
    if ( linop->isBoss() && (!me) ) { 
      linop->InverterExit();
    }
    return 0;
  }

  if ( linop->isBoss() && !me ) {
      printf("bfmbase::MCR_PolyPrec k=0 residual %le rsq %le\n",cp,rsq);
      fflush(stdout);
  }
  //  c = Mprec(p,mp,tmp,DaggerNo);
  //  d = Mprec(mp,Ap,tmp,DaggerYes);  //Ap

  SolverMatrix(p,Ap);

  c = linop->inner_real(p,Ap);
  d = linop->norm(Ap);
  linop->axpy(Ar,Ap,Ap,0);	       //Ar
  double cp0 = c;

  int max_iter = linop->max_iter;
  //c = real( dot(r,Ar) );
  for (int k=1;k<=max_iter;k++){
    linop->iter=k;

    //c = real( dot(r,Ar) );		//c = rAr
    a = c/d;

    linop->axpy(psi,p,psi,a);		//x = x + ap
    cp = linop->axpy_norm(r,Ap,r,-a);		//r = r - aAp

  
    if ( k%1 == 0 ){
      if ( linop->isBoss() && !me ) {
	printf("bfmbase::MCR_PolyPrec: k= %d r^2= %le %le\n",k,cp,sqrt(cp/ssq));
	FILE * fp = fopen("mcr.log","a");
	fprintf(fp,"k %d resid %le\n",k,sqrt(cp/ssq));
	fclose(fp);
      }
    }

    linop->axpy(true_res_vec,r,r,0.0);
    double rr = sqrt(cp/ssq);
    linop->InverterLogIteration(k, rr,src,psi,mp,mmp,tmp);


    // Stopping condition
    if ( cp <= rsq ) { 

	if ( linop->isBoss() && (!me) ) { 
	  linop->InverterExit();
	}
	
	if ( linop->isBoss() && !me ) printf("bfmbase::MCR_PolyPrec converged in %d iterations\n",k);

	PolyMdagMprecRight(psi,tmp); // Psi is now solution of PL M psi = PL eta
	linop->axpy(psi,tmp,tmp,0.0);

	double true_residual;

	if ( SolverControl == SolverDeflatedMdagMpc ) {
	  double nn;

	  // Test that MssInv works
	  ProjectToSubspace(psi,PleftProj);     
	  PromoteFromSubspace(PleftProj,tmp);  
	  
	  linop->Mprec(tmp,mp,r,0);
	  linop->Mprec(mp,Ap ,r,1); 
	  ProjectToSubspace(Ap,PleftProj);     
	  PromoteFromSubspace(PleftProj,Ap);  

	  ApplyInverse     (PleftProj,PleftMss_proj); 
	  PromoteFromSubspace(PleftMss_proj,mp);  


	  nn=linop->norm(tmp);
	  if ( linop->isBoss() && !me ) 
	    printf("bfmbase::MCR_PolyPrec: subspace norm is %le \n",nn);

	  nn=linop->norm(Ap);
	  if ( linop->isBoss() && !me ) 
	    printf("bfmbase::MCR_PolyPrec: Mss norm is %le \n",nn);


	  nn=linop->norm(mp);
	  if ( linop->isBoss() && !me ) 
	    printf("bfmbase::MCR_PolyPrec: MssInv norm is %le \n",nn);


	  linop->axpy(mp,mp,tmp,-1.0);
	  nn=linop->norm(mp);
	  if ( linop->isBoss() && !me ) 
	    printf("bfmbase::MCR_PolyPrec: MssInv error is %le \n",nn);

	  // Compute true residual of PleftM system
	  PleftM(psi,tmp);
	  linop->axpy(tmp,tmp,src_left,-1.0);
	  true_residual = sqrt(linop->norm(tmp)/linop->norm(src_left));
	  if ( linop->isBoss() && !me ) 
	    printf("bfmbase::MCR_PolyPrec: true PleftM residual is %le \n",true_residual);

	  // This should equal MPright
	  MPright(psi,tmp);
	  linop->axpy(tmp,tmp,src_left,-1.0);
	  true_residual = sqrt(linop->norm(tmp)/linop->norm(src_left));
	  if ( linop->isBoss() && !me ) 
	    printf("bfmbase::MCR_PolyPrec: true MPright residual is %le \n",true_residual);


	  nn=linop->norm(psi);
	  if ( linop->isBoss() && !me ) 
	    printf("bfmbase::MCR_PolyPrec: psi %le \n",nn);

	  // 8 Feb 2013: BUG : I am finding P_L M P_R psi != P_L M psi
	  //                 
	  //          Should find also P_L M = M P_R
	  //
	  // Apply P_R psi
	  Pright(psi,tmp);
	  linop->axpy(psi,tmp,tmp,0.0);

	  nn=linop->norm(psi);
	  if ( linop->isBoss() && !me ) 
	    printf("bfmbase::MCR_PolyPrec: P_R psi %le \n",nn);

	  Pright(psi,tmp);
	  linop->axpy(psi,tmp,tmp,0.0);

	  nn=linop->norm(psi);
	  if ( linop->isBoss() && !me ) 
	    printf("bfmbase::MCR_PolyPrec: P_R P_R psi %le \n",nn);

	  // Redetermine the defect
	  PleftM(psi,tmp);
	  linop->axpy(tmp,tmp,src_left,-1.0);
	  true_residual = sqrt(linop->norm(tmp)/linop->norm(src_left));
	  if ( linop->isBoss() && !me ) 
	    printf("bfmbase::MCR_PolyPrec: true PleftM residual is %le \n",true_residual);

	  Pleft(psi,tmp);
	  nn=linop->norm(psi);
	  if ( linop->isBoss() && !me ) 
	    printf("bfmbase::MCR_PolyPrec: true PleftM psi is %le \n",nn);
	  
	  nn=linop->norm(tmp);
	  if ( linop->isBoss() && !me ) 
	    printf("bfmbase::MCR_PolyPrec: true Pleft PleftM psi is %le \n",nn);

	  // 1-P_R determined by Mss^-1 eta_s
	  LittleDopSolverResidual = 1.0e-9;
	  ProjectToSubspace(src_uprec,PleftProj);   
 	  ApplyInverse(PleftProj,PleftMss_proj);     // Mss^inv eta_s
	  PromoteFromSubspace(PleftMss_proj,tmp);
	  
	  double term = linop->norm(tmp);
	  if ( linop->isBoss() && !me ) 
	    printf("bfmbase::MCR_PolyPrec: Subspace soln is %le \n",term);
	  term = linop->norm(psi);
	  if ( linop->isBoss() && !me ) 
	    printf("bfmbase::MCR_PolyPrec: Complement soln is %le \n",term);

	  linop->axpy(psi,tmp,psi,1.0);

	  linop->Mprec(psi,mp,tmp,0);
	  linop->Mprec(mp,Ap ,tmp,1); 

	  linop->axpy(tmp,src_uprec,Ap,-1.0);
	  true_residual = sqrt(linop->norm(tmp)/linop->norm(src_uprec));
	  if ( linop->isBoss() && !me ) 
	    printf("bfmbase::MCR_PolyPrec: true full residual is %le \n",true_residual);

	  if ( linop->isBoss() && !me ) {
	  FILE * fp = fopen("convergence.dat","a");
	  fprintf(fp,"************************\n" );
	  fprintf(fp,"Subspace [%f,%f] power %f  order %d\n", 
		  SubspaceChebyLo,
		  SubspaceChebyHi,
		  SubspaceChebyPower,
		  SubspaceChebyOrder);
	  fprintf(fp,"Precon  [%f,%f] power %f  order %d\n",
		  PreconChebyLo,
		  PreconChebyHi,
		  PreconChebyPower,
		  PreconChebyOrder);
	  fprintf(fp,"Convergence %d iter resid %le\n",k,true_residual);
	  fprintf(fp,"************************\n" );
	  fclose(fp);
	  }
	} else { 


	  linop->Mprec(psi,mp,tmp,0);
	  linop->Mprec(mp,Ap,tmp,1); 
	  linop->axpy(tmp,src_uprec,Ap,-1.0);
	  double true_residual = sqrt(linop->norm(tmp)/linop->norm(src_uprec));

	  if ( linop->isBoss() && !me ) 
	    printf("bfmbase::MCR_PolyPrec: true unprec residual is %le \n",true_residual);


	  if ( linop->isBoss() && !me ) {
	  FILE * fp = fopen("convergence.dat","a");
	  fprintf(fp,"************************\n" );
	  fprintf(fp,"NO Subspace projection\n" );
	  fprintf(fp,"Precon  [%f,%f] power %f  order %d\n",
		  PreconChebyLo,
		  PreconChebyHi,
		  PreconChebyPower,
		  PreconChebyOrder);
	  fprintf(fp,"Convergence %d iter resid %le\n",k,true_residual);
	  fprintf(fp,"************************\n" );
	  fclose(fp);
	  }
	}

	linop->threadedFreeFermion(src);
	linop->threadedFreeFermion(src_left);
	linop->threadedFreeFermion(tmp);
	linop->threadedFreeFermion(p);
	linop->threadedFreeFermion(mp);
	linop->threadedFreeFermion(mmp);
	linop->threadedFreeFermion(r);
	linop->threadedFreeFermion(Ap);
	linop->threadedFreeFermion(Ar);
	return k;
	
    }
    ci = d; //(Ap,Ap)
   //d = real( dot(r,Ar) );

    //    d = Mprec(r,mp,tmp,DaggerNo); //(r,Ar)
    //        Mprec(mp,Ar,tmp,DaggerYes);

    SolverMatrix(r,Ar);
    d = linop->inner_real(r,Ar);
    b = d/c;
    c = d;
    linop->axpy(p,p,r,b);
    d = linop->axpy_norm(Ap,Ap,Ar,b);


  }
  if ( linop->isBoss() && !me ) printf("bfmbase::MCR_PolyPrec: CG not converged \n");
  linop->threadedFreeFermion(src);
  linop->threadedFreeFermion(tmp);
  linop->threadedFreeFermion(p);
  linop->threadedFreeFermion(mp);
  linop->threadedFreeFermion(mmp);
  linop->threadedFreeFermion(r);
  linop->threadedFreeFermion(Ap);
  linop->threadedFreeFermion(Ar);
  if ( linop->isBoss() && (!me) ) { 
    linop->InverterExit();
  }

  return -1;

}

void BfmMultiGrid::PolyMdagMprecRight(Fermion_t in,Fermion_t out)
{
  int me = linop->thread_barrier();
  if ( RightPrecondition ) { 
    if ( linop->isBoss() && !me ) {
      printf("Chebyshev right preconditioning matrix x^%f\n",alpha);fflush(stdout);
    }
    PolyMdagMprec(in,out);
  } else { 
    linop->axpy(out,in,in,0.0); 
  }
  linop->thread_barrier();
}
void BfmMultiGrid::PolyMdagMprecLeft(Fermion_t in,Fermion_t out)
{

  int me = linop->thread_barrier();

  if ( KrylovPrecondition ) {

    preclinop->thread_barrier();
    double restore_rsd  = linop->residual;
    double restore_iter=  linop->max_iter;
    preclinop->thread_barrier();
    preclinop->residual = InnerKrylovResidual;
    preclinop->max_iter = InnerKrylovIterMax;
    preclinop->thread_barrier();

    Fermion_t pin;
    Fermion_t pout;

    if ( linop->Ls < preclinop->Ls ) {
      printf("Cannot increase Ls in preconditioner");
      exit(0);
    }
    if ( linop->Ls != preclinop->Ls ) {
      if ( linop->isBoss() && !me ) printf("Reducing Ls from %d to %d\n",linop->Ls,preclinop->Ls);

      pin =preclinop->threadedAllocFermion();
      pout=preclinop->threadedAllocFermion();
      for(int s=0;s<preclinop->Ls/2;s++){
	preclinop->copy_slice(in ,s,linop->Ls,pin,s,preclinop->Ls);
	preclinop->copy_slice(in ,linop->Ls-s-1,linop->Ls,pin,preclinop->Ls-s-1,preclinop->Ls);
      }
      linop->axpy(out,in,in,99);
    } else { 
      pin = in;
      pout=out;
    }

    if ( linop != preclinop) {
      preclinop->thread_barrier();
      if ( preclinop->SPIcomms() && !me ) preclinop->comm_init();
      preclinop->thread_barrier();
    }


    preclinop->axpy(pout,pin,pin,0.0);
    preclinop->thread_barrier();

    
    //    preclinop->MCR_prec(pout,pin);

    //    preclinop->CGNE_prec(pout,pin);


    int nshift=1;
    int single=0;
    double alpha[1]     = {1.0};
    double shifts[1]    = {InnerKrylovShift};
    double mresidual[1] = {InnerKrylovResidual};
    preclinop->CGNE_prec_MdagM_multi_shift(&pout,
					   pin,
					   shifts,
					   alpha,
					   nshift,
					   mresidual,
					   single);

    preclinop->thread_barrier();
    preclinop->residual=restore_rsd;
    preclinop->max_iter=restore_iter;
    preclinop->thread_barrier();

    if ( linop->cbLs != preclinop->cbLs ) {
      for(int s=0;s<preclinop->Ls/2;s++){
	preclinop->copy_slice(pout,s,preclinop->Ls,out,s,linop->Ls);
	preclinop->copy_slice(pout,preclinop->Ls-s-1,preclinop->Ls,out,linop->Ls-s-1,linop->Ls);
      }
      preclinop->threadedFreeFermion(pin);
      preclinop->threadedFreeFermion(pout);
    }

    if ( linop != preclinop) {
      linop->thread_barrier();
      if ( linop->SPIcomms() && !me )  linop->comm_init();
      linop->thread_barrier();
    }

    return;
  } else if ( LeftPrecondition ) { 
    if ( linop->isBoss() && !me ) printf("Chebyshev left preconditioning matrix x^%f\n",alpha);fflush(stdout);
    PolyMdagMprec(in,out);
    linop->thread_barrier();
  } else { 
    linop->axpy(out,in,in,0.0); 
    linop->thread_barrier();
  }
}

void BfmMultiGrid::PolyMdagMprec(Fermion_t in,Fermion_t out)
{
  PolyMdagMprec(linop,in,out);
}

void BfmMultiGrid::DeflatedSolverMpc(Fermion_t in,Fermion_t out)
{
    
  Fermion_t tmp = linop->threadedAllocFermion();
  Fermion_t tmp2= linop->threadedAllocFermion();
  Fermion_t Mtmp= linop->threadedAllocFermion();

  exit(-1);

  PolyMdagMprecRight(in,tmp2);
  PleftM(tmp2,tmp);
  PolyMdagMprecLeft(tmp,out);
    
  linop->threadedFreeFermion(tmp);
  linop->threadedFreeFermion(tmp2);
  linop->threadedFreeFermion(Mtmp);
}

void BfmMultiGrid::DeflatedSolverMdagMpc(Fermion_t in,Fermion_t out)
{
    
  Fermion_t tmp = linop->threadedAllocFermion();
  Fermion_t tmp2= linop->threadedAllocFermion();
  Fermion_t Mtmp= linop->threadedAllocFermion();

  uint64_t t1 = GetTimeBase();

  PolyMdagMprecRight(in,tmp2);

  uint64_t t2 = GetTimeBase();

  PleftM(tmp2,tmp);

  uint64_t t3 = GetTimeBase();

  PolyMdagMprecLeft(tmp,out);

  uint64_t t4 = GetTimeBase();
  
  int me = linop->thread_barrier();
  if( linop->isBoss() && !me ) {
    printf("Deflated Poly/PleftM/Poly cyc = %ld/%ld/%ld\n",t2-t1,t3-t2,t4-t3);
  }
    
  linop->threadedFreeFermion(tmp);
  linop->threadedFreeFermion(tmp2);
  linop->threadedFreeFermion(Mtmp);
}

void BfmMultiGrid::UndeflatedSolverMpc(Fermion_t in,Fermion_t out)
{
  Fermion_t tmp = linop->threadedAllocFermion();
  Fermion_t Mtmp= linop->threadedAllocFermion();
  
  linop->Mprec(in,out,Mtmp,DaggerYes);    // Non herm: M P_{-1}(Mdag M) Mdag
  PolyMdagMprecRight(out,tmp);

  linop->Mprec(tmp,out,Mtmp,DaggerNo);  
    
  linop->threadedFreeFermion(tmp);
  linop->threadedFreeFermion(Mtmp);
};

void BfmMultiGrid::UndeflatedSolverMdagMpc(Fermion_t in,Fermion_t out)
{
  Fermion_t tmp = linop->threadedAllocFermion();
  Fermion_t Mtmp= linop->threadedAllocFermion();
  
  PolyMdagMprecRight(in,tmp);

  linop->Mprec(tmp,out,Mtmp,DaggerNo);  
  linop->Mprec(out,tmp,Mtmp,DaggerYes);    

  PolyMdagMprecLeft(tmp,out);
    
  linop->threadedFreeFermion(tmp);
  linop->threadedFreeFermion(Mtmp);
};

int BfmMultiGrid::GCR_PolyPrec(Fermion_t psi, Fermion_t src,int m) // Internal
{
  double f;
  double cp,c,a,d,b;
  int me = linop->thread_barrier();

  SolverControl = SolverUndeflatedMpc;
  ChebyshevInit(PreconChebyLo,
		PreconChebyHi,
		-1.0,
		PreconChebyOrder);

  if ( linop->isBoss() && (!me) ){
    linop->InverterEnter();
  }
 
  Fermion_t p   =linop->threadedAllocFermion();
  Fermion_t tmp =linop->threadedAllocFermion();
  Fermion_t Mtmp=linop->threadedAllocFermion();
  Fermion_t Ap  =linop->threadedAllocFermion();
  Fermion_t Ar  =linop->threadedAllocFermion();
  Fermion_t r   =linop->threadedAllocFermion();
  
  //Initial residual computation & set up
  double guess = linop->norm(psi);

  SolverMatrix(psi,Ar);
  linop->axpy(r, Ar, src,-1.0);     //r = b - Ax ; b->src
  linop->copy(p,r);                 //p0 = r0
 
  a =linop->norm(p);
  cp=linop->norm(r);
 
  double ssq =  linop->norm(src);
  double residual=linop->residual;
  double rsq =  residual* residual*ssq;
 
  //Check if guess is really REALLY good :)
  if ( cp <= rsq ){
    
        if ( linop->verbose && linop->isBoss() && !me )    {
	  printf("bfmbase::GCR_prec k=0 converged - suspiciously nice guess %le %le\n",cp,rsq);
        }
        linop->threadedFreeFermion(tmp);
        linop->threadedFreeFermion(Mtmp);
        linop->threadedFreeFermion(p);
        linop->threadedFreeFermion(Ar);
        linop->threadedFreeFermion(Ap);
        linop->threadedFreeFermion(r);
        if ( linop->isBoss() && (!me) ) {
            linop->InverterExit();
        }
        return 0;
  }

  if ( linop->verbose && linop->isBoss() && !me )
    printf("bfmbase::GCR_prec k=0 residual %le rsq %le\n",cp,rsq);

 
  Fermion_t allAp[m-1];
  Fermion_t allp[m-1];
  double ApAp[m-1];
  for(int i=0;i<m-1;i++) {
    allAp[i]   = linop->threadedAllocFermion();
    allp[i]    = linop->threadedAllocFermion();
  } 

  SolverMatrix(p,Ap);

  std::complex<double> alpha;
  std::complex<double> c_rAp  = linop->inner(r,Ap);          //c = (r,Ap)
  double               c_ApAp = linop->inner_real(Ap,Ap);
 
  linop->copy(allAp[0],Ap);
  linop->copy(allp[0],p);
  ApAp[0] = c_ApAp;
  
  for (int k=1; k<=linop->max_iter; k++) {
 
    std::complex<double> beta[k];
    linop->iter=k;

    alpha = c_rAp/c_ApAp;                                      //a = (r,Ap)/(Ap,Ap) 

    linop->caxpy(psi,p,psi,real(alpha),imag(alpha));	       //x = x + ap
    linop->caxpy(r,Ap,r,-real(alpha),-imag(alpha));	       //r = r - aAp
    cp = linop->norm(r);
  
    SolverMatrix(r,Ar);
    
    linop->copy(Ap,Ar) ; // Ap=Ar         // Why

    linop->copy(p,r);    // p=r
    for(int i=0; i<k%m; i++) {
      beta[i] = linop->inner(Ar,allAp[i])/ApAp[i];                        //alpha = (Ar,Ap)/(Ap,Ap)
      linop->caxpy(Ap,allAp[i],Ap,-real(beta[i]),-imag(beta[i]));         //Ap = Ar + B Ap
      linop->caxpy(p,allp[i],p,-real(beta[i]),-imag(beta[i]));            //p  =  r + B p
    }
 
    c_ApAp = linop->norm(Ap);                    //d = (Ap,Ap)
    c_rAp  = linop->inner(r,Ap);                 //c = (r,Ap)
  
                    
    if( k%m < m-1 ) {      // we dont need to store when we need to restart
      linop->copy(allAp[k%m],Ap);
      linop->copy(allp[k%m],p);
      ApAp[k%m] = c_ApAp;
    }
 
    if ( k%10 == 0){
      if (linop->isBoss() && !me ){
	printf("bfmbase::GCR_prec: k= %d r^2= %le %le \n",k,cp,sqrt(cp/ssq));
      }
    }
 
    // Stopping condition
    if ( cp <= rsq ) {
      
      if (linop->isBoss() && !me ) printf("bfmbase::GCR_prec converged in %d iterations\n",k);

      SolverMatrix(psi,Ap);

      linop->axpy(tmp,src,Ap,-1.0);
      double true_residual = sqrt(linop->norm(tmp)/linop->norm(src));
      if ( linop->isBoss() && !me )
	printf("bfmbase::GCR_prec: true residual is %le \n",true_residual);
      
      linop->threadedFreeFermion(tmp);
      linop->threadedFreeFermion(Mtmp);
      linop->threadedFreeFermion(p);
      linop->threadedFreeFermion(r);
      linop->threadedFreeFermion(Ap);
      linop->threadedFreeFermion(Ar);
      if ( linop->isBoss() && (!me) ){
	linop->InverterExit();
      }
      return k;
    }
 
 
  }

  if (linop->isBoss() && !me ) printf("bfmbase::GCR_prec: CG not converged \n");
  linop->threadedFreeFermion(tmp);
  linop->threadedFreeFermion(Mtmp);
  linop->threadedFreeFermion(p);
  linop->threadedFreeFermion(r);
  linop->threadedFreeFermion(Ap);
  linop->threadedFreeFermion(Ar);
  if ( linop->isBoss() && (!me) ) {
    linop->InverterExit();
  }
  
  return -1;

}

  /*
   * Note on  Balancing preconditioner
   *
   * Presently solver P_L M x = P_L b    ... where
   *  
   *   P_L  = [ 1  -Mbs Mss^-1] 
   *          [ 0   0         ] 
   * and  
   *      P_L M =  P_L Mbb Mbs = Mbb - Mbs Mss^-1 Msb   0
   *                   Msb Mss   0                      0
   * 
   * Nabben and Vuik SIAM J sci comp V27 N05 :
   *
   *      P_D = I - A Z (Z^TAZ)^-1  Z^T
   *
   * Define the preconditioning operator to raise the Evals of subspace to large lambda
   *      P_B = [I - A Z (Z^TAZ)^-1  Z^T + lambda Z Z^T ]
   */ 

int  BfmMultiGrid::rCG_PolyPrec(Fermion_t psi, Fermion_t src_uprec,Fermion_t resid)
{
  int steps= 10;
  linop->thread_barrier();
  int miter=linop->max_iter;
  linop->thread_barrier();
  int iters=0;
  linop->max_iter = steps;
  linop->thread_barrier();

  while(iters < miter){
    int it = CG_PolyPrec(psi,src_uprec,resid);
    iters+=steps;
    if ( it > 0 ) {
      linop->thread_barrier();
      linop->max_iter = miter;
      linop->thread_barrier();
      return iters;
    }
  }
  return -1;
}

int BfmMultiGrid::CG_PolyPrec(Fermion_t psi, Fermion_t src,Fermion_t resid) 
{
  double f;
  double rtzp,rtz,a,d,b;

  int me = linop->thread_barrier();


  // Deflated precond CG (3.6 in Saad umsi-98-97)
  KrylovPrecondition=1;

  ChebyshevInit(PreconChebyLo,
		PreconChebyHi,
		-1.0,
		PreconChebyOrder);

  /*
  if ( linop->isBoss() && (!me) ) { 
    printf("BfmMultiGrid::CG_PolyPrec: Precondition with Chebyshev x^{%f} [%f,%f] O(%d) \n",
	   PreconChebyLo,
	   PreconChebyHi,
	   -1.0,
	   PreconChebyOrder
	   );
    linop->InverterEnter();
  }
  */

  Fermion_t p   = linop->threadedAllocFermion(); 
  Fermion_t x   = psi;
  Fermion_t z   = linop->threadedAllocFermion(); 

  Fermion_t tmp = linop->threadedAllocFermion(); 
  Fermion_t mp  = linop->threadedAllocFermion(); 
  Fermion_t mmp = linop->threadedAllocFermion(); 
  Fermion_t r   = linop->threadedAllocFermion(); 
  Fermion_t mu   = linop->threadedAllocFermion(); 

  //Initial residual computation & set up
  double guess = linop->norm(psi);

  ///////////////////////////////////
  // x_{-1} = psi
  ///////////////////////////////////
  linop->Mprec(x,mp,tmp,DaggerNo);
  linop->Mprec(mp,mmp,tmp,DaggerYes);
  linop->axpy (r, mmp, src,-1.0);        // r_{-1} = src - A x

  ///////////////////////////////////
  // Choose x_0 such that 
  // x_0 = guess +  (A_ss^inv) r_s
  // 
  // W^T (src - A x_0) = src_s - A guess_s - r_s  
  //                   = src_s - (A guess)_s - src_s  + (A guess)_s 
  //                   = 0 
  ///////////////////////////////////
  ProjectToSubspace(r,PleftProj);     
  if(linop->isBoss()&&!me){ 
    printf("Inverting little dirac op\n");fflush(stdout);
  }
  ApplyInverse(PleftProj,PleftMss_proj); // Ass^{-1} r_s
  if(linop->isBoss()&&!me){ 
    printf("Inverted little dirac op\n");fflush(stdout);
  }
  PromoteFromSubspace(PleftMss_proj,mp);  

  linop->axpy(x,x,mp,1.0);

  linop->Mprec(x,mp,tmp,DaggerNo);
  linop->Mprec(mp,mmp,tmp,DaggerYes);
  linop->axpy (r, mmp, src,-1.0);    // Recomputes r=src-x0

  ProjectToSubspace(r,PleftProj);     
  double nv=norm_vec(PleftProj);
  if(linop->isBoss()&&!me){ 
    printf("Residual after projection= %le\n",nv);
  }

  // Compute z = M x
  PolyMdagMprecLeft(r,z);
  rtzp =real(linop->inner(r,z));

  ///////////////////////////////////////
  // Solve for Mss mu = P A z and set p = z-mu
  ///////////////////////////////////////
  linop->Mprec(z,mp,tmp,DaggerNo);
  linop->Mprec(mp,mmp,tmp,DaggerYes); //Az

  ProjectToSubspace(mmp,PleftProj);     
  ApplyInverse(PleftProj,PleftMss_proj); // Mss^{-1} Az
  PromoteFromSubspace(PleftMss_proj,mu);  

  linop->axpy (p, mu, z,-1.0);

  double ssq =  linop->norm(src);
  double rsq =  linop->residual* linop->residual*ssq;

  if ( linop->isBoss() && !me ) 
    printf("BfmMultiGrid::CG_PolyPrec k=0 residual %le rsq %le\n",rtzp,rsq);

  uint64_t t_start=GetTimeBase();
  uint64_t t1;
  uint64_t t_Mss=0;
  uint64_t t_proj=0;
  for (int k=1;k<=linop->max_iter;k++){

    // Test hack -- deliberately perturb matrix
#if 0
    if ( k&0x4 ){ 
     ChebyshevInit(PreconChebyLo,
		PreconChebyHi,
		-1.0,
		PreconChebyOrder);
    } else {
     ChebyshevInit(2*PreconChebyLo,
		  PreconChebyHi,
		   -1.0,
		   2*PreconChebyOrder);
    }
#endif
    
    rtz=rtzp;
    d        =linop->Mprec(p,mp,tmp,0,1);// Dag no
              linop->Mprec(mp,mmp,tmp,1);// Dag yes
    a = rtz/d;

    //    if(linop->isBoss() && !me ) { printf("CG poly 1\n"); fflush(stdout);}
    //    linop->norm(x);

    linop->axpy(x,p,x,a);
    linop->axpy(r,mmp,r,-a);
    linop->axpy(resid,r,r,0.0);

    //    if(linop->isBoss() && !me ) { printf("CG poly 2\n"); fflush(stdout);}
    //    linop->norm(x);

    // Compute z = M x
    PolyMdagMprecLeft(r,z);

    //    if(linop->isBoss() && !me ) { printf("CG poly 3\n"); fflush(stdout);}
    //    linop->norm(x);

    rtzp =real(linop->inner(r,z));
    b = rtzp/rtz;

    //    if(linop->isBoss() && !me ) { printf("CG poly 4\n"); fflush(stdout);}
    //    linop->norm(x);

    // WAW mu = W A z
    linop->Mprec(z,mp,tmp,DaggerNo);
    linop->Mprec(mp,mmp,tmp,DaggerYes);

    //    if(linop->isBoss() && !me ) { printf("CG poly 5\n"); fflush(stdout);}
    //    linop->norm(x);

    t1=GetTimeBase();
    ProjectToSubspace(mmp,PleftProj);     // WAz
    t_proj+=GetTimeBase()-t1;

    //    if(linop->isBoss() && !me ) { printf("CG poly 6\n"); fflush(stdout);}
    //    linop->norm(x);
    t1=GetTimeBase();
    ApplyInverse(PleftProj,PleftMss_proj); // Mss^{-1} in_s       // (WtAW)^[-1] Az
    t_Mss+=GetTimeBase()-t1;

    //    if(linop->isBoss() && !me ) { printf("CG poly 7\n"); fflush(stdout);}
    //    linop->norm(x);

    t1=GetTimeBase();
    PromoteFromSubspace(PleftMss_proj,mu);     // mu = Wt  (Wt A W)^{-1} WAz
    t_proj+=GetTimeBase()-t1;

    linop->axpy(p,p,z,b);
    linop->axpy(p,mu,p,-1.0);

    //    if(linop->isBoss() && !me ) { printf("CG poly 8\n"); fflush(stdout);}
    //    linop->norm(x);

    {
      double rr=sqrt(rtzp/ssq);
      if (linop->isBoss() && !me ){
	printf("bfmbase::CG_PolyPrec: k= %d residual = %le \n",k,rr);
      }
      //  linop->InverterLogIteration(k, rr,src,x,mp,mmp,tmp);
    }

    // Stopping condition
    if ( rtzp <= rsq ) { 
	
	t1=GetTimeBase();
      if ( linop->isBoss() && !me ) {
	printf("BfmMultiGrid::CG_PolyPrec converged in %d iterations %f s\n",
					   k,1.0e-9*(t1-t_start)/1.6);
	printf("BfmMultiGrid::CG_PolyPrec MssInv %f s\n",1.0e-9*(t_Mss)/1.6);
	printf("BfmMultiGrid::CG_PolyPrec Proj   %f s\n",1.0e-9*(t_proj)/1.6);
	printf("BfmMultiGrid::CG_PolyPrec other  %f s\n",1.0e-9*(t1-t_start-t_Mss-t_proj)/1.6);
      }
	
	linop->Mprec(x,mp,tmp,0);
	linop->Mprec(mp,mmp,tmp,1); 
	linop->axpy(tmp,src,mmp,-1.0);
	
	double  mpnorm = sqrt(linop->norm(mp));
	double mmpnorm = sqrt(linop->norm(mmp));
	double psinorm = sqrt(linop->norm(x));
	double srcnorm = sqrt(linop->norm(src));
	double tmpnorm = sqrt(linop->norm(tmp));
	double true_residual = tmpnorm/srcnorm;
	if ( linop->isBoss() && !me ) {
	  printf("BfmMultiGrid::CG_PolyPrec: true residual is %le, solution %le, source %le \n",true_residual,psinorm,srcnorm);
	  printf("BfmMultiGrid::CG_PolyPrec: target residual was %le \n",linop->residual);
	  printf("BfmMultiGrid::CG_PolyPrec: mp %le, mmp %le\n",mpnorm,mmpnorm);
	}
	
	linop->threadedFreeFermion(tmp);
	linop->threadedFreeFermion(p);
	linop->threadedFreeFermion(z);
	linop->threadedFreeFermion(mu);
	linop->threadedFreeFermion(mp);
	linop->threadedFreeFermion(mmp);
	linop->threadedFreeFermion(r);
	if ( linop->isBoss() && (!me) ) { 
	  linop->InverterExit();
	}
	return k;
    }

  }
  if ( linop->isBoss() && !me ) printf("BfmMultiGrid::CG_PolyPrec: CG not converged \n");
  linop->threadedFreeFermion(tmp);
  linop->threadedFreeFermion(p);
  linop->threadedFreeFermion(mp);
  linop->threadedFreeFermion(mmp);
  linop->threadedFreeFermion(z);
  linop->threadedFreeFermion(mu);
  linop->threadedFreeFermion(r);
  if ( linop->isBoss() && (!me) ) { 
    linop->InverterExit();
  }

  return -1;
}
//
// Flexible Conjugate Gradients [Notay 2000 SIAM Vol 22 No 4 Jsci Comput]
// Combine with Deflation by Saad
//
int BfmMultiGrid::fCG_PolyPrec(Fermion_t psi, Fermion_t src) 
{
  double f;
  double a,d,b, r_dot_z,r_dot_z_p;

  int me = linop->thread_barrier();

  KrylovPrecondition=1;

  // Deflated precond CG (3.6 in Saad umsi-98-97)
  ChebyshevInit(PreconChebyLo,
		PreconChebyHi,
		-1.0,
		PreconChebyOrder);

  if ( linop->isBoss() && (!me) ) { 
    printf("BfmMultiGrid::FlexCG_PolyPrec: Precondition with Chebyshev x^{%f} [%f,%f] O(%d) \n",
	   PreconChebyLo,
	   PreconChebyHi,
	   -1.0,
	   PreconChebyOrder
	   );
    linop->InverterEnter();
  }

  Fermion_t x   = psi;
  Fermion_t tmp = linop->threadedAllocFermion(); 
  Fermion_t mp  = linop->threadedAllocFermion(); 
  Fermion_t mmp = linop->threadedAllocFermion(); 
  Fermion_t r   = linop->threadedAllocFermion(); 
  Fermion_t mu  = linop->threadedAllocFermion(); 
  Fermion_t z   = linop->threadedAllocFermion(); 

  const int mmax = 20;
  int m;
  Fermion_t  p[mmax];
  Fermion_t Ap[mmax];
  double   pAp[mmax];
  for(m=0;m<mmax;m++){
    p [m] = linop->threadedAllocFermion();
    Ap[m] = linop->threadedAllocFermion();
  }
  //Initial residual computation & set up
  double guess = linop->norm(psi);

  ///////////////////////////////////
  // x_{-1} = psi
  ///////////////////////////////////
  linop->Mprec(x,mp,tmp,DaggerNo);
  linop->Mprec(mp,Ap[0],tmp,DaggerYes);
  linop->axpy (r, Ap[0], src,-1.0);        // r_{-1} = src - A x

  ///////////////////////////////////
  // Choose x_0 such that 
  // x_0 = guess +  (A_ss^inv) r_s
  // 
  // W^T (src - A x_0) = src_s - A guess_s - r_s  
  //                   = src_s - (A guess)_s - src_s  + (A guess)_s 
  //                   = 0 
  ///////////////////////////////////
  ProjectToSubspace(r,PleftProj);     
  ApplyInverse(PleftProj,PleftMss_proj); // Ass^{-1} r_s
  PromoteFromSubspace(PleftMss_proj,mp);  

  linop->axpy(x,x,mp,1.0);
  linop->axpy(p[0],x,x,0.0);

  linop->Mprec(p[0],mp,tmp,DaggerNo);
  linop->Mprec(mp,Ap[0],tmp,DaggerYes);
  linop->axpy (r, Ap[0], src,-1.0);    // Recomputes r=src-x0

  ProjectToSubspace(r,PleftProj);     
  double nv=norm_vec(PleftProj);
  if(linop->isBoss()&&!me){ 
    printf("Residual after projection= %le\n",nv);
  }

  // Compute z = M x
  PolyMdagMprecLeft(r,z);
  // Fixme: also could run a Krylov solver
  r_dot_z =real(linop->inner(r,z));
  
  ///////////////////////////////////////
  // Solve for Mss mu = P A z and set p = z-mu
  ///////////////////////////////////////
  linop->Mprec(z,mp,tmp,DaggerNo);
  linop->Mprec(mp,mmp,tmp,DaggerYes); //Az

  ProjectToSubspace(mmp,PleftProj);     
  ApplyInverse(PleftProj,PleftMss_proj); // Mss^{-1} Az
  PromoteFromSubspace(PleftMss_proj,mu);  

  linop->axpy (p[0], mu, z,-1.0);

  double ssq =  linop->norm(src);
  double rsq =  linop->residual* linop->residual*ssq;

  if ( linop->isBoss() && !me ) 
    printf("BfmMultiGrid::FlexCG_PolyPrec k=0 residual %le rsq %le\n",r_dot_z,rsq);

  uint64_t t_start = GetTimeBase();
  uint64_t t_Mss=0;
  uint64_t t_proj=0;
  uint64_t t_ortho=0;
  uint64_t t_prec=0;
  uint64_t t1;

  for (int k=0;k<=linop->max_iter;k++){

    /*
    if ( k&0x4 ){ 
     ChebyshevInit(PreconChebyLo,
		PreconChebyHi,
		-1.0,
		PreconChebyOrder);
    } else {
     ChebyshevInit(2*PreconChebyLo,
		  PreconChebyHi,
		   -1.0,
		   2*PreconChebyOrder);
    }
    */
    int peri_k  = k % mmax;
    int peri_kp = (k+1) % mmax;
    
    pAp[peri_k]        =linop->Mprec(p[peri_k],mp,tmp,0,1);// Dag no
                        linop->Mprec(mp,Ap[peri_k],tmp,1);// Dag yes

    a = r_dot_z/pAp[peri_k];
    linop->axpy(x, p[peri_k],x,a);
    linop->axpy(r,Ap[peri_k],r,-a);

    // Compute z = M x;  This is w_i = B(r_i) in Notay
    fflush(stdout);

    t1=GetTimeBase();
    linop->thread_barrier();
    PolyMdagMprecLeft(r,z);
    linop->thread_barrier();
    t_prec+=GetTimeBase()-t1;

    double r_dot_z_p = r_dot_z;
    r_dot_z =real(linop->inner(r,z));


    // Solve WAW mu = W A z
    linop->Mprec(z,mp,tmp,DaggerNo);
    linop->Mprec(mp,mmp,tmp,DaggerYes);

    t1=GetTimeBase();
    ProjectToSubspace(mmp,PleftProj);     
    t_proj+=GetTimeBase()-t1;

    t1=GetTimeBase();
    ApplyInverse(PleftProj,PleftMss_proj); // Mss^{-1} in_s
    t_Mss+=GetTimeBase()-t1;

    t1=GetTimeBase();
    PromoteFromSubspace(PleftMss_proj,mu);  
    t_proj+=GetTimeBase()-t1;

    // New search direction p[+1] = z-mu + Sum beta_prev p[prev]
    linop->axpy(p[peri_kp],mu,z,-1.0);     // Deflation bit
    
    // peri_k == mmax-1 => peri_kp = 0.
    // Should orthogonalise 1 fields 
#if 1
    double piApk;
    double beta;
    int northog;

    //   
    northog     = (peri_kp==0)?1:peri_kp; // This is the fCG(mmax) algorithm
    //    northog     = 1;                      // This is the IPCG algorithm
    //    northog     = (k>mmax-1)?(mmax-1):k;   // This is the fCG-Tr(mmax-1) algorithm

    if  ( linop->isBoss() && !me ) printf("BfmMultiGrid::FlexCG_PolyPrec iteration %d : orthogonalising to last %d vectors\n",k,northog);

    t1=GetTimeBase();
    for(int back=0; back < northog; back++){

      int peri_back = (k-back)%mmax;

      piApk= real(linop->inner(Ap[peri_back],p[peri_kp]));

      beta = -piApk/pAp[peri_back];

      linop->axpy(p[peri_kp],p[peri_back],p[peri_kp],beta);
      linop->thread_barrier();
      if ( linop->isBoss() && !me ) printf("BfmMultiGrid::FlexCG_PolyPrec iteration %d : Beta[%d] = %le\n",k,back,beta);
      fflush(stdout);
      linop->thread_barrier();

    }
    t_ortho+=GetTimeBase() - t1;
#else 
    double beta = r_dot_z/r_dot_z_p;
    linop->axpy(p[peri_kp],p[peri_k],p[peri_kp],beta);
#endif
    // Stopping condition

    if ( linop->isBoss() && !me ) printf("BfmMultiGrid::FlexCG_PolyPrec iteration %d residual %le\n",k,sqrt(r_dot_z/ssq));
    if(1){
      double nn = linop->norm(r);
      if( linop->isBoss() && !me ) printf("resid = %le\n",nn);
      for(int s=0;s<linop->cbLs;s++){
	nn = linop->norm(r,s);
	if( linop->isBoss() && !me ) printf("resid[s] = %le\n",nn);
      }
      fflush(stdout);
      linop->thread_barrier();
    }

    if ( r_dot_z <= rsq ) { 
      
      t1=GetTimeBase();
      
      if ( linop->isBoss() && !me ) {
	printf("BfmMultiGrid::FlexCG_PolyPrec converged in %d iterations %f s\n",
					   k,1.0e-9*(t1-t_start)/1.6);
	printf("BfmMultiGrid::FlexCG_PolyPrec MssInv %f s\n",1.0e-9*(t_Mss)/1.6);
	printf("BfmMultiGrid::FlexCG_PolyPrec Proj   %f s\n",1.0e-9*(t_proj)/1.6);
	printf("BfmMultiGrid::FlexCG_PolyPrec Ortho  %f s\n",1.0e-9*(t_ortho)/1.6);
	printf("BfmMultiGrid::FlexCG_PolyPrec Prec  %f s\n",1.0e-9*(t_prec)/1.6);
	printf("BfmMultiGrid::FlexCG_PolyPrec other  %f s\n",1.0e-9*(t1-t_start-t_Mss-t_proj-t_ortho-t_prec)/1.6);
      }
	fflush(stdout);
	linop->thread_barrier();
	
	linop->Mprec(x,mp,tmp,0);
	if ( linop->isBoss() && !me ) printf("BfmMultiGrid::FlexCG_PolyPrec Mprec ");
	fflush(stdout);
	linop->thread_barrier();
	linop->Mprec(mp,mmp,tmp,1); 
	if ( linop->isBoss() && !me ) printf("BfmMultiGrid::FlexCG_PolyPrec Mprec ");
	fflush(stdout);
	linop->thread_barrier();
	linop->axpy(tmp,src,mmp,-1.0);
	if ( linop->isBoss() && !me ) printf("BfmMultiGrid::FlexCG_PolyPrec axpy ");
	fflush(stdout);
	linop->thread_barrier();
	
	double  mpnorm = sqrt(linop->norm(mp));
	double mmpnorm = sqrt(linop->norm(mmp));
	double psinorm = sqrt(linop->norm(x));
	double srcnorm = sqrt(linop->norm(src));
	double tmpnorm = sqrt(linop->norm(tmp));
	double true_residual = tmpnorm/srcnorm;
	if ( linop->isBoss() && !me ) {
	  printf("BfmMultiGrid::FlexCG_PolyPrec: true residual is %le, solution %le, source %le \n",true_residual,psinorm,srcnorm);
	  printf("BfmMultiGrid::FlexCG_PolyPrec: target residual was %le \n",linop->residual);
	  printf("BfmMultiGrid::FlexCG_PolyPrec: mp %le, mmp %le\n",mpnorm,mmpnorm);
	  fflush(stdout);
	}
	linop->thread_barrier();
	linop->threadedFreeFermion(tmp);
	linop->threadedFreeFermion(z);
	linop->threadedFreeFermion(mu);
	linop->threadedFreeFermion(mp);
	linop->threadedFreeFermion(mmp);
	linop->threadedFreeFermion(r);
	for(int m=0;m<mmax;m++){
	  linop->threadedFreeFermion(p[m]);
	  linop->threadedFreeFermion(Ap[m]);
	}
	if ( linop->isBoss() && (!me) ) { 
	  linop->InverterExit();
	}
	return k;
    }

  }
  if ( linop->isBoss() && !me ) printf("BfmMultiGrid::FlexCG_PolyPrec: CG not converged \n");
  linop->threadedFreeFermion(tmp);
  linop->threadedFreeFermion(mp);
  linop->threadedFreeFermion(mmp);
  linop->threadedFreeFermion(z);
  linop->threadedFreeFermion(mu);
  linop->threadedFreeFermion(r);
  for(int m=0;m<mmax;m++){
   linop->threadedFreeFermion(p[m]);
   linop->threadedFreeFermion(Ap[m]);
  }

  if ( linop->isBoss() && (!me) ) { 
    linop->InverterExit();
  }

  return -1;
}

void BfmMultiGrid::LdopDeflationBasisInit(int Nbasis)
{
  LdopDeflationBasisSize = Nbasis;
  int N = LdopDeflationBasisSize;
  int Ns= LocalNsubspace;
  LdopDeflationBasis.resize(N);
  LdopDeflationAv.resize(N);
  for(int v=0;v<N;v++){
    LdopDeflationBasis[v].resize(LocalNsubspace);
    LdopDeflationAv[v].resize(LocalNsubspace);
  }
  
  LdopDeflationMatrix.resize(N*N);
  LdopDeflationInverseMatrix.resize(N*N);

  // Compute vectors [Chebyshev]
  ChebyshevInit(CoarseSubspaceChebyLo,
		CoarseSubspaceChebyHi,
		CoarseSubspaceChebyPower,
		CoarseSubspaceChebyOrder);

  QDPIO::cout << "Building deflation space for little dirac op" <<endl;
  //  std::normal_distribution<double> dist;
  //  std::mt19937 eng;
  for(int v=0;v<N;v++){ 
    // Gaussian noise???
    std::vector<std::complex<double> > src(LocalNsubspace);

    for(int i=0;i<LocalNsubspace;i++){
      src[i] = std::complex<double>(drand48()-0.5,drand48()-0.5);
      src[i] = src[i]*(1.0/LocalNsubspace);
    }

    QDPIO::cout << "Vector "<< v<<endl;

#pragma omp parallel 
    {
#pragma omp for 
      for(int t=0;t<linop->nthread;t++) {
	//	PolyLdop(src,LdopDeflationBasis[v]);
	for(int i=0;i<3;i++){
	  linop->thread_barrier();
	  ApplyInverseCG(src,LdopDeflationBasis[v]);
	  double nn=norm_vec(LdopDeflationBasis[v]);
	  double scale=1.0/sqrt(nn);
	  axpby(LdopDeflationBasis[v],LdopDeflationBasis[v],LdopDeflationBasis[v],scale,0.0);
	  axpy(src,LdopDeflationBasis[v],LdopDeflationBasis[v],0.0);
	}
      }
    }
  }


  QDPIO::cout << "Orthogonalising subspace" <<endl;
  // Orthogonalise subspace
  // Remove earlier (normalised) vectors
#pragma omp parallel 
  {
#pragma omp for 
    for(int t=0;t<linop->nthread;t++) {

      std::complex<double> c;
      double n,a;

      for(int v=0;v<N;v++){
	for(int u=0;u<v;u++){
	  linop->thread_barrier();

	  //Inner product & remove component
	  c = innerProduct(LdopDeflationBasis[u],LdopDeflationBasis[v]);
	  zaxpy(LdopDeflationBasis[v],LdopDeflationBasis[u],LdopDeflationBasis[v],-c);
	}
	// Normalise this vector
	n = norm_vec(LdopDeflationBasis[v]);
	a = 1.0/sqrt(n);
	axpby(LdopDeflationBasis[v],LdopDeflationBasis[v],LdopDeflationBasis[v],a,0.0);
      }
    }
  }
  QDPIO::cout << "Computing matrix elements" <<endl;
#pragma omp parallel 
  {
#pragma omp for 
    for(int t=0;t<linop->nthread;t++) {
      // Compute matrix elements
      for(int v=0;v<N;v++){
	linop->thread_barrier();
	ApplyThread(LdopDeflationBasis[v],LdopDeflationAv[v]);
	for(int u=0;u<N;u++){
	  std::complex<double> val          = innerProduct(LdopDeflationBasis[u],LdopDeflationAv[v]);
	  LdopDeflationMatrix[N*u+v]        =val;
	  LdopDeflationInverseMatrix[N*u+v] =val;
	}
      }
    }
  }

#if 1
  QDPIO::cout << "Calling out to LAPACK for matrix inversion" <<endl;
  // Precompute inverse matrix rep on subspace
  // This could be made very efficient by LDU, and change basis as only N zaxpy's required to apply
  LapackHermitianInvert(N,(double *)&LdopDeflationInverseMatrix[0]);
#endif

  // Check the inverse
#if 1
  for(int i=0;i<N;i++){
  for(int j=0;j<N;j++){
    std::complex<double> dot=0.0;
    for(int k=0;k<N;k++){
      dot += LdopDeflationInverseMatrix[N*i+k]* LdopDeflationMatrix[k*N+j];
    }
    std::complex<double> expect;
    if(i==j) expect=1.0;
    else     expect=0.0;

    if ( abs(expect - dot) > 1.0e-10 ) { 
      printf("Oops inverse test failed\n");
      exit(0);
    }
  }}
#endif

  DeflKrylovProj.resize(N);
  DeflKrylovMss.resize(N);


}
void BfmMultiGrid::LdopDeflationProject(std::vector<std::complex<double> > &localsubspace,
					std::vector<std::complex<double> > &deflvec)
{
  if( localsubspace.size() != LocalNsubspace ) { 
    exit(0);
  }
  if( deflvec.size() != LdopDeflationBasisSize ) { 
    exit(0);
  }
  int nwork,thrlen,throff,me;

  nwork = LdopDeflationBasisSize;
  linop->thread_work(nwork,me,thrlen,throff);
  uint64_t t0=GetTimeBase();
  for(int v=throff;v<throff+thrlen;v++){
    qpx_inner((double *)&deflvec[v],LocalNsubspace,(double *)&LdopDeflationBasis[v][0],(double *)&localsubspace[0]);
  }
  linop->thread_barrier();
  uint64_t t1=GetTimeBase();
  linop->comm_gsum((double *)&deflvec[0],2*LdopDeflationBasisSize);
  //  if ( me == 0 ) {
  //    QMP_sum_double_array((double *)&deflvec[0],2*LdopDeflationBasisSize);
  //  }
  uint64_t t2=GetTimeBase();
  linop->thread_barrier();
  static int printed=0;
  if ( !printed ) { 
    if (linop->isBoss() && (!me) ) { 
      printf("DeflationProject: cycle breakdown  %ld %ld\n",t1-t0,t2-t1);
      printed=1;
    }
  }
  // for(int v=0;v<LdopDeflationBasisSize;v++){
  //   std::complex<double> cc = innerProduct(LdopDeflationBasis[v],localsubspace);
  //   if ( me == 0 ) { 
  //     printf("Projection Check: %le,%le vs %le,%le\n",
  // 	     real(deflvec[v]),imag(deflvec[v]),
  // 	     real(cc),imag(cc)
  // 	     );
  //     fflush(stdout);
  //   }
  // }
  // linop->thread_barrier();
}
void BfmMultiGrid::LdopDeflationPromote(std::vector<std::complex<double> >&deflvec,
					std::vector<std::complex<double> >&localsubspace)
{
  if( localsubspace.size() != LocalNsubspace ) { 
    exit(0);
  }
  if( deflvec.size() != LdopDeflationBasisSize ) { 
    exit(0);
  }
  int me,thrlen,throff;
#ifdef USE_XLC_OPTIMISED_CODE
  linop->thread_work(LocalNsubspace,me,thrlen,throff);
  for(int v=0;v<LdopDeflationBasisSize;v++){
    if(v==0){
      qpx_axpby(thrlen,(double *)&localsubspace[throff],(double *)&LdopDeflationBasis[v][throff],(double *)&LdopDeflationBasis[v][throff],0.0,0.0);
    }
    qpx_zaxpy(thrlen,(double *)&localsubspace[throff],(double *)&LdopDeflationBasis[v][throff],(double *)&localsubspace[throff],(double *)&deflvec[v]);
  }
  linop->thread_barrier();
#else
  zeroOut(localsubspace);
  for(int v=0;v<LdopDeflationBasisSize;v++){
    zaxpy(localsubspace,LdopDeflationBasis[v],localsubspace,deflvec[v]);
  }
#endif 
}
void BfmMultiGrid::LdopDeflationMatPromote(std::vector<std::complex<double> >&deflvec,
					   std::vector<std::complex<double> >&localsubspace)
{
  if( localsubspace.size() != LocalNsubspace ) { 
    exit(0);
  }
  if( deflvec.size() != LdopDeflationBasisSize ) { 
    exit(0);
  }
  int me,thrlen,throff;
#ifdef USE_XLC_OPTIMISED_CODE
  linop->thread_work(LocalNsubspace,me,thrlen,throff);
  for(int v=0;v<LdopDeflationBasisSize;v++){
    if(v==0){
      qpx_axpby(thrlen,(double *)&localsubspace[throff],(double *)&LdopDeflationAv[v][throff],(double *)&LdopDeflationAv[v][throff],0.0,0.0);
    }
    qpx_zaxpy(thrlen,(double *)&localsubspace[throff],(double *)&LdopDeflationAv[v][throff],(double *)&localsubspace[throff],(double *)&deflvec[v]);
  }
  linop->thread_barrier();
#else
  zeroOut(localsubspace);
  for(int v=0;v<LdopDeflationBasisSize;v++){
    zaxpy(localsubspace,LdopDeflationAv[v],localsubspace,deflvec[v]);
  }
#endif 
}


void BfmMultiGrid::LdopDeflationMatrixInverseMult(std::vector<std::complex<double> > &in,std::vector<std::complex<double> >&out)
{
  int me = linop->thread_barrier();
  if(in.size() != LdopDeflationBasisSize ) exit(0);
  if(out.size()!= LdopDeflationBasisSize ) exit(0);
  if (me==0){
#if 1
    double coeff=0.0;
    myzgemv(LdopDeflationBasisSize,(double *)&LdopDeflationInverseMatrix[0],(double *)&in[0],coeff,(double *)&out[0]);
#else
    for(int i=0;i<LdopDeflationBasisSize;i++){
      out[i] = 0.0;
      int o = i*LdopDeflationBasisSize;
      for(int j=0;j<LdopDeflationBasisSize;j++){
	out[i] += LdopDeflationInverseMatrix[o+j]*in[j];
      }
    }
#endif
  }
  linop->thread_barrier();
}
void BfmMultiGrid::LdopDeflationMatrixMult(std::vector<std::complex<double> > &in,std::vector<std::complex<double> >&out)
{
  int me = linop->thread_barrier();
  if(in.size() != LdopDeflationBasisSize ) exit(0);
  if(out.size() != LdopDeflationBasisSize ) exit(0);
  if ( me == 0 ) { 
#if 1
    double coeff=0.0;
    myzgemv(LdopDeflationBasisSize,(double *)&LdopDeflationMatrix[0],(double *)&in[0],coeff,(double *)&out[0]);
#else
    for(int i=0;i<LdopDeflationBasisSize;i++){
      out[i] = 0.0;
      int o = i*LdopDeflationBasisSize;
      for(int j=0;j<LdopDeflationBasisSize;j++){
	out[i] += LdopDeflationMatrix[o+j]*in[j];
      }
    }
#endif

  }
  linop->thread_barrier();
}

void BfmMultiGrid::PolyLdop(std::vector<std::complex<double> > &in,std::vector<std::complex<double> > &out)
{

  int N = LocalNsubspace;

  if(in.size()  != N ) exit(0);
  if(out.size() != N ) exit(0);

  std::vector<std::complex<double> > y(N);
  std::vector<std::complex<double> > Tnm(N);
  std::vector<std::complex<double> > Tn (N);
  std::vector<std::complex<double> > Tnp(N);
  
  int me = linop->thread_barrier();

  double xscale = 2.0/(hi-lo);
  double mscale = -(hi+lo)/(hi-lo);
    
  axpy(Tnm,in,in,0.0);              // Tnm=T0=in

  ApplyThread(in,y);

  axpby(Tn ,y,in,xscale,mscale);    // TN=T1 = (xscale MdagM + mscale)in 
  axpby(out,Tnm,Tn,coeffs[0]/2.0,coeffs[1]); // sum = .5c[0]T0 + c[1]T1
    
  for(int i=2;i<order;i++){
      
    ApplyThread(Tn,y);

    axpby(y,y,Tn,xscale,mscale);    // y Tn = [xscale MdagM+mscale] Tn
    axpby(Tnp,y,Tnm,2.0,-1.0);      // Tnp=2yTn - Tnm

    axpy(Tnm,Tn,Tn,0.0); 
    axpy(Tn,Tnp,Tnp,0.0);
    axpy(out,Tn,out,coeffs[i]);//Accumulate

  }
};


  /*
   * Compared to Tang-2009:  P=Pleft. P^T = PRight Q=MssInv. 
   * Script A = SolverMatrix 
   * Script P = Preconditioner
   *
   * Deflation methods considered
   *      -- Solve P A x = P b        [ like Luscher ]
   * DEF-1        M P A x = M P b     [i.e. left precon]
   * DEF-2        P^T M A x = P^T M b
   * ADEF-1       Preconditioner = M P + Q      [ Q + M + M A Q]
   * ADEF-2       Preconditioner = P^T M + Q
   * BNN          Preconditioner = P^T M P + Q
   * BNN2         Preconditioner = M P + P^TM +Q - M P A M 
   * 
   * Implement ADEF-2
   *
   * Vstart = P^Tx + Qb
   * M1 = P^TM + Q
   * M2=M3=1
   * Vout = x
   */

void BfmMultiGrid::PcgM(Fermion_t in,Fermion_t out,Fermion_t tmp)
{
  // Optional convert to single precision in this????
  int me = PcgPrecLinop->thread_barrier();

  Fermion_t pin,pout;
  if ( (linop->precision() == sizeof(double))
     &&(PcgPrecLinop->precision() == sizeof(float))
       ) {
    pin =PcgPrecLinop->threadedAllocFermion();
    pout=PcgPrecLinop->threadedAllocFermion();
    PcgPrecLinop->precisionChange(in,pin,DoubleToSingle,1);
    PcgPrecLinop->fill(out,0.0);
  } else { 
    pin=in;
    pout=out;
  }



  double restore_rsd  = PcgPrecLinop->residual;
  double restore_iter=  PcgPrecLinop->max_iter;
  PcgPrecLinop->thread_barrier();
  PcgPrecLinop->residual = InnerKrylovResidual;
  PcgPrecLinop->max_iter = InnerKrylovIterMax;
  PcgPrecLinop->thread_barrier();

  int nshift=1;
  int single=0;
  double shift = PcgShift;
  double alpha[1]     = {1.0};
  double shifts[1]    = {shift};
  double mresidual[1] = {InnerKrylovResidual};

  PcgPrecLinop->CGNE_prec_MdagM_multi_shift(&pout,
				     pin,
				     shifts,
				     alpha,
				     nshift,
				     mresidual,
				     single);
  PcgPrecLinop->thread_barrier();
  PcgPrecLinop->residual=restore_rsd;
  PcgPrecLinop->max_iter=restore_iter;
  PcgPrecLinop->thread_barrier();

  if ( (linop->precision() == sizeof(double))
     &&(PcgPrecLinop->precision() == sizeof(float))
       ) {
    PcgPrecLinop->precisionChange(pou,out,SingleToDouble,1);
    PcgPrecLinop->threadedFreeFermion(pin);
    PcgPrecLinop->threadedFreeFermion(pout);
  }

  return;
}

void BfmMultiGrid::PcgM1(Fermion_t in, Fermion_t out,Fermion_t tmp,Fermion_t mp)
{
  uint64_t t1;
  int me = linop->thread_barrier();

  Fermion_t Mtmp= linop->threadedAllocFermion();
  Fermion_t Min = linop->threadedAllocFermion();

  switch(PcgType) { 

  case PcgPrec:
  case PcgDef1:
  case PcgDef2:
    //M
    PcgM(in,out,tmp);

    break;
  case PcgMssDef:

    ProjectToSubspace(in,PleftProj);     
    ApplyInverse(PleftProj,PleftMss_proj); // Ass^{-1} r_s
    PromoteFromSubspace(PleftMss_proj,tmp);  

    linop->axpy(out,in,tmp,1.0);
    break;

  case PcgAD:
    // M+Q
    PcgM(in,out,tmp);

    ProjectToSubspace(in,PleftProj);     
    ApplyInverse(PleftProj,PleftMss_proj); // Ass^{-1} r_s
    PromoteFromSubspace(PleftMss_proj,tmp);  


    linop->axpy(out,out,tmp,1.0);
    break;

  case PcgADef1:
    //MP  + Q
    Pleft(in,mp);
    PcgM(mp,out,tmp);
    
    ProjectToSubspace(in,PleftProj);     
    ApplyInverse(PleftProj,PleftMss_proj); // Ass^{-1} r_s
    PromoteFromSubspace(PleftMss_proj,tmp);  

    linop->axpy(out,out,tmp,1.0);
    break;
  case PcgAdef2:

    // [PTM+Q] in = [1 - Q A] M in + Q in = Min + Q [ in -A Min]

    PcgM(in,Min,tmp);
    linop->Mprec(Min,tmp,Mtmp,DaggerNo);
    linop->Mprec(tmp,out,Mtmp,DaggerYes);  // out  = A Min
    linop->axpy(tmp,out,in,-1.0);          // tmp  = in - A Min

    ProjectToSubspace(tmp,PleftProj);     
    ApplyInverse(PleftProj,PleftMss_proj); // Ass^{-1} [in - A Min]_s
    PromoteFromSubspace(PleftMss_proj,tmp);// tmp = Q[in - A Min]  

    linop->axpy(out,Min,tmp,1.0); // Min+tmp
    

    // PT M + Q
    /*
    PcgM(in,tmp,out);
    Pright(tmp,out);

    ProjectToSubspace(in,PleftProj);     
    ApplyInverse(PleftProj,PleftMss_proj); // Ass^{-1} r_s
    PromoteFromSubspace(PleftMss_proj,tmp);  

    linop->axpy(out,out,tmp,1.0);
    */
    break;

  case PcgBNN:
    // PT M P +Q 
    Pleft(in,out);
    PcgM(out,mp,tmp);
    Pright(mp,out);
    
    ProjectToSubspace(in,PleftProj);     
    ApplyInverse(PleftProj,PleftMss_proj); // Ass^{-1} r_s
    PromoteFromSubspace(PleftMss_proj,tmp);  

    linop->axpy(out,out,tmp,1.0);
    break;
  default :
    exit(0);
  }
  linop->threadedFreeFermion(Mtmp);
  linop->threadedFreeFermion(Min);
}
void BfmMultiGrid::PcgM2(Fermion_t in, Fermion_t out)
{
  switch(PcgType) { 
  case PcgPrec:
  case PcgAD:
  case PcgDef1:
  case PcgADef1:
  case PcgAdef2:
  case PcgBNN:
  case PcgMssDef:
    linop->axpy(out,in,in,0.0);
    break;
  case PcgDef2:
    Pright(in,out);
    break;
  default :
    exit(0);
  }
}

double BfmMultiGrid::PcgM3(Fermion_t p, Fermion_t mp,Fermion_t mmp, Fermion_t tmp)
{
  double d;
  uint64_t t1;
 
  t1=GetTimeBase();
  int me = linop->thread_barrier();

  switch(PcgType) { 
  case PcgPrec:
  case PcgAD:
  case PcgDef2:
  case PcgADef1:
  case PcgAdef2:
  case PcgBNN:
  case PcgMssDef:
    d=linop->Mprec(p,mp,tmp,0,1);// Dag no
      linop->Mprec(mp,mmp,tmp,1);// Dag yes
    break;
  case PcgDef1:
    d=linop->Mprec(p,mmp,tmp,0,1);// Dag no
      linop->Mprec(mmp,mp,tmp,1);// Dag yes
    Pleft(mp,mmp);
    d=real(linop->inner(p,mmp));
    break;
  default :
    exit(0);
  }
  return d;
}

void BfmMultiGrid::PcgVstart(Fermion_t x, Fermion_t src,Fermion_t r,
			     Fermion_t mp, Fermion_t mmp, Fermion_t tmp)
{
  uint64_t t1;
  int me = linop->thread_barrier();
  switch(PcgType) { 
  case PcgPrec:
  case PcgMssDef:
  case PcgAD:
  case PcgDef1:
  case PcgADef1:
  case PcgBNN:

    break;

  case PcgDef2:
  case PcgAdef2:

    ///////////////////////////////////
    // Choose x_0 such that 
    // x_0 = guess +  (A_ss^inv) r_s = guess + Ass_inv [src -Aguess]
    //                               = [1 - Ass_inv A] Guess + Assinv src
    //                               = P^T guess + Assinv src 
    //                               = Vstart  [Tang notation]
    // This gives:
    // W^T (src - A x_0) = src_s - A guess_s - r_s
    //                   = src_s - (A guess)_s - src_s  + (A guess)_s 
    //                   = 0 
    ///////////////////////////////////
    linop->Mprec(x,mp,tmp,DaggerNo);
    linop->Mprec(mp,mmp,tmp,DaggerYes);

    linop->axpy (r, mmp, src,-1.0);        // r_{-1} = src - A x

    ProjectToSubspace(r,PleftProj);     
    ApplyInverse(PleftProj,PleftMss_proj); // Ass^{-1} r_s
    PromoteFromSubspace(PleftMss_proj,mp);  

    linop->axpy(x,x,mp,1.0);
    break;
  default :
    exit(0);
  }
}

void BfmMultiGrid::PcgVout  (Fermion_t in, Fermion_t out,Fermion_t src,Fermion_t tmp)
{
  switch(PcgType) { 
  case PcgDef1:
    //Qb + PT x
    ProjectToSubspace(src,PleftProj);     
    ApplyInverse(PleftProj,PleftMss_proj); // Ass^{-1} r_s
    PromoteFromSubspace(PleftMss_proj,tmp);  
    
    Pright(in,out);
    
    linop->axpy(out,tmp,out,1.0);
    break;

  case PcgPrec:
  case PcgAD:
  case PcgMssDef:
  case PcgDef2:
  case PcgADef1:
  case PcgAdef2:
  case PcgBNN:
    linop->axpy(out,in,in,0.0);
    break;
  default :
    exit(0);
  }
}
void BfmMultiGrid::PcgReport (double MssTolOuter,double MssTolInner)
{
  int me = linop->thread_barrier();

  // Timing in cycles
  PcgMssInv=0; 
  PcgMprec=0;
  PcgProj=0;
  PcgProm=0;


  if ( linop->isBoss() && (!me) ) { 
    printf("BfmMultiGrid::Pcg ******************** PARAMS ***************\n");
    switch(PcgType){
    case PcgDef1:
      printf("BfmMultiGrid::Pcg algorithm DEF1\n");
      break;
    case PcgPrec:
      printf("BfmMultiGrid::Pcg algorithm Preconditioned CG\n");
      break;
    case PcgAD:
      printf("BfmMultiGrid::Pcg algorithm AD\n");
      break;
    case PcgDef2:
      printf("BfmMultiGrid::Pcg algorithm DEF2\n");
      break;
    case PcgADef1:
      printf("BfmMultiGrid::Pcg algorithm A-DEF1\n");
      break;
    case PcgAdef2:
      printf("BfmMultiGrid::Pcg algorithm A-DEF2\n");
      break;
    case PcgBNN:
      printf("BfmMultiGrid::Pcg algorithm BNN\n");
      break;
    case PcgMssDef:
      printf("BfmMultiGrid::Pcg algorithm MssDef [equiv AD with M=1]\n");
      break;
    default:
      exit(0);
    }
    printf("BfmMultiGrid::Pcg Ldop inner solver residual %le \n",MssTolInner);
    printf("BfmMultiGrid::Pcg Ldop outer solver residual %le \n",MssTolOuter);
    printf("BfmMultiGrid::Pcg InnerKrylov residual %le \n",InnerKrylovResidual);
    printf("BfmMultiGrid::Pcg InnerKrylov iters    %d \n",InnerKrylovIterMax);
    printf("BfmMultiGrid::Pcg Preconditioner shift %le \n",PcgShift);

  }
}

int BfmMultiGrid::Pcg(Fermion_t psi, Fermion_t src,Fermion_t resid,double MssTolOuter,double MssTolInner) 
{
  double f;
  double rtzp,rtz,a,d,b;

  int me = linop->thread_barrier();
  PcgReport(MssTolOuter,MssTolInner);

  if ( linop->isBoss() && (!me) ) { 
    linop->InverterEnter();
  }


  Fermion_t x   = psi;
  Fermion_t p   = linop->threadedAllocFermion(); 
  Fermion_t z   = linop->threadedAllocFermion(); 

  Fermion_t tmp = linop->threadedAllocFermion(); 
  Fermion_t mp  = linop->threadedAllocFermion(); 
  Fermion_t mmp = linop->threadedAllocFermion(); 
  Fermion_t r   = linop->threadedAllocFermion(); 
  Fermion_t mu   = linop->threadedAllocFermion(); 

  //Initial residual computation & set up
  LittleDopSolverResidual=MssTolOuter;

  double guess = linop->norm(psi);

  //////////////////////////
  // x0 = Vstart -- possibly modify guess
  //////////////////////////
  PcgVstart(x,src,r,mp,mmp,tmp);

  // r0 = b -A x0
  linop->Mprec(x,mp,tmp,DaggerNo);
  linop->Mprec(mp,mmp,tmp,DaggerYes);
  linop->axpy (r, mmp, src,-1.0);    // Recomputes r=src-x0

  ProjectToSubspace(r,PleftProj);     
  double nv=norm_vec(PleftProj);
  if(linop->isBoss()&&!me){ 
    printf("BfmMultiGrid::Pcg subspace residual for x0 = %le\n",nv);
  }

  //////////////////////////////////
  // Compute z = M1 x
  //////////////////////////////////
  PcgM1(r,z,tmp,mp);
  rtzp =real(linop->inner(r,z));

  ///////////////////////////////////////
  // Solve for Mss mu = P A z and set p = z-mu
  // Def2: p = 1 - Q Az = Pright z 
  // Other algos M2 is trivial
  ///////////////////////////////////////
  PcgM2(z,p);

  double ssq =  linop->norm(src);
  double rsq =  linop->residual* linop->residual*ssq;


  if ( linop->isBoss() && !me ) 
    printf("BfmMultiGrid::Pcg k=0 residual %le rsq %le\n",rtzp,rsq);

  uint64_t t_start=GetTimeBase();
  uint64_t t1;
  for (int k=1;k<=linop->max_iter;k++){
    
    LittleDopSolverResidual=MssTolInner;

    rtz=rtzp;
    d= PcgM3(p,mp,mmp,tmp);
    a = rtz/d;
    
    linop->axpy(x,p,x,a);
    linop->axpy(r,mmp,r,-a);

    // Compute z = M x
    PcgM1(r,z,tmp,mp);

    rtzp =real(linop->inner(r,z));
    b = rtzp/rtz;

    PcgM2(z,mu);
    linop->axpy(p,p,mu,b);

    {
      double rr=sqrt(rtzp/ssq);
      if (linop->isBoss() && !me ){
	printf("bfmbase::Pcg: k= %d residual = %le \n",k,rr);
      }
      linop->InverterLogIteration(k, rr,src,x,mp,mmp,tmp);
    }

    // Stopping condition
    if ( rtzp <= rsq ) { 
	
      LittleDopSolverResidual=MssTolOuter;

      t1=GetTimeBase();
      if ( linop->isBoss() && !me ) {
	printf("BfmMultiGrid::Pcg converged in %d iterations %f s\n",
					   k,1.0e-9*(t1-t_start)/1.6);
	printf("BfmMultiGrid::Pcg Mprec  %f s\n",1.0e-9*(PcgMprec)/1.6);
	printf("BfmMultiGrid::Pcg MssInv %f s\n",1.0e-9*(PcgMssInv)/1.6);
	printf("BfmMultiGrid::Pcg Proj   %f s\n",1.0e-9*(PcgProj)/1.6);
	printf("BfmMultiGrid::Pcg Prom   %f s\n",1.0e-9*(PcgProm)/1.6);
	printf("BfmMultiGrid::Pcg other  %f s\n",1.0e-9*(t1-PcgMprec-PcgMssInv-PcgProj-PcgProm-t_start)/1.6);
      }
	
	linop->Mprec(x,mp,tmp,0);
	linop->Mprec(mp,mmp,tmp,1); 
	linop->axpy(tmp,src,mmp,-1.0);
	
	double  mpnorm = sqrt(linop->norm(mp));
	double mmpnorm = sqrt(linop->norm(mmp));
	double psinorm = sqrt(linop->norm(x));
	double srcnorm = sqrt(linop->norm(src));
	double tmpnorm = sqrt(linop->norm(tmp));
	double true_residual = tmpnorm/srcnorm;
	if ( linop->isBoss() && !me ) {
	  printf("BfmMultiGrid::Pcg: true residual is %le, solution %le, source %le \n",true_residual,psinorm,srcnorm);
	  printf("BfmMultiGrid::Pcg: target residual was %le \n",linop->residual);
	}
	
	linop->threadedFreeFermion(tmp);
	linop->threadedFreeFermion(p);
	linop->threadedFreeFermion(z);
	linop->threadedFreeFermion(mu);
	linop->threadedFreeFermion(mp);
	linop->threadedFreeFermion(mmp);
	linop->threadedFreeFermion(r);
	if ( linop->isBoss() && (!me) ) { 
	  linop->InverterExit();
	}
	return k;
    }

  }
  if ( linop->isBoss() && !me ) printf("BfmMultiGrid::Pcg: CG not converged \n");
  linop->threadedFreeFermion(tmp);
  linop->threadedFreeFermion(p);
  linop->threadedFreeFermion(mp);
  linop->threadedFreeFermion(mmp);
  linop->threadedFreeFermion(z);
  linop->threadedFreeFermion(mu);
  linop->threadedFreeFermion(r);
  if ( linop->isBoss() && (!me) ) { 
    linop->InverterExit();
  }

  return -1;
}
