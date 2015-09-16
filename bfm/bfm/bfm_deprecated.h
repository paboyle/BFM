
//Preserving useful bits of code which I don't want to maintain
//Just in case I need to resurrect them
#ifdef BFM_DEPRECATED
#ifdef BFM_DEPRECATED
  void Apply(std::vector<std::complex<cFloat> > &,std::vector<std::complex<cFloat> > & );
  void GetStencilData(int _idx,int _ballidx,std::vector<std::complex<double> >& nbr_data,std::vector<std::complex<double> > &my_data,int rev);

  /////////////////////////////////
  // Deflated solver controls
  /////////////////////////////////
  int KrylovPrecondition;
  int LeftPrecondition;
  int RightPrecondition;
  int order;
  double alpha;
  double hi,lo;
  double *coeffs;

  void PickBlock   (Fermion_t was, Fermion_t is,int b);
  void PickQuadrant(Fermion_t was, Fermion_t is,int q);
  void ComputeLittleMatrix (void);
  template<class Float> void PolyMdagMprec(bfm_qdp<Float>    *lop,Fermion_t in,Fermion_t out);
  ////////////// Benchmarking support /////////////
  void * BenchmarkL2Buffer;
  void BenchmarkL2Prepare(void);
  void BenchmarkL2(void);
  void BenchmarkMultiply(int Numtest);

  void HaloTest(std::vector<std::complex<cFloat> >&my_data);

  void PreconditionWithLinop(bfm_t *pop) { preclinop = pop; };

  void ChebyshevEnd(void) { delete[] coeffs; }
  void ChebyshevInit(double _lo,double _hi,double _alpha,int _order);

  void OneMinusPleft (Fermion_t in,Fermion_t out);
  void OneMinusPright(Fermion_t in,Fermion_t out);
  void PleftM (Fermion_t in,Fermion_t out);
  void MPright(Fermion_t in,Fermion_t out);

  void PolyMdagMprecRight(Fermion_t in,Fermion_t out);
  void PolyMdagMprecLeft (Fermion_t in,Fermion_t out);
  void PolyMdagMprec(Fermion_t in,Fermion_t out);
  void PolyLdop(std::vector<std::complex<double> > & in,std::vector<std::complex<double> > & out);
  void SolverMatrix(Fermion_t in,Fermion_t out);
  void DeflatedSolverMdagMpc(Fermion_t in,Fermion_t out);
  void UndeflatedSolverMdagMpc(Fermion_t in,Fermion_t out);
  void DeflatedSolverMpc(Fermion_t in,Fermion_t out);
  void UndeflatedSolverMpc(Fermion_t in,Fermion_t out);
  int  fCG_PolyPrec(Fermion_t psi, Fermion_t src_uprec);
  int  rCG_PolyPrec(Fermion_t psi, Fermion_t src_uprec,Fermion_t resid);
  int  CG_PolyPrec(Fermion_t psi, Fermion_t src_uprec,Fermion_t resid);
  int  MCR_PolyPrec(Fermion_t psi, Fermion_t src_uprec,Fermion_t true_res_vec);
  int  GCR_PolyPrec(Fermion_t psi, Fermion_t src_uprec,int m);
#endif

#ifdef BFM_DEPRECATED
void mem_benchmark_prepare(void);
void mem_benchmark_end(void);
void mem_benchmark(void);
#endif

static void * mem_benchmark_buffer;
void mem_benchmark_prepare(void)
{
  int bsize = 64*1024*1024;
  mem_benchmark_buffer = memalign(bsize,bsize);
}
void mem_benchmark_end(void)
{
  free(mem_benchmark_buffer);
}
void mem_benchmark(void)
{
  int bsize = 64*1024*1024;
  uint64_t base_addr = (uint64_t)mem_benchmark_buffer;
  double * dp=(double *)base_addr;
  uint64_t t1,t2;

  int threads = omp_get_max_threads();

  printf("OpenMP using %d threads and unlocked 64MB\n",threads);

  for( int mb=4;mb<64;mb+=4){
    
      for(int i=0;i<3;i++){
	t1=GetTimeBase();
	
#pragma omp parallel
	{
	  double d=0;
#pragma omp for
	  for(int o=0;o<mb*1024*1024;o+=32){
	    d+=*( (volatile double *) (base_addr+o)); // Load one double from each L1 line @ 32byte granularity
	  }
	}
	t2=GetTimeBase();
      }
      printf("%d MB  %ld cycles %le MB/s\n",mb,t2-t1,mb*1600.*1024.*1024./(t2-t1));
  }
  
}

template<class Float>
void BfmHDCG<cFloat>::RelaxSubspace(bfm_qdp<Float> *rop)
{
  int Ls=N5;
  multi1d<LatticeFermion> gauss(Ls);

  Fermion_t sol = rop->allocFermion();
  Fermion_t src = rop->allocFermion();
  subspace_d   = (Fermion_t *) malloc(2*Nvec*sizeof(Fermion_t));

  if ( rop->SPIcomms() ) rop->comm_init();

  if ( SubspaceRelaxChebyshev ) {
    QDPIO::cout << "Chebyshev: Building subspace using chebyshev filter order "<<SubspaceChebyOrder << " x^"
		<< SubspaceChebyPower<<"  ["<<SubspaceChebyLo<<","<<SubspaceChebyHi<<"]" <<endl;
  } else {
    QDPIO::cout << "Rational: Building subspace using rational  filter 1/(x+lambda)^2 for lambda= "
		<<SubspaceRationalLo << " tol " << SubspaceRationalResidual<<endl;
  }

#define VSKIP (1)

  for(int v=0;v<Nvec;v+=VSKIP){

    for(int s=0;s<Ls;s++) gauss[s]=zero;
    gaussian(gauss[0]);
    gaussian(gauss[Ls-1]);
    gauss[0]   = chiralProjectPlus(gauss[0]);
    gauss[Ls-1]= chiralProjectMinus(gauss[Ls-1]);

    rop->importFermion(gauss,src,1);

    int Nstep=1;
    for(int i=0;i<Nstep;i++){
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
	rop->axpy(src,sol,sol,0.0);
      }

    }
    }

    for(int i=0;i<VSKIP;i++){

      for(int s=0;s<Ls;s++) gauss[s]=zero;

      int vi=v+i;
      subspace_d[vi] = linop_d->allocFermion();
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
	int number_full_depth = 20;
	if ( (SubspaceSurfaceDepth > Ls/2) || (vi > Nvec-number_full_depth) ) { 
	  s_min[vi]=0; s_max[vi]=Ls-1;
	} else { 
	  s_min[vi]=Ls-SubspaceSurfaceDepth; s_max[vi]=SubspaceSurfaceDepth-1;
	  for(int s=s_max[vi]+1;s<s_min[vi];s++) gauss[s]=zero;
	}

      }
      for(int s=0;s<Ls;s++) { 
	QDPIO::cout << "Vec["<<vi<<"][s=" <<s<<"] == "<<norm2(gauss[s])<<endl;
      }
      linop_d->importFermion(gauss,subspace_d[vi],1);
    }

  }

  rop->freeFermion(src);
  rop->freeFermion(sol);

  if ( rop->isBoss() ) {
    printf("Relax Subspace Orthogonalising\n"); fflush(stdout);
  }

  OrthogonaliseSubspace();

  if ( linop_d->SPIcomms() ) linop_d->comm_init();
  
}

/*
void BfmHDCG<cFloat>::AddToSubspace(Fermion_t vec)
{	

  int Ls=N5;
  int vi=Nvec-1;// Replace last vec

#pragma omp parallel 
    {
#pragma omp for 
      for(int i=0;i<linop_d->threads;i++) {
	linop_d->axpy(subspace_d[vi],vec,vec,0.0);
      }
    }

    s_min[vi]=0; s_max[vi]=Ls-1;
    OrthogonaliseSubspace();
    ComputeLittleMatrixColored();  // Expensive recompute matrix
}
*/

template<class Float> 
void BfmHDCG<cFloat>::PolyMdagMprec(bfm_qdp<Float> *lop,Fermion_t in,Fermion_t out,  double *coeffs )
{
    Fermion_t Mtmp= lop->threadedAllocFermion();
    Fermion_t Tnm = lop->threadedAllocFermion();
    Fermion_t Tn  = lop->threadedAllocFermion();
    Fermion_t Tnp = lop->threadedAllocFermion();
    Fermion_t y   = lop->threadedAllocFermion();
    Fermion_t tmp=Tnp;

    int me = lop->thread_barrier();

    double xscale = 2.0/(hi-lo);
    double mscale = -(hi+lo)/(hi-lo);
    
    lop->axpy(Tnm,in,in,0.0);              // Tnm=T0=in

    lop->Mprec(in ,tmp,Mtmp,DaggerNo);  
    lop->Mprec(tmp,y  ,Mtmp,DaggerYes);    
    lop->axpby(Tn ,y,in,xscale,mscale);    // TN=T1 = (xscale MdagM + mscale)in 

    lop->axpby(out,Tnm,Tn,coeffs[0]/2.0,coeffs[1]); // sum = .5c[0]T0 + c[1]T1
    Fermion_t swizzle;


    // Optimisation :
    // y can overlap Tnp. Save cache footprint
    // MdagM + 3 axpy's reduces overhead
    // This should be more efficient than M_IRS CG
    for(int i=2;i<order;i++){

      lop->Mprec(Tn,tmp,Mtmp,DaggerNo);  
      lop->Mprec(tmp, y,Mtmp,DaggerYes);    

      lop->axpby(y,y,Tn,xscale,mscale);    // y Tn = [xscale MdagM+mscale] Tn
      lop->axpby(Tnp,y,Tnm,2.0,-1.0);      // Tnp=2yTn - Tnm
                                           // = 2 xscale MdagM Tn + 2 mscale Tn - Tnm
      lop->axpy(out,Tnp,out,coeffs[i]);    //Accumulate
      
      // Redundant, can pointer swizzle
      swizzle=Tnm;
      Tnm=Tn;
      Tn=Tnp;
      Tnp=Tnm;
    }
    lop->threadedFreeFermion(y);
    lop->threadedFreeFermion(Mtmp);
    lop->threadedFreeFermion(Tnm);
    lop->threadedFreeFermion(Tn);
    lop->threadedFreeFermion(Tnp);
};


void BfmHDCG<cFloat>::SolverMatrix(Fermion_t in,Fermion_t out)
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
void BfmHDCG<cFloat>::HaloTest(std::vector<std::complex<double> >&my_data)
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
void BfmHDCG<cFloat>::GetStencilData(int _idx,int _ballidx,
				  std::vector<std::complex<cFloat> >& nbr_data,
				  std::vector<std::complex<cFloat> >&my_data,int rev)
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
void BfmHDCG<cFloat>::Apply(std::vector<std::complex<double> > &in,std::vector<std::complex<double> > &out )
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
      
      myzgemv(Nvec,(double *)&Asparse[b*Nball*Nvec*Nvec+mu*Nvec],(double *)&nbr_data[0],coeff,(double *)&out[b*Nvec]);


    }
  }
}
void BfmHDCG<cFloat>::PickBlock (Fermion_t was,Fermion_t is,int b)
{
#pragma omp parallel 
    {
#pragma omp for 
      for(int i=0;i<linop_d->threads;i++) {
	linop_d->fill(is,0.0);
	linop_d->block_pick (is,was,is,b,&block_id[1][0]);
      }
    }
    exit(0);
}

void BfmHDCG<cFloat>::PickQuadrant(Fermion_t was, Fermion_t is,int q)
{
  
  double nn;
#pragma omp parallel 
    {
#pragma omp for 
      for(int i=0;i<linop_d->threads;i++) {
	linop_d->fill(is,0.0);
	linop_d->block_pick (is,was,is,q,&quadrant_id[1][0]);
	nn=linop_d->norm(is);
      }
    }
    if(linop_d->isBoss()) { 
      printf("Quadrant %d picked norm is %le\n",q,nn);
    }
}
void BfmHDCG<cFloat>::BenchmarkL2Prepare(void)
{
}
void BfmHDCG<cFloat>::BenchmarkL2(void)
{
}

void BfmHDCG<cFloat>::BenchmarkMultiply(int NumTest)
{
  double *Ap= (double *)&Asparse[0];
  double *oo= (double *)&Krylov_Atmp[0];
  double *ii= (double *)&Krylov_p[0];

  double flops = Nvec*Nvec*8*NumTest;
  uint64_t t1 = GetTimeBase();
  for (int i=0;i<NumTest;i++){
    myzgemv(Nvec,Ap,ii,0.0,oo);
  }
  uint64_t t2 = GetTimeBase();
  if ( linop_d->isBoss() ) {
    printf("Bench MatMul [%d]: t2-t1 cycles %ld, %f Mflop/s %f us\n",NumTest,t2-t1,flops*1600.0/(t2-t1),(t2-t1)/1600.0);
  }
}
void BfmHDCG<cFloat>::ComputeLittleMatrix (void)
{
  exit(0);
  Fermion_t phi_t = linop_d->allocFermion();
  Fermion_t tmp_t = linop_d->allocFermion();
  Fermion_t mmp   = linop_d->allocFermion();
  Fermion_t mp    = linop_d->allocFermion();
	
  for(int i=0;i<Nvec;i++){
    for(int b=0;b<Nblock;b++){            // Loop over global blocks

      // Apply DdagD to phi_i[b] -> mmp
#pragma omp parallel 
	{
#pragma omp for 
	  for(int t=0;t<linop_d->nthread;t++) {
	    PickBlock(subspace_d[i],phi_t,b);
	    linop_d->Mprec(phi_t,mp,tmp_t,0);
	    linop_d->Mprec(mp,mmp,tmp_t,1); 
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
	      Asparse[xblock*Nball*Nvec*Nvec+n*Nvec*Nvec+Nvec*j+i] = proj[xblock*Nvec+j];
	    }	    
	  }

	}
    }
  }
  linop_d->freeFermion( phi_t ); 
  linop_d->freeFermion( tmp_t ); 
  linop_d->freeFermion( mmp   ); 
  linop_d->freeFermion( mp    );
}

void BfmHDCG<cFloat>::ChebyshevInit(double _lo,double _hi,double _alpha,int _order)
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

void BfmHDCG<cFloat>::OneMinusPleft (Fermion_t in,Fermion_t out)
{
  Pleft(in,out);
  linop_d->axpy(out,out,in,-1.0);
}
void BfmHDCG<cFloat>::OneMinusPright(Fermion_t in,Fermion_t out)
{
  Pright(in,out);
  linop_d->axpy(out,out,in,-1.0);
}
void BfmHDCG<cFloat>::MPright (Fermion_t in,Fermion_t out)
{
  int me = linop_d->thread_barrier();
  Fermion_t tmp  = linop_d->threadedAllocFermion();
  Fermion_t Mtmp  = linop_d->threadedAllocFermion();

  Pright(in,out);
  linop_d->Mprec(out,tmp,Mtmp,0);
  linop_d->Mprec(tmp,out,Mtmp,1); 

  linop_d->threadedFreeFermion(Mtmp);
  linop_d->threadedFreeFermion(tmp);
}
void BfmHDCG<cFloat>::PleftM (Fermion_t in,Fermion_t out)
{
  int me = linop_d->thread_barrier();
  Fermion_t tmp  = linop_d->threadedAllocFermion();
  Fermion_t Mtmp  = linop_d->threadedAllocFermion();

  uint64_t t1=GetTimeBase();
  linop_d->Mprec(in,out,Mtmp,0);
  linop_d->Mprec(out,tmp,Mtmp,1); 
  uint64_t t2=GetTimeBase();
  Pleft(tmp,out);
  uint64_t t3=GetTimeBase();

  if( linop_d->isBoss() && !me ) {
    printf("MdagM/Pleft cyc = %ld/%ld\n",t2-t1,t3-t2);
  }

  linop_d->threadedFreeFermion(Mtmp);
  linop_d->threadedFreeFermion(tmp);
}


int BfmHDCG<cFloat>::MCR_PolyPrec(Fermion_t psi, Fermion_t src_uprec,Fermion_t true_res_vec) // Internal
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
  int me = linop_d->thread_barrier();

  if ( linop_d->isBoss() && (!me) ) { 
    FILE * fp = fopen("mcr.log","w"); // truncate log file
    fclose(fp);
    linop_d->InverterEnter();
  }

  Fermion_t src     = linop_d->threadedAllocFermion();  //(Pleft)src
  Fermion_t src_left= linop_d->threadedAllocFermion();  //(Pleft)src
  Fermion_t p   = linop_d->threadedAllocFermion(); 
  Fermion_t tmp = linop_d->threadedAllocFermion(); 
  Fermion_t mp  = linop_d->threadedAllocFermion(); 
  Fermion_t mmp = linop_d->threadedAllocFermion(); 
  Fermion_t r   = linop_d->threadedAllocFermion(); 
  Fermion_t Ap  = linop_d->threadedAllocFermion(); 
  Fermion_t Ar  = linop_d->threadedAllocFermion(); 

  //Initial residual computation & set up
  double guess = linop_d->norm(psi);
  double usrc  = linop_d->norm(src_uprec);
  if ( linop_d->isBoss() && !me ) {
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
  linop_d->axpy (r, Ar, src,-1.0);
  linop_d->axpy (p, Ar, src,-1.0);

  a =linop_d->norm(p);
  cp=linop_d->norm(r);

  double ssq =  linop_d->norm(src);
  double residual=linop_d->residual;
  double rsq =  residual* residual*ssq;

  if ( linop_d->isBoss() && !me ) {
      printf("bfmbase::MCR_PolyPrec guess %le \n",guess);
      printf("bfmbase::MCR_PolyPrec ssq %le rsq %le\n",ssq,rsq);
      printf("bfmbase::MCR_PolyPrec a %le cp %le\n",a,cp);
      fflush(stdout);
  }

  //Check if guess is really REALLY good :)
  if ( cp <= rsq ) {
    linop_d->threadedFreeFermion(src);
    linop_d->threadedFreeFermion(tmp);
    linop_d->threadedFreeFermion(p);
    linop_d->threadedFreeFermion(mp);
    linop_d->threadedFreeFermion(r);
    linop_d->threadedFreeFermion(Ap);
    linop_d->threadedFreeFermion(Ar);
    if ( linop_d->isBoss() && (!me) ) { 
      linop_d->InverterExit();
    }
    return 0;
  }

  if ( linop_d->isBoss() && !me ) {
      printf("bfmbase::MCR_PolyPrec k=0 residual %le rsq %le\n",cp,rsq);
      fflush(stdout);
  }
  //  c = Mprec(p,mp,tmp,DaggerNo);
  //  d = Mprec(mp,Ap,tmp,DaggerYes);  //Ap

  SolverMatrix(p,Ap);

  c = linop_d->inner_real(p,Ap);
  d = linop_d->norm(Ap);
  linop_d->axpy(Ar,Ap,Ap,0);	       //Ar
  double cp0 = c;

  int max_iter = linop_d->max_iter;
  //c = real( dot(r,Ar) );
  for (int k=1;k<=max_iter;k++){
    linop_d->iter=k;

    //c = real( dot(r,Ar) );		//c = rAr
    a = c/d;

    linop_d->axpy(psi,p,psi,a);		//x = x + ap
    cp = linop_d->axpy_norm(r,Ap,r,-a);		//r = r - aAp

  
    if ( k%1 == 0 ){
      if ( linop_d->isBoss() && !me ) {
	printf("bfmbase::MCR_PolyPrec: k= %d r^2= %le %le\n",k,cp,sqrt(cp/ssq));
	FILE * fp = fopen("mcr.log","a");
	fprintf(fp,"k %d resid %le\n",k,sqrt(cp/ssq));
	fclose(fp);
      }
    }

    linop_d->axpy(true_res_vec,r,r,0.0);
    double rr = sqrt(cp/ssq);
    linop_d->InverterLogIteration(k, rr,src,psi,mp,mmp,tmp);


    // Stopping condition
    if ( cp <= rsq ) { 

	if ( linop_d->isBoss() && (!me) ) { 
	  linop_d->InverterExit();
	}
	
	if ( linop_d->isBoss() && !me ) printf("bfmbase::MCR_PolyPrec converged in %d iterations\n",k);

	PolyMdagMprecRight(psi,tmp); // Psi is now solution of PL M psi = PL eta
	linop_d->axpy(psi,tmp,tmp,0.0);

	double true_residual;

	if ( SolverControl == SolverDeflatedMdagMpc ) {
	  double nn;

	  // Test that MssInv works
	  ProjectToSubspace(psi,PleftProj);     
	  PromoteFromSubspace(PleftProj,tmp);  
	  
	  linop_d->Mprec(tmp,mp,r,0);
	  linop_d->Mprec(mp,Ap ,r,1); 
	  ProjectToSubspace(Ap,PleftProj);     
	  PromoteFromSubspace(PleftProj,Ap);  

	  ApplyInverse     (PleftProj,PleftMss_proj); 
	  PromoteFromSubspace(PleftMss_proj,mp);  


	  nn=linop_d->norm(tmp);
	  if ( linop_d->isBoss() && !me ) 
	    printf("bfmbase::MCR_PolyPrec: subspace norm is %le \n",nn);

	  nn=linop_d->norm(Ap);
	  if ( linop_d->isBoss() && !me ) 
	    printf("bfmbase::MCR_PolyPrec: Mss norm is %le \n",nn);


	  nn=linop_d->norm(mp);
	  if ( linop_d->isBoss() && !me ) 
	    printf("bfmbase::MCR_PolyPrec: MssInv norm is %le \n",nn);


	  linop_d->axpy(mp,mp,tmp,-1.0);
	  nn=linop_d->norm(mp);
	  if ( linop_d->isBoss() && !me ) 
	    printf("bfmbase::MCR_PolyPrec: MssInv error is %le \n",nn);

	  // Compute true residual of PleftM system
	  PleftM(psi,tmp);
	  linop_d->axpy(tmp,tmp,src_left,-1.0);
	  true_residual = sqrt(linop_d->norm(tmp)/linop_d->norm(src_left));
	  if ( linop_d->isBoss() && !me ) 
	    printf("bfmbase::MCR_PolyPrec: true PleftM residual is %le \n",true_residual);

	  // This should equal MPright
	  MPright(psi,tmp);
	  linop_d->axpy(tmp,tmp,src_left,-1.0);
	  true_residual = sqrt(linop_d->norm(tmp)/linop_d->norm(src_left));
	  if ( linop_d->isBoss() && !me ) 
	    printf("bfmbase::MCR_PolyPrec: true MPright residual is %le \n",true_residual);


	  nn=linop_d->norm(psi);
	  if ( linop_d->isBoss() && !me ) 
	    printf("bfmbase::MCR_PolyPrec: psi %le \n",nn);

	  // 8 Feb 2013: BUG : I am finding P_L M P_R psi != P_L M psi
	  //                 
	  //          Should find also P_L M = M P_R
	  //
	  // Apply P_R psi
	  Pright(psi,tmp);
	  linop_d->axpy(psi,tmp,tmp,0.0);

	  nn=linop_d->norm(psi);
	  if ( linop_d->isBoss() && !me ) 
	    printf("bfmbase::MCR_PolyPrec: P_R psi %le \n",nn);

	  Pright(psi,tmp);
	  linop_d->axpy(psi,tmp,tmp,0.0);

	  nn=linop_d->norm(psi);
	  if ( linop_d->isBoss() && !me ) 
	    printf("bfmbase::MCR_PolyPrec: P_R P_R psi %le \n",nn);

	  // Redetermine the defect
	  PleftM(psi,tmp);
	  linop_d->axpy(tmp,tmp,src_left,-1.0);
	  true_residual = sqrt(linop_d->norm(tmp)/linop_d->norm(src_left));
	  if ( linop_d->isBoss() && !me ) 
	    printf("bfmbase::MCR_PolyPrec: true PleftM residual is %le \n",true_residual);

	  Pleft(psi,tmp);
	  nn=linop_d->norm(psi);
	  if ( linop_d->isBoss() && !me ) 
	    printf("bfmbase::MCR_PolyPrec: true PleftM psi is %le \n",nn);
	  
	  nn=linop_d->norm(tmp);
	  if ( linop_d->isBoss() && !me ) 
	    printf("bfmbase::MCR_PolyPrec: true Pleft PleftM psi is %le \n",nn);

	  // 1-P_R determined by Mss^-1 eta_s
	  LittleDopSolverResidual = 1.0e-9;
	  ProjectToSubspace(src_uprec,PleftProj);   
 	  ApplyInverse(PleftProj,PleftMss_proj);     // Mss^inv eta_s
	  PromoteFromSubspace(PleftMss_proj,tmp);
	  
	  double term = linop_d->norm(tmp);
	  if ( linop_d->isBoss() && !me ) 
	    printf("bfmbase::MCR_PolyPrec: Subspace soln is %le \n",term);
	  term = linop_d->norm(psi);
	  if ( linop_d->isBoss() && !me ) 
	    printf("bfmbase::MCR_PolyPrec: Complement soln is %le \n",term);

	  linop_d->axpy(psi,tmp,psi,1.0);

	  linop_d->Mprec(psi,mp,tmp,0);
	  linop_d->Mprec(mp,Ap ,tmp,1); 

	  linop_d->axpy(tmp,src_uprec,Ap,-1.0);
	  true_residual = sqrt(linop_d->norm(tmp)/linop_d->norm(src_uprec));
	  if ( linop_d->isBoss() && !me ) 
	    printf("bfmbase::MCR_PolyPrec: true full residual is %le \n",true_residual);

	  if ( linop_d->isBoss() && !me ) {
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


	  linop_d->Mprec(psi,mp,tmp,0);
	  linop_d->Mprec(mp,Ap,tmp,1); 
	  linop_d->axpy(tmp,src_uprec,Ap,-1.0);
	  double true_residual = sqrt(linop_d->norm(tmp)/linop_d->norm(src_uprec));

	  if ( linop_d->isBoss() && !me ) 
	    printf("bfmbase::MCR_PolyPrec: true unprec residual is %le \n",true_residual);


	  if ( linop_d->isBoss() && !me ) {
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

	linop_d->threadedFreeFermion(src);
	linop_d->threadedFreeFermion(src_left);
	linop_d->threadedFreeFermion(tmp);
	linop_d->threadedFreeFermion(p);
	linop_d->threadedFreeFermion(mp);
	linop_d->threadedFreeFermion(mmp);
	linop_d->threadedFreeFermion(r);
	linop_d->threadedFreeFermion(Ap);
	linop_d->threadedFreeFermion(Ar);
	return k;
	
    }
    ci = d; //(Ap,Ap)
   //d = real( dot(r,Ar) );

    //    d = Mprec(r,mp,tmp,DaggerNo); //(r,Ar)
    //        Mprec(mp,Ar,tmp,DaggerYes);

    SolverMatrix(r,Ar);
    d = linop_d->inner_real(r,Ar);
    b = d/c;
    c = d;
    linop_d->axpy(p,p,r,b);
    d = linop_d->axpy_norm(Ap,Ap,Ar,b);


  }
  if ( linop_d->isBoss() && !me ) printf("bfmbase::MCR_PolyPrec: CG not converged \n");
  linop_d->threadedFreeFermion(src);
  linop_d->threadedFreeFermion(tmp);
  linop_d->threadedFreeFermion(p);
  linop_d->threadedFreeFermion(mp);
  linop_d->threadedFreeFermion(mmp);
  linop_d->threadedFreeFermion(r);
  linop_d->threadedFreeFermion(Ap);
  linop_d->threadedFreeFermion(Ar);
  if ( linop_d->isBoss() && (!me) ) { 
    linop_d->InverterExit();
  }

  return -1;

}

void BfmHDCG<cFloat>::PolyMdagMprecRight(Fermion_t in,Fermion_t out)
{
  int me = linop_d->thread_barrier();
  if ( RightPrecondition ) { 
    if ( linop_d->isBoss() && !me ) {
      printf("Chebyshev right preconditioning matrix x^%f\n",alpha);fflush(stdout);
    }
    PolyMdagMprec(in,out);
  } else { 
    linop_d->axpy(out,in,in,0.0); 
  }
  linop_d->thread_barrier();
}
void BfmHDCG<cFloat>::PolyMdagMprecLeft(Fermion_t in,Fermion_t out)
{

  int me = linop_d->thread_barrier();

  if ( KrylovPrecondition ) {

    linop_d->thread_barrier();
    double restore_rsd  = linop_d->residual;
    double restore_iter=  linop_d->max_iter;
    linop_d->thread_barrier();
    linop_d->residual = PreconditionerKrylovResidual;
    linop_d->max_iter = PreconditionerKrylovIterMax;
    linop_d->thread_barrier();

    Fermion_t pin;
    Fermion_t pout;

    if ( linop_d->Ls < linop_d->Ls ) {
      printf("Cannot increase Ls in preconditioner");
      exit(0);
    }
    if ( linop_d->Ls != linop_d->Ls ) {
      if ( linop_d->isBoss() && !me ) printf("Reducing Ls from %d to %d\n",linop_d->Ls,linop_d->Ls);

      pin =linop_d->threadedAllocFermion();
      pout=linop_d->threadedAllocFermion();
      for(int s=0;s<linop_d->Ls/2;s++){
	linop_d->copy_slice(in ,s,linop_d->Ls,pin,s,linop_d->Ls);
	linop_d->copy_slice(in ,linop_d->Ls-s-1,linop_d->Ls,pin,linop_d->Ls-s-1,linop_d->Ls);
      }
      linop_d->axpy(out,in,in,0.0);
      linop_d->scale(out,1.0e-2);
    } else { 
      pin = in;
      pout=out;
    }

    if ( linop != preclinop) {
      linop_d->thread_barrier();
      if ( linop_d->SPIcomms() && !me ) linop_d->comm_init();
      linop_d->thread_barrier();
    }


    linop_d->axpy(pout,pin,pin,0.0);
    linop_d->thread_barrier();

    
    //    linop_d->MCR_prec(pout,pin);
    //    linop_d->CGNE_prec(pout,pin);


    int nshift=1;
    int single=0;
    double alpha[1]     = {1.0};
    double shifts[1]    = {PreconditionerKrylovShift};
    double mresidual[1] = {PreconditionerKrylovResidual};
    linop_d->CGNE_prec_MdagM_multi_shift(&pout,
					   pin,
					   shifts,
					   alpha,
					   nshift,
					   mresidual,
					   single);

    linop_d->thread_barrier();
    linop_d->residual=restore_rsd;
    linop_d->max_iter=restore_iter;
    linop_d->thread_barrier();

    if ( linop_d->cbLs != linop_d->cbLs ) {
      for(int s=0;s<linop_d->Ls/2;s++){
	linop_d->copy_slice(pout,s,linop_d->Ls,out,s,linop_d->Ls);
	linop_d->copy_slice(pout,linop_d->Ls-s-1,linop_d->Ls,out,linop_d->Ls-s-1,linop_d->Ls);
      }
      linop_d->threadedFreeFermion(pin);
      linop_d->threadedFreeFermion(pout);
    }

    if ( linop != preclinop) {
      linop_d->thread_barrier();
      if ( linop_d->SPIcomms() && !me )  linop_d->comm_init();
      linop_d->thread_barrier();
    }

    return;
  } else if ( LeftPrecondition ) { 
    if ( linop_d->isBoss() && !me ) printf("Chebyshev left preconditioning matrix x^%f\n",alpha);fflush(stdout);
    PolyMdagMprec(in,out);
    linop_d->thread_barrier();
  } else { 
    linop_d->axpy(out,in,in,0.0); 
    linop_d->thread_barrier();
  }
}

void BfmHDCG<cFloat>::PolyMdagMprec(Fermion_t in,Fermion_t out)
{
  PolyMdagMprec(linop,in,out);
}

void BfmHDCG<cFloat>::DeflatedSolverMpc(Fermion_t in,Fermion_t out)
{
    
  Fermion_t tmp = linop_d->threadedAllocFermion();
  Fermion_t tmp2= linop_d->threadedAllocFermion();
  Fermion_t Mtmp= linop_d->threadedAllocFermion();

  exit(-1);

  PolyMdagMprecRight(in,tmp2);
  PleftM(tmp2,tmp);
  PolyMdagMprecLeft(tmp,out);
    
  linop_d->threadedFreeFermion(tmp);
  linop_d->threadedFreeFermion(tmp2);
  linop_d->threadedFreeFermion(Mtmp);
}

void BfmHDCG<cFloat>::DeflatedSolverMdagMpc(Fermion_t in,Fermion_t out)
{
    
  Fermion_t tmp = linop_d->threadedAllocFermion();
  Fermion_t tmp2= linop_d->threadedAllocFermion();
  Fermion_t Mtmp= linop_d->threadedAllocFermion();

  uint64_t t1 = GetTimeBase();

  PolyMdagMprecRight(in,tmp2);

  uint64_t t2 = GetTimeBase();

  PleftM(tmp2,tmp);

  uint64_t t3 = GetTimeBase();

  PolyMdagMprecLeft(tmp,out);

  uint64_t t4 = GetTimeBase();
  
  int me = linop_d->thread_barrier();
  if( linop_d->isBoss() && !me ) {
    printf("Deflated Poly/PleftM/Poly cyc = %ld/%ld/%ld\n",t2-t1,t3-t2,t4-t3);
  }
    
  linop_d->threadedFreeFermion(tmp);
  linop_d->threadedFreeFermion(tmp2);
  linop_d->threadedFreeFermion(Mtmp);
}

void BfmHDCG<cFloat>::UndeflatedSolverMpc(Fermion_t in,Fermion_t out)
{
  Fermion_t tmp = linop_d->threadedAllocFermion();
  Fermion_t Mtmp= linop_d->threadedAllocFermion();
  
  linop_d->Mprec(in,out,Mtmp,DaggerYes);    // Non herm: M P_{-1}(Mdag M) Mdag
  PolyMdagMprecRight(out,tmp);

  linop_d->Mprec(tmp,out,Mtmp,DaggerNo);  
    
  linop_d->threadedFreeFermion(tmp);
  linop_d->threadedFreeFermion(Mtmp);
};

void BfmHDCG<cFloat>::UndeflatedSolverMdagMpc(Fermion_t in,Fermion_t out)
{
  Fermion_t tmp = linop_d->threadedAllocFermion();
  Fermion_t Mtmp= linop_d->threadedAllocFermion();
  
  PolyMdagMprecRight(in,tmp);

  linop_d->Mprec(tmp,out,Mtmp,DaggerNo);  
  linop_d->Mprec(out,tmp,Mtmp,DaggerYes);    

  PolyMdagMprecLeft(tmp,out);
    
  linop_d->threadedFreeFermion(tmp);
  linop_d->threadedFreeFermion(Mtmp);
};

int BfmHDCG<cFloat>::GCR_PolyPrec(Fermion_t psi, Fermion_t src,int m) // Internal
{
  double f;
  double cp,c,a,d,b;
  int me = linop_d->thread_barrier();

  SolverControl = SolverUndeflatedMpc;
  ChebyshevInit(PreconChebyLo,
		PreconChebyHi,
		-1.0,
		PreconChebyOrder);

  if ( linop_d->isBoss() && (!me) ){
    linop_d->InverterEnter();
  }
 
  Fermion_t p   =linop_d->threadedAllocFermion();
  Fermion_t tmp =linop_d->threadedAllocFermion();
  Fermion_t Mtmp=linop_d->threadedAllocFermion();
  Fermion_t Ap  =linop_d->threadedAllocFermion();
  Fermion_t Ar  =linop_d->threadedAllocFermion();
  Fermion_t r   =linop_d->threadedAllocFermion();
  
  //Initial residual computation & set up
  double guess = linop_d->norm(psi);

  SolverMatrix(psi,Ar);
  linop_d->axpy(r, Ar, src,-1.0);     //r = b - Ax ; b->src
  linop_d->copy(p,r);                 //p0 = r0
 
  a =linop_d->norm(p);
  cp=linop_d->norm(r);
 
  double ssq =  linop_d->norm(src);
  double residual=linop_d->residual;
  double rsq =  residual* residual*ssq;
 
  //Check if guess is really REALLY good :)
  if ( cp <= rsq ){
    
        if ( linop_d->verbose && linop_d->isBoss() && !me )    {
	  printf("bfmbase::GCR_prec k=0 converged - suspiciously nice guess %le %le\n",cp,rsq);
        }
        linop_d->threadedFreeFermion(tmp);
        linop_d->threadedFreeFermion(Mtmp);
        linop_d->threadedFreeFermion(p);
        linop_d->threadedFreeFermion(Ar);
        linop_d->threadedFreeFermion(Ap);
        linop_d->threadedFreeFermion(r);
        if ( linop_d->isBoss() && (!me) ) {
            linop_d->InverterExit();
        }
        return 0;
  }

  if ( linop_d->verbose && linop_d->isBoss() && !me )
    printf("bfmbase::GCR_prec k=0 residual %le rsq %le\n",cp,rsq);

 
  Fermion_t allAp[m-1];
  Fermion_t allp[m-1];
  double ApAp[m-1];
  for(int i=0;i<m-1;i++) {
    allAp[i]   = linop_d->threadedAllocFermion();
    allp[i]    = linop_d->threadedAllocFermion();
  } 

  SolverMatrix(p,Ap);

  std::complex<double> alpha;
  std::complex<double> c_rAp  = linop_d->inner(r,Ap);          //c = (r,Ap)
  double               c_ApAp = linop_d->inner_real(Ap,Ap);
 
  linop_d->copy(allAp[0],Ap);
  linop_d->copy(allp[0],p);
  ApAp[0] = c_ApAp;
  
  for (int k=1; k<=linop_d->max_iter; k++) {
 
    std::complex<double> beta[k];
    linop_d->iter=k;

    alpha = c_rAp/c_ApAp;                                      //a = (r,Ap)/(Ap,Ap) 

    linop_d->caxpy(psi,p,psi,real(alpha),imag(alpha));	       //x = x + ap
    linop_d->caxpy(r,Ap,r,-real(alpha),-imag(alpha));	       //r = r - aAp
    cp = linop_d->norm(r);
  
    SolverMatrix(r,Ar);
    
    linop_d->copy(Ap,Ar) ; // Ap=Ar         // Why

    linop_d->copy(p,r);    // p=r
    for(int i=0; i<k%m; i++) {
      beta[i] = linop_d->inner(Ar,allAp[i])/ApAp[i];                        //alpha = (Ar,Ap)/(Ap,Ap)
      linop_d->caxpy(Ap,allAp[i],Ap,-real(beta[i]),-imag(beta[i]));         //Ap = Ar + B Ap
      linop_d->caxpy(p,allp[i],p,-real(beta[i]),-imag(beta[i]));            //p  =  r + B p
    }
 
    c_ApAp = linop_d->norm(Ap);                    //d = (Ap,Ap)
    c_rAp  = linop_d->inner(r,Ap);                 //c = (r,Ap)
  
                    
    if( k%m < m-1 ) {      // we dont need to store when we need to restart
      linop_d->copy(allAp[k%m],Ap);
      linop_d->copy(allp[k%m],p);
      ApAp[k%m] = c_ApAp;
    }
 
    if ( k%10 == 0){
      if (linop_d->isBoss() && !me ){
	printf("bfmbase::GCR_prec: k= %d r^2= %le %le \n",k,cp,sqrt(cp/ssq));
      }
    }
 
    // Stopping condition
    if ( cp <= rsq ) {
      
      if (linop_d->isBoss() && !me ) printf("bfmbase::GCR_prec converged in %d iterations\n",k);

      SolverMatrix(psi,Ap);

      linop_d->axpy(tmp,src,Ap,-1.0);
      double true_residual = sqrt(linop_d->norm(tmp)/linop_d->norm(src));
      if ( linop_d->isBoss() && !me )
	printf("bfmbase::GCR_prec: true residual is %le \n",true_residual);
      
      linop_d->threadedFreeFermion(tmp);
      linop_d->threadedFreeFermion(Mtmp);
      linop_d->threadedFreeFermion(p);
      linop_d->threadedFreeFermion(r);
      linop_d->threadedFreeFermion(Ap);
      linop_d->threadedFreeFermion(Ar);
      if ( linop_d->isBoss() && (!me) ){
	linop_d->InverterExit();
      }
      return k;
    }
 
 
  }

  if (linop_d->isBoss() && !me ) printf("bfmbase::GCR_prec: CG not converged \n");
  linop_d->threadedFreeFermion(tmp);
  linop_d->threadedFreeFermion(Mtmp);
  linop_d->threadedFreeFermion(p);
  linop_d->threadedFreeFermion(r);
  linop_d->threadedFreeFermion(Ap);
  linop_d->threadedFreeFermion(Ar);
  if ( linop_d->isBoss() && (!me) ) {
    linop_d->InverterExit();
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

int  BfmHDCG<cFloat>::rCG_PolyPrec(Fermion_t psi, Fermion_t src_uprec,Fermion_t resid)
{
  int steps= 10;
  linop_d->thread_barrier();
  int miter=linop_d->max_iter;
  linop_d->thread_barrier();
  int iters=0;
  linop_d->max_iter = steps;
  linop_d->thread_barrier();

  while(iters < miter){
    int it = CG_PolyPrec(psi,src_uprec,resid);
    iters+=steps;
    if ( it > 0 ) {
      linop_d->thread_barrier();
      linop_d->max_iter = miter;
      linop_d->thread_barrier();
      return iters;
    }
  }
  return -1;
}

int BfmHDCG<cFloat>::CG_PolyPrec(Fermion_t psi, Fermion_t src,Fermion_t resid) 
{
  double f;
  double rtzp,rtz,a,d,b;

  int me = linop_d->thread_barrier();


  // Deflated precond CG (3.6 in Saad umsi-98-97)
  KrylovPrecondition=1;

  ChebyshevInit(PreconChebyLo,
		PreconChebyHi,
		-1.0,
		PreconChebyOrder);

  /*
  if ( linop_d->isBoss() && (!me) ) { 
    printf("BfmHDCG<cFloat>::CG_PolyPrec: Precondition with Chebyshev x^{%f} [%f,%f] O(%d) \n",
	   PreconChebyLo,
	   PreconChebyHi,
	   -1.0,
	   PreconChebyOrder
	   );
    linop_d->InverterEnter();
  }
  */

  Fermion_t p   = linop_d->threadedAllocFermion(); 
  Fermion_t x   = psi;
  Fermion_t z   = linop_d->threadedAllocFermion(); 

  Fermion_t tmp = linop_d->threadedAllocFermion(); 
  Fermion_t mp  = linop_d->threadedAllocFermion(); 
  Fermion_t mmp = linop_d->threadedAllocFermion(); 
  Fermion_t r   = linop_d->threadedAllocFermion(); 
  Fermion_t mu   = linop_d->threadedAllocFermion(); 

  //Initial residual computation & set up
  double guess = linop_d->norm(psi);

  ///////////////////////////////////
  // x_{-1} = psi
  ///////////////////////////////////
  linop_d->Mprec(x,mp,tmp,DaggerNo);
  linop_d->Mprec(mp,mmp,tmp,DaggerYes);
  linop_d->axpy (r, mmp, src,-1.0);        // r_{-1} = src - A x

  ///////////////////////////////////
  // Choose x_0 such that 
  // x_0 = guess +  (A_ss^inv) r_s
  // 
  // W^T (src - A x_0) = src_s - A guess_s - r_s  
  //                   = src_s - (A guess)_s - src_s  + (A guess)_s 
  //                   = 0 
  ///////////////////////////////////
  ProjectToSubspace(r,PleftProj);     
  if(linop_d->isBoss()&&!me){ 
    printf("Inverting little dirac op\n");fflush(stdout);
  }
  ApplyInverse(PleftProj,PleftMss_proj); // Ass^{-1} r_s
  if(linop_d->isBoss()&&!me){ 
    printf("Inverted little dirac op\n");fflush(stdout);
  }
  PromoteFromSubspace(PleftMss_proj,mp);  

  linop_d->axpy(x,x,mp,1.0);

  linop_d->Mprec(x,mp,tmp,DaggerNo);
  linop_d->Mprec(mp,mmp,tmp,DaggerYes);
  linop_d->axpy (r, mmp, src,-1.0);    // Recomputes r=src-x0

  ProjectToSubspace(r,PleftProj);     
  double nv=norm_vec(PleftProj);
  if(linop_d->isBoss()&&!me){ 
    printf("Residual after projection= %le\n",nv);
  }

  // Compute z = M x
  PolyMdagMprecLeft(r,z);
  rtzp =real(linop_d->inner(r,z));

  ///////////////////////////////////////
  // Solve for Mss mu = P A z and set p = z-mu
  ///////////////////////////////////////
  linop_d->Mprec(z,mp,tmp,DaggerNo);
  linop_d->Mprec(mp,mmp,tmp,DaggerYes); //Az

  ProjectToSubspace(mmp,PleftProj);     
  ApplyInverse(PleftProj,PleftMss_proj); // Mss^{-1} Az
  PromoteFromSubspace(PleftMss_proj,mu);  

  linop_d->axpy (p, mu, z,-1.0);

  double ssq =  linop_d->norm(src);
  double rsq =  linop_d->residual* linop_d->residual*ssq;

  if ( linop_d->isBoss() && !me ) 
    printf("BfmHDCG<cFloat>::CG_PolyPrec k=0 residual %le rsq %le\n",rtzp,rsq);

  uint64_t t_start=GetTimeBase();
  uint64_t t1;
  uint64_t t_Mss=0;
  uint64_t t_proj=0;
  for (int k=1;k<=linop_d->max_iter;k++){

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
    d        =linop_d->Mprec(p,mp,tmp,0,1);// Dag no
              linop_d->Mprec(mp,mmp,tmp,1);// Dag yes
    a = rtz/d;

    //    if(linop_d->isBoss() && !me ) { printf("CG poly 1\n"); fflush(stdout);}
    //    linop_d->norm(x);

    linop_d->axpy(x,p,x,a);
    linop_d->axpy(r,mmp,r,-a);
    linop_d->axpy(resid,r,r,0.0);

    //    if(linop_d->isBoss() && !me ) { printf("CG poly 2\n"); fflush(stdout);}
    //    linop_d->norm(x);

    // Compute z = M x
    PolyMdagMprecLeft(r,z);

    //    if(linop_d->isBoss() && !me ) { printf("CG poly 3\n"); fflush(stdout);}
    //    linop_d->norm(x);

    rtzp =real(linop_d->inner(r,z));
    b = rtzp/rtz;

    //    if(linop_d->isBoss() && !me ) { printf("CG poly 4\n"); fflush(stdout);}
    //    linop_d->norm(x);

    // WAW mu = W A z
    linop_d->Mprec(z,mp,tmp,DaggerNo);
    linop_d->Mprec(mp,mmp,tmp,DaggerYes);

    //    if(linop_d->isBoss() && !me ) { printf("CG poly 5\n"); fflush(stdout);}
    //    linop_d->norm(x);

    t1=GetTimeBase();
    ProjectToSubspace(mmp,PleftProj);     // WAz
    t_proj+=GetTimeBase()-t1;

    //    if(linop_d->isBoss() && !me ) { printf("CG poly 6\n"); fflush(stdout);}
    //    linop_d->norm(x);
    t1=GetTimeBase();
    ApplyInverse(PleftProj,PleftMss_proj); // Mss^{-1} in_s       // (WtAW)^[-1] Az
    t_Mss+=GetTimeBase()-t1;

    //    if(linop_d->isBoss() && !me ) { printf("CG poly 7\n"); fflush(stdout);}
    //    linop_d->norm(x);

    t1=GetTimeBase();
    PromoteFromSubspace(PleftMss_proj,mu);     // mu = Wt  (Wt A W)^{-1} WAz
    t_proj+=GetTimeBase()-t1;

    linop_d->axpy(p,p,z,b);
    linop_d->axpy(p,mu,p,-1.0);

    //    if(linop_d->isBoss() && !me ) { printf("CG poly 8\n"); fflush(stdout);}
    //    linop_d->norm(x);

    {
      double rr=sqrt(rtzp/ssq);
      if (linop_d->isBoss() && !me ){
	printf("bfmbase::CG_PolyPrec: k= %d residual = %le \n",k,rr);
      }
      //  linop_d->InverterLogIteration(k, rr,src,x,mp,mmp,tmp);
    }

    // Stopping condition
    if ( rtzp <= rsq ) { 
	
	t1=GetTimeBase();
      if ( linop_d->isBoss() && !me ) {
	printf("BfmHDCG<cFloat>::CG_PolyPrec converged in %d iterations %f s\n",
					   k,1.0e-9*(t1-t_start)/1.6);
	printf("BfmHDCG<cFloat>::CG_PolyPrec MssInv %f s\n",1.0e-9*(t_Mss)/1.6);
	printf("BfmHDCG<cFloat>::CG_PolyPrec Proj   %f s\n",1.0e-9*(t_proj)/1.6);
	printf("BfmHDCG<cFloat>::CG_PolyPrec other  %f s\n",1.0e-9*(t1-t_start-t_Mss-t_proj)/1.6);
      }
	
	linop_d->Mprec(x,mp,tmp,0);
	linop_d->Mprec(mp,mmp,tmp,1); 
	linop_d->axpy(tmp,src,mmp,-1.0);
	
	double  mpnorm = sqrt(linop_d->norm(mp));
	double mmpnorm = sqrt(linop_d->norm(mmp));
	double psinorm = sqrt(linop_d->norm(x));
	double srcnorm = sqrt(linop_d->norm(src));
	double tmpnorm = sqrt(linop_d->norm(tmp));
	double true_residual = tmpnorm/srcnorm;
	if ( linop_d->isBoss() && !me ) {
	  printf("BfmHDCG<cFloat>::CG_PolyPrec: true residual is %le, solution %le, source %le \n",true_residual,psinorm,srcnorm);
	  printf("BfmHDCG<cFloat>::CG_PolyPrec: target residual was %le \n",linop_d->residual);
	  printf("BfmHDCG<cFloat>::CG_PolyPrec: mp %le, mmp %le\n",mpnorm,mmpnorm);
	}
	
	linop_d->threadedFreeFermion(tmp);
	linop_d->threadedFreeFermion(p);
	linop_d->threadedFreeFermion(z);
	linop_d->threadedFreeFermion(mu);
	linop_d->threadedFreeFermion(mp);
	linop_d->threadedFreeFermion(mmp);
	linop_d->threadedFreeFermion(r);
	if ( linop_d->isBoss() && (!me) ) { 
	  linop_d->InverterExit();
	}
	return k;
    }

  }
  if ( linop_d->isBoss() && !me ) printf("BfmHDCG<cFloat>::CG_PolyPrec: CG not converged \n");
  linop_d->threadedFreeFermion(tmp);
  linop_d->threadedFreeFermion(p);
  linop_d->threadedFreeFermion(mp);
  linop_d->threadedFreeFermion(mmp);
  linop_d->threadedFreeFermion(z);
  linop_d->threadedFreeFermion(mu);
  linop_d->threadedFreeFermion(r);
  if ( linop_d->isBoss() && (!me) ) { 
    linop_d->InverterExit();
  }

  return -1;
}
//
// Flexible Conjugate Gradients [Notay 2000 SIAM Vol 22 No 4 Jsci Comput]
// Combine with Deflation by Saad
//
int BfmHDCG<cFloat>::fCG_PolyPrec(Fermion_t psi, Fermion_t src) 
{
  double f;
  double a,d,b, r_dot_z,r_dot_z_p;

  int me = linop_d->thread_barrier();

  KrylovPrecondition=1;

  // Deflated precond CG (3.6 in Saad umsi-98-97)
  ChebyshevInit(PreconChebyLo,
		PreconChebyHi,
		-1.0,
		PreconChebyOrder);

  if ( linop_d->isBoss() && (!me) ) { 
    printf("BfmHDCG<cFloat>::FlexCG_PolyPrec: Precondition with Chebyshev x^{%f} [%f,%f] O(%d) \n",
	   PreconChebyLo,
	   PreconChebyHi,
	   -1.0,
	   PreconChebyOrder
	   );
    linop_d->InverterEnter();
  }

  Fermion_t x   = psi;
  Fermion_t tmp = linop_d->threadedAllocFermion(); 
  Fermion_t mp  = linop_d->threadedAllocFermion(); 
  Fermion_t mmp = linop_d->threadedAllocFermion(); 
  Fermion_t r   = linop_d->threadedAllocFermion(); 
  Fermion_t mu  = linop_d->threadedAllocFermion(); 
  Fermion_t z   = linop_d->threadedAllocFermion(); 

  const int mmax = 20;
  int m;
  Fermion_t  p[mmax];
  Fermion_t Ap[mmax];
  double   pAp[mmax];
  for(m=0;m<mmax;m++){
    p [m] = linop_d->threadedAllocFermion();
    Ap[m] = linop_d->threadedAllocFermion();
  }
  //Initial residual computation & set up
  double guess = linop_d->norm(psi);

  ///////////////////////////////////
  // x_{-1} = psi
  ///////////////////////////////////
  linop_d->Mprec(x,mp,tmp,DaggerNo);
  linop_d->Mprec(mp,Ap[0],tmp,DaggerYes);
  linop_d->axpy (r, Ap[0], src,-1.0);        // r_{-1} = src - A x

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

  linop_d->axpy(x,x,mp,1.0);
  linop_d->axpy(p[0],x,x,0.0);

  linop_d->Mprec(p[0],mp,tmp,DaggerNo);
  linop_d->Mprec(mp,Ap[0],tmp,DaggerYes);
  linop_d->axpy (r, Ap[0], src,-1.0);    // Recomputes r=src-x0

  ProjectToSubspace(r,PleftProj);     
  double nv=norm_vec(PleftProj);
  if(linop_d->isBoss()&&!me){ 
    printf("Residual after projection= %le\n",nv);
  }

  // Compute z = M x
  PolyMdagMprecLeft(r,z);
  // Fixme: also could run a Krylov solver
  r_dot_z =real(linop_d->inner(r,z));
  
  ///////////////////////////////////////
  // Solve for Mss mu = P A z and set p = z-mu
  ///////////////////////////////////////
  linop_d->Mprec(z,mp,tmp,DaggerNo);
  linop_d->Mprec(mp,mmp,tmp,DaggerYes); //Az

  ProjectToSubspace(mmp,PleftProj);     
  ApplyInverse(PleftProj,PleftMss_proj); // Mss^{-1} Az
  PromoteFromSubspace(PleftMss_proj,mu);  

  linop_d->axpy (p[0], mu, z,-1.0);

  double ssq =  linop_d->norm(src);
  double rsq =  linop_d->residual* linop_d->residual*ssq;

  if ( linop_d->isBoss() && !me ) 
    printf("BfmHDCG<cFloat>::FlexCG_PolyPrec k=0 residual %le rsq %le\n",r_dot_z,rsq);

  uint64_t t_start = GetTimeBase();
  uint64_t t_Mss=0;
  uint64_t t_proj=0;
  uint64_t t_ortho=0;
  uint64_t t_prec=0;
  uint64_t t1;

  for (int k=0;k<=linop_d->max_iter;k++){

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
    
    pAp[peri_k]        =linop_d->Mprec(p[peri_k],mp,tmp,0,1);// Dag no
                        linop_d->Mprec(mp,Ap[peri_k],tmp,1);// Dag yes

    a = r_dot_z/pAp[peri_k];
    linop_d->axpy(x, p[peri_k],x,a);
    linop_d->axpy(r,Ap[peri_k],r,-a);

    // Compute z = M x;  This is w_i = B(r_i) in Notay
    t1=GetTimeBase();
    linop_d->thread_barrier();
    PolyMdagMprecLeft(r,z);
    linop_d->thread_barrier();
    t_prec+=GetTimeBase()-t1;

    double r_dot_z_p = r_dot_z;
    r_dot_z =real(linop_d->inner(r,z));

    // Solve WAW mu = W A z
    linop_d->Mprec(z,mp,tmp,DaggerNo);
    linop_d->Mprec(mp,mmp,tmp,DaggerYes);

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
    linop_d->axpy(p[peri_kp],mu,z,-1.0);     // Deflation bit
    
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

    if  ( linop_d->isBoss() && !me ) printf("BfmHDCG<cFloat>::FlexCG_PolyPrec iteration %d : orthogonalising to last %d vectors\n",k,northog);

    t1=GetTimeBase();
    for(int back=0; back < northog; back++){

      int peri_back = (k-back)%mmax;

      piApk= real(linop_d->inner(Ap[peri_back],p[peri_kp]));

      beta = -piApk/pAp[peri_back];

      linop_d->axpy(p[peri_kp],p[peri_back],p[peri_kp],beta);
      linop_d->thread_barrier();
      if ( linop_d->isBoss() && !me ) printf("BfmHDCG<cFloat>::FlexCG_PolyPrec iteration %d : Beta[%d] = %le\n",k,back,beta);
      fflush(stdout);
      linop_d->thread_barrier();

    }
    t_ortho+=GetTimeBase() - t1;
#else 
    double beta = r_dot_z/r_dot_z_p;
    linop_d->axpy(p[peri_kp],p[peri_k],p[peri_kp],beta);
#endif
    // Stopping condition

    if ( linop_d->isBoss() && !me ) printf("BfmHDCG<cFloat>::FlexCG_PolyPrec iteration %d residual %le\n",k,sqrt(r_dot_z/ssq));
    if(1){
      double nn = linop_d->norm(r);
      if( linop_d->isBoss() && !me ) printf("resid = %le\n",nn);
      for(int s=0;s<linop_d->cbLs;s++){
	nn = linop_d->norm(r,s);
	if( linop_d->isBoss() && !me ) printf("resid[s] = %le\n",nn);
      }
      fflush(stdout);
      linop_d->thread_barrier();
    }

    if ( r_dot_z <= rsq ) { 
      
      t1=GetTimeBase();
      
      if ( linop_d->isBoss() && !me ) {
	printf("BfmHDCG<cFloat>::FlexCG_PolyPrec converged in %d iterations %f s\n",
					   k,1.0e-9*(t1-t_start)/1.6);
	printf("BfmHDCG<cFloat>::FlexCG_PolyPrec MssInv %f s\n",1.0e-9*(t_Mss)/1.6);
	printf("BfmHDCG<cFloat>::FlexCG_PolyPrec Proj   %f s\n",1.0e-9*(t_proj)/1.6);
	printf("BfmHDCG<cFloat>::FlexCG_PolyPrec Ortho  %f s\n",1.0e-9*(t_ortho)/1.6);
	printf("BfmHDCG<cFloat>::FlexCG_PolyPrec Prec  %f s\n",1.0e-9*(t_prec)/1.6);
	printf("BfmHDCG<cFloat>::FlexCG_PolyPrec other  %f s\n",1.0e-9*(t1-t_start-t_Mss-t_proj-t_ortho-t_prec)/1.6);
      }
	fflush(stdout);
	linop_d->thread_barrier();
	
	linop_d->Mprec(x,mp,tmp,0);
	if ( linop_d->isBoss() && !me ) printf("BfmHDCG<cFloat>::FlexCG_PolyPrec Mprec ");
	fflush(stdout);
	linop_d->thread_barrier();
	linop_d->Mprec(mp,mmp,tmp,1); 
	if ( linop_d->isBoss() && !me ) printf("BfmHDCG<cFloat>::FlexCG_PolyPrec Mprec ");
	fflush(stdout);
	linop_d->thread_barrier();
	linop_d->axpy(tmp,src,mmp,-1.0);
	if ( linop_d->isBoss() && !me ) printf("BfmHDCG<cFloat>::FlexCG_PolyPrec axpy ");
	fflush(stdout);
	linop_d->thread_barrier();
	
	double  mpnorm = sqrt(linop_d->norm(mp));
	double mmpnorm = sqrt(linop_d->norm(mmp));
	double psinorm = sqrt(linop_d->norm(x));
	double srcnorm = sqrt(linop_d->norm(src));
	double tmpnorm = sqrt(linop_d->norm(tmp));
	double true_residual = tmpnorm/srcnorm;
	if ( linop_d->isBoss() && !me ) {
	  printf("BfmHDCG<cFloat>::FlexCG_PolyPrec: true residual is %le, solution %le, source %le \n",true_residual,psinorm,srcnorm);
	  printf("BfmHDCG<cFloat>::FlexCG_PolyPrec: target residual was %le \n",linop_d->residual);
	  printf("BfmHDCG<cFloat>::FlexCG_PolyPrec: mp %le, mmp %le\n",mpnorm,mmpnorm);
	  fflush(stdout);
	}
	linop_d->thread_barrier();
	linop_d->threadedFreeFermion(tmp);
	linop_d->threadedFreeFermion(z);
	linop_d->threadedFreeFermion(mu);
	linop_d->threadedFreeFermion(mp);
	linop_d->threadedFreeFermion(mmp);
	linop_d->threadedFreeFermion(r);
	for(int m=0;m<mmax;m++){
	  linop_d->threadedFreeFermion(p[m]);
	  linop_d->threadedFreeFermion(Ap[m]);
	}
	if ( linop_d->isBoss() && (!me) ) { 
	  linop_d->InverterExit();
	}
	return k;
    }

  }
  if ( linop_d->isBoss() && !me ) printf("BfmHDCG<cFloat>::FlexCG_PolyPrec: CG not converged \n");
  linop_d->threadedFreeFermion(tmp);
  linop_d->threadedFreeFermion(mp);
  linop_d->threadedFreeFermion(mmp);
  linop_d->threadedFreeFermion(z);
  linop_d->threadedFreeFermion(mu);
  linop_d->threadedFreeFermion(r);
  for(int m=0;m<mmax;m++){
   linop_d->threadedFreeFermion(p[m]);
   linop_d->threadedFreeFermion(Ap[m]);
  }

  if ( linop_d->isBoss() && (!me) ) { 
    linop_d->InverterExit();
  }

  return -1;
}

void BfmHDCG<cFloat>::PolyLdop(std::vector<std::complex<double> > &in,std::vector<std::complex<double> > &out)
{

  int N = LocalNsubspace;

  if(in.size()  != N ) exit(0);
  if(out.size() != N ) exit(0);

  std::vector<std::complex<double> > y(N);
  std::vector<std::complex<double> > Tnm(N);
  std::vector<std::complex<double> > Tn (N);
  std::vector<std::complex<double> > Tnp(N);
  
  int me = linop_d->thread_barrier();

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

template void BfmHDCG<float>::PolyMdagMprec<float>(bfm_qdp<float> *lop,Fermion_t in,Fermion_t out);
template void BfmHDCG<float>::PolyMdagMprec<double>(bfm_qdp<double> *lop,Fermion_t in,Fermion_t out);
template void BfmHDCG<double>::PolyMdagMprec<float>(bfm_qdp<float> *lop,Fermion_t in,Fermion_t out);
template void BfmHDCG<double>::PolyMdagMprec<double>(bfm_qdp<double> *lop,Fermion_t in,Fermion_t out);

    /*
    
    if ( Smoother ) {

      linop_f->thread_barrier();
      if ( linop_smooth->SPIcomms() && !me )  linop_smooth->comm_init();
      linop_f->thread_barrier();
      PcgMlinop<float>(out,tmp,Mtmp,linop_smooth); // Invert smaller system
      linop_f->thread_barrier();
      if ( linop_f->SPIcomms() && !me )  linop_f->comm_init();
      linop_f->thread_barrier();

      linop_f->PromIncreaseLs (linop_smooth->BlockFifthDimensionSize,
			       linop_smooth->BlockFifthDimensionWeight,
			       tmp,Min,Mtmp); // Promote Ls

      // Hack/test :: Compute true residual
      // res  = (MdagM+shift) Min - in_f
      linop_f->Mprec(Min,tmp,Mtmp,DaggerNo);
      linop_f->Mprec(tmp,out,Mtmp,DaggerYes);  
      linop_f->axpy(tmp,Min,out,PreconditionerKrylovShift); 
      linop_f->axpy(tmp,out,mgtmp,-1.0);        
      for(int s=0;s<linop_f->cbLs;s++){
	double nr,ni,ns;
	nr = linop_f->norm(tmp,s);
	ni = linop_f->norm(mgtmp,s);
	ns = linop_f->norm(Min,s);
	linop_f->ThreadBossMessage("[%d] resid %le source %le sol %le\n",s,nr,ni,ns);
      }

    } 
    */
/*
template<class cFloat> 
template<class Float> 
void BfmHDCG<cFloat>::PcgMDDlinop(Fermion_t in,Fermion_t out,Fermion_t tmp,bfm_qdp<Float> *lop,bfm_qdp<Float> *lop_wall)
{
// Optional convert to single precision in this????
  int me = lop->thread_barrier();
  double restore_rsd  = lop->residual;
  double restore_iter=  lop->max_iter;
  lop->thread_barrier();
  lop->residual = PreconditionerKrylovResidual;
  lop->max_iter = PreconditionerKrylovIterMax;
  lop->thread_barrier();
  int Ls = lop->cbLs;
  double shift = PreconditionerKrylovShift;
  double lshift= PreconditionerKrylovShift;

  uint64_t t1 = GetTimeBase();

  // Need "out" to become approximate inverse of "in".
  // Apply more iterations to surface projected vector
  
  lop->ThreadBossLog("Domain solver on surface\n");

  // r = P_surf in
  //
  // P_surface [Mdag Psurface M+lambda] Psurface g = r
  //
  // r = in - [MdagM+lambda] g  
  // 
  // out = [MdagM+lambda]^inv r + g
  // 
  // => MdagM+lambda out = in - [MdagM+lambda] g + [MdagM+lambda] g  = in.
  //
  Fermion_t r   = lop->threadedAllocFermion();
  Fermion_t g   = lop->threadedAllocFermion();
  Fermion_t mp  = lop->threadedAllocFermion();
  Fermion_t mmp = lop->threadedAllocFermion();
  double r_n,s_n,ss_n;

  lop->thread_barrier();
  lop_wall->max_iter             = 16;
  lop->thread_barrier();

  ///////////////////////////////////////////
  // restrict source to surface.. 
  ///////////////////////////////////////////
  int depth=1;
  lop->axpy(r,in,in,0.0);
  for(int s=depth;s<Ls-depth;s++){
    lop->axpby_ssp(r,0.0,r,0.0,r,s,s); // Zero middle
  }


  // Printing info
  if ( 1 ) { 
    ss_n = lop->norm(in);
    lop->ThreadBossLog("Source: %le \n",ss_n);
    
    s_n = lop->norm(r);
    lop->ThreadBossLog("Surface source: %le \n",s_n);
  }

  lop->ProjReduceLs(lop_wall->BlockFifthDimensionSize,
		    lop_wall->BlockFifthDimensionWeight,
		    r,tmp,mp); // Reduce Ls

  lop_wall->thread_barrier();
  if ( lop_wall->SPIcomms() && !me )  lop_wall->comm_init();
  lop_wall->thread_barrier();

  // Solve the small Ls system for surface modes
  std::vector<double> dummy;
  lop_wall->CGNE_single_shift(mp,tmp,lshift,0,dummy);

  lop_wall->thread_barrier();
  if ( lop->SPIcomms() && !me )  lop->comm_init();
  lop_wall->thread_barrier();

  lop->PromIncreaseLs (lop_wall->BlockFifthDimensionSize,
		       lop_wall->BlockFifthDimensionWeight,
		       mp,g,mmp); // Promote Ls

  if ( 1 ) { 

    lop->Mprec(g,mp,tmp,0);
    lop->Mprec(mp,mmp,tmp,1);
    lop->axpy(mmp,g,mmp,lshift);

    lop->axpy(r,mmp,r,-1.0); 
    r_n = lop->norm(r);
    lop->ThreadBossMessage("Surface solution residual: %le \n",sqrt(r_n/s_n));

    if ( 1 ) { 

      for(int s=0;s<Ls;s++){
	double nn = lop->norm(r,s);
	lop->ThreadBossLog("surface resid[s] = %le\n",nn);
      }

    }
    // Recomputes residual
    lop->axpy(r,mmp,in,-1.0); 

  } else { 

    lop->axpy(r,in,in,1.0);
    for(int s=0;s<depth;s++){
      lop->axpby_ssp(r,0.0,r,0.0,r,s,s); // Zero wall
      lop->axpby_ssp(r,0.0,r,0.0,r,Ls-1-s,Ls-1-s); // Zero wall
    }

  }


  if ( 0 ){
    s_n = lop->norm(r);
    lop->ThreadBossLog("Bulk source: %le \n",s_n);
  }

  lop->thread_barrier();
  lop->max_iter             = 8;
  lop->thread_barrier();

  std::vector<double> dummy;
  lop->CGNE_single_shift(out,r,shift,0,dummy);

  if ( 0 ) {
    lop->Mprec(out,mp,tmp,0);
    lop->Mprec(mp,mmp,tmp,1);
    lop->axpy(mmp,out,mmp,shift);
    lop->axpy(mmp,mmp,r,-1.0); 
    r_n = lop->norm(mmp);

    lop->ThreadBossLog("Bulk solution residual: %le \n",sqrt(r_n/s_n));

  }
  // Add together the inverses
  lop->axpy(out,g,out,1.0);

  if ( 1 ) { 
    lop->Mprec(out,mp,tmp,0);
    lop->Mprec(mp,mmp,tmp,1);
    lop->axpy(mmp,out,mmp,shift);
    lop->axpy(r,mmp,in,-1.0); 
    
    r_n = lop->norm(r);
    lop->ThreadBossLog("Combined solution residual: %le \n",sqrt(r_n/ss_n));
  }

  lop->threadedFreeFermion(r);
  lop->threadedFreeFermion(g);
  lop->threadedFreeFermion(mmp);
  lop->threadedFreeFermion(mp);

  if(!me) PcgMprec+=GetTimeBase()-t1;

  lop->thread_barrier();
  lop->residual=restore_rsd;
  lop->max_iter=restore_iter;
  lop->thread_barrier();
  return;
}
*/

#endif


  } else { 
    this->Error("Oops not implemented G5D_MprecDeriv\n");
    exit(-1);
  }
}
