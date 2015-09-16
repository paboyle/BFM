#include "littleDiracOp.h"
#include <iostream>
#include <iomanip>


BfmLittleDiracOperator::BfmLittleDiracOperator(int _N5,
					       int _Nvec,
					       int _BlockSize[5],
					       Fermion_t *_subspace,
					       bfm_qdp<double> * _linop) 
  {
    N5=_N5;
    Nvec=_Nvec;
    subspace=_subspace;
    for(int mu=0;mu<5;mu++){
      BlockSize[mu]=_BlockSize[mu];
    }
    linop=_linop;
    QDPIO::cout << "Building subspace from " << Nvec<<" vectors"<<endl;

    // Block layout
    for(int mu=0;mu<4;mu++){
      BlockSize[mu] = _BlockSize[mu];
      glatt[mu]     = Layout::lattSize()[mu];
      global_nb[mu] = glatt[mu]/BlockSize[mu];
      if ( global_nb[mu]*BlockSize[mu] !=  Layout::lattSize()[mu] ) {
	QDPIO::cout<<"global BlockSize err"<<endl;
	QMP_barrier();
	exit(0);
      }
      local_nb[mu] = QDP::Layout::subgridLattSize()[mu]/BlockSize[mu];
      if ( local_nb[mu]*BlockSize[mu] !=  QDP::Layout::subgridLattSize()[mu] ) {
	QDPIO::cout<<"local BlockSize err "<< mu << " " << local_nb[mu] << " "<< BlockSize[mu] << " "  <<endl;
	QMP_barrier();
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
    Nsubspace=Nblock*Nvec;

    // Dimensions of subspace blocks (locally)
    LocalNblock4= local_nb[0]*local_nb[1]*local_nb[2]*local_nb[3];
    LocalNblock = NblockS*LocalNblock4;
    LocalNsubspace=LocalNblock*Nvec;

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

    QDPIO::cout << "BlockSize "
		<< BlockSize[0] << "x"
		<< BlockSize[1] << "x"
		<< BlockSize[2] << "x"
		<< BlockSize[3] << "x"
		<< BlockSize[4] << endl;

    QDPIO::cout << "Making integer tables"<<endl;
    QMP_barrier();

    Nball=3*3*3*3*3;// Actually 32 of these are in fact zero
    Asparse.resize(LocalNblock);
    Asparse_check.resize(LocalNblock);

    for(int b=0;b<LocalNblock;b++){   // One per 5d-block
      Asparse[b].resize(Nball);       // One per neigbour
      Asparse_check[b].resize(Nball);       // One per neigbour
      for(int n=0;n<Nball;n++){
	Asparse[b][n].resize(Nvec*Nvec);
	Asparse_check[b][n].resize(Nvec*Nvec);
      }
    }    

    block_id.resize(2);
    local_block_id.resize(2);
    //    ball_cb.resize (2);
    //    ball_id.resize (2);
    //    ball_id_cb.resize (2);
    for(int cb=0;cb<2;cb++){
      //      ball_cb[cb].resize (cbvol);
      //      ball_id[cb].resize(cbvol);
      //      ball_id_cb[cb].resize(cbvol);
      block_id[cb].resize(cbvol);
      local_block_id[cb].resize(cbvol);
    }

    std::vector<std::vector< std::vector<int> > > coordinate(5);
    for(int mu=0;mu<5;mu++){
      coordinate[mu].resize(2);
      coordinate[mu][0].resize(cbvol); 
      coordinate[mu][1].resize(cbvol);
    }
    
    QDPIO::cout << "Mapping blockids"<<endl;
    QMP_barrier();
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

      for(int mu=0;mu<5;mu++){
	coordinate[mu][cb][cbsite] = x[mu]+ncoor[mu]*slatt[mu]; //Global coordinate
	b[mu]=coordinate[mu][cb][cbsite]/BlockSize[mu];         //Global block
	lb[mu]=x[mu]/BlockSize[mu];                         //Local block
      }
      // Needed for constructing little dirac op by picking subvectors
      block_id[cb][cbsite] = b[0]+global_nb[0]*(b[1]+global_nb[1]*(b[2]+global_nb[2]*(b[3]+global_nb[3]*b[4])));
      //      ball_cb[cb][cbsite]  = ((b[0]>>1)+(b[1]>>1)+(b[2]>>1)+(b[3]>>1)+(b[4]>>1))&0x1;
      //      ball_id[cb][cbsite]  = (b[0]&1) + 2*(b[1]&1) + 4*(b[2]&1) + 8*(b[3]&1) + 16*(b[4]&1);
      //      ball_id_cb[cb][cbsite]= ball_id[cb][cbsite] + ball_cb[cb][cbsite]*32;

      local_block_id[cb][cbsite]=lb[0]+local_nb[0]*(lb[1]+local_nb[1]*(lb[2]+local_nb[2]*(lb[3]+local_nb[3]*lb[4])));
    
    }
    QDPIO::cout << "Orthogonalising"<<endl;
    OrthogonaliseSubspace();

    QDPIO::cout << "Halo Init"<<endl;
    HaloInit();

  }

int BfmLittleDiracOperator::GetLocalBlockIdx(int d[5])
{
  return d[0] + local_nb[0]*(d[1]+local_nb[1]*(d[2]+local_nb[2]*(d[3]+local_nb[3]*d[4])));
}
void BfmLittleDiracOperator::GlobalToLocalBlock(int gb[5],int &proc,int lb[5])
{
  multi1d<int> gcoor(4);
  for(int mu=0;mu<5;mu++)   lb[mu]=gb[mu]%local_nb[mu];
  for(int mu=0;mu<4;mu++)   gcoor[mu] = gb[mu]*BlockSize[mu];
  proc = QDP::Layout::nodeNumber(gcoor);  
}
int BfmLittleDiracOperator::ReverseDirection(int mu)
{
  int d[5];
  int rd[5];
  d[0] = mu % 3; mu = mu/3;
  d[1] = mu % 3; mu = mu/3;
  d[2] = mu % 3; mu = mu/3;
  d[3] = mu % 3; mu = mu/3;
  d[4] = mu % 3; 
  for(int dir=0;dir<5;dir++){
    rd[dir] = 2-d[dir];
  }
  mu = rd[4]+3*rd[3]+9*rd[2]+27*rd[1]+81*rd[0];
  return mu;
}

void BfmLittleDiracOperator::GetDelta(int delta[5],int _ballidx,int rev)
{
  int d[5];
  int ballidx=0;
  for(d[0]=-1;d[0]<=1;d[0]++){
    for(d[1]=-1;d[1]<=1;d[1]++){
      for(d[2]=-1;d[2]<=1;d[2]++){
	for(d[3]=-1;d[3]<=1;d[3]++){
	  for(d[4]=-1;d[4]<=1;d[4]++){
	    if ( _ballidx == ballidx ) {
	      if( rev ) for(int mu=0;mu<5;mu++) delta[mu] = -d[mu];
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

int BfmLittleDiracOperator::NeighbourBlockId(int myblock,int _mu, int &me,
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
void BfmLittleDiracOperator::HaloEnd(void)
{
  int Nmsg = tproc.size();
  for(int m=0;m<Nmsg;m++) {
    QMP_free_msghandle(xmit_msghandle[m]);
    QMP_free_msghandle(recv_msghandle[m]);
    QMP_free_msgmem(xmit_msgmem[m]);
    QMP_free_msgmem(recv_msgmem[m]);
  }
}
void BfmLittleDiracOperator::HaloInit(void)
{
  tproc.resize(0);
  fproc.resize(0);
  sendbufs.resize(Nball);
  recvbufs.resize(Nball);
 
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
 
      // NB Longer term can build up face table precomputed
      // using this code sequence and fill the arrays 
      // more efficiently
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
 	rbuf_idx[idx*Nball+mu] = recvbufs[fidx].size();
 	for(int i=0;i<Nvec;i++){
	  recvbufs[fidx].push_back(0);
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
 	sbuf_idx[idx*Nball+mu] = sendbufs[tidx].size();
	sbuf_gather[tidx].push_back(fromblock);
	for(int i=0;i<Nvec;i++){
	  sendbufs[tidx].push_back(0);
	}      
      }
    } //mu
  } // idx

  if ( tproc.size() != fproc.size() ) {
    QDP_error_exit("Halo exchange logic bomb");
  }
 
  // QMP init
  int Nmsg = tproc.size();
  xmit_msgmem.resize(Nmsg);
  recv_msgmem.resize(Nmsg);
  xmit_msghandle.resize(Nmsg);
  recv_msghandle.resize(Nmsg);

  for(int m=0;m<Nmsg;m++){
    xmit_msgmem[m] = QMP_declare_msgmem((void *) &sendbufs[m][0],2*sizeof(double)*sendbufs[m].size()) ;  // send buf
    recv_msgmem[m] = QMP_declare_msgmem((void *) &recvbufs[m][0],2*sizeof(double)*recvbufs[m].size()) ; // receive buf
    
    xmit_msghandle[m] = QMP_declare_send_to(xmit_msgmem[m],tproc[m],0);
    recv_msghandle[m] = QMP_declare_receive_from(recv_msgmem[m],fproc[m],0);
   }
}
void BfmLittleDiracOperator::HaloGetStencilData(int _idx,int _ballidx,
						std::vector<std::complex<double> >& nbr_data,
						std::vector<std::complex<double> > &my_data,
						int rev)
{
  // Locate the neighbour
  int me, from, to;
  int fromblock;
  int mu;

  if ( rev ) mu =ReverseDirection(_ballidx);
  else       mu =_ballidx;

  int bid = rbuf_id [_idx*Nball+mu] ;
  int bidx= rbuf_idx [_idx*Nball+mu] ;

  if (bid==-1){ 
    for(int i=0;i<Nvec;i++){
      nbr_data[i] = my_data[i+bidx];
   }  
  } else { 
    for(int i=0;i<Nvec;i++){
      nbr_data[i] = recvbufs[bid][i+bidx];
    }  
  }
}
void BfmLittleDiracOperator::HaloExchange(std::vector<std::complex<double> >&my_data)
{
  int Nmsg = tproc.size();
 
  //  FlightLog("Filling buffers");
#pragma omp parallel
  {
#pragma omp for
  for(int tidx=0;tidx<sbuf_gather.size();tidx++){
    for(int j=0;j<sbuf_gather[tidx].size();j++){
      for(int i=0;i<Nvec;i++){
	sendbufs[tidx][j*Nvec+i] = my_data[i+sbuf_gather[tidx][j]*Nvec];
      }
    }
  }
  }

  // Do the transfer
  //  FlightLog("Initiating comms");
  for(int m=0;m<Nmsg;m++)  QMP_start(xmit_msghandle[m]);
  for(int m=0;m<Nmsg;m++)  QMP_start(recv_msghandle[m]);
  for(int m=0;m<Nmsg;m++)  QMP_wait (xmit_msghandle[m]);
  for(int m=0;m<Nmsg;m++)  QMP_wait (recv_msghandle[m]);
  //  FlightLog("Completed");
}
void BfmLittleDiracOperator::HaloTest(std::vector<std::complex<double> >&my_data)
{
  QDPIO::cout << "Halo Init"<<endl;
  HaloExchange(my_data);
  QDPIO::cout << "Halo Stencil cross check"<<endl;

  std::vector<std::complex<double> > halo_data(Nvec);
  std::vector<std::complex<double> > ref_data(Nvec);
  for(int idx=0;idx<LocalNblock;idx++){
    for(int mu=0;mu<Nball;mu++){
      int rev=0;
      GetStencilData(idx,mu,ref_data,my_data,rev);
      HaloGetStencilData(idx,mu,halo_data,my_data,rev);
      for(int i=0;i<Nvec;i++){
	if ( halo_data[i]!=ref_data[i] ) {
	  QDPIO::cout << "Stencil Mismatch " << idx<< " mu " << mu << " g "<< real(ref_data[i]) << " b "<< real(halo_data[i]) <<endl;
	}
      }
    }
  }
  QDPIO::cout << "Check complete\n" <<endl;
}

void BfmLittleDiracOperator::GetStencilData(int _idx,int _ballidx,std::vector<std::complex<double> >& nbr_data,std::vector<std::complex<double> >&my_data,int rev)
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
#ifdef REGRESS_TO_SLOW_IMPLEMENTATION
    QMP_barrier();
    if( to == 0 ) 
    {
      fprintf(stdout,"Sending %le from %d:%d to %d \n", real(my_data[fromblock*Nvec]),me,fromblock,to);
      fflush(stdout);
    }
    QMP_barrier();
#endif

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

#ifdef REGRESS_TO_SLOW_IMPLEMENTATION
    QMP_barrier();
    if( me == 0 ) 
      {
      fprintf(stdout,"Received %le from %d:%d to %d \n", real(nbr_data[0]),from,fromblock,me);
      fflush(stdout);
v    }
    QMP_barrier();
#endif
  }

}

void BfmLittleDiracOperator::OrthogonaliseSubspace(void)
{
  // Remove earlier (normalised) vectors
  std::vector<std::complex<double> > bip(LocalNblock);
  std::vector<double> bn(LocalNblock);
  for(int v=0;v<Nvec;v++){
    for(int u=0;u<v;u++){
      //Inner product & remove component
      linop->block_inner(subspace[u],subspace[v],LocalNblock,(double *)&bip[0],&local_block_id[Odd][0]);
      linop->block_zaxpy(subspace[v],subspace[u],subspace[v],-1,(double *)&bip[0],&local_block_id[Odd][0]);
    }
      // Normalise this vector
    linop->block_norm(subspace[v],LocalNblock,(double *)&bn[0],&local_block_id[Odd][0]);
    linop->block_normalise(subspace[v],LocalNblock,(double *)&bn[0],&local_block_id[Odd][0]);
  }
  
  // Orthogonality check is excessive. Remove once confident of code
  for(int v=0;v<Nvec;v++){
  for(int u=0;u<Nvec;u++){
    
    linop->block_inner(subspace[u],subspace[v],LocalNblock,(double *)&bip[0],&local_block_id[Odd][0]);

      double c = (u==v);

      for(int i=0;i<bip.size();i++) { 
	if ( abs(bip[i]-c) > 1.0e-6 ) {  
	  cout << i<<" "<< u<< " " << v << " " << bip[i] << " : " 
	       <<"Doesnt look orthogonal and should be"<<endl;
	}
      }
  }}
  QMP_barrier();
}
void BfmLittleDiracOperator::PickBlock (Fermion_t was,Fermion_t is,int b)
{
  Fermion_t blank = linop->allocFermion();
  linop->master_fill(blank,0.0);
  fflush(stdout);
  linop->block_pick(was,blank,is,b,&block_id[1][0]);//pick on global blockid
  linop->freeFermion(blank);
}

void BfmLittleDiracOperator::PickBallCb(Fermion_t was, Fermion_t is,int ball_id, int ball_cb)
{
  int id_cb = ball_id+ball_cb*32;
  Fermion_t blank = linop->allocFermion();
  linop->master_fill(blank,0.0);
  //  linop->block_pick(was,blank,is,id_cb,&ball_id_cb[1][0]);//pick on global blockid                                                                       
  linop->freeFermion(blank);
}
void BfmLittleDiracOperator::ComputeLittleMatrixColored (void)
{
  Fermion_t zero_t = linop->allocFermion();
  Fermion_t phi_t = linop->allocFermion();
  Fermion_t tmp_t = linop->allocFermion();
  Fermion_t mmp   = linop->allocFermion();
  Fermion_t mp    = linop->allocFermion();

  std::vector<std::complex<double> > phases(LocalNblock);

  QDPIO::cout << "Compute Colored"<<endl;
  linop->master_fill(zero_t,0.0);

  for(int i=0;i<Nvec;i++){

    FlightLog("Computing for vector "<< i<<" ");
    std::vector< std::vector<std::complex<double> > > proj(Nball);

    for(int b=0;b<Nball;b++){  // Loop over momenta (Nball)

      //      QDPIO::cout << "Compute Colored mu "<<b<<endl;
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
      
      for(lb[0]=0;lb[0]<local_nb[0];lb[0]++){
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

      ///////////////////////////////////////////////////////
      // We apply the matrix Nball times using these phases
      ///////////////////////////////////////////////////////
      linop->block_zaxpy(phi_t,subspace[i],zero_t,1.0,(double *)&phases[0],&local_block_id[Odd][0]);
#pragma omp parallel 
      {
	omp_set_num_threads(linop->nthread);
#pragma omp for 
	for(int t=0;t<linop->nthread;t++) {
	  linop->Mprec(phi_t,mp,tmp_t,0);
	  linop->Mprec(mp,mmp,tmp_t,1); 
	}
      }

      ///////////////////////////////////////////////////////
      // Project this into subspace
      ///////////////////////////////////////////////////////
      ProjectToSubspace(mmp,proj[b]);
      for(int lbidx=0;lbidx<LocalNblock;lbidx++) { 
	for(int j=0;j<Nvec;j++){
	  proj[b][lbidx*Nvec+j] = proj[b][lbidx*Nvec+j]*conj(phases[lbidx]);
	}
      }
    }

    QDPIO::cout << "Completed ball momenta"<<endl;
#if 0
    int zero_mom = 1+3+9+27+81;
    int max_mom = 242;
    QDPIO::cout << "XXX Zero mom for idx 0 = "<<real(proj[zero_mom][0])<<" "<<imag(proj[zero_mom][0])<<endl;
    {
      double tmp_mom[5];
      double delta[5];
      double pdotx=0;
      for(int mu=0;mu<5;mu++){
	tmp_mom[mu] = 2*M_PI/global_nb[mu];
      }
      std::complex<double> sum=0;
      for(delta[0]=-1;delta[0]<=1;delta[0]++){
      for(delta[1]=-1;delta[1]<=1;delta[1]++){
      for(delta[2]=-1;delta[2]<=1;delta[2]++){
      for(delta[3]=-1;delta[3]<=1;delta[3]++){
      for(delta[4]=-1;delta[4]<=1;delta[4]++){
	double pdotx=0;
	for(int mu=0;mu<5;mu++) pdotx+=tmp_mom[mu]*delta[mu]; // 11111
	std::complex<double> eip (cos(pdotx),sin(pdotx));
	std::complex<double> check=0;

	int sft[5];
	for(int mu=0;mu<5;mu++) sft[mu]=delta[mu]+1;
	
	int lex = sft[4]+3*sft[3]+9*sft[2]+27*sft[1]+81*sft[0];

	sum = sum + eip*Asparse[0][lex][0];
	
      }}}}}
      std::complex<double> check=proj[max_mom][0];
      QDPIO::cout <<"XXXX check "<<real(sum)<<" "<<imag(sum) << " vs " << real(check)<<" "<<imag(check)<<endl;
    }
    for(int b=0;b<Nball;b++){
      // Could check specific momenta (-1,-1,-1,-1)???
      QDPIO::cout << "XXX mom for idx "<<b<<" = "<<real(proj[b][0])<<" "<<imag(proj[b][0])<<endl;
      QDPIO::cout <<"XXX Asparse "<<b<<" : "<< real(Asparse[0][b][0])<< " " << imag(Asparse[0][b][0])<<endl;
    }
#endif

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
      for(int lbidx=0;lbidx<LocalNblock;lbidx++){
	for(int j=0;j<Nvec;j++){

	for(int mu=0;mu<5;mu++){

	  double pmu=2*M_PI/global_nb[mu];      
	  std::complex<double> pha(cos(pmu),-sin(pmu));
	  std::vector<std::complex<double> > FT(Nball);
	  std::complex<double> a = pha+1.0;
	  std::complex<double> b = (pha-1.0)*(pha-1.0);
	  std::complex<double> ab=a*b;

	  ///////////////////////////////////////////////////////////////////////////////
	  //  pha* 1 pha                 pha/a -pha        pha^2/a
	  //  1    1  1     -> inv ->    -pha  (1+pha^2)   -pha           / [ 1-pha]^2     where a=1+pha
 	  //  pha  1 pha*               pha^2/a -pha       pha/a
	  ///////////////////////////////////////////////////////////////////////////////
	  std::vector< std::vector< std::complex<double> > > imat; // Inverse matrix for this mu.
	  imat.resize(3);
	  for(int i=0;i<3;i++) imat[i].resize(3);

	  imat[0][0] = pha/ab; 	   imat[0][1] =-pha/b;	                  imat[0][2] = pha*pha/ab;
	  imat[1][0] =-pha/b;	   imat[1][1] = (1.0+pha*pha)/b;	  imat[1][2] =-pha/b;
	  imat[2][0] = pha*pha/ab; imat[2][1] =-pha/b;         	          imat[2][2] = pha/ab;
	  
	  int d[5];
	  int nd[5];

	  for(d[0]=0;d[0]<=2;d[0]++){
	  for(d[1]=0;d[1]<=2;d[1]++){
	  for(d[2]=0;d[2]<=2;d[2]++){
	  for(d[3]=0;d[3]<=2;d[3]++){
	  for(d[4]=0;d[4]<=2;d[4]++){

	    int  d_mu;
	    d_mu =   d[4]+3* d[3]+9* d[2]+27* d[1]+81* d[0];
	    for(int i=0;i<5;i++) nd[i]=d[i];
	    FT[d_mu] = 0;
	    
	    for(nd[mu]=0;nd[mu]<=2;nd[mu]++){
	      int nd_mu;
	      nd_mu = nd[4]+3*nd[3]+9*nd[2]+27*nd[1]+81*nd[0];
	      
	      FT[d_mu]+= imat[d[mu]][nd[mu]] * proj[nd_mu][lbidx*Nvec+j];
	      
	    }
	  }}}}}

	  for(int b=0;b<Nball;b++) proj[b][lbidx*Nvec+j] = FT[b];

	}//5-Directions FT
	////////////////////////
	// Copy back solution of system of eqns
	////////////////////////
	for(int b=0;b<Nball;b++)  Asparse[lbidx][b][Nvec*j+i] = proj[b][lbidx*Nvec+j];
	}//j
      }// Block
    }// Scope
  }//Nvec
  linop->freeFermion( zero_t ); 
  linop->freeFermion( phi_t ); 
  linop->freeFermion( tmp_t ); 
  linop->freeFermion( mmp   ); 
  linop->freeFermion( mp    );

}
void BfmLittleDiracOperator::ComputeLittleMatrix (void)
{
  Fermion_t phi_t = linop->allocFermion();
  Fermion_t tmp_t = linop->allocFermion();
  Fermion_t mmp   = linop->allocFermion();
  Fermion_t mp    = linop->allocFermion();
	
  for(int i=0;i<Nvec;i++){
    FlightLog("Computing for vector "<< i<<" ");
    for(int b=0;b<Nblock;b++){            // Loop over global blocks

      // Apply DdagD to phi_i[b] -> mmp
	PickBlock(subspace[i],phi_t,b);
#pragma omp parallel 
	{
	  omp_set_num_threads(linop->nthread);
#pragma omp for 
	  for(int t=0;t<linop->nthread;t++) {
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
void BfmLittleDiracOperator::ProjectToSubspace(Fermion_t vec, std::vector<std::complex<double> > &proj)
{
  std::vector<std::complex<double> > bip(LocalNblock);
  proj.resize(LocalNsubspace);

  for(int v_j=0;v_j<Nvec;v_j++){

     linop->block_inner(subspace[v_j],vec,LocalNblock,(double *)&bip[0],&local_block_id[Odd][0]);

    for(int b_j=0;b_j<LocalNblock;b_j++){
      int j = SubspaceIndex(b_j,v_j);
      proj[j] = bip[b_j];
    }
  }
}
void BfmLittleDiracOperator::PromoteFromSubspace(std::vector<std::complex<double> > &v,Fermion_t prom)
{
  linop->master_fill(prom,0.0);
  std::vector<std::complex<double> > coeff(LocalNblock);
  for(int v_i=0;v_i<Nvec;v_i++){
    for(int b_i=0;b_i<LocalNblock;b_i++){
      int i=SubspaceIndex(b_i,v_i);
      coeff[b_i]=v[i];
    }
    linop->block_zaxpy(prom,subspace[v_i],prom,1.0,(double *)&coeff[0],&local_block_id[Odd][0]);
  }
}

void myzgemv(int N,double * __restrict A, double * __restrict in,double c, double * __restrict out);

void BfmLittleDiracOperator::Apply(std::vector<std::complex<double> > &in,std::vector<std::complex<double> > &out )
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




void BfmLittleDiracOperator::ApplyInverse(std::vector<std::complex<double> > &v,std::vector<std::complex<double> > &vinv)
{
  // Conjugate gradient on A. A is hermitian posdef.
  std::vector<std::complex<double> >  p(LocalNsubspace);
  std::vector<std::complex<double> > Ap(LocalNsubspace);
  std::vector<std::complex<double> >  r(LocalNsubspace);

  for(int i=0;i<LocalNsubspace;i++){
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
  rsq = 1.0e-16*ssq;
 
  FlightLog("Source^2 = " << ssq);
  FlightLog("Target^2 = " << rsq);

  for(int k=1;k<10000;k++){
    c=cp;
    //  FlightLog("ApplyThread");
    ApplyThread(p,Ap);
    //    FlightLog("Linalg");
    d=innerProductReal(p,Ap);
    a=c/d;
#pragma omp parallel 
	{
#pragma omp for 
	  for(int i=0;i<LocalNsubspace;i++){
	    r[i]=r[i]-a*Ap[i];
	  }
	}
  
    cp=norm_vec(r);
    b=cp/c;
#pragma omp parallel 
	{
#pragma omp for 
	  for(int i=0;i<LocalNsubspace;i++){
	    vinv[i]=vinv[i]+a*p[i];
	    p[i]=b*p[i]+r[i];
	  }
	}
	//    FlightLog(" LittleDiracInversion ["<<k<<"] r^2="<<cp <<" ");
    if(cp<rsq) {
      FlightLog("LittleDiracInversion converged after "<<k<<" iterations ");
      Apply(vinv,Ap);
      for(int i=0;i<LocalNsubspace;i++){
	Ap[i]=Ap[i]-v[i];
      }
      FlightLog("LittleDiracInversion true residual "<<norm_vec(Ap));
      return;
    }
  }
  QDP_error_exit("Little Dirac Matrix inversion failed\n");
}
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

void BfmLittleDiracOperator::ApplyThread(std::vector<std::complex<double> > &in,std::vector<std::complex<double> > &out )
{
  // This is tricky to thread because comms is interleaved.                       
  // Easiest way to thread is to serially                
  // fill out a set of vectors for the neighbours, and then    
  // apply the matrix in parallel.                  
 
  //  FlightLog("Halo exchange");
  HaloExchange(in);
  //FlightLog("Stencil multiply");
  for(int mu=0;mu<Nball;mu++){
#pragma omp parallel
    {
#pragma omp for
      for(int b=0;b<LocalNblock;b++){ // Comms hard to thread so serialise
	std::vector<std::complex<double> > nbr_data(Nvec);
	int nodag=0;
	HaloGetStencilData(b,mu,nbr_data,in,nodag);
	double coeff=1.0;
	if ( mu==0)     coeff=0.0;
	myzgemv(Nvec,(double *)&Asparse[b][mu][0],(double *)&nbr_data[0],coeff,(double *)&out[b*Nvec]);
      }
    }
  }
  //  FlightLog("Done all");
}
