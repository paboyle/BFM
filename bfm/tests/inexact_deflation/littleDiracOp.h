#ifndef _LITTLE_DOP_H_

#include <qdp.h>
#include <qmp.h>
#include <bfm.h>
#include <bfm_qdp.h>
#include <bfm_qdp_g5d.h>
#include <bfm_qdp_dwf.h>
#include <chroma.h>
#include <vector>
#include <complex>

class BfmLittleDiracOperator { 
public:


  int BlockSize[5];

  int NblockS;
  int Nvec;
  int Nball;
  int N5;

  int LocalVol;
  int local_nb[5];
  int LocalNblock4;
  int LocalNblock;
  int LocalNsubspace;

  int glatt[5];
  int slatt[5];
  int ncoor[5];

  int global_nb[5];
  int Nblock4;
  int Nblock;
  int Nsubspace;

  bfm_qdp<double> *linop;
  Fermion_t *subspace;

  // Sparse matrix implementation [block][mu][Nvec x Nvec]
  std::vector< std::vector< std::vector< std::complex<double> > > > Asparse;
  std::vector< std::vector< std::vector< std::complex<double> > > > Asparse_check;
  std::vector<std::vector<int> > block_id;
  std::vector<std::vector<int> > local_block_id;
  //  std::vector<std::vector<int> > ball_cb;
  //  std::vector<std::vector<int> > ball_id;
  //  std::vector<std::vector<int> > ball_id_cb;

  //////////////////////////////////////////////
  // Faster little dirac op implementation
  //////////////////////////////////////////////
  std::vector< int > tproc; // Processor we must send to and from
  std::vector< int > fproc;
  std::vector< std::vector<std::complex<double> > > sendbufs; // data we must send and receive
  std::vector< std::vector<std::complex<double> > > recvbufs;
  std::vector<QMP_msgmem_t> xmit_msgmem;  // QMP data structures for the sending
  std::vector<QMP_msgmem_t> recv_msgmem;
  std::vector<QMP_msghandle_t> xmit_msghandle;
  std::vector<QMP_msghandle_t> recv_msghandle;

  std::vector< int > sbuf_id; // Tables mapping local block + mu combination to buffer and offset
  std::vector< int > sbuf_idx;
  std::vector< std::vector<int> > sbuf_gather;
  std::vector< int > rbuf_id;
  std::vector< int > rbuf_idx;

  void HaloExchange(std::vector<std::complex<double> >&my_data);
  void HaloInit(void);
  void HaloEnd(void);
  void HaloTest(std::vector<std::complex<double> >&my_data);
  void HaloGetStencilData(int _idx,int _ballidx,std::vector<std::complex<double> >& nbr_data,std::vector<std::complex<double> > &my_data,int rev);
  //////////////////////////////////////////////
  
  double innerProductReal(std::vector<std::complex<double> > &v1,
			  std::vector<std::complex<double> > &v2){
    double ir = 0;
    for(int i=0;i<v1.size();i++){
      ir = ir + real(v1[i])*real(v2[i]) + imag(v1[i])*imag(v2[i]);
    }
    QMP_sum_double_array(&ir,1);
    return ir;
  }
  double norm_vec(std::vector<std::complex<double> > &v) {
    double nrm=0.0;
    for(int i=0;i<v.size();i++){
      nrm += real(v[i])*real(v[i])+imag(v[i])*imag(v[i]);
    }
    QMP_sum_double_array(&nrm,1);
    return nrm;
  }


 BfmLittleDiracOperator(int _N5,
		        int _Nvec,
			int _BlockSize[5],
			Fermion_t *_subspace,
			bfm_qdp<double> * _linop) ;
  ~BfmLittleDiracOperator() {
    QDPIO::cout << "Halo End"<<endl;
    HaloEnd();
  };

  int NeighbourBlockId(int myblock,int _mu, int &me,int &from,int &xblock,int &to,int rev);
  int NeighbourGlobalBlockId(int localBlock,int mu);
  int GetLocalBlockIdx(int delta[5]);
  void GetDelta(int delta[5],int mu,int rev);
  int ReverseDirection(int mu);
  void GlobalToLocalBlock(int gb[5],int &proc,int lb[5]);
  void GetStencilData(int _idx,int _ballidx,std::vector<std::complex<double> >& nbr_data,std::vector<std::complex<double> > &my_data,int rev);
  int SubspaceIndex(int _block,int _vec) {    return _block*Nvec+ _vec ; }
  void OrthogonaliseSubspace(void);
  void PickBlock (Fermion_t was, Fermion_t is,int b);
  void PickBallCb(Fermion_t was, Fermion_t is,int ball_id, int ball_cb);
  void ComputeLittleMatrix (void);
  void ComputeLittleMatrixColored (void);
  void ProjectToSubspace(Fermion_t, std::vector<std::complex<double> > &proj);
  void PromoteFromSubspace(std::vector<std::complex<double> > &v, Fermion_t prom);
  void Apply(std::vector<std::complex<double> > &,std::vector<std::complex<double> > & );
  void ApplyThread(std::vector<std::complex<double> > &,std::vector<std::complex<double> > & );
  void ApplyInverse(std::vector<std::complex<double> > &,std::vector<std::complex<double> >&);
  int SubspaceDimension(void){return Nsubspace;};
  int SubspaceLocalDimension(void){return LocalNsubspace;};

  void Pleft (Fermion_t in,Fermion_t out);
  void Pright(Fermion_t in,Fermion_t out);
  void OneMinusPleft (Fermion_t in,Fermion_t out);
  void OneMinusPright(Fermion_t in,Fermion_t out);
  void PleftM (Fermion_t in,Fermion_t out);
  void SolvePleftM(Fermion_t src,Fermion_t sol);


  // Tests: project a linear combination of subspace vecs to subspace & repromote  //        should be compete (identity).
  //
  //        phi^dag DdagD phi = |Dphi|^2 with phi a subspace vector
  //        should be equal to Project/Apply/Promote + inner product
  // 
  //        Subspace inverse of D acting on a vector in subspace should give identity

};




#endif
