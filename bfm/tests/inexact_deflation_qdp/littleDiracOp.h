#ifndef _LITTLE_DOP_H_

#include <qdp.h>
#include <bfm.h>
#include <bfm_qdp.h>
#include <bfm_qdp_g5d.h>
#include <bfm_qdp_dwf.h>
#include <chroma.h>
#include <vector>
#include <complex>

class DeflationBlockSet : public SetFunc 
{

public:
  multi1d<int> BlockSize;
  multi1d<int> nb;

  DeflationBlockSet(multi1d<int> _BlockSize) : BlockSize(_BlockSize), nb(5)
  {
    QDPIO::cout<<"Entering DeflationBlockSet "<<endl;
    for(int mu=0;mu<Nd;mu++){
      nb[mu] = Layout::lattSize()[mu]/BlockSize[mu];
      if ( nb[mu]*BlockSize[mu] !=  Layout::lattSize()[mu] ) QDP_error_exit("BlockSize err");
      QDPIO::cout<< Layout::lattSize()[mu] << " " << BlockSize[mu] << " " << nb[mu]<<endl;
    }
    QDPIO::cout << "Made SetFunc"<<endl;
  } 

 public: int operator() (const multi1d<int>& coordinate) const
  {
    multi1d<int> b(4);
    for(int mu=0;mu<4;mu++){
      b[mu] = coordinate[mu]/BlockSize[mu];
    }
    int block = b[0]+nb[0]*(b[1]+nb[1]*(b[2]+nb[2]*b[3]));
    if(0){
    printf("MAP %d %d %d %d -> %d\n",
	   coordinate[0],
	   coordinate[1],
	   coordinate[2],
	   coordinate[3],
	   block
	   );
    fflush(stdout);
    }
    if ( block > numSubsets() ) { 
      printf("Bother said Pooh\n");
      fflush(stdout);
    }
    return block;
  }
  int numSubsets() const {return nb[0]*nb[1]*nb[2]*nb[3];}
};


class LittleDiracOperator { 
public:
  multi1d<int> BlockSize;
  int Nblock4;
  int NblockS;
  int Nblock;
  int Nvec;
  int Nsubspace;
  int Nball;
  int NmatrixNonZero;
  int N5;
  multi1d<multi1d<LatticeFermion> > subspace;
  Set Blocks;

  std::vector< std::vector< std::complex<double> > >  A;

  // Sparse matrix implementation [block][mu][Nvec x Nvec]
  std::vector< std::vector< std::vector< std::complex<double> > > > Asparse;
  std::vector< std::vector< int > > nbr;
  
  double innerProductReal(std::vector<std::complex<double> > &v1,
			  std::vector<std::complex<double> > &v2){
    double ir = 0;
    for(int i=0;i<v1.size();i++){
      ir = ir + real(v1[i])*real(v2[i]) + imag(v1[i])*imag(v2[i]);
    }
    return ir;
  }
  double norm_vec(std::vector<std::complex<double> > &v) {
    double nrm=0.0;
    for(int i=0;i<v.size();i++){
      nrm += real(v[i])*real(v[i])+imag(v[i])*imag(v[i]);
    }
    return nrm;
  }
 LittleDiracOperator(multi1d<int> _BlockSize,multi1d<multi1d<LatticeFermion> > &_subspace) : 
  N5(_subspace[0].size()),BlockSize(_BlockSize), subspace(_subspace)
  {
    QDPIO::cout << "Constructor entered"<<endl;
    Nvec   = subspace.size();
    QDPIO::cout << " building subspace from " << Nvec<<" vectors"<<endl;
    Blocks.make(DeflationBlockSet(BlockSize));
    QDPIO::cout << "Made block object"<<endl;
    // Dimensions of subspace blocks
    Nblock4= Blocks.numSubsets();
    QDPIO::cout << Nblock4<< " 4d blocks"<<endl;
    NblockS= N5/BlockSize[4];  
    if ( NblockS*BlockSize[4]!=N5 ) {
      QDPIO::cout << "BlockSize error in fifth tdim" << N5 << " " << NblockS << " " << BlockSize[4] << endl;
      exit(0);
      QMP_barrier();
      QDP_error_exit("BlockSize err");
    }
    QDPIO::cout << NblockS<< " 5fth dim partitions"<<endl;
    Nblock =NblockS*Nblock4;
    QDPIO::cout << Nblock << " 5d blocks"<<endl;
    Nsubspace=Nblock*Nvec;

    QDPIO::cout << "Resizing matrices"<<endl;
    A.resize(Nsubspace);
    for(int i=0;i<A.size();i++){
      A[i].resize(Nsubspace);
    }

    QDPIO::cout << "Sparse version"<<endl;
    Nball=3*3*3*3*3;// Actually 32 of these are in fact zero
    NmatrixNonZero = Nball*Nblock;
    Asparse.resize(Nblock);
    nbr.resize(Nblock);
    for(int b=0;b<Nblock;b++){   // One per 5d-block
      Asparse[b].resize(Nball); // One per neigbour
      nbr[b].resize(Nball);
      for(int n=0;n<Nball;n++){
	Asparse[b][n].resize(Nvec*Nvec);
      }
    }    
    
    QDPIO::cout << "Resized sparse"<<endl;
    QMP_barrier();

    multi1d<int> grid(5);
    multi1d<int> c(5);
    multi1d<int> delta(5);
    DeflationBlockSet B (BlockSize);

    for(int mu=0;mu<4;mu++) grid[mu] = Layout::lattSize()[mu];
    grid[4]=N5;

    for( c[0]=0;c[0]<grid[0] ;c[0]+=BlockSize[0]){
    for( c[1]=0;c[1]<grid[1] ;c[1]+=BlockSize[1]){
    for( c[2]=0;c[2]<grid[2] ;c[2]+=BlockSize[2]){
    for( c[3]=0;c[3]<grid[3] ;c[3]+=BlockSize[3]){
    for( c[4]=0;c[4]<grid[4] ;c[4]+=BlockSize[4]){

	int ballidx=0;
	for(delta[0]=-1;delta[0]<=1;delta[0]++){
	for(delta[1]=-1;delta[1]<=1;delta[1]++){
	for(delta[2]=-1;delta[2]<=1;delta[2]++){
	for(delta[3]=-1;delta[3]<=1;delta[3]++){
	for(delta[4]=-1;delta[4]<=1;delta[4]++){
	  
	  multi1d<int> neigh=c;
	  for(int mu=0;mu<5;mu++){
	    neigh[mu]=c[mu]+delta[mu]*BlockSize[mu];
	    if ( neigh[mu]<0 ) neigh[mu]+=grid[mu];
	    if ( neigh[mu]>=grid[mu] ) neigh[mu]-=grid[mu];
	  }

	  int me  = B(c);
	  int him = B(neigh);

	  me = me+Nblock4*c[4]/BlockSize[4];
	  him=him+Nblock4*neigh[4]/BlockSize[4];

	  nbr[me][ballidx]=him;
	  ballidx++;
	}}}}}
    }}}}}

    QDPIO::cout << "Orthogonalising"<<endl;
    OrthogonaliseSubspace();
  }
  ~LittleDiracOperator() {};

  int SubspaceIndex(int _block4,int _blockS,int _vec) {    return _block4*Nvec + Nblock4*_blockS*Nvec +  _vec ; }
  void OrthogonaliseSubspace(void);
  void PickBlock (const LatticeFermion &was, LatticeFermion &is,int b);
  void PickBlock (const multi1d<LatticeFermion> &was, multi1d<LatticeFermion> &is,int b,int s);
  void ComputeLittleMatrix (    bfm_qdp<double> &dwf);
  void ProjectToSubspace(multi1d<LatticeFermion> &vec, std::vector<std::complex<double> > &proj);
  void PromoteFromSubspace( std::vector<std::complex<double> > &v, multi1d<LatticeFermion> &prom);
  void Apply(std::vector<std::complex<double> > &,std::vector<std::complex<double> > & );
  void ApplySparse(std::vector<std::complex<double> > &,std::vector<std::complex<double> > & );
  void ApplyInverse(std::vector<std::complex<double> > &,std::vector<std::complex<double> >&);
  int SubspaceDimension(void){return Nsubspace;};

  // Tests: project a linear combination of subspace vecs to subspace & repromote  //        should be compete (identity).
  //
  //        phi^dag DdagD phi = |Dphi|^2 with phi a subspace vector
  //        should be equal to Project/Apply/Promote + inner product
  // 
  //        Subspace inverse of D acting on a vector in subspace should give identity

};

#endif
