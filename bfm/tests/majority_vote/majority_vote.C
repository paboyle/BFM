#include <unistd.h>
//#include <chroma.h>
#include <qdp.h>
//#include <actions/ferm/invert/syssolver_linop_cg_array.h>
#include <bfm.h>
#include <bfm_qdp.h>

typedef double Float;
typedef commQMPbagel1<double> bfm_t;
//
// Majority vote incarnation overrides the global sum  
// to instead broadcast result on node zero and do bit identical compare
//
// Any nodes that dribble a bit will print "Oops, I did it again!".  
//
//
class majorityVote : public bfm_t {
public:
  void IdentifyNode(void) { 
  };
  virtual void comm_gsum(double &val) {
    double local = val;
    QMP_broadcast(&val,sizeof(val));
    if ( local != val ) printf("Oops, I did it again!\n");
    IdentifyNode();
  };

  void randGauge(void)
  {
  for(int dir=0;dir<4;dir++){
  int x[4];
  int Ndircoco=36; /*4 directions stored only*/
  int Ncoco=9;

  for ( x[3]=0;x[3]<this->node_latt[3];x[3]++ ) { 
  for ( x[2]=0;x[2]<this->node_latt[2];x[2]++ ) { 
  for ( x[1]=0;x[1]<this->node_latt[1];x[1]++ ) { 
  for ( x[0]=0;x[0]<this->node_latt[0];x[0]++ ) { 
     
      for ( int coco=0;coco<9;coco++ ) { 
      for ( int reim=0;reim<2;reim++ ) { 

        Float * bagel = this->u;
	int bbase = dir*9;
        int bidx = this->bagel_idx(x,reim,coco+bbase,Ndircoco,0);
        bagel[bidx] = drand48();

      }}

  }}}}
  }}
  void randFermion(Fermion_t x_t){

    double *x = (double *) x_t;
    int nspinco = 12;
    double *xx;
    double nrm = 0.0;

    for(int site=0;site<simd_cbvol;site++){
      for(int s=0;s<Ls;s++) {
     
	xx = &x[s*nspinco*nsimd*2   + site*nspinco*nsimd*2*Ls];

	for(int spinco=0;spinco<nspinco;spinco++){
	  for(int simd_el=0;simd_el<nsimd;simd_el++){
	    for(int reim=0;reim<2;reim++){
	      int idx = reim + simd_el*2 + spinco*nsimd*2;
	      xx[idx] = (drand48()-0.5)*2; //[-1,1] uniform
	    }}}
      }
    }  

  };
};


typedef LatticeFermion T;
typedef multi1d<LatticeFermion> T5;
typedef multi1d<LatticeColorMatrix> U;

//#define Printf if ( QMP_is_primary_node() ) printf
#define Printf printf

int main (int argc,char **argv )
{
  QDP_initialize(&argc,&argv);


  /********************************************************
   * Command line parsing
   ********************************************************
   */
  if ( argc != 6 ) { 
   Printf("Usage: %s lx ly lz lt Ls\n All must be even\n",argv[0]);
   Printf("argc is %d\n",argc);
   for ( int i=0;i<argc;i++){
      Printf("%d %s\n",i,argv[i]);
   }
      exit
(-1);

  }

  /********************************************************
   * Setup QDP
   ********************************************************
   */
  multi1d<int> nrow(Nd);
  nrow[0] = atoi(argv[1]);
  nrow[1] = atoi(argv[2]);
  nrow[2] = atoi(argv[3]);
  nrow[3] = atoi(argv[4]);

  Layout::setLattSize(nrow);
  Layout::create();

  int lx = QDP::Layout::subgridLattSize()[0];
  int ly = QDP::Layout::subgridLattSize()[1];
  int lz = QDP::Layout::subgridLattSize()[2];
  int lt = QDP::Layout::subgridLattSize()[3];
  int Ls = atoi(argv[5]);


  Real M5 = 1.8;
  Real mass = 0.00;
  Real mq = mass;


  /********************************************************
   * Setup DWF operator
   ********************************************************
   */
  bfmarg dwfa;
  majorityVote  dwf;

  dwfa.node_latt[0]  = lx;
  dwfa.node_latt[1]  = ly;
  dwfa.node_latt[2]  = lz;
  dwfa.node_latt[3]  = lt;

  multi1d<int> procs = QDP::Layout::logicalSize();
  Printf("%d dim machine\n\t", procs.size());
  for(int mu=0;mu<4;mu++){
    Printf("%d ", procs[mu]);
    //    if ( procs[mu]>1 ) {
    //      dwfa.local_comm[mu] = 0;
    //    } else { 
      dwfa.local_comm[mu] = 1;
    //    }
  }
  Printf("\nLocal comm = ");
  for(int mu=0;mu<4;mu++){
    Printf("%d ", dwfa.local_comm[mu]);
  }
  Printf("\n");
  
  dwfa.precon_5d = 1;
  dwfa.Ls   = Ls;
  dwfa.M5   = toDouble(M5);
  dwfa.mass = toDouble(mq);
  dwfa.Csw  = 0.0;
  dwfa.max_iter = 1000;
  dwfa.residual = 1.e-8;
  Printf("Initialising bfm operator\n");
  dwf.init(dwfa);

  /********************************************************
   * Bagel internal single checkerboard vectors
   ********************************************************
   */
   Fermion_t psi_h[2];
   psi_h[0] = dwf.allocFermion();
   psi_h[1] = dwf.allocFermion();
   Fermion_t chi_h[2];
   chi_h[0] = dwf.allocFermion();
   chi_h[1] = dwf.allocFermion();

  dwf.randGauge();

  /*Import this checkerboard of source field to bagel*/
  Printf("Importing psi field cb %d\n",1);
  dwf.randFermion(psi_h[0]);
  dwf.randFermion(psi_h[1]);
 
  // Fill the other checkerboard of result with noise
  dwf.randFermion(chi_h[0]);
  dwf.randFermion(chi_h[1]);

#define NITER 100
  for(int i=0;i<NITER;i++) {
    dwf.axpy(chi_h[0],psi_h[0],psi_h[0],0.0)  ;
    dwf.axpy(chi_h[1],psi_h[0],psi_h[0],0.0)  ;
    dwf.CGNE(chi_h,psi_h);
  }

  Printf("Done\n"); 

}








