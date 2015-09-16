#include <chroma.h>
#include <actions/ferm/invert/syssolver_linop_cg_array.h>
#include <bfm.h>
#include <bfm_qdp.h>
#include <bfm_qdp_dwf.h>

typedef bfm bfm_t;
//typedef commQMPbagel1<double> bfm_t;

using namespace Chroma;

typedef LatticeFermion T;
typedef multi1d<LatticeFermion> T5;
typedef multi1d<LatticeColorMatrix> U;

//#define Printf if ( QMP_is_primary_node() ) printf
#define Printf printf

int main (int argc,char **argv )
{
  Chroma::initialize(&argc,&argv);


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

  /********************************************************
   * Setup DWF operator
   ********************************************************
   */
  bfmarg dwfa;
  bfm_t  dwf;
  dwfa.solver=DWF;
  dwfa.node_latt[0]  = lx;
  dwfa.node_latt[1]  = ly;
  dwfa.node_latt[2]  = lz;
  dwfa.node_latt[3]  = lt;

  multi1d<int> procs = QDP::Layout::logicalSize();
  Printf("%d dim machine\n\t", procs.size());
  for(int mu=0;mu<4;mu++){
    Printf("%d ", procs[mu]);
    if ( procs[mu]>1 ) {
      dwfa.local_comm[mu] = 0;
    } else { 
      dwfa.local_comm[mu] = 1;
    }
  }
  Printf("\nLocal comm = ");
  for(int mu=0;mu<4;mu++){
    Printf("%d ", dwfa.local_comm[mu]);
  }
  Printf("\n");
  
  multi1d<int> ncoor = QDP::Layout::nodeCoord();

  Real M5(1.8);
  Real mq(0.1);


  /********************************************************
   * Gaussian gauge field
   ********************************************************
   */

  multi1d<LatticeColorMatrix>  u(Nd); HotSt(u);

  /********************************************************
   * Gaussian wall source and result vectors
   ********************************************************
   */
  LatticePropagator psi;
  LatticePropagator chi = zero;
  LatticePropagator chi_qdp = zero;

  gaussian(psi);
  psi = where(Layout::latticeCoordinate(Nd-1) == 0, psi,LatticePropagator(zero)); 

  /********************************************************
   * Import gauge field to BAGEL
   ********************************************************
   */

  Real residual[1]; residual[0] = 1.0e-8;
  int  max_iter[1]; max_iter[0] = 1000;
  dwf_CG<double,double>(chi,psi,u,"propBagel",Ls,false,mq,M5,1,residual,max_iter);

  /********************************************************
   * QDP Linop
   ********************************************************
   */
   multi1d<int> bcs(Nd);
   bcs[0] = bcs[1] = bcs[2] = bcs[3] = 1;

   Handle< FermBC<T,U,U> > fbc(new SimpleFermBC< T, U, U >(bcs));
   Handle<CreateFermState<T,U,U> > cfs( new CreateSimpleFermState<T,U,U>(fbc));
   EvenOddPrecDWFermActArray  S_f(cfs, M5, mq, Ls);
   Handle< FermState<T,U,U> > fs( S_f.createState(u) );

   Handle< EvenOddPrecLinearOperatorArray<T,U,U> > M(S_f.precLinOp(fs,mq));
   Handle< LinearOperatorArray<T> > HM(S_f.precLinOp(fs,mq));

  // Check the result
  SysSolverCGParams invParam;
  invParam.RsdCG        = 1.0e-8;
  invParam.RsdCGRestart = 1.0e-8;
  invParam.MaxCG        = 1000;
  invParam.MaxCGRestart = 1000;

  Handle<LinOpSystemSolverArray<T> > Hsolver = 
    new LinOpSysSolverCGArray<T>(HM,invParam);
  Handle<SystemSolverArray<T> > solver =
    new PrecFermAct5DQprop<T,U,U>(M,Hsolver);

  int n_count;
  XMLFileWriter xml_qdp("QDPprop.xml");
  dwf_quarkProp4(chi_qdp,
		 xml_qdp,
		 psi,
		 0,Nd-1,
		 solver,
		 fs,
		 mq,
		 n_count);

  Double n2;
  n2 = norm2(chi-chi_qdp);
  QDPIO::cout << "|| Bagel - QDP || = "<< n2 << endl;
  

  n2 = norm2(chi);
  QDPIO::cout << "|| Bagel || = "<< n2 << endl;
  
  n2 = norm2(chi_qdp);
  QDPIO::cout << "|| QDP || = "<< n2 << endl;

  Printf("Done\n"); 

}








