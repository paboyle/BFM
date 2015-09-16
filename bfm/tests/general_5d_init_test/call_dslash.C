#include <chroma.h>
#include <bfm.h>

using namespace Chroma;
typedef multi1d<LatticeColorMatrix> U;

#define Printf if ( QMP_is_primary_node() ) printf

typedef LatticeFermion T;
typedef multi1d<LatticeFermion> T5;
typedef multi1d<LatticeColorMatrix> U;

#include <bfmcommqmp.h>
#include <bfmcommspi.h>
typedef bfmcommQMP<double> bfm_qmp;

int main (int argc,char **argv )
{
  Chroma::initialize(&argc,&argv);


  /********************************************************
   * Command line parsing
   ********************************************************
   */
#undef COMMANDLINE
#ifdef COMMANDLINE
  if ( argc != 6 ) { 
   Printf("Usage: %s lx ly lz lt Ls\n All must be even\n",argv[0]);
   Printf("argc is %d\n",argc);
   for ( int i=0;i<argc;i++){
      Printf("%d %s\n",i,argv[i]);
   }
      exit(-1);
  }
#endif

  /********************************************************
   * Setup QDP
   ********************************************************
   */
  multi1d<int> nrow(Nd);
#ifdef COMMANDLINE
  nrow[0] = atoi(argv[1]);
  nrow[1] = atoi(argv[2]);
  nrow[2] = atoi(argv[3]);
  nrow[3] = atoi(argv[4]);
  int Ls = atoi(argv[5]);
#else
  nrow[0] = 4;
  nrow[1] = 4;
  nrow[2] = 4;
  nrow[3] = 4;
  int Ls  = 8;
#endif
  Layout::setLattSize(nrow);
  Layout::create();

  int lx = QDP::Layout::subgridLattSize()[0];
  int ly = QDP::Layout::subgridLattSize()[1];
  int lz = QDP::Layout::subgridLattSize()[2];
  int lt = QDP::Layout::subgridLattSize()[3];

  /********************************************************
   * Setup DWF operator
   ********************************************************
   */
  bfmarg  dwfa;
  dwfa.node_latt[0]  = lx;
  dwfa.node_latt[1]  = ly;
  dwfa.node_latt[2]  = lz;
  dwfa.node_latt[3]  = lt;
  multi1d<int> procs = QDP::Layout::logicalSize();
  Printf("%d dim machine\n\t", procs.size());
  for(int mu=0;mu<4;mu++){
    if ( procs[mu]>1 ) {
      dwfa.local_comm[mu] = 0;
    } else { 
      dwfa.local_comm[mu] = 1;
    }
  }

  // Make up a random gauge field.
  multi1d<LatticeColorMatrix> u(Nd); HotSt(u);

  Real M5(1.8);
  Real mq(0.1);

  dwfa.mass = toDouble(mq);
  dwfa.M5   = toDouble(M5);
  dwfa.Csw = 0.0;
  dwfa.precon_5d=0;



  Printf("Kernel : Hw\n");
  Printf("Rep    : Continued Fraction \n");
  Printf("Approx : Zolotarev\n");
  dwfa.Ls = Ls;
  dwfa.solver= HwContFracZolo;
  bfm_qmp dwf_qmp;
  dwf_qmp.eps = 0.05;
  dwf_qmp.init(dwfa);
  dwf_qmp.end();


  Printf("Kernel : Hw\n");
  Printf("Rep    : Partial Fraction \n");
  Printf("Approx : Zolotarev\n");
  dwfa.Ls = Ls+1;
  dwfa.solver= HwPartFracZolo;
  dwf_qmp.init(dwfa);
  dwf_qmp.end();

  //		 HwPartFracTanh, // Never been looked at!
  //             HwContFracTanh, // Never been looked at!

  Printf("Kernel : Hw\n");
  Printf("Rep    : Cayley \n");
  Printf("Approx : Zolotarev\n");
  dwfa.Ls = Ls;
  dwfa.solver= HwCayleyZolo;
  dwf_qmp.init(dwfa);
  dwf_qmp.end();

  Printf("Kernel : H_T\n");
  Printf("Rep    : Cayley \n");
  Printf("Approx : Zolotarev\n");
  dwfa.Ls = Ls;
  dwfa.solver= HtCayleyZolo;
  dwf_qmp.init(dwfa);
  dwf_qmp.end();


  Printf("Kernel : H_w\n");
  Printf("Rep    : Cayley \n");
  Printf("Approx : Tanh\n");
  dwfa.Ls = Ls;
  dwfa.solver= HwCayleyTanh;
  dwf_qmp.init(dwfa);
  dwf_qmp.end();


  Printf("Kernel : H_t\n");
  Printf("Rep    : Cayley \n");
  Printf("Approx : Tanh\n");
  dwfa.Ls = Ls;
  dwfa.solver= HtCayleyTanh;
  dwf_qmp.init(dwfa);
  dwf_qmp.end();


  Printf("Kernel : H_m\n");
  Printf("Rep    : Cayley \n");
  Printf("Approx : Tanh\n");
  dwfa.Ls = Ls;
  dwfa.solver= HmCayleyTanh;
  dwf_qmp.init(dwfa);
  dwf_qmp.end();



}








