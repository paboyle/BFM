#include <chroma.h>
#include <bfm.h>
#include <bfm_qdp.h>

#include <bfm_su3.h>

using namespace Chroma;

typedef LatticeFermion T;
typedef multi1d<LatticeFermion> T5;
typedef multi1d<LatticeColorMatrix> U;

#define Printf if ( QMP_is_primary_node() ) printf

int main (int argc,char **argv )
{
  Chroma::initialize(&argc,&argv);

  /********************************************************
   * Setup QDP
   ********************************************************
   */
  multi1d<int> nrow(Nd);
  nrow[0] = 8;
  nrow[1] = 8;
  nrow[2] = 8;
  nrow[3] = 8;
  int Ls  = 4;

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
  bfmarg dwfa;
  bfmsu3<float> dwf;
  dwfa.solver = DWF;
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

#undef LEX_DATA
  // Make up a random gauge field.
  multi1d<LatticeColorMatrix> u(Nd); 
#ifndef LEX_DATA
  HotSt(u);
#else
  for(int mu=0;mu<4;mu++) u[mu] = Real(1.0);
  //  pokeColor(u[3],Real(1.0),0,2);
#endif
  dwfa.Ls = Ls;
  dwfa.mass = 0.0;
  dwfa.Csw = 0.0;
  dwfa.M5  = 0.0;
  dwfa.precon_5d=1;

  Printf("Initialising bfm operator\n");
  dwf.init(dwfa);

  /********************************************************
   * Gaussian source and result vectors
   ********************************************************
   */
  Matrix_t u_t[4];
  for(int i=0;i<4;i++) u_t[i] = dwf.allocMatrix();
  for(int i=0;i<4;i++) dwf.importMatrix(u[i],u_t[i]);

  QDPIO::cout << "Shifting"<<endl;
  Matrix_t r = dwf.allocMatrix();
  dwf.CovariantShift(r,u_t[0],u_t[1],0,0);

  LatticeColorMatrix check = u[0] * SHIFT(u[1],FORWARD,0);
  LatticeColorMatrix bagel;
  dwf.exportMatrix(bagel,r);

  QDPIO::cout<<"Norm2 diff = "<< norm2(check-bagel)<<endl;

  exit(0);
}





