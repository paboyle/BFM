#include <chroma.h>
#include <bfm.h>

using namespace Chroma;

#define Printf if ( QMP_is_primary_node() ) printf
#define Printf printf
typedef LatticeFermion T;
typedef multi1d<LatticeColorMatrix> U;

int main (int argc,char **argv )
{
  Chroma::initialize(&argc,&argv);


  /********************************************************
   * Command line parsing
   ********************************************************
   */
  if ( argc != 6 ) { 
   printf("Usage: %s lx ly lz lt Ls\n All must be even\n",argv[0]);
   printf("argc is %d\n",argc);
    for ( int i=0;i<argc;i++)
      printf("%d %s\n",i,argv[i]);
    exit(-1);
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
  int Ls = atoi(argv[5]);

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
  dwfa.solver = WilsonFermion;
  bfm  dwf;

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

  dwfa.Ls = Ls;
  dwfa.mass = 0.0;
  dwfa.M5  = 0.0;
  dwfa.Csw = 0.0;
  dwfa.precon_5d=0;

  Printf("Initialising bfm operator\n");
  dwf.init(dwfa);


  /********************************************************
   * Gaussian source and result vectors
   ********************************************************
   */
  multi1d<LatticeFermion> psi(Ls);
  multi1d<LatticeFermion> chi(Ls);
  multi1d<LatticeFermion> chi_qdp(Ls);
  multi1d<LatticeFermion> psi_exp(Ls);

  for(int s=0;s<Ls;s++) gaussian(psi[s]);

  /********************************************************
   * QDP Dslash
   ********************************************************
   */
#ifdef CHROMA_3
  multi1d<int> bcs(Nd); bcs[0] = bcs[1] = bcs[2] = bcs[3] = 1;
  Handle<FermBC<T,U,U> > fbc(new SimpleFermBC< T, U, U >(bcs));
  Handle<CreateFermState<T,U,U> > cfs( new CreateSimpleFermState<T,U,U>(fbc));
  Handle<FermState<T,U,U> > fs ((*cfs)(u));
  QDPWilsonDslash D(fs);
#else
  QDPWilsonDslash D(u);
#endif

  /********************************************************
   * Bagel internal single checkerboard vectors
   ********************************************************
   */
  Fermion_t psi_h = dwf.allocFermion();
  Fermion_t chi_h = dwf.allocFermion();

  /********************************************************
   * Import gauge field to BAGEL
   ********************************************************
   */
  dwf.importGauge(u);


  LatticeFermion noise; gaussian(noise);

  // Fill the other checkerboard of result with noise
  for(int s=0;s<Ls;s++) chi[s] = zero;
  for(int s=0;s<Ls;s++) chi_qdp[s] = chi[s];

  for(int dag=0;dag<2;dag++){

   /*Import this checkerboard of source field to bagel*/
   // cb is cb of result, 1-cb is cb of input field
    dwf.importFermion(psi,psi_h,1);
    dwf.dslash(psi_h,
               chi_h, 
 	       0,dag);
    dwf.exportFermion(chi,chi_h,0);

    dwf.importFermion(psi,psi_h,0);
    dwf.dslash(psi_h,
               chi_h, 
	       1,dag);
    dwf.exportFermion(chi,chi_h,1);
 
    // Check the result
    PlusMinus pm;
    if ( dag ) pm = MINUS;
    else pm = PLUS;

    for(int s=0;s<Ls;s++){
      D.apply(chi_qdp[s], psi[s], pm, 0);
      D.apply(chi_qdp[s], psi[s], pm, 1);
    }

    Double n2 = 0.0;
    for(int s=0;s<Ls;s++) { 
      Double sn2 = norm2(chi[s]-chi_qdp[s]);
      n2+=sn2;
      cout << "|| Bagel - QDP ||_s = "<< sn2 << endl;
    }
    cout << "|| Bagel - QDP || = "<< n2 << endl;

    n2 = 0.0;
    for(int s=0;s<Ls;s++) { 
      Double sn2 = norm2(chi[s]);
      n2+=sn2;
    }
    cout << "|| Bagel || = "<< n2 << endl;

    n2 = 0.0;
    for(int s=0;s<Ls;s++) {
      Double sn2 = norm2(chi_qdp[s]);
      n2+=sn2;
    }
    cout << "|| QDP || = "<< n2 << endl;
  }
  printf("Done\n"); 

}







