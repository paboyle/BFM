#include <chroma.h>
#include <bfm.h>

using namespace Chroma;
typedef LatticeFermion T;
typedef multi1d<LatticeColorMatrix> U;

//#define Printf if ( QMP_is_primary_node() ) printf
#define Printf  printf

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
    for ( int i=0;i<argc;i++)
      Printf("%d %s\n",i,argv[i]);
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
  bfm  dwf;
  dwfa.solver=DWFrb4d;
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
  multi1d<LatticeColorMatrix> u(Nd);  HotSt(u);

  Real M5(1.8);
  Real mq(0.1);

  dwfa.Ls = Ls;
  dwfa.mass = toDouble(mq);
  dwfa.M5   = toDouble(M5);
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

  for(int s=0;s<Ls;s++) gaussian(psi[s]);


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

#undef DEBUG
  // cb is cb of result, 1-cb is cb of input field
  for(int cb=0;cb<2;cb++){

    /*Import this checkerboard of source field to bagel*/
    Printf("Importing psi field cb %d\n",cb);
    dwf.importFermion(psi,psi_h,cb);

    // Fill the other checkerboard of result with noise
    for(int s=0;s<Ls;s++) chi[s] = noise;
    for(int s=0;s<Ls;s++) chi_qdp[s] = chi[s];
 
    for(int dag=0;dag<2;dag++){

      Printf("Checking cb=%d dag=%d\n",cb,dag); 

      dwf.Mooee(psi_h,
                chi_h, 
		dag);
      dwf.exportFermion(chi,chi_h,cb);
 
      // Check the result
      PlusMinus pm;
      if ( dag ) pm = MINUS;
      else pm = PLUS;

      if ( cb == 0 )
        (*M).evenEvenLinOp(chi_qdp,psi,pm);
      else 
        (*M).oddOddLinOp(chi_qdp,psi,pm);
      Double n2 = 0.0;
      for(int s=0;s<Ls;s++) { 
        Double sn2 = norm2(chi[s]-chi_qdp[s]);
	n2+=sn2;
#ifdef DEBUG
        cout << "|| Bagel - QDP || ["<<s<<"] = "<< sn2 << endl;
#endif
      }
      cout << "|| Bagel - QDP || = "<< n2 << endl;

      n2 = 0.0;
      for(int s=0;s<Ls;s++) { 
	Double sn2 = norm2(chi[s]);
        n2+=sn2;
#ifdef DEBUG
        cout << "|| Bagel || ["<<s<<"] = "<< sn2 << endl;
#endif
      }
      cout << "|| Bagel || = "<< n2 << endl;

      n2 = 0.0;
      for(int s=0;s<Ls;s++) {
	Double sn2 = norm2(chi_qdp[s]);
        n2+=sn2;
#ifdef DEBUG
        cout << "|| QDP || ["<<s<<"] = "<< sn2 << endl;
#endif
      }
      cout << "|| QDP || = "<< n2 << endl;
    }
  }
  Printf("Done\n"); 

}







