#include <chroma.h>
#include <bfm.h>
#include <bfm_qdp.h>

using namespace Chroma;

typedef LatticeFermion T;
typedef multi1d<LatticeFermion> T5;
typedef multi1d<LatticeColorMatrix> U;

#define Printf if ( QMP_is_primary_node() ) printf
//#define Printf printf
#define mybfm bfm_internal
template <class Float>
void compare_result (LatticeFermion A, LatticeFermion B, mybfm<Float> & dwf);
int main (int argc,char **argv )
{
  Chroma::initialize(&argc,&argv);


  /********************************************************
   * Command line parsing
   ********************************************************
   */

  /********************************************************
   * Setup QDP
   ********************************************************
   */
  multi1d<int> nrow(Nd);
  nrow[0] = 32;
  nrow[1] = 32;
  nrow[2] = 32;
  nrow[3] = 64;
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
  mybfm<float> dwf;
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
  multi1d<LatticeFermion> psi(Ls);
  multi1d<LatticeFermion> chi(Ls);
  multi1d<LatticeFermion> chi_qdp(Ls);
#ifndef LEX_DATA
    for(int s=0;s<Ls;s++) gaussian(psi[s]);
#else
  LatticeInteger lexcoord = Layout::latticeCoordinate(0)*1000
                          + Layout::latticeCoordinate(1)*100
                          + Layout::latticeCoordinate(2)*10
                          + Layout::latticeCoordinate(3);
  LatticeColorVector cv = zero;
  pokeColor(cv,lexcoord,0);
  for(int s=0;s<Ls;s++) {
    psi[s]=zero;
    pokeSpin(psi[s],cv,0);
  }
#endif

  /********************************************************
   * QDP Dslash
   ********************************************************
   */
  multi1d<int> bcs(Nd); bcs[0] = bcs[1] = bcs[2] = bcs[3] = 1;
  Handle<FermBC<T,U,U> > fbc(new SimpleFermBC< T, U, U >(bcs));
  Handle<CreateFermState<T,U,U> > cfs( new CreateSimpleFermState<T,U,U>(fbc));
  Handle<FermState<T,U,U> > fs ((*cfs)(u));
  QDPWilsonDslash D(fs);

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

  // cb is cb of result, 1-cb is cb of input field

  /*Import this checkerboard of source field to bagel*/
  // Fill the other checkerboard of result with noise
 
  for(int dag=0;dag<1;dag++){

    for(int s=0;s<Ls;s++) chi[s] = zero;
    for(int s=0;s<Ls;s++) chi_qdp[s] = chi[s];

    for(int cb=0;cb<2;cb++) {

      dwf.importFermion(psi,psi_h,1-cb);
      dwf.importFermion(chi,chi_h,cb);
      StopWatch sw;
      sw.reset();
      sw.start();
      dwf.dslash_generic(psi_h,
			 chi_h, 
			 chi_h, 
			 cb,dag,0,0,0); // Omit dperp for this test
      sw.stop();
      double flops = 1320.*lx*ly*lz*lt*Ls/2;	
      double swt = sw.getTimeInSeconds();
      Printf("%le Mflops %le, %le\n",flops/swt*1.0e-6,flops,swt);
      dwf.exportFermion(chi,chi_h,cb);
    } 

    // Check the result
    PlusMinus pm;
    if ( dag ) pm = MINUS;
    else pm = PLUS;

    for(int s=0;s<Ls;s++) {
      D.apply(chi_qdp[s], psi[s], pm, 0);
      D.apply(chi_qdp[s], psi[s], pm, 1);
    }

    Double n2 = 0.0;
    for(int s=0;s<Ls;s++) { 
      Double sn2 = norm2(chi[s]-chi_qdp[s]);
      n2+=sn2;
      QDPIO::cout <<"  ||diff|| ["<<s<<"] = "<< sn2<<endl;
    }
    QDPIO::cout << "|| Bagel - QDP || = "<< n2 << endl;

    n2 = 0.0;
    for(int s=0;s<Ls;s++) { 
      Double sn2 = norm2(chi[s]);
      n2+=sn2;
    }
    QDPIO::cout << "|| Bagel || = "<< n2 << endl;

    n2 = 0.0;
    for(int s=0;s<Ls;s++) {
      Double sn2 = norm2(chi_qdp[s]);
      n2+=sn2;
    }
    QDPIO::cout << "|| QDP || = "<< n2 << endl;
    for (int s=0;s<1;s++){
      QDPIO::cout << "comparing s = " << s << endl;
      compare_result(chi_qdp[s],chi[s],dwf);
    }
    exit(0);
  }
  Printf("Done\n"); 
  exit(0);
}

template <class Float>
void compare_result (LatticeFermion A, LatticeFermion B, mybfm<Float> & dwf)
{

  multi1d<int> x(4);

  int Nspinco=12;

  Printf("Result vectors \n");
  for ( x[0]=0; x[0]<dwf.node_latt[0];x[0]++ ) { 
    for ( x[1]=0; x[1]<dwf.node_latt[1];x[1]++ ) { 
      for ( x[2]=0; x[2]<dwf.node_latt[2];x[2]++ ) { 
	for ( x[3]=0; x[3]<dwf.node_latt[3];x[3]++ ) { 

    int xx[4];
    xx[0]=x[0];
    xx[1]=x[1];
    xx[2]=x[2];
    xx[3]=x[3];
    for ( int co=0;co<3;co++ ) { 
    for ( int sp=0;sp<4;sp++ ) { 

        int spco = co+3*sp;

	int reim=0;
        Fermion ferm = peekSite(A,x); 
        ColorVector cv = peekSpin(ferm,sp);
        Complex ca = peekColor(cv,co);

        ferm = peekSite(B,x); 
        cv = peekSpin(ferm,sp);
        Complex cb = peekColor(cv,co);

	if ( toDouble(norm2(cb-ca)) > 0.0001 ) { 
	  Printf("site %d %d %d %d co %d sp %d\n",x[0],x[1],x[2],x[3],co,sp);
	  Printf("%le %le\n",toDouble(real(ca)),toDouble(real(cb)));
      	  Printf("%le %le\n",toDouble(imag(ca)),toDouble(imag(cb)));
	} 


    }}


  }}}}

    QMP_barrier();

}





