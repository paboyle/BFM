#include <chroma.h>

typedef double Float;

#define Printf printf

#include <bfm.h>

typedef LatticeFermion T;
typedef multi1d<LatticeColorMatrix> U;
typedef bfm bfm_t;

using namespace Chroma;
void compare_result (LatticeFermion A, LatticeFermion B, bfm_t & dwf);

int main (int argc,char **argv )
{
  Chroma::initialize(&argc,&argv);

  if ( argc != 5 ) { 
   Printf("Usage: %s lx ly lz lt \n All must be even\n",argv[0]);
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
  bfm_t  dwf;
  dwfa.solver = WilsonNN;
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
  // Make up a gaussian source and a zero result vector
  multi1d<LatticeColorMatrix> u(Nd); 
  HotSt(u);
#if 0
  for(int mu=0;mu<4;mu++){
    u[mu]=1.0;
  }
#endif

  LatticeFermion psi;  
  LatticeFermion chi;
  LatticeFermion chi_qdp;
  gaussian(psi);

#if 0
  LatticeInteger lexcoord = Layout::latticeCoordinate(0)*1000
                          + Layout::latticeCoordinate(1)*100
                          + Layout::latticeCoordinate(2)*10
                          + Layout::latticeCoordinate(3);
  LatticeColorVector cv = zero;
  pokeColor(cv,lexcoord,0);
  psi=zero;
  pokeSpin(psi,cv,0);
#endif

  dwfa.Ls = 1;
  dwfa.mass = 0.0;
  dwfa.Csw = 0.0;

  Printf("Initialising next nearest neigbour bfm operator\n");
  dwf.init(dwfa);

  Printf("bfm operator using vector length %d\n",dwf.simd());
  Fermion_t psi_h = dwf.allocFermion();
  Fermion_t chi_h = dwf.allocFermion();
  dwf.importGaugeNN(u);


  // cb is cb of result and input field for Eguchi Kawamoto term
  for(int dag=0;dag<1;dag++){
    // Fill the other checkerboard.
    chi = zero;
    chi_qdp = zero;
    for(int cb=0;cb<2;cb++){

      /*Import this checkerboard of QDP fields to bagel*/
      Printf("Importing psi field cb %d\n",cb);


      dwf.importFermion(psi,psi_h,cb);
      dwf.dslash(psi_h,
		 chi_h, 
		 cb,dag);
      dwf.exportFermion(chi,chi_h,cb);

    }
    
    Printf("Checking dag=%d\n",dag); 

    //QDP++ implementation
    LatticeFermion fm;
    LatticeFermion tmp;
    fm = zero;

    if ( dag==0 ){
      mu=0;
      fm=fm+spinReconstructDir0Minus(u[mu]*shift(u[mu]*shift(spinProjectDir0Minus(psi), FORWARD, mu),FORWARD,mu));
      fm=fm+spinReconstructDir0Plus (shift(adj(u[mu])*
					   shift(adj(u[mu])*spinProjectDir0Plus(psi),BACKWARD, mu),BACKWARD,mu));

      mu=1;
      fm=fm+spinReconstructDir1Minus(u[mu]*shift(u[mu]*shift(spinProjectDir1Minus(psi), FORWARD, mu),FORWARD,mu));
      fm=fm+spinReconstructDir1Plus (shift(adj(u[mu])*
					   shift(adj(u[mu])*spinProjectDir1Plus(psi),BACKWARD, mu),BACKWARD,mu));

      mu=2;
      fm=fm+spinReconstructDir2Minus(u[mu]*shift(u[mu]*shift(spinProjectDir2Minus(psi), FORWARD, mu),FORWARD,mu));
      fm=fm+spinReconstructDir2Plus (shift(adj(u[mu])*
					   shift(adj(u[mu])*spinProjectDir2Plus(psi),BACKWARD, mu),BACKWARD,mu));

      mu=3;
      fm=fm+spinReconstructDir3Minus(u[mu]*shift(u[mu]*shift(spinProjectDir3Minus(psi), FORWARD, mu),FORWARD,mu));
      fm=fm+spinReconstructDir3Plus (shift(adj(u[mu])*
					   shift(adj(u[mu])*spinProjectDir3Plus(psi),BACKWARD, mu),BACKWARD,mu));

    } else { 
      mu=0;
      fm=fm+spinReconstructDir0Plus(u[mu]*shift(u[mu]*shift(spinProjectDir0Plus(psi), FORWARD, mu),FORWARD,mu));
      fm=fm+spinReconstructDir0Minus (shift(adj(u[mu])*
					   shift(adj(u[mu])*spinProjectDir0Minus(psi),BACKWARD, mu),BACKWARD,mu));

      mu=1;
      fm=fm+spinReconstructDir1Plus(u[mu]*shift(u[mu]*shift(spinProjectDir1Plus(psi), FORWARD, mu),FORWARD,mu));
      fm=fm+spinReconstructDir1Minus (shift(adj(u[mu])*
					   shift(adj(u[mu])*spinProjectDir1Minus(psi),BACKWARD, mu),BACKWARD,mu));

      mu=2;
      fm=fm+spinReconstructDir2Plus(u[mu]*shift(u[mu]*shift(spinProjectDir2Plus(psi), FORWARD, mu),FORWARD,mu));
      fm=fm+spinReconstructDir2Minus (shift(adj(u[mu])*
					   shift(adj(u[mu])*spinProjectDir2Minus(psi),BACKWARD, mu),BACKWARD,mu));

      mu=3;
      fm=fm+spinReconstructDir3Plus(u[mu]*shift(u[mu]*shift(spinProjectDir3Plus(psi), FORWARD, mu),FORWARD,mu));
      fm=fm+spinReconstructDir3Minus (shift(adj(u[mu])*
					   shift(adj(u[mu])*spinProjectDir3Minus(psi),BACKWARD, mu),BACKWARD,mu));

    }

    chi_qdp=fm;
    Double n2 = norm2(chi-chi_qdp);
    cout << "|| Bagel - QDP || = "<< n2 << endl;
    n2 = norm2(chi);
    cout << "|| Bagel || = "<< n2 << endl;
    n2 = norm2(chi_qdp);
    cout << "|| QDP || = "<< n2 << endl;
  }
  Printf("Done\n"); 
  compare_result(chi,chi_qdp,dwf);

}




void compare_result (LatticeFermion A, LatticeFermion B, bfm & dwf)
{

  multi1d<int> x(4);

  int Nspinco=12;

  printf("Result vectors \n");
  for ( x[0]=0; x[0]<dwf.node_latt[0];x[0]++ ) { 
    for ( x[1]=0; x[1]<dwf.node_latt[1];x[1]++ ) { 
      for ( x[2]=0; x[2]<dwf.node_latt[2];x[2]++ ) { 
	//  for ( x[0]=0; x[0]<1;x[0]++ ) { 
	//  for ( x[1]=0; x[1]<1;x[1]++ ) { 
	//  for ( x[2]=0; x[2]<1;x[2]++ ) { 
  for ( x[3]=0; x[3]<dwf.node_latt[3];x[3]++ ) { 

    int printed=0;
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
	  if ( printed==0){
	    printf("site %d %d %d %d\n",x[0],x[1],x[2],x[3]);
	    printed=1;
	  }
	  printf("%d %d %le + i %le != %le + i %le \n",co,sp,toDouble(real(ca)),toDouble(imag(ca)),toDouble(real(cb)),toDouble(imag(cb)));
	} else if ( toDouble(norm2(cb)) > 0.0001 ) { 
	  //if ( printed==0){
	    //	    printf("site %d %d %d %d\n",x[0],x[1],x[2],x[3]);
	    //	    printed=1;
	  //	  }
	  //	  printf("OK %d %d %le + i %le != %le + i %le \n",co,sp,toDouble(real(ca)),toDouble(imag(ca)),toDouble(real(cb)),toDouble(imag(cb)));
	}


    }}


  }}}}


}


