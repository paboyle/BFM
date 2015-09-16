#include <chroma.h>

//#define Printf if ( QMP_is_primary_node() ) printf
#define Printf printf
typedef double Float;
#include <bfm.h>

using namespace Chroma;

void compare_result (LatticeFermion A, LatticeFermion B, bfm & dwf);


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
    for ( int i=0;i<argc;i++)
      Printf("%d %s\n",i,argv[i]);
    exit(-1);
  }
#endif
  /********************************************************
   * Setup QDP
   ********************************************************
   */

  multi1d<int> nrow(Nd);
  int Ls;
#ifdef COMMANDLINE
  nrow[0] = atoi(argv[1]);
  nrow[1] = atoi(argv[2]);
  nrow[2] = atoi(argv[3]);
  nrow[3] = atoi(argv[4]);
  Ls = atoi(argv[5]);
#else 
  nrow[0] = 2;
  nrow[1] = 2;
  nrow[2] = 2;
  nrow[3] = 4;
  Ls = 4;
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
  bfmarg dwfa;
  bfm  dwf;
  dwfa.solver =DWF;
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
  multi1d<LatticeColorMatrix> u(Nd); 
  u=zero;

  dwfa.Ls = Ls;
  dwfa.mass = 0.01;
  dwfa.Csw = 0.0;
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

  for(int s=0;s<Ls;s++) gaussian(psi[s]);

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
 
  for(int dag=0;dag<2;dag++){

    for(int s=0;s<Ls;s++) chi[s] = zero;
    for(int s=0;s<Ls;s++) chi_qdp[s] = chi[s];

    /*Import this checkerboard of source field to bagel*/
    for(int cb=0;cb<2;cb++){ 
      dwf.importFermion(psi,psi_h,1-cb);
      dwf.importFermion(chi,chi_h,cb);

      //      dwf.Dperp_cpp(psi_h,
      //		    chi_h, 
      //		    cb,dag);

      dwf.dslash(psi_h,
      		 chi_h, 
		 cb,dag);

      dwf.exportFermion(chi,chi_h,cb);
    }
    
    Real mass = Real(0.01);
    for(int s=0;s<Ls-1;s++){
      if ( dag ) chi_qdp[s]    = 2*chiralProjectPlus(psi[s+1]);
      else       chi_qdp[s]    = 2*chiralProjectMinus(psi[s+1]);
    }
    if ( dag ) chi_qdp[Ls-1] = -0.02 *chiralProjectPlus(psi[0]);
    else       chi_qdp[Ls-1] = -0.02 *chiralProjectMinus(psi[0]);

    if (dag )  chi_qdp[0] += -0.02 *chiralProjectMinus(psi[Ls-1]);
    else       chi_qdp[0] += -0.02 *chiralProjectPlus (psi[Ls-1]);
    for(int s=1;s<Ls;s++){
      if ( dag ) chi_qdp[s]    += 2*chiralProjectMinus(psi[s-1]);
      else       chi_qdp[s]    += 2*chiralProjectPlus (psi[s-1]);
    }

    for(int s=0;s<Ls;s++) cout << " Bagel -qdp "<< s<< " "<< norm2(chi[s]-chi_qdp[s])<<endl;
    for(int s=0;s<Ls;s++) cout << " Bagel "<< s<< " "<< norm2(chi[s])<<endl;
    for(int s=0;s<Ls;s++) cout << " qdp   "<< s<< " "<< norm2(chi_qdp[s])<<endl;


    //    for(int s=0;s<Ls;s++) { 
      //      compare_result(chi[s],chi_qdp[s],dwf);
    //    }

  }
  Printf("Done\n"); 
}

void compare_result (LatticeFermion A, LatticeFermion B, bfm & dwf)
{

  multi1d<int> x(4);

  int Nspinco=12;

  Printf("Result vectors \n");
  for ( x[3]=0; x[3]<dwf.node_latt[3];x[3]++ ) { 
  for ( x[2]=0; x[2]<dwf.node_latt[2];x[2]++ ) { 
  for ( x[1]=0; x[1]<dwf.node_latt[1];x[1]++ ) { 
  for ( x[0]=0; x[0]<dwf.node_latt[0];x[0]++ ) { 

    int xx[4];
    xx[0]=x[0];
    xx[1]=x[1];
    xx[2]=x[2];
    xx[3]=x[3];
    Printf("site %d %d %d %d\n",x[0],x[1],x[2],x[3]);
    for ( int sp=0;sp<4;sp++ ) { 
    for ( int co=0;co<3;co++ ) { 

        int spco = co+3*sp;

	int reim=0;
        Fermion ferm = peekSite(A,x); 
        ColorVector cv = peekSpin(ferm,sp);
        Complex ca = peekColor(cv,co);

        ferm = peekSite(B,x); 
        cv = peekSpin(ferm,sp);
        Complex cb = peekColor(cv,co);

        Printf("%le %le\n",toDouble(real(ca)),toDouble(real(cb)));
        Printf("%le %le\n",toDouble(imag(ca)),toDouble(imag(cb)));

    }}

  }}}}

}






