#include <chroma.h>
#include <bfm.h>

using namespace Chroma;
typedef multi1d<LatticeColorMatrix> U;

#define Printf if ( QMP_is_primary_node() ) printf
//#define Printf  printf
typedef LatticeFermion T;
typedef multi1d<LatticeFermion> T5;
typedef multi1d<LatticeColorMatrix> U;
#include <bfmcommqmp.h>
#include <bfmcommspi.h>
typedef bfmcommspi<double> bfm_spi;
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
  nrow[0] = 16;
  nrow[1] = 16;
  nrow[2] = 16;
  nrow[3] = 16;
  int Ls  = 4;
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
  dwfa.solver = DWF;
  bfm_spi  dwf;
  bfm_qmp  dwf_qmp;
  //  bfm_internal<float>  dwf;

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

  dwfa.Ls = Ls;
  dwfa.mass = toDouble(mq);
  dwfa.M5   = toDouble(M5);
  dwfa.Csw = 0.0;
  dwfa.precon_5d=1;

  Printf("Initialising bfm operator\n");
  dwf.init(dwfa);
  dwf_qmp.init(dwfa);


  /********************************************************
   * Gaussian source and result vectors
   ********************************************************
   */
  multi1d<LatticeFermion> psi(Ls);
  multi1d<LatticeFermion> chi(Ls);
  multi1d<LatticeFermion> chi_qdp(Ls);

  for(int s=0;s<Ls;s++) gaussian(psi[s]);
  for(int s=0;s<Ls;s++) gaussian(chi[s]);
  for(int s=0;s<Ls;s++) gaussian(chi_qdp[s]);


  /********************************************************
   * QDP Linop
   ********************************************************
   */

   multi1d<int> bcs(Nd);
   bcs[0] = bcs[1] = bcs[2] = bcs[3] = 1;

   Handle< FermBC<T,U,U> > fbc(new SimpleFermBC< T, U, U >(bcs));
   Handle<CreateFermState<T,U,U> > cfs( new CreateSimpleFermState<T,U,U>(fbc));
   UnprecDWFermActArray  S_f(cfs, M5, mq, Ls);
   Handle< FermState<T,U,U> > fs( S_f.createState(u) );

   Handle< UnprecLinearOperatorArray<T,U,U> > M(S_f.unprecLinOp(fs,mq));

  /********************************************************
   * Bagel internal single checkerboard vectors
   ********************************************************
   */
  Fermion_t psi_h[2];
  for(int cb=0;cb<2;cb++) psi_h[cb] = dwf.allocFermion();
  Fermion_t chi_h[2];
  for(int cb=0;cb<2;cb++) chi_h[cb] = dwf.allocFermion();
  Fermion_t tmp_h = dwf.allocFermion();

  /********************************************************
   * Import gauge field to BAGEL
   ********************************************************
   */
  dwf.importGauge(u);
  dwf_qmp.importGauge(u);

  LatticeFermion noise; gaussian(noise);

#define DEBUG

  Printf("Starting to check\n"); fflush(stdout);

  // Fill the other checkerboard of result with noise
 
  for(int dag=0;dag<2;dag++){

    for(int s=0;s<Ls;s++) chi[s] = zero;
    for(int s=0;s<Ls;s++) chi_qdp[s] = chi[s];

    /*Import this checkerboard of source field to bagel*/
    for(int cb=0;cb<2;cb++) {
      double d = dwf.norm(psi_h[cb]);
      dwf.importFermion(psi,psi_h[cb],cb);
      Printf("norm[%d] = %le\n",cb,(double)d);
    }
    
    dwf.Munprec(psi_h, 
		chi_h,
  	        tmp_h, 
	        dag);

    for(int cb=0;cb<2;cb++) dwf.exportFermion(chi,chi_h[cb],cb);

    dwf_qmp.Munprec(psi_h, 
		    chi_h,
		    tmp_h, 
		    dag);

    fflush(stdout);
    QMP_barrier();

    for(int mu=0;mu<3;mu++){
      if( !dwfa.local_comm[mu] ) {
	if( mu==0){
	  printf("QMPSPI Recv miscompare mu%d -ve : %le %le\n",mu,dwf.recvbufs[mu*2][0],dwf_qmp.recvbufs[mu*2][0]);
	  printf("QMPSPI Recv miscompare mu%d +ve : %le %le\n",mu,dwf.recvbufs[mu*2+1][0],dwf_qmp.recvbufs[mu*2+1][0]);
	  printf("QMPSPI Send miscompare mu%d -ve : %le %le\n",mu,dwf.sendbufs[mu*2][0],dwf_qmp.sendbufs[mu*2][0]);
	  printf("QMPSPI Send miscompare mu%d +ve : %le %le\n",mu,dwf.sendbufs[mu*2+1][0],dwf_qmp.sendbufs[mu*2+1][0]);
	  printf("QMPSPI direction %d OK!\n",mu);
	}
      }
    }
 
    // Check the result
    PlusMinus pm;
    if ( dag ) pm = MINUS;
    else pm = PLUS;

    (*M)(chi_qdp,psi,pm);

    Double n2 = 0.0;
    for(int s=0;s<Ls;s++) { 
      Double sn2 = norm2(chi[s]-chi_qdp[s]);
      n2+=sn2;
#ifdef DEBUG
      QDPIO::cout << "|| Bagel - QDP || ["<<s<<"] = "<< sn2 << endl;
#endif
    }
    QDPIO::cout << "|| Bagel - QDP || = "<< n2 << endl;
    
    n2 = 0.0;
    for(int s=0;s<Ls;s++) { 
      Double sn2 = norm2(chi[s]);
      n2+=sn2;
#ifdef DEBUG
      QDPIO::cout << "|| Bagel || ["<<s<<"] = "<< sn2 << endl;
#endif
    }
    QDPIO::cout << "|| Bagel || = "<< n2 << endl;
    
    n2 = 0.0;
    for(int s=0;s<Ls;s++) {
      Double sn2 = norm2(chi_qdp[s]);
      n2+=sn2;
#ifdef DEBUG
      QDPIO::cout << "|| QDP || ["<<s<<"] = "<< sn2 << endl;
#endif
    }
    QDPIO::cout << "|| QDP || = "<< n2 << endl;

  }

  printf("Done\n"); 

}








