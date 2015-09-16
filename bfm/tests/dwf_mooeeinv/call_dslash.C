#include <chroma.h>
#include <bfm.h>

using namespace Chroma;

typedef LatticeFermion T;
typedef multi1d<LatticeColorMatrix> U;

//#define Printf if ( QMP_is_primary_node() ) printf
#define Printf  printf

void applyDiagInv(multi1d<LatticeFermion>& chi, 
             multi1d<LatticeFermion>& psi, 
             int cb,PlusMinus pm,
             Real WilsonMass, Real m_q, int N5);


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
  bfm  dwf;
  dwfa.solver = DWFrb4d;
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
  // cb is cb of result and input field
  for(int cb=0;cb<1;cb++){

    /*Import this checkerboard of source field to bagel*/
    printf("Importing psi field cb %d\n",cb);
    dwf.importFermion(psi,psi_h,cb);

    // Fill the other checkerboard of result with noise
    for(int s=0;s<Ls;s++) chi[s] = zero;
    for(int s=0;s<Ls;s++) chi_qdp[s] = chi[s];
 
    for(int dag=0;dag<2;dag++){

      printf("Checking cb=%d dag=%d\n",cb,dag); 

      dwf.MooeeInv(psi_h,
                   chi_h, 
	  	   dag);
      dwf.exportFermion(chi,chi_h,cb);
 
      // Check the result
      PlusMinus pm;
      if ( dag ) pm = MINUS;
      else pm = PLUS;

      (*M).evenEvenInvLinOp(chi_qdp,psi,pm);

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
  printf("Done\n"); 

}

void 
applyDiagInv(multi1d<LatticeFermion>& chi, 
             multi1d<LatticeFermion>& psi, 
             int cb,
	     PlusMinus pm,
             Real WilsonMass, Real m_q, int N5)
{
  Real InvTwoKappa=5.-WilsonMass;
  Real TwoKappa=1.0/InvTwoKappa;
  Real Kappa=0.5*TwoKappa;
  Real invDfactor=1.0/(1.0  + m_q/pow(InvTwoKappa,N5));

  if( chi.size() != N5 ) chi.resize(N5);

  Real invDTwoKappa = invDfactor*TwoKappa;

  switch (pm) { 

  case PLUS:
    {
      // Copy and scale by TwoKappa (1/M0)
      for(int s(0);s<N5;s++)
      chi[s][rb[cb]] = TwoKappa * psi[s] ;
      
      // First apply the inverse of Lm 
      Real fact(0.5*m_q*TwoKappa) ;
      for(int s(0);s<N5-1;s++){
	chi[N5-1][rb[cb]] -= fact * (chi[s] - GammaConst<Ns,Ns*Ns-1>()*chi[s])  ;
	fact *= TwoKappa ;
      }
      
      
      //Now apply the inverse of L. Forward elimination 
      for(int s(1);s<N5;s++)
	chi[s][rb[cb]] += Kappa*(chi[s-1] + GammaConst<Ns,Ns*Ns-1>()*chi[s-1]) ;
      
      //The inverse of D  now
      chi[N5-1][rb[cb]] *= invDfactor ;
      // That was easy....
      
      
      //The inverse of R. Back substitution...... Getting there! 
      for(int s(N5-2);s>-1;s--)
	chi[s][rb[cb]] += Kappa*(chi[s+1] - GammaConst<Ns,Ns*Ns-1>()*chi[s+1]) ;
      
      
      //Finally the inverse of Rm 
      LatticeFermion tt;
      fact = 0.5*m_q*TwoKappa;
      tt[rb[cb]] = fact*(chi[N5-1] + GammaConst<Ns,Ns*Ns-1>()*chi[N5-1]);
      for(int s(0);s<N5-1;s++){
	chi[s][rb[cb]] -= tt  ;
      tt[rb[cb]] *= TwoKappa ;
      }
    }
    break;
  case MINUS:
    {    
      // Copy and scale by TwoKappa (1/M0)
      for(int s(0);s<N5;s++)
	chi[s][rb[cb]] = TwoKappa * psi[s] ;
      
      // First apply the inverse of Lm 
      Real fact(0.5*m_q*TwoKappa) ;
      for(int s(0);s<N5-1;s++){
	chi[N5-1][rb[cb]] -= fact * (chi[s] + GammaConst<Ns,Ns*Ns-1>()*chi[s])  ;
	fact *= TwoKappa ;
      }
      
      //Now apply the inverse of L. Forward elimination 
      for(int s(1);s<N5;s++)
	chi[s][rb[cb]] += Kappa*(chi[s-1] - GammaConst<Ns,Ns*Ns-1>()*chi[s-1]) ;
      
      //The inverse of D  now
      chi[N5-1][rb[cb]] *= invDfactor ;
      // That was easy....
      
      //The inverse of R. Back substitution...... Getting there! 
      for(int s(N5-2);s>-1;s--)
	chi[s][rb[cb]] += Kappa*(chi[s+1] + GammaConst<Ns,Ns*Ns-1>()*chi[s+1]) ;
      
      //Finally the inverse of Rm 
      LatticeFermion tt;
      tt[rb[cb]] = (0.5*m_q*TwoKappa)*(chi[N5-1] - GammaConst<Ns,Ns*Ns-1>()*chi[N5-1]);
      for(int s(0);s<N5-1;s++){
	chi[s][rb[cb]] -= tt  ;
	tt[rb[cb]] *= TwoKappa ;
      }
    }
    break;

  }
  
}








