#include <chroma.h>

typedef double Float;

#define Printf printf

#include <bfm.h>

typedef LatticeFermion T;

typedef multi1d<LatticeColorMatrix> U;

typedef bfm bfm_t;

using namespace Chroma;

#if 1
class EguchiKawamoto {
public:

  double mass;
  bfmarg DWa;
  bfm_t  DW;

  bfmarg DWNNa;
  bfm_t  DWNN;

  void dslash_qdp(LatticeFermion & res, LatticeFermion &psi,multi1d<LatticeColorMatrix> u)
  {

    Real m(mass);
    LatticeFermion DwNN;
    LatticeFermion DwNNdag;
    LatticeFermion Dw;

    DwNN=zero;
    int mu=0;
    DwNN=DwNN+spinReconstructDir0Minus(u[mu]*shift(u[mu]*shift(spinProjectDir0Minus(psi), FORWARD, mu),FORWARD,mu));
    DwNN=DwNN+spinReconstructDir0Plus (shift(adj(u[mu])*
					     shift(adj(u[mu])*spinProjectDir0Plus(psi),BACKWARD, mu),BACKWARD,mu));
    
    mu=1;
    DwNN=DwNN+spinReconstructDir1Minus(u[mu]*shift(u[mu]*shift(spinProjectDir1Minus(psi), FORWARD, mu),FORWARD,mu));
    DwNN=DwNN+spinReconstructDir1Plus (shift(adj(u[mu])*
					     shift(adj(u[mu])*spinProjectDir1Plus(psi),BACKWARD, mu),BACKWARD,mu));
    
    mu=2;
    DwNN=DwNN+spinReconstructDir2Minus(u[mu]*shift(u[mu]*shift(spinProjectDir2Minus(psi), FORWARD, mu),FORWARD,mu));
    DwNN=DwNN+spinReconstructDir2Plus (shift(adj(u[mu])*
					     shift(adj(u[mu])*spinProjectDir2Plus(psi),BACKWARD, mu),BACKWARD,mu));
    
    mu=3;
    DwNN=DwNN+spinReconstructDir3Minus(u[mu]*shift(u[mu]*shift(spinProjectDir3Minus(psi), FORWARD, mu),FORWARD,mu));
    DwNN=DwNN+spinReconstructDir3Plus (shift(adj(u[mu])*
					     shift(adj(u[mu])*spinProjectDir3Plus(psi),BACKWARD, mu),BACKWARD,mu));
    
    DwNNdag=zero;
    mu=0;
    DwNNdag=DwNNdag+spinReconstructDir0Plus(u[mu]*shift(u[mu]*shift(spinProjectDir0Plus(psi), FORWARD, mu),FORWARD,mu));
    DwNNdag=DwNNdag+spinReconstructDir0Minus (shift(adj(u[mu])*
						    shift(adj(u[mu])*spinProjectDir0Minus(psi),BACKWARD, mu),BACKWARD,mu));
    
    mu=1;
    DwNNdag=DwNNdag+spinReconstructDir1Plus(u[mu]*shift(u[mu]*shift(spinProjectDir1Plus(psi), FORWARD, mu),FORWARD,mu));
    DwNNdag=DwNNdag+spinReconstructDir1Minus (shift(adj(u[mu])*
						    shift(adj(u[mu])*spinProjectDir1Minus(psi),BACKWARD, mu),BACKWARD,mu));
    
    mu=2;
    DwNNdag=DwNNdag+spinReconstructDir2Plus(u[mu]*shift(u[mu]*shift(spinProjectDir2Plus(psi), FORWARD, mu),FORWARD,mu));
    DwNNdag=DwNNdag+spinReconstructDir2Minus (shift(adj(u[mu])*
						    shift(adj(u[mu])*spinProjectDir2Minus(psi),BACKWARD, mu),BACKWARD,mu));
    
    mu=3;
    DwNNdag=DwNNdag+spinReconstructDir3Plus(u[mu]*shift(u[mu]*shift(spinProjectDir3Plus(psi), FORWARD, mu),FORWARD,mu));
    DwNNdag=DwNNdag+spinReconstructDir3Minus (shift(adj(u[mu])*
						    shift(adj(u[mu])*spinProjectDir3Minus(psi),BACKWARD, mu),BACKWARD,mu));

    
    Dw=zero;
    mu=0;
    Dw=Dw+spinReconstructDir0Minus(u[mu]*shift(spinProjectDir0Minus(psi), FORWARD, mu));
    Dw=Dw+spinReconstructDir0Plus (shift(adj(u[mu])*spinProjectDir0Plus(psi),BACKWARD, mu));
    
    mu=1;
    Dw=Dw+spinReconstructDir1Minus(u[mu]*shift(spinProjectDir1Minus(psi), FORWARD, mu));
    Dw=Dw+spinReconstructDir1Plus (shift(adj(u[mu])*spinProjectDir1Plus(psi),BACKWARD, mu));
    
    mu=2;
    Dw=Dw+spinReconstructDir2Minus(u[mu]*shift(spinProjectDir2Minus(psi), FORWARD, mu));
    Dw=Dw+spinReconstructDir2Plus (shift(adj(u[mu])*spinProjectDir2Plus(psi),BACKWARD, mu));
    
    mu=3;
    Dw=Dw+spinReconstructDir3Minus(u[mu]*shift(spinProjectDir3Minus(psi), FORWARD, mu));
    Dw=Dw+spinReconstructDir3Plus (shift(adj(u[mu])*spinProjectDir3Plus(psi),BACKWARD, mu));
    
// Wilson action = 4+mass -0.5 Dw 
//               = 4+mass -0.5 (1-gamma)e^ip -0.5 (1+gamma)e^-ip
//               = mass + i gamma sin pmu + 4 (1 - cos pmu)
//               = mass + i gmu pmu  - igmu/6 pmu^3 + 1/2 pmu^2  - pmu^4/24
//
// Eguchi-Kawamoto terms 
//
//    D^{NNr}     (  0.5*(DwNN+DwNNDag))     
//    D^{NNgamma} (  0.5*(DwNN-DwNNDag));
//
//    Add
//         C*(4 - 0.5 D^{NNr})
//         D*(  - 0.5 D^{NNgamma})
//
//    add extra 2^n in powers of p^n:
//
//                D ( 2 i gmu pmu - i 8/6 gmu pmu^3)
//                C ( 2 pmu^2 - 16/24 pmu^4) 
//
// Taking D=-1/8 
//        C=-1/4
//
// Gives 
//         S = mass + i gmu pmu (3/4) + pmu^4 (1/6 - 1/24 ) + h.o.
//
// SANITY CHECK:
//
// This is consistent with Eguchi-Kawamoto but they use ident=1, kappa conventions.
// 
// I.e....         S = - 1 
//                     + kappa [ 4/3 D_W   - 1/6 D_W^{NNgamma} - 2/6 D_W^{NNr}]
//
// Note that D_W : D_W^{NNgamma} : D_W{NNr}  arise as 8:-1:-2
//
// We have                                            1:-1/8:-1/4
// And are consistent.
// With D_W having coeff -1/2, our term on the identity is 4+m-1 = 3+m 
//
// Their corresp term on the identity obtained by multiplying by 3/(8kappa) 
//
// Thus, we identify kappa = 1/(8 (1+m/3) = 1/(8+2 m.4/3)
//
// i.e. 4/3 our mass matches theirs, due to a non-unit coeff (3/4) of gmu pmu in the above. 
//
    Real C = -1/4;
    Real D = -1/8;

    QDPIO::cout << "dslash_qdp " <<endl
		<< " psi " <<  norm2(psi) <<endl
		<< " Dw "  << norm2(Dw) <<endl
		<< " DwNN "  << norm2(DwNN) <<endl
		<< " DwNNdag "  << norm2(DwNNdag) << endl;

    res = (4+m+C*4)*psi - 0.5*Dw             
      - C * ( 0.25*(DwNN+DwNNdag))     
      - D * ( 0.25*(DwNN-DwNNdag));
    
  }

  void init(multi1d<LatticeColorMatrix> & u,double m){

    int lx = QDP::Layout::subgridLattSize()[0];
    int ly = QDP::Layout::subgridLattSize()[1];
    int lz = QDP::Layout::subgridLattSize()[2];
    int lt = QDP::Layout::subgridLattSize()[3];

    mass=m;

    /********************************************************
     * Setup DW operators
     ********************************************************
     */
    DWa.solver   = WilsonFermion;
    DWa.node_latt[0]  = lx;
    DWa.node_latt[1]  = ly;
    DWa.node_latt[2]  = lz;
    DWa.node_latt[3]  = lt;

    multi1d<int> procs = QDP::Layout::logicalSize();
    for(int mu=0;mu<4;mu++){
      if ( procs[mu]>1 ) {
	DWa.local_comm[mu] = 0;
      } else { 
	DWa.local_comm[mu] = 1;
      }
    }


    DWa.Ls = 1;
    DWa.mass = mass;
    DWa.Csw = 0.0;

    DWNNa = DWa;
    DWNNa.solver = WilsonNN;

    Printf("Initialising nearest neigbour DW operator\n");
    DW.init(DWa);
    Printf("Initialising next nearest neigbour DW operator\n");
    DWNN.init(DWNNa);
    
    DW.importGauge(u);
    DWNN.importGaugeNN(u);
    Printf("Imported gauge fields\n");
    fflush(stdout);

  }

  double norm(Fermion_t vec[2]){
    return DW.norm(vec[0])+DW.norm(vec[1]);
  }

  double axpy(Fermion_t z[2],Fermion_t x[2],Fermion_t y[2],double a){
    DW.axpy(z[0],x[0],y[0],a);
    DW.axpy(z[1],x[1],y[1],a);
  }
  double axpby(Fermion_t z[2],Fermion_t x[2],Fermion_t y[2],double a,double b){
    DW.axpby(z[0],x[0],y[0],a,b);
    DW.axpby(z[1],x[1],y[1],a,b);
  }

  void Munprec(Fermion_t in[2],Fermion_t out[2],int Dagger) {

    Fermion_t Dw[2];
    Dw[0] = DW.allocFermion();
    Dw[1] = DW.allocFermion();
    
    //    res = (4+mass)*psi - 0.5*Dw             
    //      + C * (4-0.25*(DwNN+DwNNDag))     
    //      + D * ( -0.25*(DwNN-DwNNDag));

    double C = -1/4;
    double D = -1/8;
    double c_ident = 4+mass+ 4*C;
    double c_Dnn   = -C*0.25 - D*0.25;
    double c_Dnndag= -C*0.25 + D*0.25;


    DW.dslash(in[1],Dw[0],0,Dagger); // Nearest neigh flips parity
    DW.dslash(in[0],Dw[1],1,Dagger);
    axpby(out,in,Dw,c_ident,-0.5);

    DWNN.dslash(in[0],Dw[0],0,Dagger); // next nearest neigh keeps parity
    DWNN.dslash(in[1],Dw[1],1,Dagger);
    axpy(out,Dw,out,c_Dnn);

    DWNN.dslash(in[0],Dw[0],0,1-Dagger);
    DWNN.dslash(in[1],Dw[1],1,1-Dagger);
    axpy(out,Dw,out,c_Dnndag);

    DW.freeFermion(Dw[0]);
    DW.freeFermion(Dw[1]);

  }

  // Solution, source
  int invert(LatticeFermion & psi, LatticeFermion &chi)
  {

    //////////////////////////////////////////////
    //Start the CG for unprec matrix
    //////////////////////////////////////////////
    Fermion_t sol[2]; // Both checkerboards
    Fermion_t src[2]; // Both checkerboards
    Fermion_t p[2]; // Both checkerboards
    Fermion_t mp[2];
    Fermion_t mmp[2];
    Fermion_t r[2];
    Fermion_t tmp;

    for(int cb=0;cb<2;cb++){
      sol[cb] = DW.allocFermion();
      src[cb] = DW.allocFermion();
      p[cb]   = DW.allocFermion();
      mp[cb]  = DW.allocFermion();
      mmp[cb] = DW.allocFermion();
      r[cb] = DW.allocFermion();
      DW.importFermion(psi,sol[cb],cb);
      DW.importFermion(chi,mmp[cb],cb);
    }
    QDPIO::cout << "Source ="<<norm(mmp)<<endl;
    QDPIO::cout << "Guess  ="<<norm(sol)<<endl;

    tmp=DW.allocFermion();
    QDPIO::cout << "Eguchi-Kawamoto invert" << endl;

    Munprec(mmp,src,DaggerYes); // CGNE, Src' = Mdag Src, and Sol = (MdagM)^{-1}Mdag Src 
    QDPIO::cout << "Eguchi-Kawamoto invert: Src' = Mdag Src = " << norm(mmp) <<endl;
    
    double ssq  =norm(src);
    double residual=1.0e-8;
    double rsq =  residual* residual*ssq;

    // Initial residual r = Src' - MdagM sol
    Munprec(sol,mp,DaggerNo);
    Munprec(mp,mmp,DaggerYes);
    QDPIO::cout << "Eguchi-Kawamoto invert: MdagM guess= " <<norm(mmp)<< endl;
    axpy(p,mmp,src,-1.0);
    axpy(r,mmp,src,-1.0);

    int max_iter =1000;
    double cp  = norm(r);
    QDPIO::cout << "Eguchi-Kawamoto invert: MdagM guess residual " << cp << endl;
           
    if ( cp <= rsq ) {
      QDPIO::cout <<"CGNE_unprec converged k=0 - nice guess "<<cp << " " << rsq<<endl;
      return 0;
    }  

    double c,a,d,b;

    for (int k=1;k<=max_iter;k++){

      QDPIO::cout << "Eguchi-Kawamoto invert: iteration " << k << endl;
      c=cp;

      Munprec(p,mp,DaggerNo);
      d = norm(mp);

      // Update solution vector
      a = c/d;
      axpy(sol,p,sol,a);
    
      // MdagMp
      Munprec(mp,mmp,DaggerYes); 
      axpy(r,mmp,r,-a);
      cp= norm(r);

      QDPIO::cout <<"Eguchi Kawamoto invert: CGNE_unprec k="<<k<<" r^2 = "<<cp<<endl;

      // Stopping condition
      if ( cp <= rsq ) { 
	QDPIO::cout <<"Eguchi Kawamoto invert: CGNE_unprec converged k="<<k<<" "<<cp << " " << rsq<<endl;
	Munprec(sol,mp,DaggerNo);
	Munprec(mp,mmp,DaggerYes); 
	axpy(r,mmp,src,-1.0);
	double true_residual = sqrt(norm(r)/ssq);
	QDPIO::cout <<"Eguchi Kawamoto invert: CGNE_unprec true residual "<<true_residual<<endl;
	DW.exportFermion(psi,sol[0],0);
	DW.exportFermion(psi,sol[1],1);
	return k;
      }

      // New (conjugate/M-orthogonal) search direction
      b = cp/c;
      axpy(p,p,r,b); //p= b*p+r;
    }
    QDPIO::cout <<"Eguchi Kawamoto invert: CGNE_unprec not converged"<<endl;
    return -1;
    
  }
};


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

  // Make up a random gauge field.
  // Make up a gaussian source and a zero result vector
  multi1d<LatticeColorMatrix> u(Nd); 
  HotSt(u);

  LatticeFermion psi;  
  LatticeFermion chi;
  LatticeFermion check;
  gaussian(chi);
  psi=zero;

  double mass=0.01;
  EguchiKawamoto EK;
  EK.init(u,mass);

  QDPIO::cout << "Testing matrix" <<  norm2(chi) << endl;
  Fermion_t chi_h[2];
  Fermion_t psi_h[2];
  for(int cb=0;cb<2;cb++){ 
    chi_h[cb] = EK.DW.allocFermion();
    psi_h[cb] = EK.DW.allocFermion();
    EK.DW.importFermion(chi,chi_h[cb],cb);
    EK.DW.importFermion(psi,psi_h[cb],cb);
  }
  EK.Munprec(chi_h,psi_h,0);
  EK.dslash_qdp(check,chi,u);
  for(int cb=0;cb<2;cb++){ 
    EK.DW.exportFermion(psi,psi_h[cb],cb);
  }
  for(int cb=0;cb<2;cb++){ 
    EK.DW.freeFermion(chi_h[cb]);
    EK.DW.freeFermion(psi_h[cb]);
  }
  QDPIO::cout << "norm2(diff) " << norm2(psi-check) <<endl;


  QDPIO::cout << "Testing inverter"<< endl;
  psi=zero;
  EK.invert(psi,chi);
  EK.dslash_qdp(check,psi,u);
  QDPIO::cout << "norm2(chi) "  << norm2(chi) <<endl;
  QDPIO::cout << "norm2(psi) "  << norm2(psi) <<endl;
  QDPIO::cout << "norm2(check) "<< norm2(check) <<endl;
  QDPIO::cout << "norm2(diff) " << norm2(check-chi) <<endl;

}

#endif
