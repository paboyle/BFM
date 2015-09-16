#include <chroma.h>
#include <actions/ferm/invert/syssolver_linop_cg_array.h>
#include <bfm.h>
#include <bfm_qdp.h>
#include <bfm_qdp_g5d.h>
#include <bfm_qdp_cayley_prec.h>
#include <nested_4d.h>

typedef bfm_dp bfm_t;

using namespace Chroma;

typedef LatticeFermion T;
typedef multi1d<LatticeFermion> T5;
typedef multi1d<LatticeColorMatrix> U;

#define Printf printf

int main (int argc,char **argv )
{
  Chroma::initialize(&argc,&argv);

  /********************************************************
   * Setup QDP
   ********************************************************
   */
  multi1d<int> nrow(Nd);
  nrow[0] = 16;
  nrow[1] = 16;
  nrow[2] = 16;
  nrow[3] = 16;
  int Ls = 8;
  double mass =0.1;
  double m5   =1.8;

  Layout::setLattSize(nrow);
  Layout::create();

  QDPIO::cout << "Setting up BAGEL" <<endl;
  bfmarg::Threads(64);
  bfmarg::Reproduce(0);
  bfmarg::ReproduceChecksum(0);
  bfmarg::ReproduceMasterCheck(0);
  bfmarg::Verbose(1);
  
  /********************************************************
   * Gaussian gauge field
   ********************************************************
   */

  multi1d<LatticeColorMatrix>  u(Nd); HotSt(u);

  /********************************************************
   * Gaussian source and result vectors
   ********************************************************
   */
  LatticeFermion psi;
  LatticeFermion chi;
  LatticeFermion psi_qdp;

  /* 
   * Chroma regression test set up
   */
  // BCs and gauge field
  multi1d<int> bcs(Nd);
  bcs[0] = bcs[1] = bcs[2] = bcs[3] = 1;
  Handle< FermBC<T,U,U> > fbc(new SimpleFermBC< T, U, U >(bcs));
  Handle<CreateFermState<T,U,U> > cfs( new CreateSimpleFermState<T,U,U>(fbc));

  //5d DWF even-odd operator
  //  EvenOddPrecDWFermActArray  DWF_5d(cfs, m5, mass, Ls);
  double scale=1.5;
  EvenOddPrecNEFFermActArrayParams mob;
  mob.OverMass=m5;
  mob.b5 = 0.5*(scale +1.0);
  mob.c5 = 0.5*(scale -1.0);
  mob.Mass = mass;
  mob.N5 = Ls;
  EvenOddPrecNEFFermActArray  DWF_5d(cfs,mob);

  Handle< FermState<T,U,U> >                     fs_5d( DWF_5d.createState(u) );
  //  Handle< EvenOddPrecLinearOperatorArray<T,U,U> >  M_5d(DWF_5d.precLinOp(fs_5d,mass));
  //  Handle< EvenOddPrecLinearOperatorArray<T,U,U> > PV_5d(DWF_5d.precLinOp(fs_5d,1.0));

  SysSolverCGParams invParam;
  invParam.RsdCG        = 1.0e-9;
  invParam.RsdCGRestart = 1.0e-9;
  invParam.MaxCG        = 10000;
  invParam.MaxCGRestart = 10000;

   GroupXML_t invparm_t;
   invparm_t.xml=std::string(
"   <InvertParam>\n"
"   <invType>CG_INVERTER</invType>\n"
"   <RsdCG>1.0e-9</RsdCG>\n"
"   <MaxCG>3000</MaxCG>\n"
"   </InvertParam>"
);
   invparm_t.id=std::string("CG_INVERTER");
   invparm_t.path=std::string("/InvertParam");

  // 4d effective overlap operator
  Handle<LinearOperator<T> >                     op4d( DWF_5d.linOp4D(fs_5d,mass,invparm_t));
  //  UnprecPPDWF4DLinOp<T,U,U>  Overlap_4d(M_5d,PV_5d,invParam);

  QDPIO::cout << "Initialising new linop" <<endl;
  {
    g5dParams parms;
    parms.ScaledShamirCayleyTanh(mass,m5,Ls,scale);

    bfmCayley4d<double> linop;
    linop.init(parms,u,1.0e-9,10000);

    /*
     * Forward matrix
     */
    gaussian(psi);
    gaussian(psi_qdp);
    gaussian(chi);

    QDPIO::cout << "*******************************" <<endl;
    QDPIO::cout << "Applying operator to " << norm2(chi)<<endl;
    QDPIO::cout << "*******************************" <<endl;
    linop.Munprec(chi,psi,DaggerNo);

#define REGRESS_DWF
#ifdef REGRESS_DWF
    QDPIO::cout << "*******************************" <<endl;
    QDPIO::cout << "Applying Chroma operator"<<endl;
    QDPIO::cout << "*******************************" <<endl;
    //    Overlap_4d(psi_qdp,chi,PLUS);
    (*op4d)(psi_qdp,chi,PLUS);
    QDPIO::cout << "*******************************" <<endl;
    QDPIO::cout << "|QDP-Bagel| = "<<norm2(psi-psi_qdp)<<endl;
    QDPIO::cout << "|QDP| = "<<norm2(psi_qdp)<<endl;
    QDPIO::cout << "|Bagel| = "<<norm2(psi)<<endl;
    QDPIO::cout << "*******************************" <<endl;
#endif

    gaussian(psi);
    gaussian(chi);
    QDPIO::cout << "*******************************" <<endl;
    QDPIO::cout << "Applying adjoint operator" <<endl;
    QDPIO::cout << "*******************************" <<endl;
    linop.Munprec(chi,psi,DaggerYes);

#ifdef REGRESS_DWF
    QDPIO::cout << "*******************************" <<endl;
    QDPIO::cout << "Applying Chroma operator"<<endl;
    QDPIO::cout << "*******************************" <<endl;
    //    Overlap_4d(psi_qdp,chi,MINUS);
    (*op4d)(psi_qdp,chi,MINUS);
    QDPIO::cout << "*******************************" <<endl;
    QDPIO::cout << "|QDP-Bagel| = "<<norm2(psi-psi_qdp)<<endl;
    QDPIO::cout << "|QDP| = "<<norm2(psi_qdp)<<endl;
    QDPIO::cout << "|Bagel| = "<<norm2(psi)<<endl;
    QDPIO::cout << "*******************************" <<endl;
#endif

    QDPIO::cout << "*******************************" <<endl;
    QDPIO::cout << "Running CGNE for effective 4d operator" <<endl;
    QDPIO::cout << "*******************************" <<endl;
    psi=zero;
    gaussian(chi);
    linop.CGNE_M(psi,chi);

    QDPIO::cout << "*******************************" <<endl;
    QDPIO::cout << "Subtracting contact term and normalising"<<endl;
    QDPIO::cout << "*******************************" <<endl;
    psi = (psi - chi) *Real(1.0/(1.0-mass));

    QDPIO::cout << "*******************************" <<endl;
    QDPIO::cout << "Performing valence 5d inversion"<<endl;
    QDPIO::cout << "*******************************" <<endl;
    Real resid[1] = {1.0e-8};
    int  iter [1] = {1000};
    Real M5=m5;
    Real Mass=mass;

    LatticeComplex junk;
    double rresid=1.0e-8;
    bfm_Cayley_CG<double>(psi_qdp,chi,u,junk,junk,junk,junk,parms,rresid,iter[0]);
    //    dwf_CG<double,double>(psi_qdp,chi,u,std::string("foo"),Ls,Mass,M5,1,resid,iter);

    QDPIO::cout << "*******************************" <<endl;
    QDPIO::cout << "|DWF4d-Contact-DWF5d| = "<<norm2(psi-psi_qdp)<<endl;
    QDPIO::cout << "|DWF5d| = "<<norm2(psi_qdp)<<endl;
    QDPIO::cout << "|DWF4d-Contact| = "<<norm2(psi)<<endl;
    QDPIO::cout << "*******************************" <<endl;

    linop.end();
  }


}








