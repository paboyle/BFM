#include <chroma.h>
#include <bfm.h>
#include <omp.h>

using namespace Chroma;
typedef multi1d<LatticeColorMatrix> U;

#define Printf if ( QMP_is_primary_node() ) printf

typedef LatticeFermion T;
typedef multi1d<LatticeFermion> T5;
typedef multi1d<LatticeColorMatrix> U;

//#include <bfmcommqmp.h>
//#include <bfmcommspi.h>

void test_solver(BfmSolver solver);
Handle< UnprecLinearOperatorArray<T,U,U> > GetLinOp(U u, int solver,int Ls, double M5, double mq,double eps_lo,double eps_hi);
double mobius_scale = 1.7;

int main (int argc,char **argv )
{
  omp_set_num_threads(1);
  Chroma::initialize(&argc,&argv);
  WilsonTypeFermActs4DEnv::registerAll(); 

  /********************************************************
   * Command line parsing
   ********************************************************
   */
#define COMMANDLINE
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
#else
  nrow[0] = 4;
  nrow[1] = 2;
  nrow[2] = 2;
  nrow[3] = 4;
#endif

  Layout::setLattSize(nrow);
  Layout::create();

  
  test_solver(HtCayleyTanh);
  test_solver(HwCayleyTanh);
  test_solver(HmCayleyTanh);

  test_solver(HtCayleyZolo);
  test_solver(HwCayleyZolo);
  test_solver(HwContFracZolo);
  test_solver(HwPartFracZolo);

}

void test_solver(BfmSolver solver)
{
  int Ls  = 8;

  multi1d<LatticeColorMatrix> u(4);
  HotSt(u);

  double M5=1.8;
  double mq=0.0;

  int lx = QDP::Layout::subgridLattSize()[0];
  int ly = QDP::Layout::subgridLattSize()[1];
  int lz = QDP::Layout::subgridLattSize()[2];
  int lt = QDP::Layout::subgridLattSize()[3];

  multi1d<int> procs = QDP::Layout::logicalSize();
  /********************************************************
   * Setup DWF operator
   ********************************************************
   */
  bfmarg  dwfa;
  dwfa.node_latt[0]  = lx;
  dwfa.node_latt[1]  = ly;
  dwfa.node_latt[2]  = lz;
  dwfa.node_latt[3]  = lt;
  for(int mu=0;mu<4;mu++){
    if ( procs[mu]>1 ) {
      dwfa.local_comm[mu] = 0;
    } else { 
      dwfa.local_comm[mu] = 1;
    }
  }

  bfm dwf_qmp;

  dwfa.mass = toDouble(mq);
  dwfa.M5   = toDouble(M5);
  dwfa.Csw = 0.0;
  dwfa.precon_5d=0;

  dwfa.zolo_lo   = 0.0;
  dwfa.zolo_hi   = 1.0;
  Printf("************************************\n");
  if ( solver == HtCayleyTanh ) { 
    Printf("Testing HtCayleyTanh aka DWF\n");
    dwfa.Ls = Ls;
  } else if ( solver == HmCayleyTanh ) { 
    Printf("Testing HmCayleyTanh Moebius\n");
    dwfa.mobius_scale = mobius_scale;
    dwfa.Ls = Ls;
  } else if ( solver == HwCayleyTanh ) { 
    Printf("Testing HwCayleyTanh aka Borici\n");
    dwfa.Ls = Ls;
  } else if ( solver == HtCayleyZolo ) { 
    Printf("Testing HtCayleyZolo\n");
    dwfa.zolo_lo   = 0.05;
    dwfa.zolo_hi   = 1.8;
    dwfa.Ls = Ls;
  } else if ( solver == HwCayleyZolo ) { 
    Printf("Testing HwCayleyZolo aka Chiu\n");
    dwfa.Ls = Ls;
    dwfa.zolo_lo   = 0.05;
    dwfa.zolo_hi   = 7.0;
  } else if ( solver == HwPartFracZolo ) { 
    Printf("Testing HwPartFracZolo aka KEK\n");
    dwfa.Ls = Ls = 9;
    dwfa.zolo_lo   = 0.05;
    dwfa.zolo_hi   = 7.0;
    Printf("Partial fraction: setting Ls =%d\n",Ls);
  } else if ( solver == HwContFracZolo ) { 
    Printf("Testing HwContFracZolo \n");
    dwfa.Ls = Ls = 9;
    Printf("Cont fraction: setting Ls =%d\n",Ls);
    dwfa.zolo_lo   = 0.05;
    dwfa.zolo_hi   = 7.0;
  }
  Printf("************************************\n");

  dwfa.solver= solver;

  Printf("Initialising bagel\n"); fflush(stdout);
  dwf_qmp.init(dwfa);
  Printf("Importing gauge fieldp\n"); fflush(stdout);
  dwf_qmp.importGauge(u);

  Fermion_t psi_t[2];
  Fermion_t chi_t[2];

  multi1d<LatticeFermion>  check(Ls);
  multi1d<LatticeFermion> result(Ls);
  multi1d<LatticeFermion>  psi(Ls);
  for(int s=0;s<dwfa.Ls;s++) {
    if ( s==0) gaussian(psi[s]);
    //    else psi[s]=zero;
    else gaussian(psi[s]);
  }

  for(int cb=0;cb<2;cb++){
    psi_t [cb] = dwf_qmp.allocFermion();
    chi_t [cb] = dwf_qmp.allocFermion();
    dwf_qmp.importFermion(psi,psi_t[cb],cb);
  }

  Printf("Getting Chroma Linop\n"); fflush(stdout);
  Handle< UnprecLinearOperatorArray<T,U,U> > linop =GetLinOp(u, 
							     dwf_qmp.solver,
							     dwf_qmp.Ls, 
							     dwf_qmp.M5, 
							     dwf_qmp.mass,
							     dwf_qmp.zolo_lo,
							     dwf_qmp.zolo_hi
							     );
  Printf("Got Chroma Linop\n"); fflush(stdout);

  for (int dag=0;dag<2;dag++){ 
    Printf("Calling Munprec\n"); fflush(stdout);
    dwf_qmp.G5D_Munprec(psi_t,chi_t,dag);
    Printf("Called Munprec\n"); fflush(stdout);
    double d=dwf_qmp.norm(chi_t[0])+dwf_qmp.norm(chi_t[1]);
    Printf("Norm G5D %le\n",d);fflush(stdout);

    dwf_qmp.exportFermion(result,chi_t[0],0);
    dwf_qmp.exportFermion(result,chi_t[1],1);

    Printf("Calling Chroma Linop\n"); fflush(stdout);

    if ( dag ) 
      (*linop)(check,psi,MINUS);
    else
      (*linop)(check,psi,PLUS);

    Printf("Called Chroma Linop\n"); fflush(stdout);

    Real rmc=0.0;
    Real r=0.0;
    Real c=0.0;

    for(int s=0;s<Ls;s++){
      QDPIO::cout << " s " << s << "  diff " << norm2(result[s]-check[s]) << " bagel " 
		  << norm2(result[s]) << " chroma " << norm2(check[s]) <<endl;

      rmc+=norm2(result[s]-check[s]);
      r+=norm2(result[s]);
      c+=norm2(check[s]);
    }

    QDPIO::cout << "Norm diff "<<  rmc << " " << r << " " << c << endl;
    
    if ( toDouble(rmc)/toDouble(r) > 1.0e-8  ) { 
      Printf ("LOOKS BAD!!!\n");
      exit(-1);
    }
  }
  for(int cb=0;cb<2;cb++){
    dwf_qmp.freeFermion(psi_t[cb]);
    dwf_qmp.freeFermion(chi_t[cb]);
  }
  dwf_qmp.end();
  

}

Handle< UnprecLinearOperatorArray<T,U,U> > GetLinOp(U u, int solver,int Ls,double M5_d, double mq_d,
						    double eps_l, double eps_h)
{
  Real M5(M5_d);
  Real mq(mq_d);
  Real eps_lo(eps_l);
  Real eps_hi(eps_h);
   multi1d<int> bcs(Nd);
   bcs[0] = bcs[1] = bcs[2] = bcs[3] = 1;

   Handle< FermBC<T,U,U> > fbc(new SimpleFermBC< T, U, U >(bcs));
   Handle<CreateFermState<T,U,U> > cfs( new CreateSimpleFermState<T,U,U>(fbc));

   if ( solver == HtCayleyTanh ) {
     UnprecDWFermActArray  S_f(cfs, M5, mq, Ls);
     Handle< FermState<T,U,U> > fs( S_f.createState(u) );
     Handle< UnprecLinearOperatorArray<T,U,U> > M(S_f.unprecLinOp(fs,mq));
     return M;
   }
   if ( solver == HwCayleyTanh ) {
     Real b5 = 1.0;
     Real c5 = 1.0;
     UnprecNEFFermActArray  S_f(cfs, M5,b5,c5, mq, Ls);
     Handle< FermState<T,U,U> > fs( S_f.createState(u) );
     Handle< UnprecLinearOperatorArray<T,U,U> > M(S_f.unprecLinOp(fs,mq));
     return M;
   }
   if ( solver == HmCayleyTanh ) {
     Real b5 = 0.5*(mobius_scale +1.0);
     Real c5 = 0.5*(mobius_scale -1.0);
     UnprecNEFFermActArray  S_f(cfs, M5,b5,c5, mq, Ls);
     Handle< FermState<T,U,U> > fs( S_f.createState(u) );
     Handle< UnprecLinearOperatorArray<T,U,U> > M(S_f.unprecLinOp(fs,mq));
     return M;
   }
   if ( solver == HwCayleyZolo ) {
     UnprecZoloNEFFermActArrayParams params;
     params.OverMass=M5;
     params.Mass=mq;
     params.b5=1.0;
     params.c5=1.0;
     params.N5=Ls;
     params.approximation_type = COEFF_TYPE_ZOLOTAREV;
     params.ApproxMin=eps_lo;
     params.ApproxMax=eps_hi;
     UnprecZoloNEFFermActArray  S_f(cfs, params);
     Handle< FermState<T,U,U> > fs( S_f.createState(u) );
     Handle< UnprecLinearOperatorArray<T,U,U> > M(S_f.unprecLinOp(fs,mq));
     return M;
   }
   if ( solver == HtCayleyZolo ) {
     UnprecZoloNEFFermActArrayParams params;
     params.OverMass=M5;
     params.Mass=mq;
     params.b5=1.0;
     params.c5=0.0;
     params.N5=Ls;
     params.approximation_type = COEFF_TYPE_ZOLOTAREV;
     params.ApproxMin=eps_lo;
     params.ApproxMax=eps_hi;
     UnprecZoloNEFFermActArray  S_f(cfs, params);
     Handle< FermState<T,U,U> > fs( S_f.createState(u) );
     Handle< UnprecLinearOperatorArray<T,U,U> > M(S_f.unprecLinOp(fs,mq));
     return M;
   }
   if ( solver == HwPartFracZolo ) {
     if ( Ls%2 == 0 ) { 
       Printf("Ls is not odd\n");
       exit(-1);
     }
     UnprecOvExtFermActArrayParams param;
     param.OverMass=M5; 
     param.Mass=mq;
     param.RatPolyDeg = Ls;
     param.ApproxMin =eps_lo;
     param.ApproxMax =eps_hi;
     param.b5 =1.0;
     param.c5 =1.0;
     param.approximation_type=COEFF_TYPE_ZOLOTAREV;
     //     param.approximation_type=COEFF_TYPE_TANH_UNSCALED;
     //     param.approximation_type=COEFF_TYPE_TANH;
     param.tuning_strategy_xml=
"<TuningStrategy><Name>OVEXT_CONSTANT_STRATEGY</Name></TuningStrategy>\n";
     UnprecOvExtFermActArray S_f(cfs,param);
     Handle< FermState<T,U,U> > fs( S_f.createState(u) );
     Handle< UnprecLinearOperatorArray<T,U,U> > M(S_f.linOp(fs));
     return M;
   }
   if ( solver == HwContFracZolo ) {
     UnprecOvlapContFrac5DFermActParams param;
     param.Mass=mq; // How is M5 set? Wilson mass In AuxFermAct
     param.ApproxMin=eps_lo;
     param.ApproxMax=eps_hi;
     param.approximation_type=COEFF_TYPE_ZOLOTAREV;
     param.RatPolyDeg=Ls;
     // The following is why I think Chroma made some directional errors:
     param.AuxFermAct= std::string(
"<AuxFermAct>\n"
"  <FermAct>UNPRECONDITIONED_WILSON</FermAct>\n"
"  <Mass>-1.8</Mass>\n"
"  <b5>1</b5>\n"
"  <c5>0</c5>\n"
"  <MaxCG>1000</MaxCG>\n"
"  <RsdCG>1.0e-7</RsdCG>\n"
"  <FermionBC>\n"
"      <FermBC>SIMPLE_FERMBC</FermBC>\n"
"      <boundary>1 1 1 1</boundary>\n"
"   </FermionBC> \n"
"</AuxFermAct>"
);
     param.AuxFermActGrp= std::string("");
     UnprecOvlapContFrac5DFermActArray S_f(fbc,param);
     Handle< FermState<T,U,U> > fs( S_f.createState(u) );
     Handle< UnprecLinearOperatorArray<T,U,U> > M(S_f.linOp(fs));
     return M;
   }
   exit(0);
}









