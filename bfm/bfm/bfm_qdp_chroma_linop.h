#ifndef _BFM_QDP_CHROMA_LINOP_H_
#define _BFM_QDP_CHROMA_LINOP_H_

#include <chroma.h>
#include <actions/ferm/invert/syssolver_linop_cg_array.h>
#include <actions/ferm/invert/syssolver_linop_aggregate.h>
using namespace Chroma;
typedef multi1d<LatticeColorMatrix> U;
typedef LatticeFermion T4;
typedef multi1d<LatticeFermion> T5;

// Get Chroma equivalents for regression purposes for DeltaL, Solver, linear operator
Handle< SystemSolver<LatticeFermion> >      GetSolver(multi1d<LatticeColorMatrix> u, g5dParams &parms);
Handle< LinearOperatorArray<T4> >  GetLinOp(multi1d<LatticeColorMatrix> u, g5dParams &parms);
Handle< LinearOperator<T4> >                 GetDeltaL(multi1d<LatticeColorMatrix> u, g5dParams &parms);


Handle< LinearOperator<T4> >  GetDeltaL(multi1d<LatticeColorMatrix> u, g5dParams &parms)
{
  Real M5(parms.M5);
  Real mq(parms.mass);
  Real eps_lo(parms.zolo_lo);
  Real eps_hi(parms.zolo_hi);
  Real scale(parms.mobius_scale);
  int Ls = parms.Ls;
   multi1d<int> bcs(Nd);
   bcs[0] = bcs[1] = bcs[2] = bcs[3] = 1;
   Handle< FermBC<T4,U,U> > fbc(new SimpleFermBC< T4, U, U >(bcs));
   Handle<CreateFermState<T4,U,U> > cfs( new CreateSimpleFermState<T4,U,U>(fbc));


   GroupXML_t invparm;
   invparm.xml=std::string(
"   <InvertParam>\n"
"   <invType>CG_INVERTER</invType>\n"
"   <RsdCG>1.0e-9</RsdCG>\n"
"   <MaxCG>3000</MaxCG>\n"
"   </InvertParam>"
);
   invparm.id=std::string("CG_INVERTER");
   invparm.path=std::string("/InvertParam");

   if ( parms.solver == HtCayleyTanh ) {
     UnprecDWFermActArray  S_f(cfs, M5, mq, Ls);
     Handle< FermState<T4,U,U> > fs( S_f.createState(u) );
     Handle< LinearOperator<T4> >  M(S_f.DeltaLs(fs,invparm));
     return  M;
   }
   if ( parms.solver == HwCayleyTanh ) {
     Real b5 = 1.0;
     Real c5 = 1.0;
     UnprecNEFFermActArray  S_f(cfs, M5,b5,c5, mq, Ls);
     Handle< FermState<T4,U,U> > fs( S_f.createState(u) );
     Handle< LinearOperator<T4> >  M(S_f.DeltaLs(fs,invparm));
     return  M;
   }
   if ( parms.solver == HmCayleyTanh ) {
     Real b5 = 0.5*(scale +1.0);
     Real c5 = 0.5*(scale -1.0);
     UnprecNEFFermActArray  S_f(cfs, M5,b5,c5, mq, Ls);
     Handle< FermState<T4,U,U> > fs( S_f.createState(u) );
     Handle< LinearOperator<T4> >  M(S_f.DeltaLs(fs,invparm));
     return  M;
   }
   if ( parms.solver == HwCayleyZolo ) {
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
     Handle< FermState<T4,U,U> > fs( S_f.createState(u) );
     Handle< LinearOperator<T4> >  M(S_f.DeltaLs(fs,invparm));
     return  M;
   }
   if ( parms.solver == HtCayleyZolo ) {
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
     Handle< FermState<T4,U,U> > fs( S_f.createState(u) );
     Handle< LinearOperator<T4> >  M(S_f.DeltaLs(fs,invparm));
     return M;
   }
   if ( parms.solver == HwPartFracZolo ) {
     if ( Ls%2 == 0 ) { 
       printf("Ls is not odd\n");
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
     Handle< FermState<T4,U,U> > fs( S_f.createState(u) );
     Handle< LinearOperator<T4> >  M(S_f.DeltaLs(fs,invparm));
     return M;
   }
   if ( parms.solver == HwContFracZolo ) {
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
"  <RsdCG>1.0e-9</RsdCG>\n"
"  <FermionBC>\n"
"      <FermBC>SIMPLE_FERMBC</FermBC>\n"
"      <boundary>1 1 1 1</boundary>\n"
"   </FermionBC> \n"
"</AuxFermAct>"
);
     param.AuxFermActGrp= std::string("");
     UnprecOvlapContFrac5DFermActArray S_f(fbc,param);
     Handle< FermState<T4,U,U> > fs( S_f.createState(u) );
     Handle< LinearOperator<T4> >  M(S_f.DeltaLs(fs,invparm));
     return  M;
   }
   exit(0);
}

Handle< LinearOperatorArray<T4> >  GetLinOp(multi1d<LatticeColorMatrix> u, g5dParams &parms)
{
  Real M5(parms.M5);
  Real mq(parms.mass);
  Real eps_lo(parms.zolo_lo);
  Real eps_hi(parms.zolo_hi);
  Real scale(parms.mobius_scale);
  int Ls = parms.Ls;
   multi1d<int> bcs(Nd);
   bcs[0] = bcs[1] = bcs[2] = bcs[3] = 1;
   Handle< FermBC<T4,U,U> > fbc(new SimpleFermBC< T4, U, U >(bcs));
   Handle<CreateFermState<T4,U,U> > cfs( new CreateSimpleFermState<T4,U,U>(fbc));


   GroupXML_t invparm;
   invparm.xml=std::string(
"   <InvertParam>\n"
"   <invType>CG_INVERTER</invType>\n"
"   <RsdCG>1.0e-9</RsdCG>\n"
"   <MaxCG>3000</MaxCG>\n"
"   </InvertParam>"
);
   invparm.id=std::string("CG_INVERTER");
   invparm.path=std::string("/InvertParam");

   if ( parms.solver == HtCayleyTanh ) {
     UnprecDWFermActArray  S_f(cfs, M5, mq, Ls);
     Handle< FermState<T4,U,U> > fs( S_f.createState(u) );
     Handle< LinearOperatorArray<T4> > M(S_f.unprecLinOp(fs,mq));
     return  M;
   }
   if ( parms.solver == HwCayleyTanh ) {
     Real b5 = 1.0;
     Real c5 = 1.0;
     UnprecNEFFermActArray  S_f(cfs, M5,b5,c5, mq, Ls);
     Handle< FermState<T4,U,U> > fs( S_f.createState(u) );
     Handle< LinearOperatorArray<T4> > M(S_f.unprecLinOp(fs,mq));
     return  M;
   }
   if ( parms.solver == HmCayleyTanh ) {
     QDPIO::cout << "Mobius Cayley Tanh Linop for scale "<<scale << endl;
     QDPIO::cout << "Mobius Cayley Tanh Linop for M5    "<<M5 << endl;
     QDPIO::cout << "Mobius Cayley Tanh Linop for mq    "<<mq << endl;
     QDPIO::cout << "Mobius Cayley Tanh Linop for Ls    "<<Ls << endl;
     Real b5 = 0.5*(scale +1.0);
     Real c5 = 0.5*(scale -1.0);
     UnprecNEFFermActArray  S_f(cfs, M5,b5,c5, mq, Ls);
     Handle< FermState<T4,U,U> > fs( S_f.createState(u) );
     Handle< LinearOperatorArray<T4> > M(S_f.unprecLinOp(fs,mq));
     return  M;
   }
   if ( parms.solver == HwCayleyZolo ) {
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
     Handle< FermState<T4,U,U> > fs( S_f.createState(u) );
     Handle< LinearOperatorArray<T4> > M(S_f.unprecLinOp(fs,mq));
     return  M;
   }
   if ( parms.solver == HtCayleyZolo ) {
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
     Handle< FermState<T4,U,U> > fs( S_f.createState(u) );
     Handle< LinearOperatorArray<T4> > M(S_f.unprecLinOp(fs,mq));
     return M;
   }
   if ( parms.solver == HwPartFracZolo ) {
     if ( Ls%2 == 0 ) { 
       printf("Ls is not odd\n");
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
     Handle< FermState<T4,U,U> > fs( S_f.createState(u) );
     Handle< LinearOperatorArray<T4> > M(S_f.linOp(fs));
     return M;
   }
   if ( parms.solver == HwContFracZolo ) {
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
"  <RsdCG>1.0e-9</RsdCG>\n"
"  <FermionBC>\n"
"      <FermBC>SIMPLE_FERMBC</FermBC>\n"
"      <boundary>1 1 1 1</boundary>\n"
"   </FermionBC> \n"
"</AuxFermAct>"
);
     param.AuxFermActGrp= std::string("");
     UnprecOvlapContFrac5DFermActArray S_f(fbc,param);
     Handle< FermState<T4,U,U> > fs( S_f.createState(u) );
     Handle< LinearOperatorArray<T4> > M(S_f.linOp(fs));
     return  M;
   }
   exit(0);
}
Handle< SystemSolver<LatticeFermion> > GetSolver(multi1d<LatticeColorMatrix> u, g5dParams &parms)
{
  Real M5(parms.M5);
  Real mq(parms.mass);
  Real eps_lo(parms.zolo_lo);
  Real eps_hi(parms.zolo_hi);
  Real scale(parms.mobius_scale);
  int Ls = parms.Ls;
   multi1d<int> bcs(Nd);
   bcs[0] = bcs[1] = bcs[2] = bcs[3] = 1;
   Handle< FermBC<T4,U,U> > fbc(new SimpleFermBC< T4, U, U >(bcs));
   Handle<CreateFermState<T4,U,U> > cfs( new CreateSimpleFermState<T4,U,U>(fbc));


   GroupXML_t invparm;
   invparm.xml=std::string(
"   <InvertParam>\n"
"   <invType>CG_INVERTER</invType>\n"
"   <RsdCG>1.0e-10</RsdCG>\n"
"   <MaxCG>3000</MaxCG>\n"
"   </InvertParam>"
);
   invparm.id=std::string("CG_INVERTER");
   invparm.path=std::string("/InvertParam");

   if ( parms.solver == HtCayleyTanh ) {
     UnprecDWFermActArray  S_f(cfs, M5, mq, Ls);
     QDPIO::cout << "GetSolver: DWF 4d prec "<<endl;
     QDPIO::cout << "GetSolver: M5 "<<M5<<endl;
     QDPIO::cout << "GetSolver: mq "<<mq<<endl;
     QDPIO::cout << "GetSolver: Ls "<<Ls<<endl;
     Handle< FermState<T4,U,U> > fs( S_f.createState(u) );
     Handle< LinearOperatorArray<T4> > M(S_f.unprecLinOp(fs,mq));
     return  S_f.qprop(fs,invparm);
   }
   if ( parms.solver == HwCayleyTanh ) {
     Real b5 = 1.0;
     Real c5 = 1.0;
     UnprecNEFFermActArray  S_f(cfs, M5,b5,c5, mq, Ls);
     Handle< FermState<T4,U,U> > fs( S_f.createState(u) );
     Handle< LinearOperatorArray<T4> > M(S_f.unprecLinOp(fs,mq));
     return S_f.qprop(fs,invparm);
   }
   if ( parms.solver == HmCayleyTanh ) {
     Real b5 = 0.5*(scale +1.0);
     Real c5 = 0.5*(scale -1.0);
     QDPIO::cout << "GetSolver: b5 "<<b5<<endl;
     QDPIO::cout << "GetSolver: c5 "<<c5<<endl;
     QDPIO::cout << "GetSolver: M5 "<<M5<<endl;
     QDPIO::cout << "GetSolver: mq "<<mq<<endl;
     QDPIO::cout << "GetSolver: Ls "<<Ls<<endl;
     UnprecNEFFermActArray  S_f(cfs, M5,b5,c5, mq, Ls);
     Handle< FermState<T4,U,U> > fs( S_f.createState(u) );
     Handle< LinearOperatorArray<T4> > M(S_f.unprecLinOp(fs,mq));
     return S_f.qprop(fs,invparm);
   }
   if ( parms.solver == HwCayleyZolo ) {
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
     Handle< FermState<T4,U,U> > fs( S_f.createState(u) );
     Handle< LinearOperatorArray<T4> > M(S_f.unprecLinOp(fs,mq));
     return S_f.qprop(fs,invparm);
   }
   if ( parms.solver == HtCayleyZolo ) {
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
     Handle< FermState<T4,U,U> > fs( S_f.createState(u) );
     Handle< LinearOperatorArray<T4> > M(S_f.unprecLinOp(fs,mq));
     return S_f.qprop(fs,invparm);
   }
   if ( parms.solver == HwPartFracZolo ) {
     if ( Ls%2 == 0 ) { 
       printf("Ls is not odd\n");
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
     Handle< FermState<T4,U,U> > fs( S_f.createState(u) );
     Handle< LinearOperatorArray<T4> > M(S_f.linOp(fs));
     return S_f.qprop(fs,invparm);
   }
   if ( parms.solver == HwContFracZolo ) {
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
"  <RsdCG>1.0e-10</RsdCG>\n"
"  <FermionBC>\n"
"      <FermBC>SIMPLE_FERMBC</FermBC>\n"
"      <boundary>1 1 1 1</boundary>\n"
"   </FermionBC> \n"
"</AuxFermAct>"
);
     param.AuxFermActGrp= std::string("");
     UnprecOvlapContFrac5DFermActArray S_f(fbc,param);
     Handle< FermState<T4,U,U> > fs( S_f.createState(u) );
     Handle< LinearOperatorArray<T4> > M(S_f.linOp(fs));
     return S_f.qprop(fs,invparm);
   }
   exit(0);
}
#endif
