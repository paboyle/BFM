#include <chroma.h>
#include <bfm.h>
#include <bfm_qdp.h>
#include <bfm_qdp_g5d.h>
#include <math.h>

#include <bfm_qdp_chroma_linop.h>

#define Printf if ( QMP_is_primary_node() ) printf

#include <bfmcommqmp.h>
#include <bfmcommspi.h>
typedef bfm_dp bfm_qmp;

void test_solver(BfmSolver solver);

int main (int argc,char **argv )
{
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

  for(int i=0;i<1;i++){
    test_solver(HtCayleyTanh);
    test_solver(HwCayleyTanh);
    test_solver(HmCayleyTanh);
    test_solver(HtCayleyZolo);
    test_solver(HwCayleyZolo);
    test_solver(HwContFracZolo);
    test_solver(HwPartFracZolo);
  }
}

void test_solver(BfmSolver solver)
{
  g5dParams parms;

  int Ls=8;
  double M5=1.8;
  double mq=0.1;
  double wilson_lo = 0.05;
  double wilson_hi = 6.8;
  double shamir_lo = 0.025;
  double shamir_hi = 1.7;
  double ht_scale=1.7;
  double hw_scale=1.0;

  if ( solver == HtCayleyTanh ) { 
    Printf("Testing HtCayleyTanh aka DWF\n");
    parms.ShamirCayleyTanh(mq,M5,Ls);
  } else if ( solver == HmCayleyTanh ) { 
    parms.ScaledShamirCayleyTanh(mq,M5,Ls,ht_scale);
    Printf("Testing HmCayleyTanh Moebius\n");
  } else if ( solver == HwCayleyTanh ) { 
    parms.WilsonCayleyTanh(mq,M5,Ls,hw_scale);
    Printf("Testing HwCayleyTanh aka Borici\n");
  } else if ( solver == HtCayleyZolo ) { 
    parms.ShamirCayleyZolo(mq,M5,Ls,shamir_lo,shamir_hi);
    Printf("Testing HtCayleyZolo\n");
  } else if ( solver == HwCayleyZolo ) { 
    parms.WilsonCayleyZolo(mq,M5,Ls,wilson_lo,wilson_hi);
    Printf("Testing HwCayleyZolo aka Chiu\n");
  } else if ( solver == HwPartFracZolo ) { 
    Printf("Testing HwPartFracZolo aka KEK\n");
    Ls = 9;
    parms.WilsonPartFracZolo(mq,M5,Ls,wilson_lo,wilson_hi);
    Printf("Partial fraction: setting Ls =%d\n",Ls);
  } else if ( solver == HwContFracZolo ) { 
    Printf("Testing HwContFracZolo \n");
    Ls = 9;
    parms.WilsonContFracZolo(mq,M5,Ls,wilson_lo,wilson_hi);
  } else { 
    Printf("Unknown testcase \n");
    exit(-1);
  }

  multi1d<LatticeColorMatrix> u(4);
  HotSt(u);
  //  u=zero;
  

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

  bfm_qmp dwf_qmp;

  dwfa.mass = parms.mass;
  dwfa.M5   = parms.M5;
  dwfa.Csw = 0.0;
  dwfa.precon_5d=0;
  dwfa.Ls = parms.Ls;
  dwfa.solver= parms.solver;
  dwfa.zolo_lo   = parms.zolo_lo;
  dwfa.zolo_hi   = parms.zolo_hi;
  dwfa.mobius_scale = parms.mobius_scale;

  dwf_qmp.init(dwfa);
  dwf_qmp.importGauge(u);

  Fermion_t psi_t[2];
  Fermion_t chi_t[2];
  Fermion_t eo_t[2];
  Fermion_t ee_t[2];

  multi1d<LatticeFermion> check(Ls);
  multi1d<LatticeFermion> result(Ls);
  multi1d<LatticeFermion> ee(Ls);
  multi1d<LatticeFermion> eetmp(Ls);
  multi1d<LatticeFermion> eo(Ls);
  multi1d<LatticeFermion> mprec(Ls);
  multi1d<LatticeFermion> psi(Ls);
  for(int s=0;s<dwfa.Ls;s++) {
    if ( s==0) gaussian(psi[s]);
    else gaussian(psi[s]);
  }

  for(int cb=0;cb<2;cb++){
    psi_t [cb] = dwf_qmp.allocFermion();
    chi_t [cb] = dwf_qmp.allocFermion();
    eo_t  [cb] = dwf_qmp.allocFermion();
    ee_t  [cb] = dwf_qmp.allocFermion();
    dwf_qmp.importFermion(psi,psi_t[cb],cb);
  }
  

  Printf("Getting Chroma Linop\n"); fflush(stdout);
  Handle< LinearOperatorArray<LatticeFermion> > linop =GetLinOp(u, parms);


  Printf("Checking unprec and prec match\n"); fflush(stdout);
  for (int dag=0;dag<2;dag++){ 

    Printf("Calling Chroma Linop\n"); fflush(stdout);
    if ( dag ) 
      (*linop)(check,psi,MINUS);
    else
      (*linop)(check,psi,PLUS);
    Printf("Called Chroma Linop\n"); fflush(stdout);
    QDPIO::cout << "Chroma "<<norm2(check)<<endl;

    dwf_qmp.importFermion(psi,psi_t[0],0);
    dwf_qmp.importFermion(psi,psi_t[1],1);

    for (int cb=0;cb<2;cb++){ 
      dwf_qmp.G5D_Meo  (psi_t[1-cb],eo_t[cb],cb,dag);
      dwf_qmp.G5D_Mooee(psi_t[cb],ee_t[cb],dag);    
    }
    dwf_qmp.G5D_Munprec(psi_t,chi_t,dag);

    for(int cb=0;cb<2;cb++){
      dwf_qmp.exportFermion(psi,psi_t[cb],cb);
      dwf_qmp.exportFermion(result,chi_t[cb],cb);
      dwf_qmp.exportFermion(ee,ee_t[cb],cb);
      dwf_qmp.exportFermion(eo,eo_t[cb],cb);
    }

    mprec = ee+eo;

    QDPIO::cout <<"Mprec   "<< norm2(mprec) <<endl;
    QDPIO::cout <<"Munprec "<< norm2(result) <<endl;
    QDPIO::cout <<"N2diff of Mprec and Munprec is       "<< norm2(result-mprec) <<endl;
    QDPIO::cout <<"N2diff of the Munprec and Chroma is  "<< norm2(check-result) <<endl;
    QDPIO::cout <<"N2diff of the Mprec and Chroma is  "<< norm2(check-mprec) <<endl;


    Printf("Checking Mee Mee^{-1} is unity, dag = %d\n",dag); fflush(stdout);
    for(int cb=0;cb<2;cb++){
#if 1
      dwf_qmp.G5D_Mooee(psi_t[cb],ee_t[0],dag);    
      dwf_qmp.G5D_MooeeInv(ee_t[0],ee_t[1],dag);    
      dwf_qmp.exportFermion(ee,ee_t[1],cb);
#else
      dwf_qmp.G5D_Mooee(psi_t[cb],ee_t[0],dag);    
      dwf_qmp.G5D_MooeeInv(psi_t[cb],ee_t[1],dag);    
      dwf_qmp.exportFermion(ee,ee_t[0],cb);
      dwf_qmp.exportFermion(eetmp,ee_t[1],cb);
#endif
    }
    multi1d<LatticeFermion> check(psi);
    QDPIO::cout << "norm of psi_t" << dwf_qmp.norm(psi_t[0]) <<endl;
    QDPIO::cout << "norm of psi_t" << dwf_qmp.norm(psi_t[1]) <<endl;

    QDPIO::cout << "Inverse check "<< norm2(ee-check)<< endl;
    for(int s=0;s<Ls;s++){
      QDPIO::cout << "Inverse check ["<<s<<"] "<< norm2(ee[s]-check[s])<< " " 
      		  << norm2(ee[s]) << " " 
		  << norm2(check[s]) << endl;

    }
  }


  for(int cb=0;cb<2;cb++){
    dwf_qmp.freeFermion(psi_t[cb]);
    dwf_qmp.freeFermion(ee_t[cb]);
    dwf_qmp.freeFermion(eo_t[cb]);
    dwf_qmp.freeFermion(chi_t[cb]);
  }
  dwf_qmp.end();
  

}

