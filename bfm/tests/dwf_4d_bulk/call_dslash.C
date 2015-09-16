#include <chroma.h>
#include <bfm.h>
#include <bfm_qdp.h>
#include <bfm_wrapper.h>
//#include <bfm_qdp_dwf.h>
//#include <bfm_qdp_g5d.h>
//#include <bfm_qdp_chroma_linop.h>

//#include <actions/ferm/invert/syssolver_linop_cg_array.h>
typedef bfm_dp bfm_t;
 
using namespace Chroma;

typedef LatticeFermion T;
typedef multi1d<LatticeFermion> T5;
typedef multi1d<LatticeColorMatrix> U;

#define Printf if ( QMP_is_primary_node() ) printf

void Tinv(LatticeFermion &sol, 
	  LatticeFermion &src, 
	  multi1d<LatticeColorMatrix> u,
	  double mass,
	  double M5,int inv=1);
double mobius_scale;

int main (int argc,char **argv )
{
  Chroma::initialize(&argc,&argv);

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


  Layout::setLattSize(nrow);
  Layout::create();


  multi1d<LatticeColorMatrix>  u(Nd); 
  if (0)
  {
    ArchivGauge_t Header ;
    readArchiv(Header,u,argv[6]);
  } else {
    HotSt(u);
  }
  double mass =0.0;
  double m5   =1.4;

  QDPIO::cout << "Setting up BAGEL for dimensions "<<
    nrow[0]<<" " <<
    nrow[1]<<"x" <<
    nrow[2]<<"x" <<
    nrow[3]<<" Ls=" <<
    Ls<<endl;
  bfmarg::Threads(1);
  bfmarg::Reproduce(0);
  bfmarg::ReproduceChecksum(0);
  bfmarg::ReproduceMasterCheck(0);
  bfmarg::Verbose(2);
  

  LatticeFermion src;
  LatticeFermion dwf_prop;
  LatticeFermion ovlap_surface;
  multi1d<LatticeFermion> bulk_dwf(Ls);

  
    gaussian(src);

  // Obtain surface DWF propagator
  QDPIO::cout << "Getting surface propagator"<<endl;

  multi1d<LatticeFermion> psi(Ls);
  multi1d<LatticeFermion> chi(Ls);
  {
    BfmWrapperParams BfmWP;

    BfmWP.BfmMatrix = BfmMat_M;
    BfmWP.BfmInverter=BfmInv_CG; 
    BfmWP.BfmPrecision=Bfm64bit;
    BfmWP.MaxIter=10000;
    BfmWP.RsdTarget.resize(1);
    BfmWP.RsdTarget[0]=1.0e-18;
    BfmWP.Delta=1.0e-3;

    mobius_scale=1.2;
    BfmWP.BAP.ScaledShamirCayleyTanh(mass,m5,Ls,mobius_scale);

    BfmWrapper BfmWrap(BfmWP);
    BfmWrap.links = u;
    
    for(int s=0;s<Ls;s++)psi[s]=zero;
    for(int s=0;s<Ls;s++)chi[s]=zero;
    psi[0]   = chiralProjectPlus(src); 
    psi[Ls-1]= chiralProjectMinus(src);
    BfmWrap.bfmInvert5d<double>(chi,psi);
    dwf_prop = chiralProjectMinus(chi[0]);
    dwf_prop+= chiralProjectPlus(chi[Ls-1]);

  }


  // Add contact term to get effective overlap prop
  QDPIO::cout << "Adding contact term and normalising"<<endl;
  ovlap_surface = dwf_prop* (1.0-mass) + src;

  //Now try to build the full 5d solution out of 4d operations
  QDPIO::cout << "Reconstructing the bulk field"<<endl;
  bulk_dwf[0] = dwf_prop;

  Gamma G5(15);
  LatticeFermion bulk_source =
          chiralProjectPlus(ovlap_surface)
    -mass*chiralProjectMinus(ovlap_surface)
    -G5*src;

  bulk_source = bulk_source/(1.0-mass);

  //Apply transfer matrix repeatedly into bulk
  Tinv(bulk_dwf[Ls-1],bulk_source,u,mass,m5,1);

  // A Little test
  Tinv(bulk_dwf[1],bulk_dwf[Ls-1],u,mass,m5,0);
  QDPIO::cout<< "T Tinv = " << norm2(bulk_dwf[1]-bulk_source)<<endl;

  for(int s=Ls-2;s>0;s--){
    Tinv(bulk_dwf[s],bulk_dwf[s+1],u,mass,m5,1);
  }

  // Unwrap basis permutation on the result vector
  QDPIO::cout << "Applying permutation matrix"<<endl;
  multi1d<LatticeFermion> dwf(Ls);
  for(int s=0;s<Ls;s++){
    int sp = (s+1)%Ls;
    dwf[s]=chiralProjectMinus(bulk_dwf[s]) + chiralProjectPlus(bulk_dwf[sp]);
  }
  ///////////////////////////////////////////////////
  // Compare this to DWF prop from Bagel
  ///////////////////////////////////////////////////
  for(int s=0;s<Ls;s++) { 
    Double sn2 = norm2(chi[s]-dwf[s]);
    QDPIO::cout << "|| Bagel - QDP || ["<<s<<"] = "<< sn2 << endl;
    sn2 = norm2(chiralProjectPlus(chi[s]));
    QDPIO::cout << "|| QDP+ || ["<<s<<"] = "<< sn2 << endl;
    sn2 = norm2(chiralProjectPlus(dwf[s]));
    QDPIO::cout << "|| dwf+ || ["<<s<<"] = "<< sn2 << endl;
    sn2 = norm2(chiralProjectMinus(chi[s]));
    QDPIO::cout << "|| QDP- || ["<<s<<"] = "<< sn2 << endl;
    sn2 = norm2(chiralProjectMinus(dwf[s]));
    QDPIO::cout << "|| dwf- || ["<<s<<"] = "<< sn2 << endl;
  }

  // J5q formed as follows:
  //	psi = chiralProjectPlus(psi_5[Ls/2-1]);
  //    psi += chiralProjectMinus(psi_5[Ls/2]);
  //	PJ5q+= localInnerProduct(psi,psi);
  //
  // In 5d overlap field basis this is
  //  Q = Pminus bulk_dwf[Ls/2]
  //    + Pplus  bulk_dwf[Ls/2] (remember we start at zero)
  // 
  // This is just (T^{Ls/2}+T^{-Ls/2})^{-1} Dov(m)
  // So, taking bulk_dwf[Ls/2] and multiplying by
  //
  // T^{Ls/2}+T^{-Ls/2} Psi should restore the Dov(m)
  QDPIO::cout << "Reconstructing the overlap surface"<<endl;
  
  LatticeFermion dov_reconstruct=zero;
  LatticeFermion midpoint = bulk_dwf[Ls/2];

  LatticeFermion tmp=midpoint;
  LatticeFermion tmp1;
  for(int s=0;s<Ls/2;s++){
    Tinv(tmp1,tmp,u,mass,m5,1);
    tmp=tmp1;
  }
  dov_reconstruct = tmp;

  tmp=midpoint;
  for(int s=0;s<Ls/2;s++){
    Tinv(tmp1,tmp,u,mass,m5,0);
    tmp=tmp1;
  }
  dov_reconstruct += tmp;
  
  tmp=dov_reconstruct - ovlap_surface;
  QDPIO::cout << "Reconstruction of the surface prop error"<<norm2(tmp)<<endl;
  QDPIO::cout << "Reconstruction of the surface prop dov_reconstruc"<<norm2(dov_reconstruct)<<endl;
  QDPIO::cout << "Reconstruction of the surface prop ovlap_surface"<<norm2(ovlap_surface)<<endl;
  
}

void Tinv(LatticeFermion &sol, 
	  LatticeFermion &src, 
	  multi1d<LatticeColorMatrix> u,
	  double mass,
	  double M5, 
	  int inv)
{
  // Set up BAGEL object
  double residual=1.0e-18;
  int max_iter   = 10000;

  BfmSolver solver = DWFTransfer;
  if ( inv )solver = DWFTransferInv;

  int lx = QDP::Layout::subgridLattSize()[0];
  int ly = QDP::Layout::subgridLattSize()[1];
  int lz = QDP::Layout::subgridLattSize()[2];
  int lt = QDP::Layout::subgridLattSize()[3];

  multi1d<int> procs = QDP::Layout::logicalSize();

  bfmarg bfma;

  //Physics parameters
  bfma.Ls           = 1;
  bfma.M5           = toDouble(M5);
  bfma.mass         = toDouble(mass);
  bfma.precon_5d    = 0;
  bfma.mobius_scale = mobius_scale;
  bfma.solver       = solver;
  bfma.max_iter     = max_iter;
  bfma.residual     = toDouble(residual);
  bfmarg::Threads(1);

  //Geometry
  bfma.node_latt[0] = lx;
  bfma.node_latt[1] = ly;
  bfma.node_latt[2] = lz;
  bfma.node_latt[3] = lt;
  for(int mu=0;mu<4;mu++){
    if (procs[mu]>1) bfma.local_comm[mu] = 0;
    else             bfma.local_comm[mu] = 1;
  }

  QDPIO::cout << "bfm_qdp:: Initialising BAGEL-2 G5D solver for transfer matrix "<<endl;

  bfm_qdp<double> bfm; 

  bfm.init(bfma);
  QDPIO::cout << "bfm_qdp:: Importing gauge field "<<endl;
  bfm.importGauge(u);



  Fermion_t sol_t[2];
  Fermion_t src_t[2];
  for(int cb=0;cb<2;cb++){
    sol_t[cb] = bfm.allocFermion();
    src_t[cb] = bfm.allocFermion();
  QDPIO::cout << "bfm_qdp:: Importing src fermions "<<endl;
    bfm.importPhysicalFermion(src,src_t[cb],cb);
  QDPIO::cout << "bfm_qdp:: Importing sol fermions "<<endl;
    bfm.importPhysicalFermion(src,sol_t[cb],cb);
  }
  double nn = bfm.norm(src_t);
  Printf("src_norm = %le\n",nn);
  bfm.G5D_Transfer(src_t,sol_t);
  nn = bfm.norm(sol_t);
  Printf("sol_norm = %le\n",nn);
  
  for(int cb=0;cb<2;cb++){
    bfm.exportPhysicalFermion(sol,sol_t[cb],cb);
    bfm.freeFermion(sol_t[cb]);
    bfm.freeFermion(src_t[cb]);
  }
  bfm.end();

  // Could get CHROMA linop and verify now. but need 5d solution.
}







