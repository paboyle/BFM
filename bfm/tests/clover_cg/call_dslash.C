#include <chroma.h>
#include <bfm.h>
#include <bfm_qdp.h>
#include <bfm_qdp_wilson.h>

typedef bfm bfm_t;

using namespace Chroma;

typedef LatticeFermion T;
typedef multi1d<LatticeFermion> T5;
typedef multi1d<LatticeColorMatrix> U;
typedef double Float;

//#define Printf if ( QMP_is_primary_node() ) printf
#define Printf //printf
//#define NTHREAD 2
// Lazy - just make global
bfmarg clova;
bfm_t  clov;
void *thr_cg ( void * ptr);
Fermion_t psi_h[2];
Fermion_t chi_h[2];

int main (int argc,char **argv )
{
    Chroma::initialize(&argc,&argv);
 
 
    /********************************************************
     * Command line parsing
     ********************************************************
     */
    if ( argc != 6 )
    {
        Printf("Usage: %s lx ly lz lt nthreads\n All must be even\n",argv[0]);
        Printf("argc is %d\n",argc);
        for ( int i=0; i<argc; i++)
        {
            Printf("%d %s\n",i,argv[i]);
        }
        exit
        (-1);
 
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
    int NTHREAD = atoi(argv[5]);
    int Ls = 1;
 
    Layout::setLattSize(nrow);
    Layout::create();
 
    int lx = QDP::Layout::subgridLattSize()[0];
    int ly = QDP::Layout::subgridLattSize()[1];
    int lz = QDP::Layout::subgridLattSize()[2];
    int lt = QDP::Layout::subgridLattSize()[3];
 
    /********************************************************
     * Setup Clover operator
     ********************************************************
     */
    bfmarg::Threads(NTHREAD);
    bfmarg::Reproduce(0);
    bfmarg::ReproduceChecksum(0);
    bfmarg::ReproduceMasterCheck(1);
    bfmarg::Verbose(1);


    clova.solver = CloverFermion;
    clova.node_latt[0]  = lx;
    clova.node_latt[1]  = ly;
    clova.node_latt[2]  = lz;
    clova.node_latt[3]  = lt;
    clova.verbose=1;
    clova.threads = NTHREAD;
 
    multi1d<int> procs = QDP::Layout::logicalSize();
    Printf("%d dim machine\n\t", procs.size());
    for(int mu=0; mu<4; mu++)
    {
        Printf("%d ", procs[mu]);
        if ( procs[mu]>1 )
        {
            clova.local_comm[mu] = 0;
        }
        else
        {
            clova.local_comm[mu] = 1;
        }
    }
    Printf("\nLocal comm = ");
    for(int mu=0; mu<4; mu++)
    {
        Printf("%d ", clova.local_comm[mu]);
    }
    Printf("\n");

    double M5 = -0.004;
    double mass = 0.004;
    double mq = mass;
 
    clova.precon_5d = 0;
    clova.Ls   = Ls;
    clova.M5   = M5;
    clova.mass = mq;
    clova.Csw  = 0.0;
    clova.max_iter = 10000;
    clova.residual = 1.e-8;
 
    Printf("Initialising bfm operator\n");
    clov.init(clova);
 
    /********************************************************
     * Bagel internal single checkerboard vectors
     ********************************************************
     */
    psi_h[0] = clov.allocFermion();
    psi_h[1] = clov.allocFermion();
 
    chi_h[0] = clov.allocFermion();
    chi_h[1] = clov.allocFermion();
 
    multi1d<LatticeColorMatrix>  u(Nd);
    HotSt(u);
    multi1d< RScalar< Float > > clover_diag(12*lx*ly*lz*lt);
    multi1d< RComplex< Float > > clover_offdiag(30*lx*ly*lz*lt);
    clov.importGauge(u);

    // Make up a gaussian source and a zero result vector
    LatticeFermion psi;
    gaussian(psi);

    LatticeFermion chi;
    chi = zero;

    // Need some clover params
    CloverFermActParams params;
    params.clovCoeffR = Real(1.92);
    params.clovCoeffT = Real(0.57);
    params.Mass=Real(mass);


    for(int rep=0;rep<12;rep++)
    {
 
    multi1d<int> bcs(Nd);
    bcs[0] = bcs[1] = bcs[2] = bcs[3] = 1;
    Handle<FermBC<T,U,U> > fbc(new SimpleFermBC< T, U, U >(bcs));
    Handle<CreateFermState<T,U,U> > cfs( new CreateSimpleFermState<T,U,U>(fbc));
    Handle<FermState<T,U,U> > fs ((*cfs)(u));
 
    // Need a clover operator - reference
    QDPCloverTerm qdp_clov , qdp_invclov;
    qdp_clov.create(fs, params);
 
    qdp_invclov.create(fs,params,qdp_clov);  // make a copy
    qdp_invclov.choles(0);                   // invert the cb=0 part
    qdp_invclov.choles(1);                   // invert the cb=1 part
 
    //clov.importClover(qdp_clov,qdp_invclov);
    qdp_clov.packForBAGEL(clover_diag,clover_offdiag);
    clov.importClover(clover_diag,clover_offdiag, clov.A);
 
    qdp_invclov.packForBAGEL(clover_diag,clover_offdiag);
    clov.importClover(clover_diag,clover_offdiag, clov.Ainv);
 
    /*Import this checkerboard of source field to bagel*/
    Printf("Importing psi field cb \n");
    Float *psi_p = (Float *) &psi.elem(0).elem(0).elem(0).real();
    Float *chi_p = (Float *) &chi.elem(0).elem(0).elem(0).real();
 
    clov.importFermion(psi_p,psi_h[0],0);
    clov.importFermion(psi_p,psi_h[1],1);   

    pthread_t threads[NTHREAD];
 
    Printf("Creating threads\n");
    for(int t=0; t<NTHREAD; t++)
    {
        pthread_create(&threads[t],NULL,thr_cg,NULL);
    }
    Printf("Joining threads\n");
    for(int t=0; t<NTHREAD; t++)
    {
        pthread_join(threads[t],NULL);
    }

    clov.exportFermion(chi_p,chi_h[0],0);
    clov.exportFermion(chi_p,chi_h[1],1);

    LatticeFermion regress;
    UnprecCloverLinOp M(fs,params);
     
    // Check the result
    PlusMinus pm = PLUS;
    M(regress,chi,pm);

    regress = regress - psi;

    QDPIO::cout << "bfm_qdp:: QDP regression check :  |M sol - src| = " << norm2(regress) << endl;
    if ( toDouble( norm2(regress) / norm2(psi) ) > 1.0e-5 ) { 
      QDPIO::cout << "bfm_qdp:: QDP regression check : This is worryingly large - PANIC"<< endl;
      exit(-1);
    } 

    }

    Printf("Done\n");
 
}
 
void *thr_cg ( void * ptr)
{
    #define NITER 1
    Printf("CG thread\n");
    for(int i=0; i<NITER; i++)
    {
        Printf("CG thread norms %le %le \n",clov.norm(chi_h[0]),
               clov.norm(psi_h[0]));
        fflush(stdout);
        clov.axpy(chi_h[0],psi_h[0],psi_h[0],0.0)  ;
        clov.axpy(chi_h[1],psi_h[1],psi_h[1],0.0)  ;
        clov.CGNE(chi_h,psi_h);
        Printf("CG thread done\n");
        fflush(stdout);
    }
    return NULL;
}
 
 
 
 
 
 
 
 
