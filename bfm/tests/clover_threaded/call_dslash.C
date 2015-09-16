#include <chroma.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
 
typedef double Float;
 
//#define Printf if ( QMP_is_primary_node() ) printf
#define Printf //printf
//#define NTHREAD 64
 
#include <bfm.h>
#include <bfm_qdp.h>
#include <bfm_qdp_wilson.h>
 
typedef LatticeFermion T;
typedef multi1d<LatticeColorMatrix> U;
typedef bfm bfm_t;
Fermion_t psi_h;
Fermion_t chi_h;
Fermion_t tmp_h;
int cb,dag;
bfmarg clova;
bfm_t  clov;
void * thr_Mooee(void *p);
void * thr_MooeeInv(void *p);
void * thr_Meo(void *p);
void * thr_Mprec(void *p);
 
using namespace Chroma;
void compare_result (LatticeFermion A, LatticeFermion B, bfm_t & clov);

int timeval_subtract(struct timeval *result, struct timeval *t2, struct timeval *t1)
{
    long int diff = (t2->tv_usec + 1000000 * t2->tv_sec) - (t1->tv_usec + 1000000 * t1->tv_sec);
    result->tv_sec = diff / 1000000;
    result->tv_usec = diff % 1000000;

    return (diff<0);
}

void timeval_print(struct timeval *tv)
{
    char buffer[30];
    time_t curtime;

    printf("%ld.%06ld", tv->tv_sec, tv->tv_usec);
    curtime = tv->tv_sec;
    strftime(buffer, 30, "%m-%d-%Y  %T", localtime(&curtime));
    printf(" = %s.%06ld\n", buffer, tv->tv_usec);
}
 
void      testFermion  (LatticeFermion &psi,  bfm_t   & clov)
{
    multi1d<int> x(4);
 
    int Nspinco=12;
    Float *Psi_p = (Float *) &psi.elem(0).elem(0).elem(0).real();
 
    cout << "testFermion : Checking I get the CHROMA array ordering correct by comparing to peekSite\n";
 
    for ( x[3]=0; x[3]<clov.node_latt[3]; x[3]++ )
    {
        for ( x[2]=0; x[2]<clov.node_latt[2]; x[2]++ )
        {
            for ( x[1]=0; x[1]<clov.node_latt[1]; x[1]++ )
            {
                for ( x[0]=0; x[0]<clov.node_latt[0]; x[0]++ )
                {
                    int xx[4];
                    xx[0]=x[0];
                    xx[1]=x[1];
                    xx[2]=x[2];
                    xx[3]=x[3];
 
                    int print = 0;
 
                    for ( int co=0; co<3; co++ )
                    {
                        for ( int sp=0; sp<4; sp++ )
                        {
 
                            int spco = co+3*sp;
 
                            int reim=0;
                            int cidx = clov.chroma_idx(xx,reim,spco,Nspinco);
                            Fermion ferm = peekSite(psi,x);
                            ColorVector cv = peekSpin(ferm,sp);
                            Complex cp = peekColor(cv,co);
                            if ( Psi_p[cidx] != toDouble(real(cp) ) )
                            {
                                print = 1;
                            }
 
                            reim=1;
                            cidx = clov.chroma_idx(xx,reim,spco,Nspinco);
                            ferm = peekSite(psi,x);
                            cv = peekSpin(ferm,sp);
                            cp = peekColor(cv,co);
                            if ( Psi_p[cidx] != toDouble(imag(cp) ) )
                            {
 
                                print = 1;
 
                            }
 
                        }
                    }
 
                    if ( print )
                        cout <<x[0]<<" "<<x[1]<<" "<<x[2]<<" "<<x[3]<< endl;
 
                    for ( int co=0; co<3; co++ )
                    {
                        for ( int sp=0; sp<4; sp++ )
                        {
 
                            int spco = co+3*sp;
 
                            int reim=0;
                            int cidx = clov.chroma_idx(xx,reim,spco,Nspinco);
                            Fermion ferm = peekSite(psi,x);
                            ColorVector cv = peekSpin(ferm,sp);
                            Complex cp = peekColor(cv,co);
                            if ( print )
                            {
                                cout << "("<<Psi_p[cidx] << "," <<Psi_p[cidx+1] << ") " << cp<<endl;
                            }
 
                        }
                    }
 
                }
            }
        }
    }
 
}
 
int main (int argc,char **argv )
{
    Chroma::initialize(&argc,&argv);
 
    if ( argc != 7 )
    {
        Printf("Usage: %s lx ly lz lt nthreads testcase\n All must be even (except the last)\n",argv[0]);
        Printf("testcase 0 : Meo\n");
        Printf("testcase 1 : Mooee\n");
        Printf("testcase 2 : Mooeeinv\n");
        Printf("testcase 3 : Mprec\n");
        Printf("testcase 4 : Munprec\n");
        Printf("argc is %d\n",argc);
        for ( int i=0; i<argc; i++)
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
    int NTHREAD  = atoi(argv[5]);
    int testcase = atoi(argv[6]);

    struct timeval tvBegin, tvEnd, tvDiff;
 
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
 
    bfmarg::Threads(NTHREAD);
    bfmarg::Reproduce(1);
    bfmarg::ReproduceChecksum(0);
    bfmarg::ReproduceMasterCheck(0);
    bfmarg::Verbose(1);
 
    clova.solver = CloverFermion;
    clova.node_latt[0]  = lx;
    clova.node_latt[1]  = ly;
    clova.node_latt[2]  = lz;
    clova.node_latt[3]  = lt;
    clova.threads = NTHREAD;
 
    multi1d<int> procs = QDP::Layout::logicalSize();
    Printf("%d dim machine\n\t", procs.size());
    for(int mu=0; mu<4; mu++)
    {
        if ( procs[mu]>1 )
        {
            clova.local_comm[mu] = 0;
        }
        else
        {
            clova.local_comm[mu] = 1;
        }
    }
 
    // Make up a random gauge field.
    multi1d<LatticeColorMatrix> u(Nd);
    HotSt(u);
 
    // Make up a gaussian source and a zero result vector
    LatticeFermion psi;
    gaussian(psi);
 
 
    LatticeFermion chi;
    chi = zero;
 
    clova.Ls = 1;
    clova.mass = 0.0;
    clova.Csw = 0.0;
 
    Printf("Initialising bfm operator\n");
    clov.init(clova);
 
    clov.importGauge(u);
    multi1d< RScalar< Float > > clover_diag(12*lx*ly*lz*lt);
    multi1d< RComplex< Float > > clover_offdiag(30*lx*ly*lz*lt);
 
    //testFermion  (psi,clov);
 
    // Need some clover params
    CloverFermActParams params;
    params.clovCoeffR = Real(1.92);
    params.clovCoeffT = Real(0.57);
    params.Mass=Real(0.0);
 
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
 
    tmp_h = clov.allocFermion();
    LatticeFermion chi_qdp;    
 
    EvenOddPrecCloverLinOp D(fs,params);
 
 
    LatticeFermion noise;
    LatticeFermion tmp;
    gaussian(noise);
 
    Printf("bfm operator using vector length %d\n",clov.simd());
    psi_h = clov.allocFermion();
    chi_h = clov.allocFermion();
 
    // cb is cb of result, 1-cb is cb of input field
    for(cb=0; cb<2; cb++)
    {
 
        /*Import this checkerboard of QDP fields to bagel*/
        Printf("Importing psi field cb %d\n",cb);
 
        // Fill the other checkerboard.
        chi = zero;
        chi_qdp = zero;
 
        Float *psi_p = (Float *) &psi.elem(0).elem(0).elem(0).real();
        Float *chi_p = (Float *) &chi.elem(0).elem(0).elem(0).real();
        clov.importFermion(psi_p,psi_h,1-cb);
 
        for(dag=0; dag<2; dag++)
        {
 
            Printf("Checking cb=%d dag=%d\n",cb,dag);
 
            // Check the result
            PlusMinus pm;
            if ( dag ) pm = MINUS;
            else pm = PLUS;
 
            pthread_t threads[NTHREAD];
 
            switch(testcase)
            {
            case 0:
 
                Printf("testcase 0 : Meo\n");

                gettimeofday(&tvBegin, NULL);
 
                for(int t=0; t<NTHREAD; t++)
                {
                    pthread_create(&threads[t],NULL,thr_Meo,NULL);
                }
                for(int t=0; t<NTHREAD; t++)
                {
                    pthread_join(threads[t],NULL);
                }

                gettimeofday(&tvEnd, NULL);
  
                timeval_subtract(&tvDiff, &tvEnd, &tvBegin);
                Printf("bfm : %ld.%06ld\n", tvDiff.tv_sec, tvDiff.tv_usec);

 
                gettimeofday(&tvBegin, NULL);
 
                if(cb==0)
                    D.evenOddLinOp(chi_qdp, psi, pm);
                else
                    D.oddEvenLinOp(chi_qdp, psi, pm);

                gettimeofday(&tvEnd, NULL);
  
                timeval_subtract(&tvDiff, &tvEnd, &tvBegin);
                Printf("chroma : %ld.%06ld\n", tvDiff.tv_sec, tvDiff.tv_usec);
                Printf("chroma flops : %ld\n",D.nFlops());

                clov.exportFermion(chi_p,chi_h,cb);
 
                break;
 
            case 1:
 
                Printf("testcase 1 : Mooee\n");

                gettimeofday(&tvBegin, NULL); 
 
                for(int t=0; t<NTHREAD; t++)
                {
                    pthread_create(&threads[t],NULL,thr_Mooee,NULL);
                }
                for(int t=0; t<NTHREAD; t++)
                {
                    pthread_join(threads[t],NULL);
                }

		gettimeofday(&tvEnd, NULL);
  
                timeval_subtract(&tvDiff, &tvEnd, &tvBegin);
                Printf("bfm : %ld.%06ld\n", tvDiff.tv_sec, tvDiff.tv_usec);

                gettimeofday(&tvBegin, NULL);
                if((1-cb)==0)
                    D.evenEvenLinOp(chi_qdp, psi, pm);
                else
                    D.oddOddLinOp(chi_qdp, psi, pm);

                gettimeofday(&tvEnd, NULL);
  
                timeval_subtract(&tvDiff, &tvEnd, &tvBegin);
                Printf("chroma : %ld.%06ld\n", tvDiff.tv_sec, tvDiff.tv_usec);
                Printf("chroma flops : %ld\n",D.nFlops());

 
                clov.exportFermion(chi_p,chi_h,1-cb);
 
                break;
 
            case 2:
 
                Printf("testcase 2 : Mooeeinv\n");
 
                gettimeofday(&tvBegin, NULL);

                for(int t=0; t<NTHREAD; t++)
                {
                    pthread_create(&threads[t],NULL,thr_MooeeInv,NULL);
                }
                for(int t=0; t<NTHREAD; t++)
                {
                    pthread_join(threads[t],NULL);
                }

		gettimeofday(&tvEnd, NULL);
  
                timeval_subtract(&tvDiff, &tvEnd, &tvBegin);
                Printf("bfm : %ld.%06ld\n", tvDiff.tv_sec, tvDiff.tv_usec);

 
                gettimeofday(&tvBegin, NULL);

                D.evenEvenInvLinOp(chi_qdp, psi, pm);

                gettimeofday(&tvEnd, NULL);
  
                timeval_subtract(&tvDiff, &tvEnd, &tvBegin);
                Printf("chroma : %ld.%06ld\n", tvDiff.tv_sec, tvDiff.tv_usec);
                Printf("chroma flops : %ld\n",D.nFlops());
 
                clov.exportFermion(chi_p,chi_h,1-cb);
 
                break;
 
            case 3:
 
                gettimeofday(&tvBegin, NULL);

                Printf("testcase 3 : Mprec\n");
 
                for(int t=0; t<NTHREAD; t++)
                {
                    pthread_create(&threads[t],NULL,thr_Mprec,NULL);
                }
                for(int t=0; t<NTHREAD; t++)
                {
                    pthread_join(threads[t],NULL);
                }

		gettimeofday(&tvEnd, NULL);
  
                timeval_subtract(&tvDiff, &tvEnd, &tvBegin);
                Printf("bfm : %ld.%06ld\n", tvDiff.tv_sec, tvDiff.tv_usec);

 
                gettimeofday(&tvBegin, NULL);

                D(chi_qdp, psi, pm);

                gettimeofday(&tvEnd, NULL);
  
                timeval_subtract(&tvDiff, &tvEnd, &tvBegin);
                Printf("chroma : %ld.%06ld\n", tvDiff.tv_sec, tvDiff.tv_usec);
                Printf("chroma flops : %ld\n",D.nFlops());
 
                clov.exportFermion(chi_p,chi_h,1-cb);
 
                break;
 
            }
 
            Double n2 = norm2(chi-chi_qdp);
 
            if(toBool( n2 > Double(1.0e-6) ) && clov.node_latt[0] < 5)
                compare_result(chi,chi_qdp,clov);
 
            cout << "|| Bagel - QDP || = "<< n2 << endl;
 
            n2 = norm2(chi);
            cout << "|| Bagel || = "<< n2 << endl;
            n2 = norm2(chi_qdp);
            cout << "|| QDP || = "<< n2 << endl;
 
        }
    }
    Printf("Done\n");
 
}
 
void * thr_Mooee(void *p)
{
 
    int me = clov.thread_barrier();
    //Printf("thr_Mooee : %d %d threads\n",me,NTHREAD);
    clov.Mooee(psi_h,chi_h,dag,1-cb);
 
}
 
void * thr_MooeeInv(void *p)
{
 
    int me = clov.thread_barrier();
    //Printf("thr_MooeeInv : %d %d threads\n",me,NTHREAD);
    clov.MooeeInv(psi_h,chi_h,dag,1-cb);
 
}
 
 
void * thr_Meo(void *p)
{
 
    int me = clov.thread_barrier();
    //Printf("thr_Meo : %d %d threads\n",me,NTHREAD);
    clov.Meo(psi_h,chi_h,cb,dag);
 
}
 
 
void * thr_Mprec(void *p)
{
 
    int me = clov.thread_barrier();
    //Printf("thr_Mprec : %d %d threads\n",me,NTHREAD);
    clov.Mprec(psi_h,chi_h,tmp_h,dag);
 
}
 
 
void compare_result (LatticeFermion A, LatticeFermion B, bfm_t & clov)
{
 
    multi1d<int> x(4);
 
    int Nspinco=12;
 
    Printf("Result vectors \n");
    for ( x[3]=0; x[3]<clov.node_latt[3]; x[3]++ )
    {
        for ( x[2]=0; x[2]<clov.node_latt[2]; x[2]++ )
        {
            for ( x[1]=0; x[1]<clov.node_latt[1]; x[1]++ )
            {
                for ( x[0]=0; x[0]<clov.node_latt[0]; x[0]++ )
                {
 
                    int xx[4];
                    xx[0]=x[0];
                    xx[1]=x[1];
                    xx[2]=x[2];
                    xx[3]=x[3];
 
                    Printf("site %d %d %d %d\n",x[0],x[1],x[2],x[3]);
                    for ( int sp=0; sp<4; sp++ )
                    {
                        for ( int co=0; co<3; co++ )
                        {
 
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
 
                        }
                    }
 
 
                }
            }
        }
    }
 
 
}
 
 
 
 
 
 
