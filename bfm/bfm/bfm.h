/****************************************************************************/
/* Mar/2007                                                                 */
/* Peter Boyle.                                                             */
/*                                                                          */
/* bfm.h                                                                   */
/* C++ header file DWF & Wilson CG inverter with SIMD support               */
/*                                                                          */
/****************************************************************************/

#ifndef INCLUDED_DWF_PAB_H
#define INCLUDED_DWF_PAB_H

#undef PLL_HACK
#define BAGEL_THREADS (64)
#include <math.h>
#include <complex>
#include <mcheck.h>

#include <string.h>

#include "bfm_options.h"
#include <bagel_int.h>
#include <stdlib.h>
#include <vector>

#undef BGQ
#ifndef BGQ
uint64_t GetTimeBase(void) ;
#endif

#ifdef BFM_QDP_BINDING
#include <qdp.h>
#define CHROMA_3
#endif


#define BFM_REPRO
#define bzero(A,B) memset(A,0,B);

typedef double QDPdouble ;
typedef void * Fermion_t;
typedef void * Matrix_t;

typedef Fermion_t bfm_fermion [2];


enum BfmMemType   { mem_fast, mem_sendbuf, mem_slow };
enum BfmReproType { ReproIdle, ReproRecord , ReproCheck };
  /* Constant definitions */
  enum {
    ND  = 4,
    SPINOR_SIZE = 24,
    HALF_SPINOR_SIZE = 12,
    GAUGE_SIZE = 72,
    CLOVER_MAT_SIZE = 42,
    HALF_CLOVER_MAT_SIZE = 21,
    Nmu = 4,
    Ncb = 2,
    NMinusPlus = 2,
    Minus=0,
    Plus =1,
    DaggerYes=1,
    DaggerNo =0,
    SingleToDouble=1,
    DoubleToSingle=0,
    Odd=1,
    Even=0
  };


#include <bfm_thread.h>
#include <bfmActionParams.h>

enum BfmVerbosity { BfmDebug=8,BfmPerformance=4,BfmMessage = 2, BfmError=1 };

class bfmarg : public bfmActionParams {

public:
  ////////////////////////////
  // Machine information
  ////////////////////////////
   static int threads;
   int node_latt[4];  /*Node subvolume*/
   int ncoor[4];      /*location of node in 4d node grid*/
   int local_comm[4]; /*are comms periodic/trivial in this axis      */
   
   //////////////////////////////////////////////////////
   // Controls the solver iteration and performance reporting
   /////////////////////////////////////////////////////
   static void Verbose(int def) { default_verbose = def; }
   static int default_verbose;
   int verbose;
   int rb_precondition_cb;
   int time_report_iter;
   int max_iter;      // CG control
   double residual;   // CG control

   ///////////////////////////////////////////////////////////////////
   // Reproduction control. Watchfile allows a job termination request
   ///////////////////////////////////////////////////////////////////
   static char * watchfile;
   static int reproduce;
   static int reproduce_checksum;
   static int reproduce_master_check;
   
   ///////////////////////////////////////////////////
   // Infrequently changed matrix convention controls
   ///////////////////////////////////////////////////
   static int part_frac_chroma_convention;
   static int onepluskappanorm;
   static int CGdiagonalMee;
   /////////////////////////////////////////////////////
   // Accessors for the static behavioural controls
   /////////////////////////////////////////////////////
   static void UseCGdiagonalMee(int useit){ CGdiagonalMee = useit;};
   static void WatchFile(char * file) {watchfile=file;};
   static void ReproduceChecksum(int rp) { reproduce_checksum = rp ;};
   static void ReproduceMasterCheck(int rp) { reproduce_master_check = rp ;};
   static void Reproduce(int r) { reproduce = r ;};

   static void Threads  (int t) { threads =t;}
   static void OnePlusKappaNorm  (int yesno) { onepluskappanorm = yesno ; }

   bfmarg(){
     time_report_iter = -100;
     ncoor[0]=ncoor[1]=ncoor[2]=ncoor[3]=0;
     rb_precondition_cb = Odd;
     verbose=default_verbose;
    };

 };

 //////////////////////////////////////////////////////////////////////////
 // Bagel fermion matrix class
 //////////////////////////////////////////////////////////////////////////
 // Internally the code works only with a logical/data-parallel subnodes.
 // Interaction between logical subnodes only happens during comms
 // This class operates on "internal" layout vectors

 template <class Float>
 class bfmbase : public bfmarg , public ThreadModel {

  public:
   //////////////////////////////////////
   // Internal copy of gauge field
   //////////////////////////////////////
   Float *u;
   // Clover
   Float *A, *Ainv;

   /////////////////////////////////////////////////
   // Node vol / Layout information
   // Integer data tables for faces and comms buffers
   /////////////////////////////////////////////////
   integer ncb; // 2
   integer base_parity;// Worked out from ncoor
   integer node_cbvol; // Can either reduce Ls, or reduce vol by ncb
   integer cbLs;       // precon_5d controls this

   integer nsimd;          // Simd reordering of data
   integer simd_cbvol;
   integer simd_latt[ND];
   integer simd_grid[ND];
   integer simd_nbound[ND];
   integer simd_allbound;

   // Use to support next-nearest neighbour Hamber-Wu-Eguchi-Kawamoto
   integer halo_depth; 

   // Buffers
   Float *sendbufs[8];
   Float *recvbufs[8];

   // Table for assembler kernels
   //Densely packed tables, store CB info and is_halo mask
   short_integer *shift_table[Ncb];
   short_integer *shift_table_half[Ncb];
   std::vector<short_integer> lebesgue_sites;

   integer *cb_table[Ncb];
   integer *face_table[Ncb][NMinusPlus][ND];
   int      comm_offset[NMinusPlus][ND];
   double *complex_i_simd;  // Used to hold i,-i, 0 etc for use in assembler kernels

   int precision_test;

   // Output wrappers

   void ThreadBossDebug(const char *fmt,...);   // Printf only on boss processor & thread 0
   void ThreadBossPerformance(const char *fmt,...);   // Printf only on boss processor & thread 0
   void ThreadBossMessage(const char *fmt,...); // Printf only on boss processor & thread 0
   void ThreadBossLog(const char *fmt,...); // Printf only on boss processor & thread 0
   void BossDebug(const char *fmt,...);   // Printf only on boss processor
   void BossMessage(const char *fmt,...); // Printf only on boss processor
   void BossLog(const char *fmt,...); // Printf only on boss processor
   void NodeDebug(const char *fmt,...);   // Printf only on boss processor
   void NodeMessage(const char *fmt,...); // Printf only on boss processor
   void Error(const char *fmt,...);       // Printf only on boss processor

   ////////////////////////////////////////////////////////
   // Initialisation
   ////////////////////////////////////////////////////////
   void init(bfmarg &arg);
   void end (void) ;

   ////////////////////////////////////////////////////////
   // Layout
   ////////////////////////////////////////////////////////
   int simd(void) const 
   {
     return BAGEL_NSIMD;
   }
   void pointers_init(void) ;
   void pointers_end(void) ;
   void lebesgue_order(void);

   int  parity(int x[4]);
   int  psite (int x[4],integer latt[4]);
   int  interleave_site(int pm,int mu, int site);

 inline integer bagel_idx5d (int simd_site,int s,int nspinco) 
 {
   if ( this->precon_5d ) s = s/2;
   return s*nspinco*this->simd()*2 + simd_site*nspinco*this->simd()*2*this->cbLs;
 }
 inline
 integer bagel_idx5d_tmp (int lls, int cbsite,int cbs,int reim,int i,int i_size) 
 {
   // Take 4d bagel ordered site index, + s index and map to 5d bagel ordered.
   int nvec = BAGEL_NSIMD; 
   int msk  = nvec-1;
   int ssite   = cbsite/nvec;
   int simdlane= cbsite&msk;
   return reim
	+ 2*simdlane
	+ 2*nvec*i
	+ 2*nvec*i_size*cbs
	+ 2*nvec*i_size*lls*ssite;
 }
 inline
 integer bagel_idx5d_tmp (int cbsite,int cbs,int reim,int i,int i_size) 
 {
   return this->bagel_idx5d_tmp (this->cbLs,cbsite,cbs,reim, i, i_size) ;
 }

 inline integer bagel_idx5d (int x[4],int s,int reim,int i, int i_size, int half_cb)
 {
   if ( this->precon_5d ) s = s/2;
   return this->bagel_idx(x,reim,i+s*i_size,i_size*this->cbLs,half_cb);
 }
 inline 
 integer bagel_idx (int x[4],int reim,int i, int i_size, int half_cb)
 {
   int gx[4];
   int mu;

   for (mu=0;mu<4;mu++ ) gx[mu] = x[mu]/this->simd_latt[mu];

   integer simd_el = gx[0]+this->simd_grid[0]*( gx[1]
			  +this->simd_grid[1]*( gx[2]
			  +this->simd_grid[2]*  gx[3] ));
   int scoor[4]; // SIMD coordinate

   scoor[0] = x[0] % this->simd_latt[0]; /*Work out our location within the SIMD vol*/
   scoor[1] = x[1] % this->simd_latt[1];
   scoor[2] = x[2] % this->simd_latt[2];
   scoor[3] = x[3] % this->simd_latt[3];

   int scb   = this->parity(scoor);          /*Work out local checkerboard of site*/
   int ssite = this->psite(scoor,this->simd_latt); /*Work out checkerboarded site number*/

   // Real/Imaginary innermost
   int idx = reim + simd_el*2
		  + i*simd()*2
       + ssite*i_size*simd()*2;

   if ( !half_cb ) idx+= scb*this->simd()*2*i_size*this->simd_cbvol;

   return idx;
 }

   integer chroma_idx(int site[4],int reim,int i, int i_size) ; // for import
   integer cps_idx(int site[4],int reim,int i, int i_size) ; // for import


   //////////////////////////////////////////////////////////////////
   // Comms related
   //////////////////////////////////////////////////////////////////
   void gather(int result_cb, Fermion_t psi,int dag) ; 
   void gather(int mu,int result_cb, Fermion_t psi,int dag) ; 

   virtual void comm_init(void) = 0;
   virtual void comm_end(void) = 0;
   virtual void comm(int result_cb, Fermion_t psi,int dag) = 0; 
   virtual void comm_start(int result_cb, Fermion_t psi,int dag) = 0; 
   virtual void comm_complete(int result_cb, Fermion_t psi) = 0; 
   virtual void comm_gsum(double &val);
   virtual void comm_gsum(double *val,int N);
   virtual void comm_csum(uint64_t &csum) {}; // Place holder
   virtual bool isBoss();
   virtual int  SPIcomms(void) { return 0; };
   void check_is_same(double d, char *s);


   /////////////////////
   // Linear algebra
   /////////////////////
   void axpby_ssp_internal(Fermion_t out,double al, double ah,Fermion_t x,
			   double bl, double bh,Fermion_t y,int sxo,int sy);
   void axpby_ssp     (Fermion_t out,double a,Fermion_t x,double b,Fermion_t y,int sxo,int sy);
   void axpby_ssp_proj(Fermion_t out,double a,Fermion_t x,double b,Fermion_t y,int sxo,int sy,int psign);
   void ag5xpby_ssp_proj(Fermion_t out,double a,Fermion_t x,double b,Fermion_t y,int sxo,int sy,int psign);

   void ag5xpby_ssp   (Fermion_t out,double a,Fermion_t x,double b,Fermion_t y,int sxo,int sy);
   void axpbg5y_ssp   (Fermion_t out,double a,Fermion_t x,double b,Fermion_t y,int sxo,int sy);
   void ag5xpbg5y_ssp (Fermion_t out,double a,Fermion_t x,double b,Fermion_t y,int sxo,int sy);


   ////////////////////////////////////////////////////////////
   // Opaque type fermion vector interface:
   // Allocation, precision change, import export, linear operations
   ////////////////////////////////////////////////////////////

   // Block support for inexact deflation
   void      block_pick (Fermion_t out,Fermion_t match,Fermion_t nomatch,int block,int *blockids);
   void      block_norm (Fermion_t in,int nblock,double *result,int *blockids);
   void      block_normalise(Fermion_t in,int nblock,double *result,int *blockids);
   void      block_inner(Fermion_t in1,Fermion_t int2,int nblock,double *result,double *reduce,int *blockids);
   void      block_inner(Fermion_t in1,Fermion_t int2,int nblock,double *result,double *reduce,int *blockids,int smin,int smax);

   template <class cFloat>
   void      block_project(Fermion_t *basis,int nbasis,Fermion_t in,
			   int nblock,cFloat *result,
			   double *reduce,int *blockids,int smin,int smax)  __attribute__ ((optimize("-O6")));
   template <class cFloat>
   void      block_promote(Fermion_t result,Fermion_t *basis, int nbasis,
			   cFloat *zap,int *blockids,int smin,int smax) __attribute__ ((optimize("-O6")));
   template <class cFloat>
   void      block_zaxpy(Fermion_t result,Fermion_t x,Fermion_t y, double a, cFloat *zap,int *blockids);

   template <class cFloat>
   void      block_zaxpy(Fermion_t result,Fermion_t x,Fermion_t y, double a, cFloat *zap,int *blockids,int smin,int smax);

   // Linear algebra, scaling, filling etc...
   void      master_fill(Fermion_t psi,double a);
   double    norm(Fermion_t psi);
   double    local_norm(Fermion_t psi);
   void      fill(Fermion_t psi,double a);
   void      set_zero(Fermion_t psi);
   void      scale(Fermion_t psi,double a);
   void      scale(Fermion_t psi,double ar,double ai);
   double   inner_real(Fermion_t x_t, Fermion_t y_t);
   std::complex<double>   dot(Fermion_t x_t, Fermion_t y_t,int sx,int sy);
   std::complex<double>   dot(Fermion_t x_t, Fermion_t y_t);
   void zaxpy(Fermion_t r,Fermion_t x, Fermion_t y, std::complex<double> z){
     caxpy(r,x,y,real(z),imag(z));
   };
   void      caxpy(Fermion_t result, Fermion_t x, Fermion_t y,double a,double b);
   void      caxpby(Fermion_t result, Fermion_t x, Fermion_t y,double a,double b,double c, double d);
   void      conj(Fermion_t psi);
   void      copy(Fermion_t psi,Fermion_t copyme);
   void      axpy(Fermion_t result, Fermion_t x, Fermion_t y,double a);
   void      axpby(Fermion_t result, Fermion_t x, Fermion_t y,double a,double b);
   double    axpy_norm(Fermion_t result, Fermion_t x, Fermion_t y,double a);
   double    cg_psi_p (Fermion_t psi, Fermion_t p, Fermion_t r,Fermion_t mmp,double a,double b);
   double    axpby_norm(Fermion_t result, Fermion_t x, Fermion_t y,double a,double b);
   void      axpibg5x(Fermion_t r, Fermion_t x, double a,double b);
   void      merge(Float *r, Float *lo, Float *hi, int nsite);
   void four_to_five(Fermion_t input[2], Fermion_t result5d[2]);
   void five_to_four(Fermion_t input5d[2], Fermion_t result[2], int sgn);

   /////////////////////////////////////////////
   // Both parity linalg for unprec CG support
   //////////////////////////////////////////////
   double normMat(Matrix_t vec);
   double norm(Fermion_t vec,int s);
   void   UniformRandom(Fermion_t vec);
   void   copy_slice(Fermion_t in_t ,int si,int Lsi,Fermion_t out_t,int so,int Lso);

   double norm(Fermion_t vec[2]){
     return norm(vec[0])+norm(vec[1]);
   }
   double norm(Fermion_t vec[2],int s){
     return norm(vec[0],s)+norm(vec[1],s);
   }
   void axpy(Fermion_t z[2],Fermion_t x[2],Fermion_t y[2],double a){
     axpy(z[0],x[0],y[0],a);
     axpy(z[1],x[1],y[1],a);
   }
   void axpby(Fermion_t z[2],Fermion_t x[2],Fermion_t y[2],double a,double b){
     axpby(z[0],x[0],y[0],a,b);
     axpby(z[1],x[1],y[1],a,b);
   }
   void scale(Fermion_t z[2],double a){
     scale(z[0],a);
     scale(z[1],a);
   }


   // Allocation
   Matrix_t allocMatrix   (void); 
   void     zeroMatrix   (Matrix_t handle); 
   void     freeMatrix    (Matrix_t handle);
   Fermion_t allocFermion   (int mem_type=mem_slow); 
   void      freeFermion    (Fermion_t handle);
   Fermion_t threadedAllocFermion   (int mem_type=mem_slow); 
   void      threadedFreeFermion    (Fermion_t handle);
   //copy from Qi
   Fermion_t allocCompactFermion   (int mem_type=mem_slow); 
   Fermion_t threadedAllocCompactFermion   (int mem_type=mem_slow); 


   // Precision change
   int  precision(void) { return sizeof(Float); };
   void precisionChange (Fermion_t in, Fermion_t out, int control,int cb);

   // Try to precondition with reduced Ls and smarter handling of center
   // View middle as a generalised multigrid aggregate consisting of the
   // (weighted) average of a spin/color across a single site block in 5th dimension. 
   //
   // View M = Mee Meo
   //          Moe Moo
   //
   // as having defined and controlled matrix elements on these aggregates, and
   // represent M exactly on these blocks.
   // 
   // MdagM is what is inverted, but since M is faithfully represented perhaps MdagM
   // is not so bad. Can same trick be used to represent M in even-odd prec as Ldop, and
   // then MdagM is defined?
   //
   // Perhaps then two smaller ldop's are applied and much less non-local comms?
   //
   //   void ProjReduceLs (int Block, std::vector<double> &Weights,Fermion_t in, Fermion_t out, Fermion_t tmp);
   //   void PromIncreaseLs (int Block, std::vector<double> &Weights,Fermion_t in, Fermion_t out,Fermion_t tmp);

   //Import export
   //        Slow                              Fast
   // double psi   [x%2][t][z][y][x] [spin][color][realimag]
   // double gauge [x%2][t][z][y][x] [row][column][realimag]
   // Must call import gauge for EIGHT mu's.
   //
   void impexFermion  (QDPdouble *psi, Fermion_t handle, int import, int cb,int s) ;
   void      importFermion  (QDPdouble *psi, Fermion_t handle,int cb);
   void      exportFermion  (QDPdouble *psi, Fermion_t handle,int cb);
   void exportForce(Matrix_t handle, QDPdouble *force, int mu,int cb) ;
   virtual void importGauge(QDPdouble *gauge, int dir) ;
   virtual void exportGauge(QDPdouble *gauge, int dir) ;
   virtual int doubleStored(void);

  template <typename FloatEXT>  void cps_importGauge(FloatEXT *importme);
 #ifdef BFM_QDP_BINDING
   void exportForce(Matrix_t handle,multi1d<LatticeColorMatrix> &u,int cb) ;
   virtual void   importGauge(multi1d<LatticeColorMatrix> &u);
   virtual void   exportGauge(multi1d<LatticeColorMatrix> &u,multi1d<LatticeColorMatrix> &udag);
   virtual void   importGaugeNN(multi1d<LatticeColorMatrix> &u);//next nearest
   void           importFermion  (LatticeFermion &psi, Fermion_t handle,int cb,int s=0);
   void           exportFermion  (LatticeFermion &psi, Fermion_t handle,int cb,int s=0);
   void           importClover(multi1d< RScalar <REAL> > &clover_diag, multi1d< RComplex <REAL> > &clover_offdiag, Float *bagel_clover);
   // Handles DWF and G5D surface terms
   void           importPhysicalFermion(LatticeFermion &psi, Fermion_t handle,int cb);
   void           exportPhysicalFermion(LatticeFermion &psi, Fermion_t handle,int cb);
   // for DWF
   void      importFermion  (multi1d<LatticeFermion> &psi, Fermion_t handle,int cb);
   void      exportFermion  (multi1d<LatticeFermion> &psi, Fermion_t handle,int cb);
   complex<double> inner(Fermion_t x,Fermion_t y);
 #endif


   /////////////////////////////////
   // Abstract matrix routines
   /////////////////////////////////
   void Munprec(Fermion_t chi[2], 
		Fermion_t psi[2], 
		Fermion_t tmp,
		int dag) ;

   virtual double Mprec(Fermion_t chi, 
			Fermion_t psi, 
			Fermion_t tmp, 
			int dag,int donrm=0);

   virtual double MprecTilde(Fermion_t chi, 
			     Fermion_t psi, 
			     Fermion_t tmp, 
			     int dag,int donrm=0);

   //add by Jianglei
   virtual double CompactMprec(Fermion_t compact_chi, 
			       Fermion_t compact_psi, 
			       Fermion_t chi,
			       Fermion_t psi,
			       Fermion_t tmp, 
			       int dag,int donrm=0);

   //add by Jianglei
   void CompactMunprec(Fermion_t compact_chi[2], 
		       Fermion_t compact_psi[2], 
		       Fermion_t chi[2],
		       Fermion_t psi[2],
		       Fermion_t tmp,
		       int dag);

   void Meo(Fermion_t chi, 
	    Fermion_t psi,
	    int result_cb,
	    int dag);

   void Mooee(Fermion_t chi, 
	      Fermion_t psi, 
	      int dag,int cb=0) ;

   void MooeeInv(Fermion_t chi, 
		 Fermion_t psi, 
		 int dag,int cb=0) ;

   /////////////////////////////////
   // HMC support
   // Useful force terms for preconditioned odd-odd operator
   /////////////////////////////////
   void MprecDeriv(Fermion_t X, 
		   Fermion_t Y, 
		   Matrix_t  Force[2],
		   int dag,double coeff=1.0);
   void MeoDeriv(Fermion_t X, 
		 Fermion_t Y, 
		 Matrix_t  Force,
		 int cb,
		 int dag);

   void DslashDeriv(Fermion_t X, 
		      Fermion_t Y, 
		      Matrix_t  Force,
		      int cb,
		      int dag);

   ////////////////////////////////////////////////////////////////////////////////////
   // Some solver logging facilities
   ////////////////////////////////////////////////////////////////////////////////////

   // Timing controls
   // Count matrix ops
   uint64_t InversionEnteredUsec;
   uint64_t InversionMprecCount;
   uint64_t InversionDwCount;

   // Additional cost controls
   int LoggingSuppressed;
   int InverterLogTrueResidual;
   int InverterCompareExactSolution;
   Fermion_t InverterExactSolution;

   // Where to stick info
   std::string InverterLogFile;
   FILE *InverterLogFP;

   void InverterLoggerInit(void){
     InverterLogTrueResidual = 0;
     InverterCompareExactSolution=0;
     InverterLogFP=NULL;    
     LoggingSuppressed=1;
   };

   void InverterLoggingBegin(std::string);
   void InverterLoggingEnd(void);

   uint64_t Nanosecond() {return (1000*GetTimeBase())/MHz();}
   uint64_t Microsecond() {return GetTimeBase()/MHz();}
   uint64_t CpuClock() {return GetTimeBase(); }
   int MHz(void) { return 1600;}

   void InverterEnter(void);
   void InverterExit (void);
   void InverterRegisterExactSolution(Fermion_t solution);
 void InverterLogIteration(int iter, double residual, Fermion_t src, Fermion_t current_solution,Fermion_t tmp1,Fermion_t tmp2,Fermion_t mtmp);

   ////////////////////////////////////////////////////////////////////////////////////
   // Solvers
   ////////////////////////////////////////////////////////////////////////////////////
 #ifdef BFM_QDP_BINDING
   int CGNE_M(Fermion_t sol_guess[2], Fermion_t source[2], multi1d<bfm_fermion> &evecs, multi1d<double> &evals);
   int CGNE(Fermion_t sol_guess[2], Fermion_t source[2], multi1d<bfm_fermion> &evecs, multi1d<double> &evals);
   int MCR(Fermion_t sol_guess[2], Fermion_t source[2], multi1d<bfm_fermion> &evecs, multi1d<double> &evals);
 #endif
   int CGNE_M(Fermion_t sol_guess[2], Fermion_t source[2], int Nev,bfm_fermion *evecs, double *evals);
   int CGNE(Fermion_t sol_guess[2], Fermion_t source[2]  , int Nev,bfm_fermion *evecs, double *evals);
   int MCR(Fermion_t sol_guess[2], Fermion_t source[2]   , int Nev,bfm_fermion *evecs, double *evals);

   int CGNE_M(Fermion_t sol_guess[2], Fermion_t source[2]);
   int CGNE_Mdag(Fermion_t sol_guess[2], Fermion_t source[2]);
   int CGNE_prec_MdagM(Fermion_t sol_guess, Fermion_t source);
   int CGNE_prec_MdagM_multi_shift(Fermion_t sol_guess[], 
				   Fermion_t source,
				   double    shifts[],
				   double    alpha[],
				   int       nshift,
				   double mresidual[],
				   int single, // sum single result
				   int forcecontinue=0);


   int CGNE_MdagM(Fermion_t sol_guess[2], Fermion_t source[2]);
   int CGNE(Fermion_t sol_guess[2], Fermion_t source[2]);

   int CGNE_prec(Fermion_t sol_guess, Fermion_t source); // Internal
   int CGNE_prec(Fermion_t psi, Fermion_t src,Fermion_t *residuals,int nresid); // Internal
   int CGNE_single_shift(Fermion_t psi, Fermion_t src,double shift,std::vector<double> & polynomial);

   int EIG_CGNE_M(Fermion_t solution[2], Fermion_t source[2]);
   int Eig_CGNE_prec(Fermion_t psi, Fermion_t src);

   int MCR(Fermion_t sol_guess[2], Fermion_t source[2]);
   int MCR_G5R(Fermion_t sol_guess[2], Fermion_t source[2]);
   int MCR_prec(Fermion_t sol_guess, Fermion_t source); // Internal
   int MCR_single_shift(Fermion_t sol_guess, Fermion_t source,double shift); // Internal
   int MCR_prec(Fermion_t sol_guess, Fermion_t source,Fermion_t *residuals,int nresid); // Internal

   int MCR_PolyPrec(Fermion_t sol_guess, Fermion_t source); // Internal

   int MCR_prec_G5R(Fermion_t sol_guess, Fermion_t source); // Internal
   void G5R(Fermion_t input, Fermion_t &result);

   // Rudy implemented PauliVillars by fourier transform in 5th dim
   // If we do this better would NOT replicate th CG routine
   int CGNE_PV_M(Fermion_t solution[2], Fermion_t source[2],int n, int L5);
   int CGNE_PV_Mdag(Fermion_t solution[2], Fermion_t source[2],int n, int L5);
   int CGNE_PV_prec_M(Fermion_t psi, Fermion_t src,int n, int L5);
   int CGNE_PV_prec(Fermion_t psi, Fermion_t src,int n, int L5);

   int CGNE_Mdag_unprec(Fermion_t sol_guess[2], Fermion_t source[2]);
   int CGNE_MdagM_unprec(Fermion_t sol_guess[2], Fermion_t source[2]);
   int CGNE_M_unprec(Fermion_t sol_guess[2], Fermion_t source[2]);
   int CGNE_unprec(Fermion_t sol_guess[2], Fermion_t source[2]);


   ///////////////////////////////////////////
   // Implementation specific subroutines
   ///////////////////////////////////////////
   // Action dependent Dslash implementation routines
   // This does a lot of "if/else" but splitting into 
   // distinct virtual classes would also replicate code and
   // be harder to maintain (ala Chroma).
   void bgq_l1p_optimisation(int mode);

   /*
   void gather_proj(integer *sites,
		    Float     *psi,
		    integer comm_off, 
		    Float  *comm_buf,
		    integer nface,
		    int permute, 
		    int Ls,
		    int pm,
		    int mu
		    );
   */
   virtual void dslash(Fermion_t chi, 
		       Fermion_t psi, 
		       int cb,
		       int dag) ;

   virtual double dslash_generic(Fermion_t chi, 
		       Fermion_t psi_in, 
		       Fermion_t psi_out, 
		       int cb,
		       int dag,int dperp=1, int addscale=0, double a=0,double b=0) ;

   virtual double dslash_generic_serial(Fermion_t chi, 
					Fermion_t psi_in, 
					Fermion_t psi_out, 
					int cb,
					int dag,int dperp=1, int addscale=0, double a=0,double b=0) ;
   virtual double dslash_generic_sloppy(Fermion_t chi, 
					Fermion_t psi_in, 
					Fermion_t psi_out, 
					int cb,
					int dag, int addscale=0, double a=0,double b=0) ;
   // Clover term
   void   applyClover(Fermion_t psi,Fermion_t chi,int cb,int inv=0);
   void   applyCloverCpp(Fermion_t psi,Fermion_t chi,int cb,int inv=0);

   // DWF routines
   void Mooee5dprec(Fermion_t psi, 
			    Fermion_t chi, 
			    int dag);
   void MooeeInv5dprec(Fermion_t psi, 
		     Fermion_t chi, 
		     int dag);
   void Mooee4dprec(Fermion_t psi, 
		    Fermion_t chi, 
		    int dag);
   void MooeeInv4dprec(Fermion_t psi, 
		     Fermion_t chi, 
		     int dag);

   void MunprecTW(Fermion_t chi[2], 
		  Fermion_t psi[2], 
		  Fermion_t tmp,
		  int n,int L5,
		  int dag) ;

  double MprecTW(Fermion_t chi, 
		 Fermion_t psi, 
		 Fermion_t tmp, int n, int L5,
		 int dag,int donrm=0) ;

   void Mooee5dprec_TW(Fermion_t psi_t, 
				   Fermion_t chi_t, int n, int L5,
				   int dag);

   void MooeeInv5dprec_TW(Fermion_t psi_t, 
				      Fermion_t chi_t, int n, int L5, 
				      int dag);


   /////////////////////////////////////////////////////////////////
   // GeneralisedFiveDim support routines
   /////////////////////////////////////////////////////////////////
   Fermion_t g5d_tmp;
   void GeneralisedFiveDimEnd(void);
   void GeneralisedFiveDimInit(void);
   int IsGeneralisedFiveDim(void){
     if( solver == HwPartFracZolo   ) return 1;
     if( solver == HwContFracZolo   ) return 1;
     if( solver == HwPartFracTanh   ) return 1;
     if( solver == HwContFracTanh   ) return 1;
     if( solver == HwCayleyTanh      ) return 1;
     if( solver == HtCayleyTanh      ) return 1;
     if( solver == HmCayleyTanh      ) return 1;
     if( solver == HwCayleyZolo      ) return 1;
     if( solver == HtCayleyZolo      ) return 1;
     if( solver == DWFTransfer       ) return 1;
     if( solver == DWFTransferInv    ) return 1;
     return 0;
   }
   char *SolverString( BfmSolver solver )
   {
     if(solver ==DWF)           return (char *) "DWF";
     if(solver ==DWFrb4d)       return (char *) "DWFrb4d";
     if(solver ==WilsonFermion) return (char *) "WilsonFermion";
     if(solver ==CloverFermion) return (char *) "CloverFermion";
     if(solver ==WilsonTM)      return (char *) "WilsonTM";
     if(solver ==WilsonNN)      return (char *) "WilsonNN";
     ////////////////////////////////////
     // Group of generalised 5d actions
     // Kerne/Representation/approximation Combos
     ////////////////////////////////////
     if(solver == HwPartFracZolo) return (char *) "HwPartFracZolo";
     if(solver == HwContFracZolo) return (char *) "HwContFracZolo";
     if(solver == HwPartFracTanh) return (char *) "HwPartFracTanh";
     if(solver == HwContFracTanh) return (char *) "HwContFracTanh";
     if(solver == HwCayleyZolo) return (char *) "HwCayleyZolo";
     if(solver == HtCayleyZolo) return (char *) "HtCayleyZolo";
     if(solver == HwCayleyTanh) return (char *) "HwCayleyTanh";
     if(solver == HmCayleyTanh) return (char *) "HmCayleyTanh";
     if(solver == HtCayleyTanh) return (char *) "HtCayleyTanh";
     if(solver == DWFTransfer) return (char *) "DWFTransfer";
     if(solver == DWFTransferInv) return (char *) "DWFTransferInv";
     printf("Unknown solver type");
     exit(-1);
   };

   private:
   double eps;
   public:
   zolotarev_data *zdata;

   // Part frac
   std::vector<double> p; 
   std::vector<double> q;

   // Cont frac
   std::vector<double> Beta;
   std::vector<double> cc;;
   std::vector<double> cc_d;;
   std::vector<double> sqrt_cc;

   // Cayley form Moebius (tanh and zolotarev)
   std::vector<double> omega; 
   std::vector<double> bs;    // S dependent coeffs
   std::vector<double> cs;    
   std::vector<double> as;    


  // For preconditioning Cayley form
  std::vector<double> bee;    
  std::vector<double> cee;    
  std::vector<double> aee;    
  std::vector<double> beo;    
  std::vector<double> ceo;    
  std::vector<double> aeo;    
  // LDU factorisation of the eeoo matrix
  std::vector<double> lee;    
  std::vector<double> leem;    
  std::vector<double> uee;    
  std::vector<double> ueem;    
  std::vector<double> dee;    
  std::vector<double> See;
  std::vector<double> Aee;

  void G5D_gnuplot_approx(void); 

  void G5D_DW(Fermion_t psi[2],Fermion_t chi,int cb,int dag,double scale=1.0); 
  void G5D_HW(Fermion_t psi[2],Fermion_t chi,int cb,int dag); 
  void G5D_Dminus(Fermion_t psi[2],Fermion_t chi[2],int dag);
  void G5D_Dplus(Fermion_t psi[2],Fermion_t chi[2],int dag);

  void G5D_ConsvCurrents(Fermion_t sol[2]);

  void G5D_Transfer   (Fermion_t psi[2],Fermion_t chi[2]); 
  void G5D_TransferNumDen(Fermion_t psi[2],Fermion_t chi[2],int sgn,int dag);  

  void G5D_Munprec(Fermion_t psi[2], 
		   Fermion_t chi[2], 
		   int dag);

  double G5D_Mprec(Fermion_t chi, 
		   Fermion_t psi, 
		   Fermion_t tmp, 
		   int dag,int donrm=0);

  //  void G5D_MprecDeriv(Fermion_t chi, 
  //			Fermion_t psi, 
  //			Matrix_t force[2], 
  //			int dag);


  void G5D_Munprec_Cayley(Fermion_t psi[2], 
			    Fermion_t chi[2], 
			    int dag);

  void G5D_Munprec_PartFrac(Fermion_t psi[2], 
				  Fermion_t chi[2], 
				  int dag);
  void G5D_Munprec_ContFrac(Fermion_t psi[2], 
			    Fermion_t chi[2], 
			    int dag);

  // Preconditioned support
  void G5D_Meo(Fermion_t chi, 
	       Fermion_t psi,
	       int result_cb,
	       int dag);
  void G5D_MeoDeriv(Fermion_t chi, 
		    Fermion_t psi,
		    Matrix_t force,
		    int result_cb,
		    int dag);

  void G5D_Mooee(Fermion_t chi, 
		 Fermion_t psi, 
		 int dag) ;

  void G5D_MooeeInv(Fermion_t chi, 
		    Fermion_t psi, 
		    int dag) ;

  // Assembly Optimised version
  void G5D_CayleyMooeeFiveD(Fermion_t psi,
		       Fermion_t chi,
		       int dag);

  void G5D_CayleyMooeeInvFiveD(Fermion_t psi,
		       Fermion_t chi,
		       int dag);

  void G5D_CayleyMeoFiveD(Fermion_t psi,
			  Fermion_t chi,
			  int dag);


  void DperpDWFcompat_cpp (Fermion_t psi,Fermion_t chi,int cb,int dag);


  //////////////////////////////////////////////////
  // Reproducibility checking and performance monitoring
  //////////////////////////////////////////////////
  int iter;
  int mprecFlopsPerSite ( void ) ;
  int axpyBytes ( void ) ;
  int axpyNormFlops ( void ) ;
  int axpyFlops ( void ) ;
  int mprecFlops( void ) ;
  int cbSites( void ) ;
#ifdef BFM_REPRO
  int repro_entry;
  int csum_entry;
  void repro_control(int state) ;
  void repro(double d,int who) ;
  void repro(uint64_t sum) ;
  std::vector< std::vector<double> > repro_log;
  std::vector<uint64_t> csum_log;
#else
  void repro_control(int state) {};
  void repro(double d, int who) {};
  void repro(uint64_t) {};
#endif


  int repro_state;
  int CheckStopWatchFile(void);
  uint64_t  checksum(Fermion_t psi);
  uint64_t  checksumGauge(void);
  uint64_t  saved_checksum;
  void dump(Fermion_t psi, char *filename, double scale,int cb);

};

// Assembler kernels
#include <bfm_vmx.h>

// System dependent allocator abstraction
void *bfm_alloc(size_t size,int mem_type = mem_slow);
void  bfm_free (void *ptr,int length=0);

#ifdef BFM_QMP
#include <bfmcommqmp.h>
#ifndef BFM_BGQ_SPI_DSLASH
typedef bfmcommQMP<double> bfm;
typedef bfmcommQMP<double> bfm_dp;
typedef bfmcommQMP<float>  bfm_sp;
template <class Float>
class bfm_internal : public bfmcommQMP<Float> { 
};
#else
#include <bfmcommspi.h> // Override to use SPI comms inside dslash
typedef bfmcommspi<double> bfm;
typedef bfmcommspi<double> bfm_dp;
typedef bfmcommspi<float>  bfm_sp;
template <class Float>
class bfm_internal : public bfmcommspi<Float> { 
};
#endif
#endif

#ifdef BFM_FAKE
#include <bfmcommfake.h>
typedef bfmcommfake<double> bfm;
typedef bfmcommfake<double> bfm_dp;
typedef bfmcommfake<float>  bfm_sp;
template <class Float>
class bfm_internal : public bfmcommfake<Float> { 
};
#endif

#ifdef BFM_IROIRO
#include <bfmcommiroiro.h>
typedef bfmcommIroIro<double> bfm;
typedef bfmcommIroIro<double> bfm_dp;
typedef bfmcommIroIro<float>  bfm_sp;
template <class Float>
class bfm_internal : public bfmcommIroIro<Float> { 
};
#endif

#endif
