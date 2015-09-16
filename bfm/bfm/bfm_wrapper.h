#ifndef __bfm_wrapper_h__
#define __bfm_wrapper_h__

#include <string>
// Bfm Headers
#include <bfm.h>
#include <bfm_qdp.h>
#include <bfm_hdcg_wrapper.h>


  enum BfmInvType { 
    BfmInv_CG,
    BfmInv_MCR,
    BfmInv_MCRG5,
    BfmInv_HDCG
  };
  enum BfmMatType {
    BfmMat_M,
    BfmMat_MdagM_oo,
    BfmMat_CayleyDminusUnprec,
    BfmMat_Mdag
  };
  struct BfmWrapperParams {
  public:
    BfmInvType BfmInverter; 
    BfmMatType BfmMatrix;   // M, MpcdagMpc
    BfmPrecisionType BfmPrecision;
    int  MaxIter;
    multi1d<Real> RsdTarget;
    Real Delta;
    bfmActionParams BAP;
  };

  /////////////////////////////////////////////////////////////
  // Wrapper class for the BFM inverter
  /////////////////////////////////////////////////////////////
  class BfmWrapper {

  public:
    BfmWrapperParams invParam;
      // These are the links
      // They may be smeared and the BC's may be applied
    typedef LatticeFermion T4;
    typedef LatticeColorMatrix U;
    typedef multi1d<LatticeColorMatrix> Q;


    //begin karthee clover
    multi1d<RScalar<REAL> >  CloverDiag;
    multi1d<RComplex<REAL> > CloverOffDiag;
    multi1d<RScalar<REAL> >  CloverInvDiag;
    multi1d<RComplex<REAL> > CloverInvOffDiag;
    //end karthee clover

    Q links;


    // Does the inversion
    int bfmInvert4d_mixed(T4 &sol,const T4 &chi);
    int bfmInvert5d_mixed(multi1d<T4> &sol,const multi1d<T4> &chi);

    template<class Float> void SpectralRange(bfm_qdp<Float>  &bfm);

    template<class Float> 
    int bfmInvert4d(T4 &sol,const T4 &chi);

    template<class Float> int bfmDminusInvert(multi1d<T4> &sol,multi1d<T4> &src);

    template<class Float> 
    int bfmInvert5d(multi1d<T4> &sol,const multi1d<T4> &chi);

    template<class Float>     int bfmCayleyDminus(multi1d<T4> &vec);

    template<class Float> 
    int bfmMultiInvert5d(multi1d< multi1d<T4> >& psi, 
			 const multi1d<Real>& shifts, 
			 const multi1d<Real>& Residuals, 
			 const multi1d<T4>& chi);

    template<class Float> 
    int bfmInvert_mixed(bfm_qdp<Float> &bfm_dp,
			Fermion_t sol_d[2],
			Fermion_t src_d[2]);
    template<class Float> 
    int bfmInvert_mixed(bfm_qdp<Float>  &bfm_dp,
			bfm_qdp<float> &bfm_sp,
			Fermion_t sol_d[2],
			Fermion_t src_d[2]);
    template<class Float> 
    void bfmTwoFlavorForce5d(multi1d<T4> &phi,multi1d<LatticeColorMatrix> &force);

    template<class Float> 
      void bfmTwoFlavorRatioForce5d(multi1d<T4> &phi,multi1d<LatticeColorMatrix> &force,
				    Real M_pv, Real M_f);

    template<class Float> 
      void bfmMixedInvert_prec(bfm_qdp<Float>&M, Fermion_t sol,Fermion_t src);
    
    template<class Float> 
      void bfmTwoFlavorRatioForce(bfm_qdp<Float>&M, bfm_qdp<Float>&PV, Fermion_t phi,Matrix_t force[2]);

    void RationalSanityCheck(double power, int n_p, double a0, double ak[], double bk);

    template<class Float> 
    void bfmOneFlavorRatioRationalForce(multi1d<T4> &phi,
					multi1d<LatticeColorMatrix> &force,
					Real M_pv, Real M_f,
					Real nrm_f , multi1d<Real> &shifts_pv, multi1d<Real> &residue_pv,
					Real nrm_pv, multi1d<Real> &shifts_f, multi1d<Real> &residue_f
					);
    template<class Float> 
      void bfmOneFlavorRatioRationalForce(bfm_qdp<Float>&M, bfm_qdp<Float>&PV,
					  Fermion_t phi, 
					  Matrix_t force[2],
					  int n_pv, double a0_pv, double ak_pv[],double bk_pv[],double rk_pv[],
					  int n_f , double a0_f , double ak_f[] ,double bk_f[],double rk_f[]
					  );

    void RationalSanityCheck(double power, int n_p, double a0, double ak[], double bk[]);

    
    template<class Float> 
    int bfmMultiInvert4d(multi1d<T4> & psi, 
			 const multi1d<Real>& shifts, 
			 const multi1d<Real>& Residuals, 
			 const T4& chi);

    BfmWrapper(const BfmWrapperParams& invParam_) : invParam(invParam_) , links(Nd) {
      QDPIO::cout << "BfmWrapper constructor resid" << invParam.RsdTarget[0] <<endl;
    };
  protected:
    BfmWrapper() : links(Nd) {};

  };

#endif 

