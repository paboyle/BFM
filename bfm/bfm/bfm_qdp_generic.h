#ifndef _BFM_QDP_GENERIC_H_
#define _BFM_QDP_GENERIC_H_


enum inv_type_t {
    CG_PREC_M,
    CG_PREC_MDAG,
    CG_PREC_MDAGM, // done on odd parity
    CG_PREC_MDAGM_MULTI,
    CG_UNPREC_M,
    CG_UNPREC_MDAG,
    CG_UNPREC_MDAGM,
    M_UNPREC,
    MDAG_UNPREC
  } ;


template<class Float>
class bfm_qdp : public bfm_internal<Float> {
 public:
  // Used to control what the threading routine does
  inv_type_t inv_type;
  // Used for thread argument glue
  Fermion_t qdp_chi_h[2];
  Fermion_t qdp_psi_h[2];
  Fermion_t mtmp;
  // Eigenvectors
  multi1d<bfm_fermion> qdp_evecs;
  multi1d<double> qdp_evals;


  // Following used by multishift cg thread glue
  int nshift;
  Fermion_t *qdp_chi_multi_h;
  double    *shifts;
  double    *alpha;
  double    *mresidual;
  int single;
  // return the number of iters
  int iter;
};


template<class Float> void * bfm_thread_cg_routine(void * dwf_p);
template<class Float> void   bfm_spawn_cg (bfm_qdp<Float> &dwf);


#ifdef BFM_CHROMA_BINDING
#include <chroma.h>
#define BfmFermToProp FermToProp
#define BfmPropToFerm PropToFerm
#else

  void BfmFermToProp(const LatticeFermionD& a, LatticePropagatorD& b, 
		     int color_index, int spin_index);

  void BfmFermToProp(const LatticeFermionF& a, LatticePropagatorF& b, 
		  int color_index, int spin_index);

  void BfmPropToFerm(const LatticePropagatorF& b, LatticeFermionF& a, 
		     int color_index, int spin_index);

  void BfmPropToFerm(const LatticePropagatorD& b, LatticeFermionD& a, 
		     int color_index, int spin_index);
#endif

#ifdef CHROMA_3
static Set timeslice;
#else
static UnorderedSet timeslice;
#endif
static int timeslice_is_made = 0;

class TSliceFunc : public SetFunc
{
public:
  TSliceFunc(int dir): dir_decay(dir) {}

  int operator() (const multi1d<int>& coordinate) const {return coordinate[dir_decay];}
  int numSubsets() const {return Layout::lattSize()[dir_decay];}

  int dir_decay;

private:
  TSliceFunc() {}  // hide default constructor
};
Set &GetTimeslice(void);

/*******************************************************
 * What calls the routine below? It is DWF specific. Likely deprecated
 *******************************************************/
#ifdef DONT_INCLUDE_IT
template <class Float>
int bfmbase<Float>::CG4d(LatticePropagator &sol, 
			  LatticePropagator &src,
			  multi1d<LatticeColorMatrix> &U,
			  std::string output_stem)
{

  this->importGauge(U);

  Fermion_t chi_h[2];
  Fermion_t psi_h[2];
  for(int cb=0;cb<2;cb++){
    chi_h[cb] = this->allocFermion();
    psi_h[cb] = this->allocFermion();
  }

  for(int color_source = 0; color_source < Nc; ++color_source) {
    for(int spin_source = 0; spin_source < Ns; ++spin_source) {
	
      LatticeFermion psi = zero;
      LatticeFermion chi;

      // Extract a fermion source
      BfmPropToFerm(src, chi, color_source, spin_source);

      /*
       * Normalize the source in case it is really huge or small -
       * a trick to avoid overflows or underflows
       */
      Real nrm = sqrt(norm2(chi));
      if (toFloat(nrm) == 0.0) {
	psi = zero;
      } else {
	Real fact = 1.0/nrm;
	chi *= fact;

	for(int cb=0;cb<2;cb++){
	  this->importFermion(chi,chi_h[cb],cb);
	  this->importFermion(psi,psi_h[cb],cb);
	}
	int iter = this->CGNE(psi_h, chi_h);
	
	for(int cb=0;cb<2;cb++){
	  this->exportFermion(psi,psi_h[cb],cb);
	}	
        // Unnormalize the source following the inverse of the normalization above
        fact = Real(1) / fact;
        psi *= fact;

        /*
         * Move the solution to the appropriate components
         * of quark propagator.
         */
        BfmFermToProp(psi, sol, color_source, spin_source);
      }

    } /* end loop over spin_source */
  } /* end loop over color_source */

  for(int cb=0;cb<2;cb++){
    this->freeFermion(psi_h[cb]);
    this->freeFermion(chi_h[cb]);
  }
  return 0;
}
#endif













template<class Float>
void bfm_spawn_cg (bfm_qdp<Float> &dwf)
{
#pragma omp parallel 
    {
#pragma omp for 
      for(int i=0;i<dwf.nthread;i++) {
	bfm_thread_cg_routine<Float>((void *)&dwf);
      }
    }
}

template<class Float>
void * bfm_thread_cg_routine(void * arg)
{
  bfm_qdp<Float> *dwf = (bfm_qdp<Float> *)arg;
  int me = dwf->thread_barrier();

  if ( dwf->inv_type == CG_PREC_MDAGM ) {
    dwf->iter = dwf->CGNE_prec_MdagM(dwf->qdp_chi_h[1],dwf->qdp_psi_h[1]);
  } else if ( dwf->inv_type == CG_PREC_M ) {
    dwf->iter = dwf->CGNE_M(dwf->qdp_chi_h,dwf->qdp_psi_h,dwf->qdp_evecs,dwf->qdp_evals);
  } else if ( dwf->inv_type == CG_PREC_MDAG ) {
    dwf->iter = dwf->CGNE_Mdag(dwf->qdp_chi_h,dwf->qdp_psi_h);
  } else if ( dwf->inv_type == CG_PREC_MDAGM_MULTI ) {
    dwf->iter = dwf->CGNE_prec_MdagM_multi_shift(dwf->qdp_chi_multi_h,
						 dwf->qdp_psi_h[1],// odd parity
						 dwf->shifts,
						 dwf->alpha,
						 dwf->nshift,
						 dwf->mresidual,
						 dwf->single); // sum single result
  } else if ( dwf->inv_type == M_UNPREC ) {
    dwf->Munprec(dwf->qdp_psi_h,dwf->qdp_chi_h,dwf->mtmp,DaggerNo);
  } else if ( dwf->inv_type == MDAG_UNPREC ) {
    dwf->Munprec(dwf->qdp_psi_h,dwf->qdp_chi_h,dwf->mtmp,DaggerYes);
  } else if ( dwf->inv_type == CG_UNPREC_M ) {
    dwf->iter = dwf->CGNE_unprec(dwf->qdp_chi_h,dwf->qdp_psi_h);
  } else if ( (dwf->inv_type == CG_UNPREC_MDAG)
	   || (dwf->inv_type == CG_UNPREC_MDAGM) ){
    QDP_error_exit("Not implemented UNPREC cg as threaded yet");
  }
  fflush(stdout);
  dwf->thread_barrier();
  return (void *)dwf;
}

/*******************************************************************
 * This is the QDP impex glue layer
 *******************************************************************
 */

#include <qdp.h>

template <class Float>
void      bfmbase<Float>::exportForce(Matrix_t handle, multi1d<LatticeColorMatrix> &f, int cb)
{
  if ( f.size() != Nd ) exit(-1);

  for ( int mu=0;mu<Nd;mu++){

    QDPdouble * f_p;
    f_p = (QDPdouble *)&f[mu].elem(0).elem();
    exportForce(handle,f_p,mu,cb);

  }
}

// This is the base Bagel 2 import. Use QDP to do the shifting, conjugation why not....
template <class Float>
void      bfmbase<Float>::importGauge(multi1d<LatticeColorMatrix> &u)
{
  if ( u.size() != Nd ) exit(-1);
  for ( int mu=0;mu<Nd;mu++){

    QDPdouble * U_p;

    int dir = 2*mu+1;
    U_p = (QDPdouble *)&u[mu].elem(0).elem();
    importGauge(U_p,dir);

    dir = 2*mu;
    LatticeColorMatrix udag;
    udag = adj(shift(u[mu],BACKWARD,mu));

    U_p = (QDPdouble *)&udag.elem(0).elem();
    importGauge(U_p,dir);

  }
}


// This is the base Bagel 2 import. Use QDP to do the shifting, conjugation why not....
template <class Float>
void      bfmbase<Float>::exportGauge(multi1d<LatticeColorMatrix> &u,multi1d<LatticeColorMatrix> &udag)
{
  if ( u.size() != Nd ) exit(-1);
  for ( int mu=0;mu<Nd;mu++){

    QDPdouble * U_p;

    int dir = 2*mu+1;
    U_p = (QDPdouble *)&u[mu].elem(0).elem();
    exportGauge(U_p,dir);

    dir = 2*mu;
    LatticeColorMatrix utmp;
    U_p = (QDPdouble *)&utmp.elem(0).elem();
    exportGauge(U_p,dir);

    udag[mu] = shift(utmp,FORWARD,mu);

  }
}



template <class Float>
void      bfmbase<Float>::importGaugeNN(multi1d<LatticeColorMatrix> &u)
{
  if ( u.size() != Nd ) exit(-1);
  for ( int mu=0;mu<Nd;mu++){

    LatticeColorMatrix utmp;
    QDPdouble * U_p;

    int dir = 2*mu+1;
    utmp = u[mu] * shift(u[mu],FORWARD,mu);
    U_p = (QDPdouble *)&utmp.elem(0).elem();
    importGauge(U_p,dir);

    dir = 2*mu;
    utmp = adj(shift(shift(utmp,BACKWARD,mu),BACKWARD,mu));
    U_p = (QDPdouble *)&utmp.elem(0).elem();
    importGauge(U_p,dir);

  }
}


template <class Float>
void      bfmbase<Float>::importPhysicalFermion(LatticeFermion &psi, Fermion_t handle, int cb)
{

  LatticeFermion tmp=zero;
  for(int s=0;s<Ls;s++){
    importFermion(tmp,handle,cb,s);
  }

  switch(solver) { 
  case DWF:
  case DWFrb4d:
  case HtCayleyTanh: // Plain old DWF.
  case HtCayleyZolo: // Chiu Improved
    tmp = chiralProjectPlus(psi);
    importFermion(tmp,handle,cb,0);
    tmp = chiralProjectMinus(psi);
    importFermion(tmp,handle,cb,Ls-1);
    break;

  case HwCayleyZolo: // Chiu Optimal
  case HwCayleyTanh: // Scaled shamir
  case HmCayleyTanh: // Scaled shamir

    tmp = chiralProjectPlus(psi);
    importFermion(tmp,handle,cb,0);
    tmp = chiralProjectMinus(psi);
    importFermion(tmp,handle,cb,Ls-1);
    break;
  case WilsonFermion:
    // Clover, // Unimplemented
  case WilsonTM:  
  case WilsonNN: // Next Nearest neighbour Wilson term (Eguchi-Kawamoto support)
  case DWFTransfer:
  case DWFTransferInv:
    importFermion(psi,handle,cb,0);
    break;
  case HwPartFracZolo: // KEK's approach
  case HwContFracZolo: // Edwards, Kennedy et al prefer this
  case HwPartFracTanh: // Never been looked at!
  case HwContFracTanh: // Never been looked at!
      importFermion(psi,handle,cb,Ls-1);
      break;
  default:
    exit(-1);
    break;
  }
}


template <class Float>
void      bfmbase<Float>::exportPhysicalFermion(LatticeFermion &psi, Fermion_t handle, int cb)
{
  LatticeFermion tmp=zero;

  switch(solver) { 
  case DWF:
  case DWFrb4d:
  case HwCayleyZolo: // Chiu Optimal
  case HtCayleyZolo: // Chiu Improved
  case HwCayleyTanh: // Scaled shamir
  case HmCayleyTanh: // Scaled shamir
  case HtCayleyTanh: // Plain old DWF.
    exportFermion(tmp,handle,cb,0);
    psi = chiralProjectMinus(tmp);
    exportFermion(tmp,handle,cb,Ls-1);
    psi += chiralProjectPlus(tmp);
    break;
  case WilsonFermion:
    // Clover, // Unimplemented
  case WilsonTM:  
  case WilsonNN: // Next Nearest neighbour Wilson term (Eguchi-Kawamoto support)
  case DWFTransfer:
  case DWFTransferInv:
    exportFermion(psi,handle,cb,0);
  break;
  case HwPartFracZolo: // KEK's approach
  case HwContFracZolo: // Edwards, Kennedy et al prefer this
  case HwPartFracTanh: // Never been looked at!
  case HwContFracTanh: // Never been looked at!
    exportFermion(psi,handle,cb,Ls-1);
  }
}


template <class Float>
void      bfmbase<Float>::importFermion(LatticeFermion & psi, Fermion_t handle, int cb,int s)
{
// Need to worry about base parity 
  QDPdouble *psi_p = (QDPdouble *) &psi.elem(0).elem(0).elem(0).real();
  impexFermion(psi_p,handle,1,cb,s);
}
template <class Float>
void      bfmbase<Float>::exportFermion  (LatticeFermion &psi, Fermion_t handle, int cb,int s)
{
// Need to worry about base parity
  QDPdouble *psi_p = (QDPdouble *) &psi.elem(0).elem(0).elem(0).real();
  impexFermion(psi_p,handle,0,cb,s);
}

template <class Float>
void      bfmbase<Float>::importFermion(multi1d<LatticeFermion> &psi, Fermion_t handle, int cb)
{
// Need to worry about base parity 
  if ( psi.size() != Ls ) exit(-1);
  for(int s=0;s<Ls;s++){
    QDPdouble *psi_p = (QDPdouble *) &(psi[s].elem(0).elem(0).elem(0).real());
    impexFermion(psi_p,handle,1,cb,s);
  }
}
template <class Float>
void      bfmbase<Float>::exportFermion  (multi1d<LatticeFermion> &psi, Fermion_t handle, int cb)
{
// Need to worry about base parity
  if ( psi.size() != Ls ) exit(-1);
  LatticeFermion tmp;
  for(int s=0;s<Ls;s++){
    QDPdouble *psi_p = (QDPdouble *) &(psi[s].elem(0).elem(0).elem(0).real());
    impexFermion(psi_p,handle,0,cb,s);
  }
}


#endif
