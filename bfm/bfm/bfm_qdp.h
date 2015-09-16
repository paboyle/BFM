#ifndef BFM_QDP_H
#define BFM_QDP_H

#include <bfm.h>


#ifdef BFM_QDP_BINDING

/**************************************************************
 * Internal interfaces
 **************************************************************
 */
typedef LatticeFermion T4;
typedef multi1d<LatticeFermion> T5;
typedef multi1d<LatticeColorMatrix> G;

#include <bfm_qdp_generic.h>
#include <bfm_cg_mixed_prec.h> // Uses bfm_qdp datatype. Probably shouldnt
/*********************************************
 * Action interfaces
 *********************************************
 */


#include <bfm_qdp_dwf.h>
#include <bfm_qdp_wilson.h>
#include <bfm_qdp_wilson_tm.h>
//#include <bfm_qdp_dwf_gminus2.h>
#include <bfmActionParams.h>
#include <bfm_qdp_g5d.h>
#include <bfm_qdp_cayley_prec.h>                  
///////////// Wish list below /////////////////
//include <bfm_qdp_wilson_clover.h>


#endif
#endif
