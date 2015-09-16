#ifndef _BFM_CG_MIXED_PREC_H_
#define _BFM_CG_MIXED_PREC_H_

int bfm_mixed_precision_CG(  bfm_qdp<double> &bfm_double,
			     bfm_qdp<float>  &bfm_single,
			     double residual,int max_iter,
			     Fermion_t sol_d[2],
			     Fermion_t src_d[2]
			     );

#endif
