#ifndef BFM_VMX_H
#define BFM_VMX_H

extern "C" {

  //Linalg
 void vmx_vaxpy       (double*out,double *A,double *InScale, double *Add, integer len);
 void vmx_vaxpby      (double *out,double *A,double *B,double *InScale, double *Add, integer len);
 void vmx_caxpy       (double *out,double *A,double *InScale, double *Add, integer len);
 void vmx_caxpby      (double *out,double *A,double *B,double *InScale, double *Add, integer len);
 void vmx_inner  (double *InScale, double *Add, integer len, double *norm);
 void vmx_vaxpy_norm  (double *out,double *A,double *InScale, double *Add, integer len, double *norm);
 void vmx_vaxpby_norm (double *out,double *A,double *B,double *InScale, double *Add, integer len, double *norm);
 void vmx_1pG5vaxpy   (double *addout,double *A,double *InScaleProj, integer len );
 void vmx_tmass       (double *out,double *coeff_ig5,double *coeff_ident,double *in,integer len);
 void vmx_merge       (double *out,double *inlo, double *inhi,integer len);
 void vmx_vaxpby_ssp      (integer *argp);

 void vmx_vaxpy_s       (float *out,double *A,float *InScale, float *Add, integer len);
 void vmx_vaxpby_s      (float *out,double *A,double *B,float *InScale, float *Add, integer len);
 void vmx_caxpy_s       (float *out,double *A,float *InScale, float *Add, integer len);
 void vmx_caxpby_s      (float *out,double *A,double *B,float *InScale, float *Add, integer len);
 void vmx_inner_s  (float *InScale, float *Add, integer len, double *norm);
 void vmx_vaxpy_norm_s  (float *out,double *A,float *InScale, float *Add, integer len, double *norm);
 void vmx_vaxpby_norm_s (float *out,double *A,double *B,float *InScale, float *Add, integer len, double *norm);
 void vmx_1pG5vaxpy_s   (float *addout,double *A,float *InScaleProj, integer len );
 void vmx_tmass_s       (float  *out,double *coeff_ig5,double *coeff_ident,float  *in,integer len);
 void vmx_merge_s       (float *out,float *inlo, float *inhi,integer len);
 void vmx_vaxpby_ssp_s      (integer *argp);


  void vmx_cg       (integer *args);
  void vmx_cg_s     (integer *args);

  //////////////////////////////////////////
  // Dirac ops
  //////////////////////////////////////////

  //DWF incl D5
 void vmx_dwf_dp_int(integer arg_p);       //interior links only
 void vmx_dwf_dp_int_dag(integer arg_p);
 void vmx_dwf_dp_int_s(integer arg_p);
 void vmx_dwf_dp_int_dag_s(integer arg_p);

 void vmx_dwf_dp_a_int(integer arg_p);       //interior links only + a*psi
 void vmx_dwf_dp_a_int_dag(integer arg_p);
 void vmx_dwf_dp_a_int_s(integer arg_p);
 void vmx_dwf_dp_a_int_dag_s(integer arg_p);

  //Wilson(xLs) (excludes D5)
 void vmx_wil_int(integer arg_p); 
 void vmx_wil_int_dag(integer arg_p);
 void vmx_wil_int_s(integer arg_p);
 void vmx_wil_int_dag_s(integer arg_p);
 void vmx_wil_a_int(integer arg_p); 
 void vmx_wil_a_int_dag(integer arg_p);
 void vmx_wil_a_int_s(integer arg_p);
 void vmx_wil_a_int_dag_s(integer arg_p);

 // Pass 2 pieces of wilson hopping term
 void vmx_dwf_ext(integer arg_p);       //exterior links only & adds to
 void vmx_dwf_ext_dag(integer arg_p);
 void vmx_dwf_ext_s(integer arg_p);
 void vmx_dwf_ext_dag_s(integer arg_p);

 void vmx_dwf_ext_scale(integer arg_p);       //exterior links only & adds to & scales result
 void vmx_dwf_ext_scale_dag(integer arg_p);
 void vmx_dwf_ext_scale_s(integer arg_p);
 void vmx_dwf_ext_scale_dag_s(integer arg_p);

  /*
   * No comms overlap versions
   */
 void vmx_dwf_dp(integer arg_p);       //interior links only
 void vmx_dwf_dp_dag(integer arg_p);
 void vmx_dwf_dp_s(integer arg_p);
 void vmx_dwf_dp_dag_s(integer arg_p);

 void vmx_dwf_dp_as(integer arg_p);       //interior links only + a*psi
 void vmx_dwf_dp_as_dag(integer arg_p);
 void vmx_dwf_dp_as_s(integer arg_p);
 void vmx_dwf_dp_as_dag_s(integer arg_p);

 void vmx_wil(integer arg_p);       //interior links only
 void vmx_wil_dag(integer arg_p);
 void vmx_wil_s(integer arg_p);
 void vmx_wil_dag_s(integer arg_p);

 void vmx_wil_hs(integer arg_p);
 void vmx_wil_dag_hs(integer arg_p);
 void vmx_wil_as_hs(integer arg_p);
 void vmx_wil_as_dag_hs(integer arg_p);
 void vmx_merge_hs       (float *out,float *inlo, float *inhi,integer len);

 void vmx_wil_as(integer arg_p);       //interior links only + a*psi
 void vmx_wil_as_dag(integer arg_p);
 void vmx_wil_as_s(integer arg_p);
 void vmx_wil_as_dag_s(integer arg_p);

  // Force terms
 void vmx_deriv(integer arg_p); 
 void vmx_deriv_dag(integer arg_p); 
 void vmx_deriv_s(integer arg_p); 
 void vmx_deriv_dag_s(integer arg_p); 

  /*
   * Clover multiply
   */
  void vmx_clov_apply(integer arg_p);
  void vmx_clov_apply_s(integer arg_p);

  /*
   * Cayley matrix support
   */
  void vmx_cayley(integer argp);
  void vmx_cayley_s(integer argp);
  void vmx_cayley_inv(integer argp);
  void vmx_cayley_inv_s(integer argp);
  void vmx_cayley_inv_dag(integer argp);
  void vmx_cayley_inv_dag_s(integer argp);

 /**************************************************************
  * comms buffer suport
  **************************************************************
  */
 void vmx_gather_proj(integer arg_p);
 void vmx_gather_proj_s(integer arg_p);
 void vmx_gather_proj_hs(integer arg_p);
 void vmx_gather_proj_perm(integer arg_p);
 void vmx_gather_proj_perm_s(integer arg_p);
 void vmx_gather_proj_perm_hs(integer arg_p);
 void vmx_gather_proj_extr_h(integer arg_p);
 void vmx_gather_proj_extr_l(integer arg_p);
 void vmx_gather_proj_extr_h_s(integer arg_p);
 void vmx_gather_proj_extr_l_s(integer arg_p);
 void vmx_gather_proj_extr_h_hs(integer arg_p);
 void vmx_gather_proj_extr_l_hs(integer arg_p);

}
#endif 
