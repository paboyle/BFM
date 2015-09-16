
# Substitutions
buildCXX = @bagel_CXX@
buildCXXFLAGS = @bagel_CXXFLAGS@ 
buildLDFLAGS = @bagel_LDFLAGS@
buildLIBS = @bagel_LIBS@
buildLINKFLAGS = $(buildLDFLAGS) $(buildLIBS)

# Override compile and links for build system compilation
SUFFIXES: .C $(EXEEXT)

noinst_PROGRAMS = \
	write_vaxpy$(EXEEXT) \
	write_caxpy$(EXEEXT) \
	write_inner$(EXEEXT) \
	write_vaxpy_ssp$(EXEEXT) \
	write_vaxpy_norm$(EXEEXT) \
	write_vmx_vaxpy_norm$(EXEEXT) \
	write_face_scatter$(EXEEXT) \
	write_1pg5axpy$(EXEEXT) \
	qcdoc_decom$(EXEEXT) \
	qcdoc_recon$(EXEEXT) \
	bgl_recon$(EXEEXT) \
	bgl_decom$(EXEEXT) \
	vmx_dwf$(EXEEXT) \
	vmx_clov_apply$(EXEEXT) \
	vmx_cayley$(EXEEXT) \
	vmx_cayley_inv$(EXEEXT) \
	vmx_cg$(EXEEXT) \
	qpx_dwf$(EXEEXT) \
	qpx_dwf_singlepass$(EXEEXT) \
	qpx_deriv$(EXEEXT) \
	qpx_dwf_comms$(EXEEXT) \
	vmx_dperp$(EXEEXT) \
	vmx_gather$(EXEEXT) \
	vmx_gather_proj$(EXEEXT) \
	vmx_tmass$(EXEEXT) \
	vmx_merge$(EXEEXT) \
	vmx_scatter$(EXEEXT)


write_vaxpy_SOURCES = write_vaxpy.C
write_vaxpy$(EXEEXT): write_vaxpy.C

write_caxpy_SOURCES = write_caxpy.C
write_caxpy$(EXEEXT): write_caxpy.C

write_inner_SOURCES = write_inner.C
write_inner$(EXEEXT): write_inner.C

write_vaxpy_ssp_SOURCES = write_vaxpy_ssp.C
write_vaxpy_ssp$(EXEEXT): write_vaxpy_ssp.C

vmx_tmass_SOURCES = vmx_tmass.C
vmx_tmass$(EXEEXT): vmx_tmass.C

vmx_cg_SOURCES = vmx_cg.C
vmx_cg$(EXEEXT): vmx_cg.C

vmx_merge_SOURCES = vmx_merge.C
vmx_merge$(EXEEXT): vmx_merge.C

write_vaxpy_norm_SOURCES = write_vaxpy_norm.C
write_vaxpy_norm$(EXEEXT): write_vaxpy_norm.C

write_vmx_vaxpy_norm_SOURCES = write_vmx_vaxpy_norm.C
write_vmx_vaxpy_norm$(EXEEXT): write_vmx_vaxpy_norm.C

write_face_scatter_SOURCES = write_face_scatter.C
write_face_scatter$(EXEEXT): write_face_scatter.C

write_1pg5axpy_SOURCES = write_1pg5axpy.C
write_1pg5axpy$(EXEEXT): write_1pg5axpy.C

qcdoc_decom_SOURCES = qcdoc_decom.C
qcdoc_decom$(EXEEXT): qcdoc_decom.C

qcdoc_recon_SOURCES = qcdoc_recon.C
qcdoc_recon$(EXEEXT): qcdoc_recon.C

bgl_decom_SOURCES = bgl_decom.C
bgl_decom$(EXEEXT): bgl_decom.C

bgl_recon_SOURCES = bgl_recon.C
bgl_recon$(EXEEXT): bgl_recon.C

vmx_dwf_SOURCES = vmx_dwf.C
vmx_dwf$(EXEEXT): vmx_dwf.C

vmx_clov_apply_SOURCES = vmx_clov_apply.C
vmx_clov_apply$(EXEEXT): vmx_clov_apply.C

vmx_cayley_SOURCES = vmx_cayley.C
vmx_cayley$(EXEEXT): vmx_cayley.C

vmx_cayley_inv_SOURCES = vmx_cayley_inv.C
vmx_cayley_inv$(EXEEXT): vmx_cayley_inv.C

qpx_dwf_SOURCES = qpx_dwf.C
qpx_dwf$(EXEEXT): qpx_dwf.C

qpx_dwf_singlepass_SOURCES = qpx_dwf_singlepass.C
qpx_dwf_singlepass$(EXEEXT): qpx_dwf_singlepass.C

qpx_deriv_SOURCES = qpx_deriv.C
qpx_deriv$(EXEEXT): qpx_deriv.C

qpx_dwf_comms_SOURCES = qpx_dwf_comms.C
qpx_dwf_comms$(EXEEXT): qpx_dwf_comms.C

vmx_gather_SOURCES = vmx_gather.C
vmx_gather$(EXEEXT): vmx_gather.C

vmx_scatter_SOURCES = vmx_scatter.C
vmx_scatter$(EXEEXT): vmx_scatter.C

vmx_gather_proj_SOURCES = vmx_gather_proj.C
vmx_gather_proj$(EXEEXT): vmx_gather_proj.C

vmx_dperp_SOURCES = vmx_dperp.C
vmx_dperp$(EXEEXT): vmx_dperp.C

%$(EXEEXT): %.C
	$(buildCXX) $(buildCXXFLAGS) -o $@ $< $(buildLDFLAGS) $(buildLIBS)
	chmod u+x $@

