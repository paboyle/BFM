
nobase_include_HEADERS = \
	bagel_int.h \
	bfm_options.h \
	bfm.h\
	bfm_linalg.h \
	bfm_mprec.h \
	bfm_vmx.h \
	bfmActionParams.h\
	bfm_thread.h \
	bfmcommfake.h \
	bfmcommqmp.h \
	bfmcommiroiro.h \
	bfmcommspi.h \
	bfm_cg.h \
	bfm_cg_mixed_prec.h \
	bfm_cg_unprec.h \
	bfm_cg_multi.h \
	bfm_import_export.h \
	bfm_iroiro.h	\
	zolotarev.h\
	bfm_qdp.h \
	bfm_qdp_dwf.h \
	bfm_qdp_wilson.h \
	bfm_qdp_wilson_tm.h \
	bfm_qdp_generic.h \
	bfm_qdp_g5d.h \
	bfm_qdp_cayley_prec.h \
	bfm_qdp_dwf_gminus2.h \
	bfm_qdp_chroma_linop.h \
	bfm_wrapper.h \
	bfm_hdcg_wrapper.h \
	BfmHDCG.h \
	nested_4d.h

lib_LIBRARIES = libbfm.a 

libbfm_a_SOURCES = \
	blas_complex_invert.c\
	zolotarev.c \
	bfmbase.C \
	bfmdop.C\
	bfm_block.C\
	bfmdperp.C \
	bfm_memalign.C \
	BfmHDCG.C 

gen_C_sources = \
	vmx_gather_proj.C \
	vmx_gather_proj_perm.C \
	vmx_gather_proj_extr_h.C \
	vmx_gather_proj_extr_l.C \
	vmx_gather_proj_s.C \
	vmx_gather_proj_perm_s.C \
	vmx_gather_proj_extr_h_s.C \
	vmx_gather_proj_extr_l_s.C \
	vmx_cg.C \
	vmx_cg_s.C \
	vmx_vaxpy.C \
	vmx_caxpy.C \
	vmx_caxpby.C \
	vmx_inner.C \
	vmx_vaxpy_norm.C \
	vmx_vaxpby.C  \
	vmx_vaxpby_ssp.C  \
	vmx_vaxpby_norm.C \
	vmx_vaxpy_s.C \
	vmx_caxpy_s.C \
	vmx_caxpby_s.C \
	vmx_inner_s.C \
	vmx_vaxpy_norm_s.C \
	vmx_vaxpby_ssp_s.C  \
	vmx_vaxpby_s.C  \
	vmx_vaxpby_norm_s.C \
	vmx_dwf_dp_int.C \
	vmx_dwf_dp_int_dag.C \
	vmx_dwf_dp_int_s.C \
	vmx_dwf_dp_int_dag_s.C \
	vmx_dwf_dp_a_int.C \
	vmx_dwf_dp_a_int_dag.C \
	vmx_dwf_dp_a_int_s.C \
	vmx_dwf_dp_a_int_dag_s.C \
	vmx_wil_int.C \
	vmx_wil_int_dag.C \
	vmx_wil_int_s.C \
	vmx_wil_int_dag_s.C \
	vmx_wil_a_int.C \
	vmx_wil_a_int_dag.C \
	vmx_wil_a_int_s.C \
	vmx_wil_a_int_dag_s.C \
	vmx_deriv.C \
	vmx_deriv_dag.C \
	vmx_deriv_s.C \
	vmx_deriv_dag_s.C \
	vmx_dwf_ext.C  \
	vmx_dwf_ext_dag.C \
	vmx_dwf_ext_s.C \
	vmx_dwf_ext_dag_s.C \
	vmx_dwf_ext_scale.C \
	vmx_dwf_ext_scale_dag.C \
	vmx_dwf_ext_scale_s.C \
	vmx_dwf_ext_scale_dag_s.C \
	vmx_dwf_dp.C \
	vmx_dwf_dp_dag.C \
	vmx_dwf_dp_s.C \
	vmx_dwf_dp_dag_s.C \
	vmx_dwf_dp_as.C \
	vmx_dwf_dp_as_dag.C \
	vmx_dwf_dp_as_s.C \
	vmx_dwf_dp_as_dag_s.C \
	vmx_cayley.C \
	vmx_cayley_s.C \
	vmx_cayley_inv.C \
	vmx_cayley_inv_s.C \
	vmx_cayley_inv_dag.C \
	vmx_cayley_inv_dag_s.C \
	vmx_wil.C \
	vmx_wil_dag.C \
	vmx_wil_s.C \
	vmx_wil_dag_s.C \
	vmx_wil_as.C \
	vmx_wil_as_dag.C \
	vmx_wil_as_s.C \
	vmx_wil_as_dag_s.C \
	vmx_tmass.C \
	vmx_tmass_s.C \
	vmx_merge.C \
	vmx_merge_s.C \
	vmx_1pg5vaxpy.C \
	vmx_1pg5vaxpy_s.C \
	vmx_clov_apply.C\
	vmx_clov_apply_s.C


halfprec_C_sources = \
	vmx_gather_proj_hs.C \
	vmx_gather_proj_perm_hs.C \
	vmx_gather_proj_extr_h_hs.C \
	vmx_gather_proj_extr_l_hs.C \
	vmx_merge_hs.C \
	vmx_wil_as_hs.C \
	vmx_wil_as_dag_hs.C \
	vmx_wil_dag_hs.C \
	vmx_wil_hs.C 

gen_asm_sources=$(gen_C_sources:.C=_asm.S)
halfprec_asm_sources=$(halfprec_C_sources:.C=_asm.S)

extra_sources=
###############################
# Thread model
###############################
if BUILD_SPI
  extra_sources+=bfm_thread_spi.C
endif

libbfm_a_LIBADD=
xlcobj=
if BUILD_XLCQPX
#  libbfm_a_LIBADD+=qpx_zgemv.o blas_complex_invert.o
#  xlcobj+=qpx_zgemv.o blas_complex_invert.o
  libbfm_a_LIBADD+=qpx_zgemv.o 
  xlcobj+=qpx_zgemv.o 
endif

##########################################
# Comms
##########################################
if BUILD_COMM_IROIRO
 libbfm_a_SOURCES+=bfmcommiroiro.C
endif
if BUILD_COMM_QMP
 libbfm_a_SOURCES+=bfmcommqmp.C
endif
if BUILD_COMM_SPI
 libbfm_a_SOURCES+=bfmcommqmp.C
 libbfm_a_SOURCES+=bfmcommspi.C
endif
if BUILD_COMM_FAKE
 libbfm_a_SOURCES+=bfmcommfake.C
endif

######################################
# Code generation targets
######################################
if BUILD_BFM_NOARCH
 extra_sources+=$(gen_C_sources)
endif

if BUILD_BFM_NOARCHBGQ
 extra_sources+=$(gen_C_sources)
endif

if BUILD_BFM_BGQ
 extra_sources+=$(gen_asm_sources)
 extra_sources+=$(halfprec_asm_sources)
endif

if BUILD_BFM_KNC
 extra_sources+=$(gen_asm_sources)
endif

if BUILD_BFM_BGL
  extra_sources+=$(gen_asm_sources)
endif

#####################################
# CHROMA/QDP++ compile options if needed
#####################################
if BUILD_QDP

if BUILD_CHROMA_CONFIG
AM_CXX= `chroma-config --cxx`
AM_CXXFLAGS= -g `chroma-config --cxxflags` 
 extra_sources+=bfm_qdp.C
 extra_sources+=bfm_wrapper.C 
 extra_sources+=bfm_cg_mixed_prec.C 
else
AM_CXX= `qdp++-config --cxx`
AM_CXXFLAGS= -g `qdp++-config --cxxflags` 
 extra_sources+=bfm_qdp.C
 extra_sources+=bfm_wrapper.C 
 extra_sources+=bfm_cg_mixed_prec.C 
endif

endif

#########################################

# Mark the extra sources for cleaning
CLEANFILES = execution.map *~ \
	$(gen_C_sources) \
	$(gen_asm_sources) 


# Add the extra sources to the library, but not to the dist
nodist_libbfm_a_SOURCES = $(extra_sources)

# manual dependencies 
libbfm_a_DEPENDENCIES = $(extra_sources) $(xlcobj)

####################################################
# Bagel 2 assembler files
####################################################
#
# How to call BAGEL code gen routines
#
vmx_bin=$(bagelWriterDir)/qpx_dwf$(EXEEXT)
#vmx_bin=$(bagelWriterDir)/vmx_dwf$(EXEEXT)
lockcache=-l
bagelTargetName = @bagel_MUNGED_TARGET@
bagelWriterDir  = @top_builddir@/code_gen
singleFlag=-s


#
# XLC version
#
qpx_zgemv.o: @srcdir@/qpx_zgemv.c
	mpixlc -O3 -c @srcdir@/qpx_zgemv.c

#blas_complex_invert.o: @srcdir@/blas_complex_invert.c
#	mpixlc -O3 -c @srcdir@/blas_complex_invert.c -I/opt/ibmmath/essl/5.1/include/

####################################################
# Wilson first pass routines (or smart 5d overlap)
####################################################
vmx_deriv.C vmx_deriv_asm.S: $(bagelWriterDir)/qpx_deriv
	$(bagelWriterDir)/qpx_deriv -nvmx_deriv -P$(bagelTargetName) > $@

vmx_deriv_dag.C vmx_deriv_dag_asm.S: $(bagelWriterDir)/qpx_deriv
	$(bagelWriterDir)/qpx_deriv -d -nvmx_deriv_dag  -P$(bagelTargetName)> $@

vmx_deriv_s.C vmx_deriv_s_asm.S: $(bagelWriterDir)/qpx_deriv
	$(bagelWriterDir)/qpx_deriv -nvmx_deriv_s -P$(bagelTargetName) > $@

vmx_deriv_dag_s.C vmx_deriv_dag_s_asm.S: $(bagelWriterDir)/qpx_deriv
	$(bagelWriterDir)/qpx_deriv -d -nvmx_deriv_dag_s -P$(bagelTargetName) > $@

vmx_wil_int.C vmx_wil_int_asm.S: $(vmx_bin)
	$(vmx_bin) $(lockcache) -nvmx_wil_int -i -P$(bagelTargetName) > $@

vmx_wil_int_dag.C vmx_wil_int_dag_asm.S: $(vmx_bin)
	$(vmx_bin) $(lockcache) -nvmx_wil_int_dag -i -d -P$(bagelTargetName) > $@

vmx_wil_int_s.C vmx_wil_int_s_asm.S: $(vmx_bin)
	$(vmx_bin) $(lockcache) -nvmx_wil_int_s -i -P$(bagelTargetName)  $(singleFlag) > $@

vmx_wil_int_dag_s.C vmx_wil_int_dag_s_asm.S: $(vmx_bin)
	$(vmx_bin) $(lockcache) -nvmx_wil_int_dag_s -d -i -P$(bagelTargetName)  $(singleFlag) > $@

vmx_wil_a_int.C vmx_wil_a_int_asm.S: $(vmx_bin)
	$(vmx_bin) $(lockcache) -nvmx_wil_a_int -i -a -P$(bagelTargetName) > $@

vmx_wil_a_int_dag.C vmx_wil_a_int_dag_asm.S: $(vmx_bin)
	$(vmx_bin) $(lockcache) -nvmx_wil_a_int_dag -i -a -d -P$(bagelTargetName) > $@

vmx_wil_a_int_s.C vmx_wil_a_int_s_asm.S: $(vmx_bin)
	$(vmx_bin) $(lockcache) -nvmx_wil_a_int_s -i -a -P$(bagelTargetName)  $(singleFlag) > $@

vmx_wil_a_int_dag_s.C vmx_wil_a_int_dag_s_asm.S: $(vmx_bin)
	$(vmx_bin) $(lockcache) -nvmx_wil_a_int_dag_s -a -d -i -P$(bagelTargetName)  $(singleFlag) > $@

##################################################
# half precision comms buffer, single internal variants
# useful for polynomial preconditioning.
# Avoid overlapping comms/compute and
# get faster single pass dslash implementation
##################################################

vmx_wil_as_hs.C vmx_wil_as_hs_asm.S: $(vmx_bin)
	$(vmx_bin) $(lockcache) -nvmx_wil_as_hs -k -a -o -h -P$(bagelTargetName)  $(singleFlag) > $@

vmx_wil_as_dag_hs.C vmx_wil_as_dag_hs_asm.S: $(vmx_bin)
	$(vmx_bin) $(lockcache) -nvmx_wil_as_dag_hs -h -k -d -a -o  -P$(bagelTargetName)  $(singleFlag) > $@

vmx_wil_hs.C vmx_wil_hs_asm.S: $(vmx_bin)
	$(vmx_bin) $(lockcache) -nvmx_wil_hs -h -P$(bagelTargetName)  $(singleFlag) > $@

vmx_wil_dag_hs.C vmx_wil_dag_hs_asm.S: $(vmx_bin)
	$(vmx_bin) $(lockcache) -nvmx_wil_dag_hs -h -d  -P$(bagelTargetName)  $(singleFlag) > $@

#################################
# Cayley form dperp routines
#################################
vmx_cayley.C vmx_cayley_asm.S: $(bagelWriterDir)/vmx_cayley
	$(bagelWriterDir)/vmx_cayley -nvmx_cayley -P$(bagelTargetName) > $@

vmx_cayley_s.C vmx_cayley_s_asm.S: $(bagelWriterDir)/vmx_cayley
	$(bagelWriterDir)/vmx_cayley -nvmx_cayley_s -P$(bagelTargetName)  $(singleFlag) > $@

vmx_cayley_inv.C vmx_cayley_inv_asm.S: $(bagelWriterDir)/vmx_cayley_inv
	$(bagelWriterDir)/vmx_cayley_inv -nvmx_cayley_inv -P$(bagelTargetName) > $@

vmx_cayley_inv_s.C vmx_cayley_inv_s_asm.S: $(bagelWriterDir)/vmx_cayley_inv
	$(bagelWriterDir)/vmx_cayley_inv -nvmx_cayley_inv_s -P$(bagelTargetName)  $(singleFlag)> $@

vmx_cayley_inv_dag.C vmx_cayley_inv_dag_asm.S: $(bagelWriterDir)/vmx_cayley_inv
	$(bagelWriterDir)/vmx_cayley_inv -nvmx_cayley_inv_dag -d -P$(bagelTargetName) > $@

vmx_cayley_inv_dag_s.C vmx_cayley_inv_dag_s_asm.S: $(bagelWriterDir)/vmx_cayley_inv
	$(bagelWriterDir)/vmx_cayley_inv -nvmx_cayley_inv_dag_s -d -P$(bagelTargetName)  $(singleFlag)> $@


####################################################
# DWF first pass routines
####################################################

# result = int dslash to the result (including dperp)
# Is it possible to drop the no addto versions?
vmx_dwf_dp_int.C vmx_dwf_dp_int_asm.S: $(vmx_bin)
	$(vmx_bin) $(lockcache) -nvmx_dwf_dp_int -p -i -P$(bagelTargetName) > $@

vmx_dwf_dp_int_dag.C vmx_dwf_dp_int_dag_asm.S: $(vmx_bin)
	$(vmx_bin) $(lockcache) -nvmx_dwf_dp_int_dag -d -i -p -P$(bagelTargetName) > $@

vmx_dwf_dp_int_s.C vmx_dwf_dp_int_s_asm.S: $(vmx_bin)
	$(vmx_bin) $(lockcache) -nvmx_dwf_dp_int_s -p -i -P$(bagelTargetName)  $(singleFlag) > $@

vmx_dwf_dp_int_dag_s.C vmx_dwf_dp_int_dag_s_asm.S: $(vmx_bin)
	$(vmx_bin) $(lockcache) -nvmx_dwf_dp_int_dag_s -p -i -d -P$(bagelTargetName)  $(singleFlag) > $@


# Scale preexisting contents of result and then add int dlsash to the result (including dperp)
vmx_dwf_dp_a_int.C vmx_dwf_dp_a_int_asm.S: $(vmx_bin)
	$(vmx_bin) $(lockcache) -nvmx_dwf_dp_a_int -a -p -i -P$(bagelTargetName) > $@

vmx_dwf_dp_a_int_dag.C vmx_dwf_dp_a_int_dag_asm.S: $(vmx_bin)
	$(vmx_bin) $(lockcache) -nvmx_dwf_dp_a_int_dag -a -d -i -p -P$(bagelTargetName) > $@

vmx_dwf_dp_a_int_s.C vmx_dwf_dp_a_int_s_asm.S: $(vmx_bin)
	$(vmx_bin) $(lockcache) -nvmx_dwf_dp_a_int_s -a -p -i -P$(bagelTargetName)  $(singleFlag) > $@

vmx_dwf_dp_a_int_dag_s.C vmx_dwf_dp_a_int_dag_s_asm.S: $(vmx_bin)
	$(vmx_bin) $(lockcache) -nvmx_dwf_dp_a_int_dag_s -a -p -i -d -P$(bagelTargetName)  $(singleFlag) > $@


####################################################################################
# Following are deprecated I believe
################################################################
# result = dslash to the result (including dperp)
vmx_dwf_dp_as.C vmx_dwf_dp_as_asm.S: $(vmx_bin)
	$(vmx_bin) $(lockcache) -nvmx_dwf_dp_as -p -a -o  -P$(bagelTargetName) > $@

vmx_dwf_dp_as_dag.C vmx_dwf_dp_as_dag_asm.S: $(vmx_bin)
	$(vmx_bin) $(lockcache) -nvmx_dwf_dp_as_dag -d -p -a -o  -P$(bagelTargetName) > $@

vmx_dwf_dp_as_s.C vmx_dwf_dp_as_s_asm.S: $(vmx_bin)
	$(vmx_bin) $(lockcache) -nvmx_dwf_dp_as_s -p -a -o  -P$(bagelTargetName)  $(singleFlag) > $@

vmx_dwf_dp_as_dag_s.C vmx_dwf_dp_as_dag_s_asm.S: $(vmx_bin)
	$(vmx_bin) $(lockcache) -nvmx_dwf_dp_as_dag_s -p -d -a -o  -P$(bagelTargetName)  $(singleFlag) > $@

vmx_dwf_dp.C vmx_dwf_dp_asm.S: $(vmx_bin)
	$(vmx_bin) $(lockcache) -nvmx_dwf_dp -p  -P$(bagelTargetName) > $@

vmx_dwf_dp_dag.C vmx_dwf_dp_dag_asm.S: $(vmx_bin)
	$(vmx_bin) $(lockcache) -nvmx_dwf_dp_dag -d -p  -P$(bagelTargetName) > $@

vmx_dwf_dp_s.C vmx_dwf_dp_s_asm.S: $(vmx_bin)
	$(vmx_bin) $(lockcache) -nvmx_dwf_dp_s -p  -P$(bagelTargetName)  $(singleFlag) > $@

vmx_dwf_dp_dag_s.C vmx_dwf_dp_dag_s_asm.S: $(vmx_bin)
	$(vmx_bin) $(lockcache) -nvmx_dwf_dp_dag_s -p -d  -P$(bagelTargetName)  $(singleFlag) > $@



vmx_wil_as.C vmx_wil_as_asm.S: $(vmx_bin)
	$(vmx_bin) $(lockcache) -nvmx_wil_as -k -a -o  -P$(bagelTargetName) > $@

vmx_wil_as_dag.C vmx_wil_as_dag_asm.S: $(vmx_bin)
	$(vmx_bin) $(lockcache) -nvmx_wil_as_dag -k -d -a -o  -P$(bagelTargetName) > $@

vmx_wil_as_s.C vmx_wil_as_s_asm.S: $(vmx_bin)
	$(vmx_bin) $(lockcache) -nvmx_wil_as_s -k -a -o  -P$(bagelTargetName)  $(singleFlag) > $@

vmx_wil_as_dag_s.C vmx_wil_as_dag_s_asm.S: $(vmx_bin)
	$(vmx_bin) $(lockcache) -nvmx_wil_as_dag_s -k -d -a -o  -P$(bagelTargetName)  $(singleFlag) > $@

vmx_wil.C vmx_wil_asm.S: $(vmx_bin)
	$(vmx_bin) $(lockcache) -nvmx_wil  -P$(bagelTargetName) > $@

vmx_wil_dag.C vmx_wil_dag_asm.S: $(vmx_bin)
	$(vmx_bin) $(lockcache) -nvmx_wil_dag -d  -P$(bagelTargetName) > $@

vmx_wil_s.C vmx_wil_s_asm.S: $(vmx_bin)
	$(vmx_bin) $(lockcache) -nvmx_wil_s  -P$(bagelTargetName)  $(singleFlag) > $@

vmx_wil_dag_s.C vmx_wil_dag_s_asm.S: $(vmx_bin)
	$(vmx_bin) $(lockcache) -nvmx_wil_dag_s -d  -P$(bagelTargetName)  $(singleFlag) > $@



########################################################################################################
# Pass 2 for DWF and Wilson: the (halo) externally connected fields in D_W in a second pass after comms
########################################################################################################
vmx_dwf_ext.C vmx_dwf_ext_asm.S: $(vmx_bin)
	$(vmx_bin) $(lockcache) -nvmx_dwf_ext -e -a -P$(bagelTargetName) > $@

vmx_dwf_ext_dag.C vmx_dwf_ext_dag_asm.S: $(vmx_bin)
	$(vmx_bin) $(lockcache) -nvmx_dwf_ext_dag -d -e -a -P$(bagelTargetName) > $@

vmx_dwf_ext_s.C vmx_dwf_ext_s_asm.S: $(vmx_bin)
	$(vmx_bin) $(lockcache) -nvmx_dwf_ext_s -e -a -P$(bagelTargetName)  $(singleFlag) > $@

vmx_dwf_ext_dag_s.C vmx_dwf_ext_dag_s_asm.S: $(vmx_bin)
	$(vmx_bin) $(lockcache) -nvmx_dwf_ext_dag_s -e -a -d -P$(bagelTargetName)  $(singleFlag) > $@
#
# with scale. Would like to split interior exterior here. This forces full reload of
# output field to scale it, even for interior points
vmx_dwf_ext_scale.C vmx_dwf_ext_scale_asm.S: $(vmx_bin)
	$(vmx_bin) $(lockcache) -nvmx_dwf_ext_scale -e -o -a -P$(bagelTargetName) > $@

vmx_dwf_ext_scale_dag.C vmx_dwf_ext_scale_dag_asm.S: $(vmx_bin)
	$(vmx_bin) $(lockcache) -nvmx_dwf_ext_scale_dag -d -o -e -a -P$(bagelTargetName) > $@

vmx_dwf_ext_scale_s.C vmx_dwf_ext_scale_s_asm.S: $(vmx_bin)
	$(vmx_bin) $(lockcache) -nvmx_dwf_ext_scale_s -o -e -a -P$(bagelTargetName)  $(singleFlag) > $@

vmx_dwf_ext_scale_dag_s.C vmx_dwf_ext_scale_dag_s_asm.S: $(vmx_bin)
	$(vmx_bin) $(lockcache) -nvmx_dwf_ext_scale_dag_s -o -e -a -d -P$(bagelTargetName)  $(singleFlag) > $@


###########################################
# Support routines
###########################################


vmx_gather_proj.C vmx_gather_proj_asm.S: $(bagelWriterDir)/vmx_gather_proj$(EXEEXT)
	$(bagelWriterDir)/vmx_gather_proj$(EXEEXT) -nvmx_gather_proj -P$(bagelTargetName) > $@

vmx_gather_proj_perm.C vmx_gather_proj_perm_asm.S: $(bagelWriterDir)/vmx_gather_proj$(EXEEXT)
	$(bagelWriterDir)/vmx_gather_proj$(EXEEXT) -nvmx_gather_proj_perm -P$(bagelTargetName) -p > $@

vmx_gather_proj_extr_h.C vmx_gather_proj_extr_h_asm.S: $(bagelWriterDir)/vmx_gather_proj$(EXEEXT)
	$(bagelWriterDir)/vmx_gather_proj$(EXEEXT) -nvmx_gather_proj_extr_h -P$(bagelTargetName) -e -h > $@

vmx_gather_proj_extr_l.C vmx_gather_proj_extr_l_asm.S: $(bagelWriterDir)/vmx_gather_proj$(EXEEXT)
	$(bagelWriterDir)/vmx_gather_proj$(EXEEXT) -nvmx_gather_proj_extr_l -P$(bagelTargetName) -e > $@

vmx_cg.C vmx_cg_asm.S: $(bagelWriterDir)/vmx_cg$(EXEEXT)
	$(bagelWriterDir)/vmx_cg$(EXEEXT) -nvmx_cg -P$(bagelTargetName) > $@

vmx_cg_s.C vmx_cg_s_asm.S: $(bagelWriterDir)/vmx_cg$(EXEEXT)
	$(bagelWriterDir)/vmx_cg$(EXEEXT) -nvmx_cg_s -P$(bagelTargetName) $(singleFlag)> $@


vmx_vaxpy.C vmx_vaxpy_asm.S: $(bagelWriterDir)/write_vaxpy$(EXEEXT)
	$(bagelWriterDir)/write_vaxpy$(EXEEXT) -nvmx_vaxpy -P$(bagelTargetName) > $@

vmx_vaxpby.C vmx_vaxpby_asm.S: $(bagelWriterDir)/write_vaxpy$(EXEEXT)
	$(bagelWriterDir)/write_vaxpy$(EXEEXT) -b -nvmx_vaxpby -P$(bagelTargetName) > $@

vmx_vaxpby_ssp.C vmx_vaxpby_ssp_asm.S: $(bagelWriterDir)/write_vaxpy_ssp$(EXEEXT)
	$(bagelWriterDir)/write_vaxpy_ssp$(EXEEXT) -nvmx_vaxpby_ssp -P$(bagelTargetName) > $@

vmx_vaxpy_norm.C vmx_vaxpy_norm_asm.S: $(bagelWriterDir)/write_vmx_vaxpy_norm$(EXEEXT)
	$(bagelWriterDir)/write_vmx_vaxpy_norm$(EXEEXT) -nvmx_vaxpy_norm -P$(bagelTargetName) > $@

vmx_vaxpby_norm.C vmx_vaxpby_norm_asm.S: $(bagelWriterDir)/write_vmx_vaxpy_norm$(EXEEXT)
	$(bagelWriterDir)/write_vmx_vaxpy_norm$(EXEEXT) -b -nvmx_vaxpby_norm -P$(bagelTargetName) > $@

vmx_caxpy.C vmx_caxpy_asm.S: $(bagelWriterDir)/write_caxpy$(EXEEXT)
	$(bagelWriterDir)/write_caxpy$(EXEEXT) -nvmx_caxpy -P$(bagelTargetName) > $@

vmx_caxpby.C vmx_caxpby_asm.S: $(bagelWriterDir)/write_caxpy$(EXEEXT)
	$(bagelWriterDir)/write_caxpy$(EXEEXT) -b -nvmx_caxpby -P$(bagelTargetName) > $@

vmx_inner.C vmx_inner_asm.S: $(bagelWriterDir)/write_inner$(EXEEXT)
	$(bagelWriterDir)/write_inner$(EXEEXT) -nvmx_inner -P$(bagelTargetName) > $@

#
# S precision versions
#
vmx_gather_proj_s.C vmx_gather_proj_s_asm.S: $(bagelWriterDir)/vmx_gather_proj$(EXEEXT)
	$(bagelWriterDir)/vmx_gather_proj$(EXEEXT) -nvmx_gather_proj_s -P$(bagelTargetName) $(singleFlag)  > $@

vmx_gather_proj_perm_s.C vmx_gather_proj_perm_s_asm.S: $(bagelWriterDir)/vmx_gather_proj$(EXEEXT)
	$(bagelWriterDir)/vmx_gather_proj$(EXEEXT) -nvmx_gather_proj_perm_s -P$(bagelTargetName) $(singleFlag)  -p> $@

vmx_gather_proj_extr_h_s.C vmx_gather_proj_extr_h_s_asm.S: $(bagelWriterDir)/vmx_gather_proj$(EXEEXT)
	$(bagelWriterDir)/vmx_gather_proj$(EXEEXT) -nvmx_gather_proj_extr_h_s -P$(bagelTargetName) $(singleFlag)  -e -h> $@

vmx_gather_proj_extr_l_s.C vmx_gather_proj_extr_l_s_asm.S: $(bagelWriterDir)/vmx_gather_proj$(EXEEXT)
	$(bagelWriterDir)/vmx_gather_proj$(EXEEXT) -nvmx_gather_proj_extr_l_s -P$(bagelTargetName) $(singleFlag)  -e > $@

#
# half precision comms versions
#
vmx_gather_proj_hs.C vmx_gather_proj_hs_asm.S: $(bagelWriterDir)/vmx_gather_proj$(EXEEXT)
	$(bagelWriterDir)/vmx_gather_proj$(EXEEXT) -c -nvmx_gather_proj_hs -P$(bagelTargetName) $(singleFlag)  > $@

vmx_gather_proj_perm_hs.C vmx_gather_proj_perm_hs_asm.S: $(bagelWriterDir)/vmx_gather_proj$(EXEEXT)
	$(bagelWriterDir)/vmx_gather_proj$(EXEEXT) -c -nvmx_gather_proj_perm_hs -P$(bagelTargetName) $(singleFlag)  -p> $@

vmx_gather_proj_extr_h_hs.C vmx_gather_proj_extr_h_hs_asm.S: $(bagelWriterDir)/vmx_gather_proj$(EXEEXT)
	$(bagelWriterDir)/vmx_gather_proj$(EXEEXT) -c -nvmx_gather_proj_extr_h_hs -P$(bagelTargetName) $(singleFlag)  -e -h> $@

vmx_gather_proj_extr_l_hs.C vmx_gather_proj_extr_l_hs_asm.S: $(bagelWriterDir)/vmx_gather_proj$(EXEEXT)
	$(bagelWriterDir)/vmx_gather_proj$(EXEEXT) -c -nvmx_gather_proj_extr_l_hs -P$(bagelTargetName) $(singleFlag)  -e > $@

vmx_caxpy_s.C vmx_caxpy_s_asm.S: $(bagelWriterDir)/write_caxpy$(EXEEXT)
	$(bagelWriterDir)/write_caxpy$(EXEEXT) -nvmx_caxpy_s -P$(bagelTargetName)  $(singleFlag) > $@

vmx_caxpby_s.C vmx_caxpby_s_asm.S: $(bagelWriterDir)/write_caxpy$(EXEEXT)
	$(bagelWriterDir)/write_caxpy$(EXEEXT) -b -nvmx_caxpby_s -P$(bagelTargetName)  $(singleFlag) > $@

vmx_inner_s.C vmx_inner_s_asm.S: $(bagelWriterDir)/write_inner$(EXEEXT)
	$(bagelWriterDir)/write_inner$(EXEEXT) -nvmx_inner_s -P$(bagelTargetName)  $(singleFlag) > $@

vmx_vaxpy_s.C vmx_vaxpy_s_asm.S: $(bagelWriterDir)/write_vaxpy$(EXEEXT)
	$(bagelWriterDir)/write_vaxpy$(EXEEXT) -nvmx_vaxpy_s -P$(bagelTargetName)  $(singleFlag) > $@

vmx_vaxpby_s.C vmx_vaxpby_s_asm.S: $(bagelWriterDir)/write_vaxpy$(EXEEXT)
	$(bagelWriterDir)/write_vaxpy$(EXEEXT) -b -nvmx_vaxpby_s -P$(bagelTargetName)  $(singleFlag) > $@

vmx_vaxpby_ssp_s.C vmx_vaxpby_ssp_s_asm.S: $(bagelWriterDir)/write_vaxpy_ssp$(EXEEXT)
	$(bagelWriterDir)/write_vaxpy_ssp$(EXEEXT) -nvmx_vaxpby_ssp_s -P$(bagelTargetName)  $(singleFlag) > $@

vmx_vaxpy_norm_s.C vmx_vaxpy_norm_s_asm.S: $(bagelWriterDir)/write_vmx_vaxpy_norm$(EXEEXT)
	$(bagelWriterDir)/write_vmx_vaxpy_norm$(EXEEXT) -nvmx_vaxpy_norm_s -P$(bagelTargetName)  $(singleFlag) > $@

vmx_vaxpby_norm_s.C vmx_vaxpby_norm_s_asm.S: $(bagelWriterDir)/write_vmx_vaxpy_norm$(EXEEXT)
	$(bagelWriterDir)/write_vmx_vaxpy_norm$(EXEEXT) -b -nvmx_vaxpby_norm_s -P$(bagelTargetName)  $(singleFlag) > $@

vmx_1pg5vaxpy.C vmx_1pg5vaxpy_asm.S: $(bagelWriterDir)/write_1pg5axpy$(EXEEXT)
	$(bagelWriterDir)/write_1pg5axpy$(EXEEXT) -nvmx_1pG5vaxpy -P$(bagelTargetName)  > $@

vmx_1pg5vaxpy_s.C vmx_1pg5vaxpy_s_asm.S: $(bagelWriterDir)/write_1pg5axpy$(EXEEXT)
	$(bagelWriterDir)/write_1pg5axpy$(EXEEXT) -nvmx_1pG5vaxpy_s -P$(bagelTargetName) $(singleFlag) > $@

vmx_tmass.C vmx_tmass_asm.S: $(bagelWriterDir)/vmx_tmass$(EXEEXT)
	$(bagelWriterDir)/vmx_tmass$(EXEEXT) -nvmx_tmass -P$(bagelTargetName)  > $@

vmx_tmass_s.C vmx_tmass_s_asm.S: $(bagelWriterDir)/vmx_tmass$(EXEEXT)
	$(bagelWriterDir)/vmx_tmass$(EXEEXT) -nvmx_tmass_s -P$(bagelTargetName) $(singleFlag) > $@

vmx_merge.C vmx_merge_asm.S: $(bagelWriterDir)/vmx_merge$(EXEEXT)
	$(bagelWriterDir)/vmx_merge$(EXEEXT) -nvmx_merge -P$(bagelTargetName)  > $@

vmx_merge_s.C vmx_merge_s_asm.S: $(bagelWriterDir)/vmx_merge$(EXEEXT)
	$(bagelWriterDir)/vmx_merge$(EXEEXT) -nvmx_merge_s -P$(bagelTargetName) $(singleFlag) > $@

vmx_merge_hs.C vmx_merge_hs_asm.S: $(bagelWriterDir)/vmx_merge$(EXEEXT)
	$(bagelWriterDir)/vmx_merge$(EXEEXT) -h -nvmx_merge_hs -P$(bagelTargetName) $(singleFlag) > $@


vmx_clov_apply.C vmx_clov_apply_asm.S: $(bagelWriterDir)/vmx_clov_apply$(EXEEXT)
	$(bagelWriterDir)/vmx_clov_apply$(EXEEXT) -nvmx_clov_apply -P$(bagelTargetName)  > $@

vmx_clov_apply_s.C vmx_clov_apply_s_asm.S: $(bagelWriterDir)/vmx_clov_apply$(EXEEXT)
	$(bagelWriterDir)/vmx_clov_apply$(EXEEXT) -nvmx_clov_apply_s -P$(bagelTargetName) $(singleFlag) > $@




