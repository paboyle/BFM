
AM_CXXFLAGS = -I@top_srcdir@/bfm ${CHROMA_CXXFLAGS} `chroma-config --cxxflags` 

noinst_PROGRAMS = call_dslash4 testLittleDiracOp

call_dslash4_SOURCES = BfmMultiGrid.C call_dslash.C
call_dslash4_OBJS=
call_dslash4_LDFLAGS =   blas_complex_invert.o qpx_zgemv.o  -L@top_builddir@/bfm `chroma-config --ldflags`  -L/opt/ibmmath/lib64
call_dslash4_LDADD=-lbfm  -lbfm @CHROMA_LIBS@  -lessl
testLittleDiracOp_SOURCES = BfmMultiGrid.C testLittleDiracOp.C
testLittleDiracOp_LDFLAGS =blas_complex_invert.o qpx_zgemv.o   -L@top_builddir@/bfm `chroma-config --ldflags` 
testLittleDiracOp_LDADD=-lbfm  -lbfm @CHROMA_LIBS@

include ../../Make.bgq
