
AM_CXXFLAGS = -I@top_srcdir@/bfm ${CHROMA_CXXFLAGS} `chroma-config --cxxflags`

noinst_PROGRAMS = call_dslash testLittleDiracOp

call_dslash_SOURCES = littleDiracOp.C call_dslash.C 
call_dslash_LDFLAGS =   -L@top_builddir@/bfm `chroma-config --ldflags` 
call_dslash_LDADD=-lbfm  -lbfm @CHROMA_LIBS@ 
testLittleDiracOp_SOURCES = littleDiracOp.C testLittleDiracOp.C 
testLittleDiracOp_LDFLAGS =   -L@top_builddir@/bfm `chroma-config --ldflags` 
testLittleDiracOp_LDADD=-lbfm  -lbfm @CHROMA_LIBS@

include ../../Make.bgq
