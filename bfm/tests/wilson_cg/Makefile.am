

AM_CXXFLAGS = -I@top_srcdir@/bfm ${CHROMA_CXXFLAGS} `chroma-config --cxxflags`

noinst_PROGRAMS = conjugate_gradient

conjugate_gradient_SOURCES = conjugate_gradient.C
conjugate_gradient_LDFLAGS =   -L@top_builddir@/bfm `chroma-config --ldflags` 
conjugate_gradient_LDADD=-lbfm  -lbfm @CHROMA_LIBS@
include ../../Make.bgq
