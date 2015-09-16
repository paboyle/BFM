
AM_CXXFLAGS = -I@top_srcdir@/bfm 
AM_LDFLAGS = -L@top_builddir@/bfm 

if BUILD_QDP
AM_CXX= `chroma-config --cxx` 
AM_CXXFLAGS+= `chroma-config --cxxflags` 
AM_LDFLAGS+= `chroma-config --ldflags`  `chroma-config --libs`
endif

noinst_PROGRAMS = conjugate_gradient

conjugate_gradient_SOURCES = conjugate_gradient.C

LIBS = -lbfm -lm 

include ../../Make.bgq
