
AM_CXXFLAGS = -I@top_srcdir@/bfm 
AM_LDFLAGS = -L@top_builddir@/bfm 

if BUILD_QDP
AM_CXX= `chroma-config --cxx` 
AM_CXXFLAGS+= `chroma-config --cxxflags` 
AM_LDFLAGS+= `chroma-config --ldflags`  `chroma-config --libs`
endif

noinst_PROGRAMS = call_dslash  

call_dslash_SOURCES = call_dslash.C cb0dag0.C cb0dag1.C cb1dag0.C cb1dag1.C

LIBS = -lbfm -lm 

include ../../Make.bgq
