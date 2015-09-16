

AM_CXXFLAGS = -I@top_srcdir@/bfm ${CHROMA_CXXFLAGS} `chroma-config --cxxflags`

noinst_PROGRAMS = call_dslash

call_dslash_SOURCES = call_dslash.C
call_dslash_LDFLAGS =   -L@top_builddir@/bfm `chroma-config --ldflags` 
call_dslash_LDADD=-lbfm  -lbfm @CHROMA_LIBS@ -lpthread
include ../../Make.bgq
