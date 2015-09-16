

AM_CXXFLAGS = -I@top_srcdir@/bfm ${CHROMA_CXXFLAGS} `chroma-config --cxxflags`
LIBS=`chroma-config --libs` -lbfm @LIBS@
noinst_PROGRAMS = call_dslash
call_dslash_SOURCES = call_dslash.C
call_dslash_LDFLAGS = -L@top_builddir@/bfm `chroma-config --ldflags`


include ../../Make.bgq
