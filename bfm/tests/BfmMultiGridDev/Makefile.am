AM_CXXFLAGS =-I@top_srcdir@/bfm ${CHROMA_CXXFLAGS} `chroma-config --cxxflags` 

noinst_PROGRAMS = TestBfmMultiGrid

TestBfmMultiGrid_SOURCES = TestBfmMultiGrid.C
TestBfmMultiGrid_OBJS=
TestBfmMultiGrid_LDFLAGS=-L@top_builddir@/bfm `chroma-config --ldflags` 
TestBfmMultiGrid_LDADD=@CHROMA_LIBS@ -lbfm -lesslbg -lxlf90_r -lxlfmath -lxl -L/opt/ibmmath/lib64 -L/opt/ibmcmp/xlf/bg/14.1/lib64  -lbfm

include ../../Make.bgq
