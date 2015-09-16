AM_CXXFLAGS = -I@top_srcdir@/bfm  `qdp++-config --cxxflags`
AM_LDFLAGS = -L@top_builddir@/bfm `qdp++-config --ldflags`

noinst_PROGRAMS = majority_vote

majority_vote_SOURCES = majority_vote.C
LIBS = -lbfm `qdp++-config --libs`  

include ../../Make.bgq
