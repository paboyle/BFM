
LIST= bfmbase.C \
	bfmcommfake.C \
	bfmdop.C \
	bfmdperp.C \
	bfm_memalign.C \
	bfm_thread_spi.C \
	zolotarev.c \
	conjugate_gradient.C

PPCF=/bgsys/bgq/drivers/ppcfloor/

MYLIBS=-lmpich -lopa  -lmpl -lpami  -lSPI_cnk -lrt -lpthread
MYLDF=-L${PPCF}/comm/gcc/lib/  -L${PPCF}/comm/sys/lib -L${PPCF}/spi/lib
MYINC=-I. -I${PPCF}/comm/gcc/include/ -fpermissive

release: 
	cp @top_srcdir@/bfm/*.C .
	cp @top_srcdir@/bfm/*.c .
	cp @top_srcdir@/bfm/*.h .
	- cp @top_builddir@/bfm/*.C .
	- cp @top_builddir@/bfm/*.h .
	- cp @top_builddir@/bfm/*.S .
	cp @top_srcdir@/tests/ibm_qcd_qpx_diag/*.C .
	$(CXX) $(MYINC) $(CXXFLAGS) $(LDFLAGS)  $(LIST) *.S -o ibm_qcd_qpx_diag $(MYLDF) $(MYLIBS) -fopenmp
