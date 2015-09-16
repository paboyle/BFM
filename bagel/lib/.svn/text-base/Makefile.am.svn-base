INCFLAGS = -I@top_srcdir@/lib -I@top_builddir@/include

AM_CFLAGS = $(INCFLAGS)

lib_LIBRARIES = libbagel.a

nobase_include_HEADERS = \
	bagel_config.h \
	bagel_config_internal.h \
	bagel_int.h\
	alpha.h \
	alpha264.h \
	mic.h \
	powerpc.h \
	prefetch.h \
	processor.h \
	registers.h \
	unknown.h \
	usparc.h

libbagel_a_SOURCES = \
	alpha064.C \
	alpha064_single.C \
	alpha164.C \
	alpha164_single.C \
	alpha264.C \
	alpha264_single.C \
	arguments.C \
	arralloc.C \
	bgl.C \
	bgq.C \
	dump.C \
	mic.C \
	offsets.C \
	powerIII.C \
	powerIII_single.C \
	powerpc.C \
	powerpc_single.C \
	prefetch.C \
	processor.C \
	queue_complex.C \
	registers.C \
	schedule.C \
	setarch.C \
	stack.C \
	stream.C \
	unknown.C \
	usparcII.C \
	usparcII_single.C \
	externs.C 
