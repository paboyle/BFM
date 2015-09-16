/**
 * \file
 * \brief Prefetching support.
 * \author Bernhard Mendl
 */


#include "prefetch.h"
#include "processor.h"
#include <stdio.h>


extern struct processor *PROC;


#define MYDEBUG 1


void prefetch0_data(int bytes, int basereg, int offset)
{
	prefetch_data(L1, bytes, basereg, offset, 0);
}


void prefetch1_data(int bytes, int basereg, int offset)
{
	prefetch_data(L2, bytes, basereg, offset, 0);
}


void prefetch_data(
	enum PREFETCH_LEVEL level,
	int bytes,
	int basereg,
	int offset,
	int offset_closer)
{
	int i, cl;
	int unaligned, dangling, nlines, aloff, aloff_closer;
	int offhandle, offhandle_closer;

	// don't insert prefetches in C output
	if (PROC->id == UNKNOWN || PROC->id == UNKNOWN_SINGLE) {
		return;
	}

	cl = PROC->CacheLine;  // size of cache line

	unaligned = offset % cl;             // unaligned bytes at beginning
	dangling  = (bytes-unaligned) % cl;  // unaligned bytes at end
	nlines    = ((bytes-unaligned-dangling) / cl) + !(!dangling) + !(!unaligned);
	aloff = offset - unaligned;  // align offset to CL bdry

#ifdef MYDEBUG
	fprintf(stderr, "Debug: %s: off=%d, bytes=%d, pre=%d, post=%d, aloff=%d, nlines=%d\n",
		__FUNCTION__,
		offset, bytes, unaligned, dangling, aloff, nlines);
#endif

	unaligned = offset_closer % cl;             // unaligned bytes at beginning
	dangling  = (bytes-unaligned) % cl;  // unaligned bytes at end
	nlines    = ((bytes-unaligned-dangling) / cl) + !(!dangling) + !(!unaligned);
	aloff_closer = offset_closer - unaligned;  // align offset to CL bdry

#ifdef MYDEBUG
	if (level == L2L1) {
		fprintf(stderr, "Debug: %s: off=%d, bytes=%d, pre=%d, post=%d, aloff=%d, nlines=%d\n",
			__FUNCTION__,
			offset_closer, bytes, unaligned, dangling, aloff_closer, nlines);
	}
#endif

	for (i = 0; i < nlines; i++) {
		offhandle        = get_offset_handle(aloff + (cl*i), Byte);
		offhandle_closer = get_offset_handle(aloff_closer + (cl*i), Byte);

		switch (level)
		{
		case L1:
			make_inst(CACHPIPE, PREF0, offhandle, basereg);
			break;
		case L2:
			make_inst(CACHPIPE, PREF1, offhandle, basereg);
			break;
		case L2L1:
			make_inst(CACHPIPE, PREF1, offhandle, basereg);
			make_inst(CACHPIPE, PREF0, offhandle_closer, basereg);
			break;
		default:
			fprintf(stderr, "Error: Invalid prefetch level.\n");
			abort();
		}
	}
}



void prefetch_data_nt(
	enum PREFETCH_LEVEL level,
	int bytes,
	int basereg,
	int offset,
	int offset_closer)
{
	int i, cl;
	int unaligned, dangling, nlines, aloff, aloff_closer;
	int offhandle, offhandle_closer;

	// don't insert prefetches in C output
	if (PROC->id == UNKNOWN || PROC->id == UNKNOWN_SINGLE) {
		return;
	}

	cl = PROC->CacheLine;  // size of cache line

	unaligned = offset % cl;             // unaligned bytes at beginning
	dangling  = (bytes-unaligned) % cl;  // unaligned bytes at end
	nlines    = ((bytes-unaligned-dangling) / cl) + !(!dangling) + !(!unaligned);
	aloff = offset - unaligned;  // align offset to CL bdry

	unaligned = offset_closer % cl;             // unaligned bytes at beginning
	dangling  = (bytes-unaligned) % cl;  // unaligned bytes at end
	nlines    = ((bytes-unaligned-dangling) / cl) + !(!dangling) + !(!unaligned);
	aloff_closer = offset_closer - unaligned;  // align offset to CL bdry

	for (i = 0; i < nlines; i++) {
		offhandle        = get_offset_handle(aloff + (cl*i), Byte);
		offhandle_closer = get_offset_handle(aloff_closer + (cl*i), Byte);

		switch (level)
		{
		case L1:
			make_inst(CACHPIPE, PREF0NT, offhandle, basereg);
			break;
		case L2:
			make_inst(CACHPIPE, PREF1NT, offhandle, basereg);
			break;
		case L2L1:
			make_inst(CACHPIPE, PREF1NT, offhandle, basereg);
			make_inst(CACHPIPE, PREF0NT, offhandle_closer, basereg);
			break;
		default:
			fprintf(stderr, "Error: Invalid prefetch level.\n");
			abort();
		}
	}
}


void prefetch0_data_nt(int bytes, int basereg, int offset)
{
	prefetch_data_nt(L1, bytes, basereg, offset, 0);
}


void prefetch1_data_nt(int bytes, int basereg, int offset)
{
	prefetch_data_nt(L2, bytes, basereg, offset, 0);
}
