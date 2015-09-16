/**
 * \file
 * \brief Prefetching support.
 * \author Bernhard Mendl
 */


#ifndef PREFETCH_H_INCLUDED
#define PREFETCH_H_INCLUDED 1


#include <limits.h>


enum PREFETCH_LEVEL
{
	L1,   // L2 -> L1
	L2,   // mem -> L2
	L2L1  // combined prefetch: mem -> L2 -> L1
};


/// Prefetch into L1.
void prefetch0_data(int bytes, int basereg, int offset);


/// Prefetch into L2.
void prefetch1_data(int bytes, int basereg, int offset);


/**
 * \brief Prefetches \p bytes of data from <tt>offset(basereg)</tt> into cache \p level.
 *
 * NOTE: It is assumed that the address in \p basereg is aligned to a cache line boundary!
 *
 * \param[in] level    The cache level to be prefetched into.
 * \param[in] bytes    Number of Bytes to be prefetched beginning at <tt>offset(basereg)</tt>.
 * \param[in] basereg  An index of Iregs containing a pointer.
 * \param[in] offset   Offset relative to \p basereg. Note: This is not an offset _handle_, but the real offset in Bytes!
 * \param[in] offset_closer If \p level is L2L1, 2 prefetches are issued. Then this optional parameter must be set to an offset for L1.
 */
void prefetch_data(enum PREFETCH_LEVEL level, int bytes, int basereg, int offset, int offset_closer = 0);


/**
 * See prefetch_data(), but prefetches have non-temporal hint set.
 */
void prefetch_data_nt(enum PREFETCH_LEVEL level, int bytes, int basereg, int offset, int offset_closer = 0);
void prefetch0_data_nt(int bytes, int basereg, int offset);
void prefetch1_data_nt(int bytes, int basereg, int offset);


#endif  /* PREFETCH_H_INCLUDED */
