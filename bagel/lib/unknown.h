
/*
 *
 *  Copyright UKQCD Collaboration, October 2000.
 *  Written by Peter Boyle.
 *  This software is provided for NON-COMMERCIAL use only,
 *  It is provided under the GNU pubic License V2
 *  It is provided as is and is not guaranteed fit for any purpose.
 *
 */

#ifndef _UNKNOWN_C_
#define _UNKNOWN_C_


void set_unknown_type_float(void );
void set_unknown_type_double(void);
void set_unknown_bgq(void);
void set_unknown_vmx(void);
void create_unknown_instructions( struct processor *proc);
struct processor *create_unknown (void);
int start_loop_unknown(int counter);
void stop_loop_unknown(int branchno,int counter);
void check_iterations_unknown(int cntreg,int retno);
void directive_print_unknown(struct instruction *i);
void queue_cmovge_unknown(int cond,int src,int dst);
void cache_print_unknown(struct instruction *i);

#endif
