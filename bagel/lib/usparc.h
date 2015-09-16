
/*
 *
 *  Copyright UKQCD Collaboration, October 2000.
 *  Written by Peter Boyle.
 *  This software is provided for NON-COMMERCIAL use only,
 *  It is provided under the GNU pubic License V2
 *  It is provided as is and is not guaranteed fit for any purpose.
 *
 */

#ifndef _ULTRASPARC_
#define _ULTRASPARC_

/*Ultra-II double precision*/
void create_usparcII_instructions( struct processor *proc);
struct processor *create_usparcII (void);
int start_loop_usparcII(int counter);
void stop_loop_usparcII(int branchno,int counter);
void check_iterations_usparcII(int cntreg,int retno);
void directive_print_usparcII(struct instruction *i);
void queue_cmovge_usparcII(int cond,int src,int dst);
void cache_print_usparcII(struct instruction *i);

/*PPC 440 single precision*/
void create_usparcIIs_instructions( struct processor *proc);
struct processor *create_usparcIIs(void);
int start_loop_usparcIIs(int counter);
void stop_loop_usparcIIs(int branchno,int counter);
void check_iterations_usparcIIs(int cntreg,int retno);
void directive_print_usparcIIs(struct instruction *i);
void queue_cmovge_usparcIIs(int cond,int src,int dst);
void cache_print_usparcIIs(struct instruction *i);


#endif
