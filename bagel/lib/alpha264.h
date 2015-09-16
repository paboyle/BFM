
/*
 *
 *  Copyright Peter Boyle & Glasgow University, 2000.
 *  This software is provided for NON-COMMERCIAL use only,
 *  It is provided under the GNU pubic License V2
 *  It is provided as is and is not guaranteed fit for any purpose.
 *
 */

#ifndef _ALPHA_
#define _ALPHA_

/*21264 double precision*/
void create_alpha264_instructions( struct processor *proc);
struct processor *create_alpha264 (void);
int start_loop_alpha264(int counter);
void stop_loop_alpha264(int branchno,int counter);
void check_iterations_alpha_264(int cntreg,int retno);
void directive_print_alpha264(struct instruction *i);
void queue_cmovge_alpha264(int cond,int src,int dst);
void cache_print_alpha264(struct instruction *i);

/*21264 single precision*/
void create_alpha264s_instructions( struct processor *proc);
struct processor *create_alpha264 s(void);
int start_loop_alpha264s(int counter);
void stop_loop_alpha264s(int branchno,int counter);
void check_iterations_alpha_264s(int cntreg,int retno);
void directive_print_alpha264s(struct instruction *i);
void queue_cmovge_alpha264s(int cond,int src,int dst);
void cache_print_alpha264s(struct instruction *i);

/*21164 single precision*/
void create_alpha164s_instructions( struct processor *proc);
struct processor *create_alpha164 s(void);
int start_loop_alpha164s(int counter);
void stop_loop_alpha164s(int branchno,int counter);
void check_iterations_alpha_164s(int cntreg,int retno);
void directive_print_alpha164s(struct instruction *i);
void queue_cmovge_alpha164s(int cond,int src,int dst);
void cache_print_alpha164s(struct instruction *i);



#endif
