
/*
 *
 *  Copyright Peter Boyle & Glasgow University, 2000.
 *  This software is provided for NON-COMMERCIAL use only,
 *  It is provided under the GNU pubic License V2
 *  It is provided as is and is not guaranteed fit for any purpose.
 *
 */

#ifndef _POWER_PC_
#define _POWER_PC_

/*PPC 440 double precision*/
void create_ppc440_instructions( struct processor *proc);
struct processor *create_ppc440 (void);
int start_loop_ppc440(int counter);
void stop_loop_ppc440(int branchno,int counter);
void check_iterations_ppc440(int cntreg,int retno);
void directive_print_ppc440(struct instruction *i);
void queue_cmovge_ppc440(int cond,int src,int dst);
void cache_print_ppc440(struct instruction *i);
void cache_print_bgq(struct instruction *i);

/*PPC 440 single precision*/
void create_ppc440s_instructions( struct processor *proc);
struct processor *create_ppc440s(void);
int start_loop_ppc440s(int counter);
void stop_loop_ppc440s(int branchno,int counter);
void check_iterations_ppc440s(int cntreg,int retno);
void directive_print_ppc440s(struct instruction *i);
void queue_cmovge_ppc440s(int cond,int src,int dst);
void cache_print_ppc440s(struct instruction *i);


/*PowerIII double precision*/
void create_powerIII_instructions( struct processor *proc);
struct processor *create_powerIII (void);
int start_loop_powerIII(int counter);
void stop_loop_powerIII(int branchno,int counter);
void check_iterations_powerIII(int cntreg,int retno);
void directive_print_powerIII(struct instruction *i);
void queue_cmovge_powerIII(int cond,int src,int dst);
void cache_print_powerIII(struct instruction *i);

/*PowerIII single precision*/
void create_powerIIIs_instructions( struct processor *proc);
struct processor *create_powerIIIs (void);
int start_loop_powerIIIs(int counter);
void stop_loop_powerIIIs(int branchno,int counter);
void check_iterations_powerIIIs(int cntreg,int retno);
void directive_print_powerIIIs(struct instruction *i);
void queue_cmovge_powerIIIs(int cond,int src,int dst);
void cache_print_powerIIIs(struct instruction *i);

/*bgl goodies*/
struct processor *create_bgl440 ();
void create_bgl440_instructions( struct processor *proc);

/*bgq goodies*/
struct processor *create_bgq ();
void create_bgq_instructions( struct processor *proc);

void queue_cmovgt_bgq(int cond,int src,int dst);
void queue_cmovlt_bgq(int cond,int src,int dst);
void queue_cmovge_bgq(int cond,int src,int dst);
void queue_cmovle_bgq(int cond,int src,int dst);

#endif
