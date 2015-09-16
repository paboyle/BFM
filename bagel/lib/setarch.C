// TODO: a lot more MIC support to be included!!!!!!!
/*
 *
 *  Copyright UKQCD Collaboration October 2000.
 *  Written by Peter Boyle.
 *
 *  It is provided under the GNU pubic License V2
 *  It is provided as is and is not guaranteed fit for any purpose.
 *
 */

#include "processor.h"
#include "alpha.h"
#include "mic.h"
#include "powerpc.h"
#include "usparc.h"
#include "unknown.h"
#include <stdio.h>
#include <string.h>


extern struct processor *PROC;


// NOTE: must have corresponding elements to processor.h: enum procids!!!
const char *procargs[] =
{
	"alpha264",
	"alpha264s",
	"alpha164",
	"alpha164s",
	"alpha064",
	"alpha064s",
	"knc",   /* Intel MIC Knights Corner */
	"kncs",  /* Intel MIC Knights Corner, single precision */
	"ppc440",
	"ppc440s",
	"bgl",
	"bgq",
	"powerIII",
	"powerIIIs",
	"usparcII",
	"usparcIIs",
	"noarch",
	"noarchs",
	"noarch-bgq",
	"noarch-vmx"
};


void fix_short_integer(struct processor *proc);

void queue_cmovle_alpha264(int cond, int src, int dst);
void queue_cmovgt_alpha264(int cond, int src, int dst);
void queue_cmovlt_alpha264(int cond, int src, int dst);
void queue_cmovle_knc(int cond, int src, int dst);
void queue_cmovgt_knc(int cond, int src, int dst);
void queue_cmovlt_knc(int cond, int src, int dst);
void queue_cmovgt_powerIII(int cond, int src, int dst);
void queue_cmovlt_powerIII(int cond, int src, int dst);
void queue_cmovle_powerIII(int cond, int src, int dst);
void queue_cmovgt_powerIIIs(int cond, int src, int dst);
void queue_cmovlt_powerIIIs(int cond, int src, int dst);
void queue_cmovle_powerIIIs(int cond, int src, int dst);
void queue_cmovgt_ppc440(int cond, int src, int dst);
void queue_cmovlt_ppc440(int cond, int src, int dst);
void queue_cmovle_ppc440(int cond, int src, int dst);
void queue_cmovle_unknown(int cond, int src, int dst);
void queue_cmovgt_unknown(int cond, int src, int dst);
void queue_cmovlt_unknown(int cond, int src, int dst);
void queue_cmovgt_usparcII(int cond, int src, int dst);
void queue_cmovle_usparcII(int cond, int src, int dst);
void queue_cmovlt_usparcII(int cond, int src, int dst);

void select_processor_code(int procid);
void add_pragmas(void);


/* {...Set the processor type */
void set_processor_optarg(char *opt)
{
	struct processor *myproc;
	int i;

	// check whether opt is NULL! I had a SEGFAULT (when procargs was not prepared for MIC though).
	if (opt == NULL) {
		fprintf(stderr, "Error: No processor name specified!\n");
		abort();
	}

	for (i = 0; i < nprocids; i++)
	{
		if (!(strcmp(opt, procargs[i])))
		{
			myproc = create_processor(i);
			myproc->id = i;
			fix_short_integer(myproc);
			if ((i == UNKNOWN_BGQ) || (i == UNKNOWN_VMX))
			{
				myproc->id = UNKNOWN;
			}
			set_processor(myproc);
			set_label_prefix("_%s_lab%d");
			add_pragmas();
			return;
		}
	}
	for (i = 0; i < nprocids; i++)
	{
		printf("\t-P%s\n", procargs[i]);
	}
	exit(0);
}


void fix_short_integer(struct processor *proc)
{
  if ( (proc->id != PPCBGQ)  ) { // Use int == long == void * away from BGQ
    fprintf(stderr,"*** Fixing short integers to be long\n");
    proc->iset->instructions[LOADPIPE][ILD_EA_SHORT]
      =proc->iset->instructions[LOADPIPE][ILD_EA];
    proc->IS_size = proc->I_size;
  }
}


/*Processor dependent functions -> use set of function pointers*/
struct processor *create_processor(int id)
{
	switch (id)
	{
	case ALPHA264:
		return create_alpha264();
		break;
	case ALPHA164:
		return create_alpha164();
		break;
	case ALPHA064:
		return create_alpha064();
		break;
	case ALPHA264_SINGLE:
		return create_alpha264s();
		break;
	case ALPHA164_SINGLE:
		return create_alpha164s();
		break;
	case ALPHA064_SINGLE:
		return create_alpha064s();
		break;
	case KNC:  // FIXME: distinguish between single/double precision
	case KNC_SINGLE:
		return create_knc();
		break;
	case PPC440:
		return create_ppc440();
		break;
	case PPC440_SINGLE:
		return create_ppc440s();
		break;
	case POWERIII:
		return create_powerIII();
		break;
	case POWERIII_SINGLE:
		return create_powerIIIs();
		break;
	case PPC440BGL:
		return create_bgl440();
		break;
	case PPCBGQ:
		return create_bgq();
		break;
	case USPARCII:
		return create_usparcII();
		break;
	case USPARCII_SINGLE:
		return create_usparcIIs();
		break;
	case UNKNOWN_SINGLE:
		set_unknown_type_float();
		return create_unknown();
		break;
	case UNKNOWN:
		set_unknown_type_double();
		return create_unknown();
		break;
	case UNKNOWN_BGQ:
		set_unknown_type_double();
		set_unknown_bgq();
		return create_unknown();
		break;
	case UNKNOWN_VMX:
		set_unknown_type_double();
		set_unknown_vmx();
		return create_unknown();
		break;
	default:
		printf("Bad type of processor\n");
		exit(-1);
		break;
	}
}


void add_pragmas(void)
{
	struct instruction *pragmas = new instruction[npragmas];
	int i;

	pragmas[DCBT_SPACE].mnemonic = (char *)"pragma_dcbt_space";
	pragmas[DCBT_PRE].mnemonic   = (char *)"pragma_dcbt_pre";
	pragmas[DCBT_POST].mnemonic  = (char *)"pragma_dcbt_post";
	pragmas[LOAD_LIM].mnemonic   = (char *)"pragma_load_lim";
	pragmas[LOAD_SPACE].mnemonic = (char *)"pragma_load_space";
	pragmas[STORE_LIM].mnemonic  = (char *)"pragma_store_lim";
	pragmas[LS_LIM].mnemonic     = (char *)"pragma_ls_lim";
	pragmas[STORE_INORDER].mnemonic = (char *)"pragma_store_inorder";

	for (i = 0; i < npragmas; i++)
	{
		pragmas[i].narg = 1;
		pragmas[i].argtypes[0] = Cycles;
		pragmas[i].inout[0]    = inarg;
		pragmas[i].pipe = PRAGMA;
		pragmas[i].code = i;
	}
	PROC->iset->instructions[PRAGMA] = pragmas;
}


int start_loop(int counter)
{
	switch (PROC->id)
	{
	case ALPHA264:
		return start_loop_alpha264(counter);
		break;
	case ALPHA164:
		return start_loop_alpha164(counter);
		break;
	case ALPHA064:
		return start_loop_alpha064(counter);
		break;
	case ALPHA264_SINGLE:
		return start_loop_alpha264s(counter);
		break;
	case ALPHA164_SINGLE:
		return start_loop_alpha164s(counter);
		break;
	case ALPHA064_SINGLE:
		return start_loop_alpha064s(counter);
		break;
	case KNC:  // FIXME: single/double precision
	case KNC_SINGLE:
		return start_loop_knc(counter);
		break;
	case PPC440:
	case PPC440BGL:
	case PPCBGQ:
		return start_loop_ppc440(counter);
		break;
	case PPC440_SINGLE:
		return start_loop_ppc440s(counter);
		break;
	case POWERIII:
		return start_loop_powerIII(counter);
		break;
	case POWERIII_SINGLE:
		return start_loop_powerIIIs(counter);
		break;
	case USPARCII:
		return start_loop_usparcII(counter);
		break;
	case USPARCII_SINGLE:
		return start_loop_usparcIIs(counter);
		break;
	default:
		return start_loop_unknown(counter);
		break;
	}
}


void stop_loop(int branchno, int counter)
{
	switch (PROC->id)
	{
	case ALPHA264:
		stop_loop_alpha264(branchno, counter);
		break;
	case ALPHA164:
		stop_loop_alpha164(branchno, counter);
		break;
	case ALPHA064:
		stop_loop_alpha064(branchno, counter);
		break;
	case ALPHA264_SINGLE:
		stop_loop_alpha264s(branchno, counter);
		break;
	case ALPHA164_SINGLE:
		stop_loop_alpha164s(branchno, counter);
		break;
	case ALPHA064_SINGLE:
		stop_loop_alpha064s(branchno, counter);
		break;
	case KNC:  // FIXME: single/double precision
	case KNC_SINGLE:
		stop_loop_knc(branchno, counter);
		break;
	case PPC440:
	case PPC440BGL:
	case PPCBGQ:
		stop_loop_ppc440(branchno, counter);
		break;
	case PPC440_SINGLE:
		stop_loop_ppc440s(branchno, counter);
		break;
	case POWERIII:
		stop_loop_powerIII(branchno, counter);
		break;
	case POWERIII_SINGLE:
		stop_loop_powerIIIs(branchno, counter);
		break;
	case USPARCII:
		stop_loop_usparcII(branchno, counter);
		break;
	case USPARCII_SINGLE:
		stop_loop_usparcIIs(branchno, counter);
		break;
	default:
		stop_loop_unknown(branchno, counter);
		break;
	}
}


void directive_print(struct instruction *i)
{
	if (i->code == LS_BARRIER)
	{
		return;
	}

	switch (PROC->id)
	{
	case ALPHA264:
		directive_print_alpha264(i);
		break;
	case ALPHA164:
		directive_print_alpha164(i);
		break;
	case ALPHA064:
		directive_print_alpha064(i);
		break;
	case ALPHA264_SINGLE:
		directive_print_alpha264s(i);
		break;
	case ALPHA164_SINGLE:
		directive_print_alpha164s(i);
		break;
	case ALPHA064_SINGLE:
		directive_print_alpha064s(i);
		break;
	case KNC:  // FIXME: single/double precision
	case KNC_SINGLE:
		directive_print_knc(i);
		break;
	case PPC440:
	case PPC440BGL:
	case PPCBGQ:
		directive_print_ppc440(i);
		break;
	case PPC440_SINGLE:
		directive_print_ppc440s(i);
		break;
	case POWERIII:
		directive_print_powerIII(i);
		break;
	case POWERIII_SINGLE:
		directive_print_powerIIIs(i);
		break;
	case USPARCII:
		directive_print_usparcII(i);
		break;
	case USPARCII_SINGLE:
		directive_print_usparcIIs(i);
		break;
	default:
		directive_print_unknown(i);
		break;
	}
}


void queue_cmovge(int cond, int src, int dst)
{
	switch (PROC->id)
	{
	case ALPHA264:
		queue_cmovge_alpha264(cond, src, dst);
		break;
	case ALPHA164:
		queue_cmovge_alpha164(cond, src, dst);
		break;
	case ALPHA064:
		queue_cmovge_alpha064(cond, src, dst);
		break;
	case ALPHA264_SINGLE:
		queue_cmovge_alpha264s(cond, src, dst);
		break;
	case ALPHA164_SINGLE:
		queue_cmovge_alpha164s(cond, src, dst);
		break;
	case ALPHA064_SINGLE:
		queue_cmovge_alpha064s(cond, src, dst);
		break;
	case KNC:  // FIXME: single/double precision
	case KNC_SINGLE:
		queue_cmovge_knc(cond, src, dst);
		break;
	case PPC440:
	case PPC440BGL:
		queue_cmovge_ppc440(cond, src, dst);
		break;
	case PPCBGQ:
#define USE_ISEL
#ifdef USE_ISEL
		queue_cmovge_bgq(cond, src, dst);
#else
		queue_cmovge_ppc440(cond, src, dst);
#endif
		break;
	case PPC440_SINGLE:
		queue_cmovge_ppc440s(cond, src, dst);
		break;
	case POWERIII:
		queue_cmovge_powerIII(cond, src, dst);
		break;
	case POWERIII_SINGLE:
		queue_cmovge_powerIIIs(cond, src, dst);
		break;
	case USPARCII:
		queue_cmovge_usparcII(cond, src, dst);
		break;
	case USPARCII_SINGLE:
		queue_cmovge_usparcIIs(cond, src, dst);
		break;
	default:
		queue_cmovge_unknown(cond, src, dst);
		break;
	}
}


void cache_print(struct instruction *i)
{
	switch (PROC->id)
	{
	case ALPHA264:
		cache_print_alpha264(i);
		break;
	case ALPHA164:
		cache_print_alpha164(i);
		break;
	case ALPHA064:
		cache_print_alpha064(i);
		break;
	case ALPHA264_SINGLE:
		cache_print_alpha264s(i);
		break;
	case ALPHA164_SINGLE:
		cache_print_alpha164s(i);
		break;
	case ALPHA064_SINGLE:
		cache_print_alpha064s(i);
		break;
	case KNC:  // FIXME: single/double precision
	case KNC_SINGLE:
		cache_print_knc(i);
		break;
	case PPC440:
	case PPC440BGL:
		cache_print_ppc440(i);
		break;
	case PPCBGQ:
		cache_print_bgq(i);
		break;
	case PPC440_SINGLE:
		cache_print_ppc440s(i);
		break;
	case POWERIII:
		cache_print_powerIII(i);
		break;
	case POWERIII_SINGLE:
		cache_print_powerIIIs(i);
		break;
	case USPARCII:
		cache_print_usparcII(i);
		break;
	case USPARCII_SINGLE:
		cache_print_usparcIIs(i);
		break;
	default:
		break;
	}
}


void check_iterations(int cntreg, int retno)
{
	switch (PROC->id)
	{
	case ALPHA264:
		check_iterations_alpha264(cntreg, retno);
		break;
	case ALPHA164:
		check_iterations_alpha164(cntreg, retno);
		break;
	case ALPHA064:
		check_iterations_alpha064(cntreg, retno);
		break;
	case ALPHA264_SINGLE:
		check_iterations_alpha264s(cntreg, retno);
		break;
	case ALPHA164_SINGLE:
		check_iterations_alpha164s(cntreg, retno);
		break;
	case ALPHA064_SINGLE:
		check_iterations_alpha064s(cntreg, retno);
		break;
	case KNC:  // FIXME: single/double precision
	case KNC_SINGLE:
		check_iterations_knc(cntreg, retno);
		break;
	case PPC440:
	case PPC440BGL:
	case PPCBGQ:
		check_iterations_ppc440(cntreg, retno);
		break;
	case PPC440_SINGLE:
		check_iterations_ppc440s(cntreg, retno);
		break;
	case POWERIII:
		check_iterations_powerIII(cntreg, retno);
		break;
	case POWERIII_SINGLE:
		check_iterations_powerIIIs(cntreg, retno);
		break;
	case USPARCII:
		check_iterations_usparcII(cntreg, retno);
		break;
	case USPARCII_SINGLE:
		check_iterations_usparcIIs(cntreg, retno);
		break;
	default:
		check_iterations_unknown(cntreg, retno);
		break;
	}
}


/*...}*/


void queue_cmovgt(int cond, int src, int dst)
{
	switch (PROC->id)
	{
	case ALPHA064:
	case ALPHA264_SINGLE:
	case ALPHA164_SINGLE:
	case ALPHA064_SINGLE:
	case ALPHA264:
	case ALPHA164:
		queue_cmovgt_alpha264(cond, src, dst);
		break;
	case KNC:  // FIXME: single/double precision
	case KNC_SINGLE:
		queue_cmovgt_knc(cond, src, dst);
		break;
	case PPCBGQ:
#ifdef USE_ISEL
		queue_cmovgt_bgq(cond, src, dst);
		break;
#endif
	case PPC440:
	case PPC440BGL:
	case PPC440_SINGLE:
		queue_cmovgt_ppc440(cond, src, dst);
		break;
	case POWERIII:
	case POWERIII_SINGLE:
		queue_cmovgt_powerIII(cond, src, dst);
		break;
	case USPARCII:
	case USPARCII_SINGLE:
		queue_cmovgt_usparcII(cond, src, dst);
		break;
	default:
		queue_cmovgt_unknown(cond, src, dst);
		break;
	}
}


void queue_cmovle(int cond, int src, int dst)
{
	switch (PROC->id)
	{
	case ALPHA064:
	case ALPHA264_SINGLE:
	case ALPHA164_SINGLE:
	case ALPHA064_SINGLE:
	case ALPHA264:
	case ALPHA164:
		queue_cmovle_alpha264(cond, src, dst);
		break;
	case KNC:  // FIXME: single/double precision
	case KNC_SINGLE:
		queue_cmovle_knc(cond, src, dst);
		break;
	case PPCBGQ:
#ifdef USE_ISEL
		queue_cmovle_bgq(cond, src, dst);
		break;
#endif
	case PPC440:
	case PPC440BGL:
	case PPC440_SINGLE:
		queue_cmovle_ppc440(cond, src, dst);
		break;
	case POWERIII:
	case POWERIII_SINGLE:
		queue_cmovle_powerIII(cond, src, dst);
		break;
	case USPARCII:
	case USPARCII_SINGLE:
		queue_cmovle_usparcII(cond, src, dst);
		break;
	default:
		queue_cmovle_unknown(cond, src, dst);
		break;
	}

}


void queue_cmovlt(int cond,int src,int dst)
{
	switch (PROC->id)
	{
	case ALPHA064:
	case ALPHA264_SINGLE:
	case ALPHA164_SINGLE:
	case ALPHA064_SINGLE:
	case ALPHA264:
	case ALPHA164:
		queue_cmovlt_alpha264(cond, src, dst);
		break;
	case KNC:  // FIXME: single/double precision
	case KNC_SINGLE:
		queue_cmovlt_knc(cond, src, dst);
		break;
	case PPCBGQ:
#ifdef USE_ISEL
		queue_cmovlt_bgq(cond, src, dst);
		break;
#endif
	case PPC440:
	case PPC440BGL:
	case PPC440_SINGLE:
		queue_cmovlt_ppc440(cond, src, dst);
		break;
	case POWERIII:
	case POWERIII_SINGLE:
		queue_cmovlt_powerIII(cond, src, dst);
		break;
	case USPARCII:
	case USPARCII_SINGLE:
		queue_cmovlt_usparcII(cond, src, dst);
		break;
	default:
		queue_cmovlt_unknown(cond, src, dst);
		break;
	}
}
