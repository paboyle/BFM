/*
 *
 *  Copyright Peter Boyle and Glasgow University, 2000.
 *
 *  It is provided under the GNU pubic License V2
 *  It is provided as is and is not guaranteed fit for any purpose.
 *
 */

#include "processor.h"


struct processor *PROC = NULL;

extern int argcount;
extern int argcopied[];


/**
 * MIC swizzle patterns used in assembly.
 */
const char *mic_swizzles[] =  // NB: order MUST correspond to enum MIC_SWIZ!
{
	"",  // "{dcba}" or empty string
	"{cdab}",
	"{badc}",
	"{dacb}",
	"{aaaa}",
	"{bbbb}",
	"{cccc}",
	"{dddd}"
};


/*These routines should be visible internally to this file only*/


void set_processor(struct processor *proc)
{
	PROC = proc;
}

void set_abi_caller_save_all(void)
{
	PROC->Isave = 0L;
	PROC->Fsave = 0L;
}

void test_proc(void)
{
	if (PROC == NULL)
	{
		puts("Error: Processor not set");
		exit(0);
	}
}


int is_x86() {
	return PROC->id == KNC || PROC->id == KNC_SINGLE;
}


/* ...}*/

/* {... LOOP management routines*/
int nxtlab = 0;

int get_target_label(void)
{
	return(nxtlab++);
}
/*
*Start loop and end loop are in architecture files since
* inner loop on power pc uses the CTR register for better
* branch prediction
*/
/*...}*/

/* {...Insert an instruction*/
void make_inst(int PIPE, int OP, ...)
{
	struct instruction *inst = new instruction;
	va_list ap;
	int i;
	int ivar;
	int argtype;
	int ignore_alloc;

	if (!inst)
	{
		perror("Malloc failed");
		exit(0);
	}
	test_proc();

	ignore_alloc = 0;
	if ((PIPE == STORPIPE && OP == SAVE_REG_EA)
		|| (PIPE == LOADPIPE && OP == LOAD_REG_EA)
		|| (PIPE == FSTRPIPE && OP == FSAVE_REG_EA)
		|| (PIPE == FLODPIPE && OP == FLOAD_REG_EA))
	{
		ignore_alloc = 1;
	}

	*inst = PROC->iset->instructions[PIPE][OP];

	// Detect undefined instruction mnemonics
	if (inst->mnemonic == NULL && PROC->id < UNKNOWN) {
		fprintf(stderr, "Error: Undefined mnemonic. pipe=%d, op=%d.\n", PIPE, OP);
		abort();
	}


	#if 0  /* DEBUGGING */
	printf("make_inst,args: ");
	for (int k = 0; k < inst->narg; k++)
		printf("%d ", inst->args[k]);
	printf("\n");
	#endif


	va_start(ap, OP);

	//printf("Info: make_inst: pipe=%d,op=%d,narg=%d\n", PIPE, OP, inst->narg);

	for (i = 0; i < inst->narg; i++)
	{
		/* MIC hack {... */
		if (PROC->id == KNC || PROC->id == KNC_SINGLE)
		{
			/* Skip one argument if its value is already pre-set by create_*_instructions().
			 * At the moment this applies to scalar floating-point ops on MIC only: Vector
			 * mask (args[1]) is already set.
			 * This hack saves us modifying more than 150 make_inst() calls and adding
			 * have_mic() switches.
			 */
			if ((PIPE == FMACPIPE || PIPE == FADDPIPE || PIPE == FMULPIPE || PIPE == FMOVPIPE)
				&& (i == 1))  // NB: mask is args[1]
			{
				fprintf(stderr, "Info: %s: Skipped mask arg for scalar flop on MIC VPU: pipe=%d,op=%d\n",
					__FUNCTION__, PIPE, OP);
				++i;
			}

			/* Some FSTOREs use mask which is the last arg, but which is not explicitly passed to make_inst(). Therefore we must break this loop in time for those stores. */
			if (PIPE == FSTRPIPE
				&& (OP == FSAVE_EA || OP == FSSAVE_EA)
				&& i == (inst->narg - 1))
			{
				fprintf(stderr, "Info: %s: Break before MIC mask (narg[%d]): pipe=%d,op=%d\n",
					__FUNCTION__, i, PIPE, OP);
				break;
			}
		}
		/* ...} MIC hack */


		argtype = inst->argtypes[i];
		switch (argtype)
		{
		case Fregs:
			ivar = va_arg(ap, int);
			if (ivar < 0 || ivar > PROC->fregs->regfile)
			{
				printf("Error: Invalid floating register %d arg %d\n", ivar, i);
				printf("instr: pipe=%d,code=%d\n", inst->pipe, inst->code);
				exit(0);
			}
			if (!(PROC->fregs->allocated[ivar]) && !ignore_alloc)
			{
				printf("Error: unallocated freg %d used operand %d\n", ivar, i);
				exit(0);
			}
			inst->args[i] = ivar;
			break;
		case Iregs:
			ivar = va_arg(ap, int);
			if (ivar < 0 || ivar > PROC->iregs->regfile)
			{
				printf("Error: Invalid register %d\n", ivar);
				printf("instr: pipe=%d,code=%d\n", inst->pipe, inst->code);
				exit(0);
			}
			if (!(PROC->iregs->allocated[ivar]) && !ignore_alloc)
			{
				printf("Error: unallocated reg %d used\n", ivar);
				printf("instr: pipe=%d,code=%d\n", inst->pipe, inst->code);
				exit(0);
			}
			inst->args[i] = ivar;
			break;
		case SIMM:
			/* Get one number argument */
			ivar = va_arg(ap, int);
			if (abs(get_offset(ivar)) > 16000 )
			{
				puts("Error: Immediate arg");
				exit(0);
			}
			inst->args[i] = ivar;
			break;
		case UIMM:
			ivar = va_arg(ap, int);
			inst->args[i] = ivar;
			break;
		case Cycles:
			/* Get one number argument */
			ivar = va_arg(ap, int);
			inst->args[i] = ivar;
			break;
		case Label:
			/* Get one number argument */
			if (inst->pipe == DIRECTIVE)
			{
				if (inst->code == Target || inst->code == CmovTarget)
				{
					ivar = va_arg(ap, int);
					if (ivar >= nxtlab || ivar < 0)
					{
						puts("Error: Unknown Label");
						exit(0);
					}
					inst->args[i] = ivar;
				}
				else if (inst->code == COMMENT)
				{
					if (sizeof(int) != sizeof(char *)) {
						fprintf(stderr, "Warning: Assembly comments hack might fail!\n");
					}
					inst->args[0] = va_arg(ap, int);  // pretty dirty HACK: interpret (char *) as (int) here and cast back to (char *) in directive_print_knc(). Should work fine if sizeof(int)==sizeof(void *).
				}
				else
				{
					inst->mnemonic = va_arg(ap, char *);
				}
			}
			else
			{
				ivar = va_arg(ap, int);
				if (ivar >= nxtlab || ivar < 0)
				{
					printf("Error: Unknown label %d,%d:%d %d\n", PIPE, OP, ivar, nxtlab);
					exit(0);
				}
				inst->args[i] = ivar;
			}
			break;
		case MICmask:
			ivar = va_arg(ap, int);
			if (ivar >= PROC->nmicmask || ivar < 0)
			{
				printf("Error: Invalid MIC mask %d\n", ivar);
				printf("instr: pipe=%d,code=%d\n", inst->pipe, inst->code);
				exit(0);
			}
			inst->args[i] = ivar;
			break;
		case MICswiz:
			ivar = va_arg(ap, int);
			if (ivar >= nmicswiz || ivar < 0)
			{
				printf("Error: Invalid MIC swizzle %d\n", ivar);
				exit(0);
			}
			inst->args[i] = ivar;
			break;
		default:
			puts("Error: bad argtype");
			exit(0);
		}
	}
	va_end(ap);

	/* And add the instruction to the queue */
	queue_instruction(inst);
}


/* ...End instruction insertion}*/



int get_I_size(void)
{
	return (PROC->I_size);
}

/* {...QUEUEING convenience routines */
void queue_load_addr( int dest, int immhandle, int src )
{
	if ( PROC->iset->instructions[LOADPIPE][LOAD_ADDR].mnemonic == NULL ){
		make_inst(IALUPIPE,IADD_IMM,dest,src,immhandle);
	}
	else {
		make_inst(LOADPIPE,LOAD_ADDR,dest,immhandle,src);
	}
}



/**
 * Loads single float/double into reg: "MOV reg <== offhandle(basereg)"
 */
void queue_fload(int reg, int offhandle, int basereg, enum Datum t)
{
	if (t == Double)
		make_inst(FLODPIPE,FLOAD_EA,reg,offhandle,basereg);
	else
		make_inst(FLODPIPE,FSLOAD_EA,reg,offhandle,basereg);
}


/**
 * \brief Stores scalar float.
 * \param reg The source register that is to be stored.
 * \param offhandle
 * \param basereg
 * \param t
 */
void queue_fstore( int reg, int offhandle, int basereg , enum Datum t)
{
	if ( t == Double )
		make_inst(FSTRPIPE,FSAVE_EA,reg,offhandle,basereg );
	else
		make_inst(FSTRPIPE,FSSAVE_EA,reg,offhandle,basereg );
}
void queue_iload( int reg, int offhandle, int basereg )
{
	make_inst(LOADPIPE,ILD_EA,reg,offhandle,basereg );
}
void queue_iload_short(int reg,int offhandle,int basereg)
{
  // Should be BGQ only?
  fprintf(stderr,"Queueing a short load\n");
  make_inst ( LOADPIPE,ILD_EA_SHORT,reg,offhandle,basereg);
}

/**
 * Stores reg at "basereg+offhandle".
 */
void queue_istore( int reg, int offhandle, int basereg )
{
	make_inst(STORPIPE,IST_EA,reg,offhandle,basereg );
}


/**
 * Scalar integer add: dest = a + b
 */
void queue_iadd(int dest, int a, int b)
{
	if (is_x86()) {
		// we don't have non-destructive scalar int add on x86-64
		if (dest == a)      make_inst(IALUPIPE, IADD, dest, b);
		else if (dest == b) make_inst(IALUPIPE, IADD, dest, a);
		else {
			make_inst(IALUPIPE, MOV, dest, a);
			make_inst(IALUPIPE, IADD, dest, b);
		}
	}
	else make_inst(IALUPIPE, IADD, dest, a, b);
}

/**
 * Scalar integer subtract: dest = a - b
 */
void queue_isub( int dest, int a, int b )
{
	if (is_x86()) {
		if (dest == a) make_inst(IALUPIPE, ISUB, dest, b);
		else if (dest == b) {
			// no reverse sub on x86 :(
			make_inst(IALUPIPE, ISUB, dest, a);  // dest = b - a
			make_inst(IALUPIPE, INEG, dest);     // dest = -dest
		}
		else {
			make_inst(IALUPIPE, MOV, dest, a);
			make_inst(IALUPIPE, ISUB, dest, b);
		}
	}
	else make_inst(IALUPIPE, ISUB, dest, a, b);
}

/**
 * Scalar signed integer multiply: dest = a * b
 */
void queue_imul(int dest, int a, int b)
{
	if (is_x86()) {
		if (dest == a)      make_inst(IALUPIPE, IMUL, dest, b);
		else if (dest == b) make_inst(IALUPIPE, IMUL, dest, a);
		else {
			make_inst(IALUPIPE, MOV, dest, a);
			make_inst(IALUPIPE, IMUL, dest, b);
		}
	}
	else make_inst(IALUPIPE, IMUL, dest, a, b);
}

/**
 * Scalar signed integer multiply with immediate: dest = a * himm
 */
void queue_imul_imm( int dest, int a, int himm )
{
	if (is_x86()) {
		// There's a three-operand form with 32-bit immediates available.
		// See Intel 64 Architectures Software Developer's Manual Vol. 2A - IMUL

		make_inst(IALUPIPE, IMUL_IMM, dest, himm, a);

		// If immediate is larger than 32 bits, we must use the following implementation. Adjust IMUL_IMM's narg etc.!
		/*
		make_inst(IALUPIPE, MOV, dest, a);
		make_inst(IALUPIPE, IMUL_IMM, dest, himm);
		*/
	}
	else make_inst(IALUPIPE, IMUL_IMM, dest, a, himm);
}

/**
 * Integer add with immediate: dest = src + imm
 */
void queue_iadd_imm( int dest, int src, int immhandle )
{
	if (is_x86()) {
		/* x86 doesn't have 3-operand scalar int add ... */
		if (dest == src) {
			/* This is: dest += imm */
			make_inst(IALUPIPE, IADD_IMM, dest, immhandle);
		}
		else if (get_offset(immhandle) == 0) {
			/* We effectively have: dest = src! This is used at least in write_vaxpy, so we implement an efficient branch here. Otherwise the less efficient mov+add branch would be chosen... */
			make_inst(IALUPIPE, MOV, dest, src);
			// QUESTION: MOV doesn't set FLAGS register. Is that a problem?
		}
		else {
			make_inst(IALUPIPE, MOV, dest, src);
			make_inst(IALUPIPE, IADD_IMM, dest, immhandle);
		}
	}
	else make_inst(IALUPIPE, IADD_IMM, dest, src, immhandle);
}

void queue_isub_imm( int dest, int src, int immhandle )
{
  int offset    = get_offset(immhandle);
  int negoffset = -offset;
  int negoffseth= get_offset_handle (negoffset,Byte);
  queue_iadd_imm(dest,src,negoffseth);
}


/**
 * AND with immediate: dest = src & immhandle
 */
void queue_iand_imm( int dest, int src, int immhandle )
{
	if (is_x86()) {
		if (dest == src) make_inst(IALUPIPE, IAND_IMM, dest, immhandle);
		else {
			make_inst(IALUPIPE, MOV, dest, src);
			make_inst(IALUPIPE, IAND_IMM, dest, immhandle);
		}
	}
	else make_inst(IALUPIPE, IAND_IMM, dest, src, immhandle);
}


/**
 * Logical right shift: dest = src >> get_offset(shifthandle)
 */
void queue_rshift(int dest, int src, int shifthandle)
{
  if (have_mic()) {
    if (dest == src) make_inst(IALUPIPE, RSHIFT, dest, shifthandle);
    else {
      queue_mov(src, dest);
      make_inst(IALUPIPE, RSHIFT, dest, shifthandle);
    }
  } else { 
    make_inst(IALUPIPE,RSHIFT,dest,src,shifthandle);
  }
}
void queue_lshift(int dst,int src,int himm)
{
  make_inst(IALUPIPE,LSHIFT,dst,src,himm);
}
void queue_ori   (int dst,int src,int huimm)
{
  make_inst(IALUPIPE,IOR_IMM,dst,src,huimm);
}
void queue_and   (int dst,int src1,int src2)
{
  make_inst(IALUPIPE,AND,dst,src1,src2);
}
void queue_andc   (int dst,int src1,int src2)
{
  make_inst(IALUPIPE,ANDC,dst,src1,src2);
}


/**
 * Scalar floating-point mul.
 * dest=a*b
 */
void queue_fmul(int dest, int a, int b)
{
	make_inst(FMULPIPE, FMUL, dest, a, b);
}


/**
 * Scalar floating-point add: dest = a + b.
 */
void queue_fadd(int dest, int a, int b)
{
	make_inst(FADDPIPE, FADD, dest, a, b);  // On MIC we use VPU for scalar FP ops => We must broadcast scalar to vector and extract one element from result... (The same with SUB, MUL etc)
}


/**
 * dest = a - b
 */
void queue_fsub(int dest, int a, int b)
{
	make_inst(FADDPIPE, FSUB, dest, a, b);
}


/**
 * Move floating-point register.
 */
void queue_fmov(int dest,int src)
{
	make_inst(FMOVPIPE,FMOV,dest,src);
}
void queue_fneg(int dest,int src)
{
  make_inst(FMOVPIPE,FNEG,dest,src);
}
void queue_iload_imm(int dest,int imm)
{
	make_inst(LOADPIPE,LOAD_IMM,dest,imm);
}


/* ...END queueing convenience routines */


/**
 * \see enter_routine()
 * \see exit_routine()
 */
static char const *routine_name = NULL;


/**
 * \brief Writes function preamble, grabs stack, saves regs.
 * \param[in] name The routine name.
 * \param[in] nbyte The number of bytes to be grabbed on stack.
 */
void enter_routine(const char *const name, int nbyte)
{
	if (routine_name != NULL) {
		fprintf(stderr, "Error: Entering new routine without having exited the previous one.\n");
		abort();
	}

	routine_name = name;
	make_inst(DIRECTIVE, Enter_Routine, name);
	grab_stack(nbyte);
	save_regs();
}


/**
 * \brief Restores regs, frees stack, writes function exit.
 *
 * Also checks the numbers of arguments used and defined.
 */
void exit_routine()
{
	int i, count_args_copied;

	if (routine_name == NULL) {
		fprintf(stderr, "Error: Trying to exit a routine.\n");
		abort();
	}

	restore_regs();
	free_stack();
	make_inst(DIRECTIVE, Exit_Routine, routine_name);
	routine_name = NULL;

	// check usage of arguments
	count_args_copied = 0;
	for (i = 0; i < MAX_FUNC_ARG; i++) {
		if (argcopied[i])  // count number of args copied
			++count_args_copied;
	}
	if (count_args_copied < argcount) {
		fprintf(stderr, "Warning: %s: defargcount is %d, but used only %d args.\n",
			routine_name, argcount, count_args_copied);
	}
	else if (count_args_copied > argcount) {
		fprintf(stderr, "Error: %s: More args used than defined by defargcount().\n",
			routine_name);
		abort();
	}
}


/**
 * \brief Sets a label / jump target.
 */
void target(int label)
{
	make_inst(DIRECTIVE, Target, label);
}




void queue_barrier(void)
{
	make_inst(DIRECTIVE,Pipe_Flush, "Barrier");
}
void queue_mbar(void)
{
	make_inst(DIRECTIVE,LS_BARRIER, "mbar");
}

/**
 * Copies contents from reg src to reg dst. Both registers are non-vector regs.
 * NOTE: order of arguments is src, dst!!!
 */
void queue_mov(int src, int dst)
{
	if (is_x86())
		make_inst(IALUPIPE, MOV, dst, src);
	else
		make_inst(IALUPIPE, IOR, dst, src, src);
}

void pragma(int code,int val)
{
	make_inst(PRAGMA,code,val);
}

void conditional_branch_cmpzero(int branch_type, int condreg, int target  )
{
	if ( PROC->iset->instructions[IALUPIPE][IOR_PREDICATE].mnemonic !=NULL ) {
		if (is_x86()) {
			make_inst(IALUPIPE, IOR_PREDICATE, condreg, condreg);
		}
		else {
			make_inst(IALUPIPE,IOR_PREDICATE,condreg,condreg,condreg); /*Sets the predicate reg*/
		}
		make_inst(BRCHPIPE,branch_type,target);
	} else {
		make_inst(BRCHPIPE,branch_type,condreg,target);
	}
}
/* ...}*/


