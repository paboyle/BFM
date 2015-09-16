/*
 *
 *  Copyright Peter Boyle and Glasgow University 2000.
 *
 *  It is provided under the GNU pubic License V2
 *  It is provided as is and is not guaranteed fit for any purpose.
 *
 */
#include "processor.h"
extern struct processor *PROC;

/* {... ARGUMENT management routines */
int argcount = -1;
int argcopied [MAX_FUNC_ARG];
static int curarg=-1;
static int argbase;

int args_on_stack;
extern int frame_up;


int defargcount(int narg)
{
  char *str;
  int i;

  if (narg > MAX_FUNC_ARG) {
	fprintf(stderr, "Error: Too much arguments.\n");
	abort();
  }

  test_proc();
  curarg = 0;
  argcount = narg;
  argbase = PROC->ABI_argbase;
  args_on_stack = 0;  // counts number of args passed via stack
  for (i = 0; i < MAX_FUNC_ARG; i++) {
	argcopied[i] = 0;
  }

  for(i=0;i<narg;i++) {
	if (i >= PROC->MAX_ARGS_BY_REGS) {
		++args_on_stack;
		continue;
	}

    str = (char *)malloc(8);
    sprintf(str,"Arg%d",i);
    PROC->iregs->allocated[argbase+i] = 1;
    PROC->iregs->mnemonic[argbase+i] = str;
  }
  initialise_streams();
  return(0);
}

void  getarg(int ireg)
{
  int argreg;

  argcopied[curarg] = 1;
  if (curarg >= PROC->MAX_ARGS_BY_REGS)
  {
	fprintf(stderr, "Info: Loading arg %d from stack into Ireg %d.\n", curarg, ireg);
	int offhandle = def_offset(PROC->IREGS_SIZE * (curarg - PROC->MAX_ARGS_BY_REGS + 1 + 1), Byte, "ARG_OFFSET");  // extra +1 because we do a PUSHQ at the beginning which moves the SP.
	int scratch = allocate_tmp_reg(Iregs);
	queue_iload(scratch, frame_up, PROC->StackPointer);  // load "old" SP into scratch
	queue_iload(ireg, offhandle, scratch);  // load arg value from stack (relative to old SP)
	free_reg(scratch, Iregs);
  }
  else {
    argreg = curarg+argbase;
    if ( ireg > argreg && ireg < argbase + argcount ){
      printf("Ouch: register allocation over top of arguments\n");
      exit(0);
    }
    queue_mov(argreg, ireg);
    free_reg(Iregs,argreg);
  }
  curarg ++;
}
/* ...} */
