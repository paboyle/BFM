/*
 *
 *  Copyright Peter Boyle and Glasgow University, 2000.
 *
 *  It is provided under the GNU pubic License V2
 *  It is provided as is and is not guaranteed fit for any purpose.
 *
 */

#include "processor.h"

extern struct processor *PROC;
extern int args_on_stack;
int N_Ahead = 1;
int frame_down;
int frame_up;
static int frame_regbase;
static int frame_whbase;/*Prefetch bit bucket*/
int frame_membase;  // NB: referenced in unknown.C
/*
 * Some cpu's have rules about the stack pointer being aligned on a
 * 16 byte boundary. (US-II) If we ever support something with other
 * alignment constraints, easiest to just use lcm as a unified constraint.
 * Return offset from stack pointer for *LOW* end of allocated buffer
 */
// change to 64 just to make sure that MIC vectors are stored well aligned
#define STACK_ALIGN 64


/*
 * This routine returns a byte offset handle to a bias from the stack pointer
 * that gives the *low* end of the alloced space
 *
 * The registers are saved in a space *above* the allocated memory
 * If REGWINDOWS are used we cannot allocate stack portably, and raise an error
 * If no regwindows are used then we can only produce leaf routines on US-II
 *
 * We guarantee it is safe to both prefetch and write hint to the stack pointer
 *
 */
int grab_stack ( int nbyte )
{
  //int unitsize = PROC->FP_size;
  int frame;
  int bias;
  int alignbyte;


  if (PROC->id >= UNKNOWN ) {  // Note: potentially dangerous if we define new procs after UNKNOWN.
    def_off ( bias , Byte, 0) ;
    /*    fprintf(stderr,"nbyte =%d\n",nbyte );*/
    frame_membase = nbyte;
    return bias;
  }

  /*Make sure we only align frames on (STACK_ALIGN byte) bdy's*/

  alignbyte = STACK_ALIGN* ( ((nbyte+STACK_ALIGN-1)/STACK_ALIGN) );
  if ( alignbyte != nbyte ) {
    fprintf(stderr, "grab_stack_alloc\n");
    fprintf(stderr, "Warning aligning stack frame: %d bytes requested\n",nbyte);
  }


  /*
   * We could redefine the locations of mem and regs to match the reg windows
   * on US-II,III
   * Haven't done as yet, and this restricts us to leaf routines
   * on US...
   */

  if (args_on_stack ) {
    fprintf(stderr,"oops args are on stack\n");
    fflush(stderr);
    exit(-1);
  }
  
  int frame_reg_bytes = 
     (PROC->iregs->regfile * PROC->IREGS_SIZE)
    +(PROC->fregs->regfile * PROC->FREGS_SIZE);

  frame_membase = 0;
  frame_regbase = frame_membase + alignbyte;

  frame = alignbyte + frame_reg_bytes;

  def_off( frame_down_tmp, Byte, -frame ); frame_down= frame_down_tmp;
  def_off( frame_up_tmp, Byte, frame);      frame_up = frame_up_tmp;

  /*SPARC v9 has an annoying 2047 implicit offset from sp, to stack addr   */
  /* this is called bias in documentation. Hence PROC->Bias member         */
  /* I corrupt this terminology  slightly to                               */
  /* mean a bias to sp that points to the data storage on stack for routine*/
  bias = def_offset ( PROC->Bias+frame_membase, Byte, "bias" );

  if ( PROC->RegWindows ) {
    make_inst(DIRECTIVE,PUSH_WINDOW,frame_down);
    return(bias);
  }else {
    queue_load_addr(PROC->StackPointer,frame_down,PROC->StackPointer);
    return(bias);
  }
}

void free_stack ( void )
{
  if (PROC->id >= UNKNOWN ) return;

  if ( PROC->RegWindows ) {
    make_inst(DIRECTIVE,POP_WINDOW,frame_up);
  }
  else{
    queue_load_addr(PROC->StackPointer,frame_up,PROC->StackPointer);
  }
}

static int SAVE_I[32];
static int SAVE_F[32];


void save_regs(void)
{
  int i;
  char *str;

  if (PROC->id >= UNKNOWN ) return;

  if ( !PROC->RegWindows ){
    for ( i = 0 ; i <PROC->iregs->regfile ; i ++ ){
      if ( PROC->Isave & (1<<i) ) {
        str=(char *)malloc(10);
        sprintf(str,"S_I_%d",i);

	/*Add in PROC->Bias for US-II in v9 mode*/

			SAVE_I[i] = def_offset( PROC->Bias 
					       +frame_regbase 
					       +(PROC->fregs->regfile)*PROC->FREGS_SIZE
					        + i*PROC->IREGS_SIZE, Byte, str);
			make_inst(STORPIPE, SAVE_REG_EA, i, SAVE_I[i], PROC->StackPointer);
      }
    }

    for ( i = 0 ; i <PROC->fregs->regfile ; i ++ ){
      if ( PROC->Fsave & (1<<i) ) {
        str=(char *)malloc(10);
        sprintf(str,"S_F_%d",i);

	/*Add in PROC->Bias for US-II in v9 mode*/

		if (PROC->id == KNC || PROC->id == KNC_SINGLE) {
			// NB: 64 is size of MIC vector regs
		  SAVE_F[i] =
            def_offset(frame_regbase + i*64 + PROC->Bias, Byte, str);
		}
		else {
		  SAVE_F[i] =
            def_offset(frame_regbase + i*8 + PROC->Bias, Byte, str);
		}
		//printf("FSAVE_REG: %d, %d, %d\n", i, SAVE_F[i], PROC->StackPointer);  // DEBUGGING
        make_inst(FSTRPIPE,FSAVE_REG_EA,i,SAVE_F[i],PROC->StackPointer);
     }
    }
  }
}

void restore_regs(void)
{
  int i;

  if (PROC->id >= UNKNOWN ) return;

  if ( ! PROC->RegWindows ){
	// restore regs in reverse order as they were stored.
	// restore fregs first as they where saved last.
	// Order doesn't matter because I decided not to use PUSHQ/POPQ on MIC...

	for ( i = PROC->fregs->regfile - 1; i >= 0; i -- ){
      if ( PROC->Fsave & (1<<i) ) {
        make_inst(FLODPIPE,FLOAD_REG_EA,i,SAVE_F[i],PROC->StackPointer);
      }
    }

    for ( i = PROC->iregs->regfile - 1; i >= 0; i -- ){
      if ( PROC->Isave & (1<<i) ) {
		  make_inst(LOADPIPE,LOAD_REG_EA,i,SAVE_I[i],PROC->StackPointer);
      }
    }
  }
}
