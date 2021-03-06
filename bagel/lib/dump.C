/*
 *
 *  Copyright Peter Boyle and Glasgow University 2000.
 *  It is provided under the GPL as is and is not guaranteed fit for any purpose.
 *
 */

#include "processor.h"
#include <string.h>

extern struct processor *PROC;
extern struct inst_queue_entry *HEAD;
void stall_inst_print(struct instruction *i);

void arg_sprint(struct instruction *i, int j,char *buf);
void inst_print_unknown(struct instruction *i);

/*{... OUTPUT routines*/
int human_readable = 0;
void set_human_readable (int h)
{
  human_readable = h;
}

static char *label_prefix;
static char *label_fmt = (char *)"_%s_lab%d";
void set_label_fmt(char *fmt)
{
  label_fmt = fmt;
}
void set_label_prefix(const char *pfx)
{
  label_prefix = (char *)pfx;
}
char *make_label(int i)
{
  static char buf[84];
  if (label_prefix)
    sprintf(buf,label_fmt,label_prefix,i);
  else
    sprintf(buf,".lab%d",i);

  return buf;
}


/**
 * Maximum bytes of a (final) register name (eg. "%zmm31") including terminating '\0'.
 */
#define MAX_REGNAME_LEN 8


char *make_reg_name (int kind, int reg )
{
  char * str = (char *)malloc(MAX_REGNAME_LEN);
  switch(kind){
    case Iregs:
		if (PROC->id == KNC || PROC->id == KNC_SINGLE)
		{
			// MIC is x86-based architecture with ugly reg names...
			switch (reg)
			{
				// the regs of cases 0..5 comply with ABI for passing arguments on 64-bit Linux
				case 0:  sprintf(str, "%%rdi"); break;
				case 1:  sprintf(str, "%%rsi"); break;
				case 2:  sprintf(str, "%%rdx"); break;
				case 3:  sprintf(str, "%%rcx"); break;
				case 4:  sprintf(str, "%%r8");  break;
				case 5:  sprintf(str, "%%r9");  break;

				// assign the rest of "uglily named" regs. (No special order here, but must comply with definition of callee-saved regs and stack pointer in mic.C::create_knc()!!)
				case 6:  sprintf(str, "%%rax"); break;
				case 7:  sprintf(str, "%%rsp"); break;  // stack pointer
				case 8:  sprintf(str, "%%rbx"); break;
				case 9:  sprintf(str, "%%rbp"); break;

				default: sprintf(str, PROC->iregs->fmt, reg); break;
			}
		}
		else {
			sprintf(str, PROC->iregs->fmt, reg);
		}
      break;
    case Fregs:
      switch ( PROC->id ){
      case USPARCII: /* Ultrasparc nastily uses 2 single fregs as one double*/
                     /* Bleuchh.... so dfreg 31 is %f62 etc.. */
	sprintf(str,PROC->fregs->fmt,2*reg);
        break;
      default:
	sprintf(str,PROC->fregs->fmt,reg);
        break;
      }
      break;
    default:
      puts("Error: Bad register kind");
      exit(0);
      break;
  }
  return(str);
}

void dump_instruction_queue(void)
{
  struct inst_queue_entry *inst;
  struct offset *off;
  int handle;
  inst = HEAD;

  handle = 0;
  printf("/* \n");
  printf(" * BAGEL machine generated output.  \n");
  printf(" *   \n");
  printf(" * It is provided under the GNU pubic License V2  \n");
  printf(" * It is provided as is, and is not guaranteed fit for any purpose.\n");
  printf(" * BAGEL was written by Peter Boyle  \n");
  printf(" */  \n");

  if (PROC->id >= UNKNOWN){
    human_readable = 1;
    while ( (off= get_offset_struct (handle++))!=NULL ) {
      printf("#define %s %d \n",off->name, off->offst);
    }
    printf("\n");
    while (inst !=NULL){
      inst_print_unknown(inst->inst);
      inst = inst->next;
    }
    return;

  } else {
    if ( human_readable ){
      while ( (off= get_offset_struct (handle++))!=NULL ) {
        printf("#define %s %d \n",off->name, off->offst);
      }
      printf("\n");
      dump_regs(Fregs);
      dump_regs(Iregs);
    }
    while (inst !=NULL){
      inst_print(inst->inst);
      if ( human_readable && inst->inst->staller ){
        stall_inst_print(inst->inst->staller);
      }
      inst = inst->next;
    }
  }
  fflush(stdout);
}

void arg_print(struct instruction *i, int j)
{
  int atype=i->argtypes[j];

  if ( human_readable ) {
    switch( atype ) {
    case Fregs:
      if ( PROC->fregs->allocated[i->args[j]] ){
        printf("%s",PROC->fregs->mnemonic[i->args[j]]);
      }
      else{
        printf("%s",make_reg_name(Fregs,i->args[j]));
      }
      break;

    case Iregs:
      if ( PROC->iregs->allocated[i->args[j]] ){
        printf("%s",PROC->iregs->mnemonic[i->args[j]]);
      }
      else{
        printf("%s",make_reg_name(Iregs,i->args[j]));
      }
      break;
    case SIMM:
    case UIMM:
      printf("%s",get_offset_name(i->args[j]));break;
    case Label:
      printf("%s", make_label(i->args[j]));break;
    case MICmask:
		if (i->pipe == FMACPIPE || i->pipe == FADDPIPE || i->pipe == FMULPIPE || i->pipe == FMOVPIPE
			|| i->pipe == FSTRPIPE || i->pipe == SIMDMICPIPE || i->pipe == SIMDMICSWIZPIPE)
		{
			// The mask reg %k0 cannot be used as a write-mask.
			if (i->args[j] == 0) fprintf(stderr, "Warning: Applying MIC mask %%k0 enables ALL elements.\n");
			else printf("{%%k%d}", i->args[j]);
		}
		else {
			printf("%%k%d", i->args[j]);
		}
		break;
	case MICswiz:
	    printf("%s", mic_swizzles[i->args[j]]);
	    break;
    default:
      fprintf(stderr, "Error: Bad instruction\n");break;
    }
  }
  /* non-human readable... */
  else {
      switch( atype ) {
      case Fregs:
        printf("%s",make_reg_name(Fregs,i->args[j]));break;
      case Iregs:
        printf("%s",make_reg_name(Iregs,i->args[j]));break;
      case SIMM:
        if (is_x86()
			&& 
	    (i->pipe == IALUPIPE 
	     || (i->pipe == LOADPIPE && i->code == LOAD_ADDR) 
	     || (i->pipe == LOADPIPE && i->code == LOAD_IMM) 
	     || (i->pipe == SIMDPERMUTEPIPE && i->code == VPERMCOARSE)))
		{
			printf("$%d", get_offset(i->args[j]));   // immediate numbers
		}
        else {
			printf("%4d",  get_offset(i->args[j]));  // mem offsets
		}
        break;
      case UIMM:
        if (is_x86()
	    && (i->pipe == IALUPIPE 
		|| (i->pipe == LOADPIPE && i->code == LOAD_ADDR) 
		|| (i->pipe == LOADPIPE && i->code == LOAD_IMM) 
		|| (i->pipe == SIMDPERMUTEPIPE && i->code == VPERMCOARSE)))
		{
			printf("$%u", get_offset(i->args[j]));   // immediate numbers
		}
        else {
			printf("0x%4x",  get_offset(i->args[j]));  // mem offsets
		}
        break;
      case Label:
        printf("%s", make_label(i->args[j]));break;
      case MICmask:
		if (i->pipe == FMACPIPE || i->pipe == FADDPIPE || i->pipe == FMULPIPE || i->pipe == FMOVPIPE
			|| i->pipe == FSTRPIPE || i->pipe == SIMDMICPIPE || i->pipe == SIMDMICSWIZPIPE)
		{
			if (i->args[j] == 0) fprintf(stderr, "Warning: Applying MIC mask %%k0 enables ALL elements.\n");
			else printf("{%%k%d}", i->args[j]);  // mask reg used as write-mask
		}
		else {
			printf("%%k%d", i->args[j]);  // used as "normal" register
		}
		break;
	  case MICswiz:
	    printf("%s", mic_swizzles[i->args[j]]);
	    break;
      default:
        fprintf(stderr, "Error: Bad instruction\n");break;
      }
  }

}

void inst_print(struct instruction *i)
{
  int j,p;




  /*
   * We've taken the knowledge of effective addressing out of the arg types to
   * make dependency checking logic simpler.
   * Consequently we need to check here whether the opcode demands construction
   * of n(reg) or n,reg arguments
   */
  if ( i-> mnemonic == NULL ) {
    fprintf(stderr,"Bad operation pipe %d code %d\n",i->pipe,i->code);
    exit(-1);
  }


  if ( i-> pipe == PRAGMA ) return;

  if ( i->code == Enter_Routine && i->pipe == DIRECTIVE) {
    set_label_prefix(i->mnemonic);
  }

  switch ( i-> pipe ){
  case DIRECTIVE:
    directive_print(i); /*These are so architecture specific I give up generality*/
    break;
  case CACHPIPE:
    cache_print ( i); /*These are so architecture specific I give up generality*/
    break;
  case LOADPIPE:
    if ( i->code == LOAD_IMM ) {
	  printf("\t%s ",i->mnemonic);
	  if ( PROC->switchargs[i->pipe] ){
		 arg_print(i,1);
		 printf(",");
		 arg_print(i,0);
	  }
	  else{
		 arg_print(i,0);
		 printf(",");
		 arg_print(i,1);
	  }
	  printf("\n");
      break;
    } else {
		if ((PROC->id == KNC || PROC->id == KNC_SINGLE) && i->code == LOAD_ADDR) {
		  printf("\t%s ", i->mnemonic);
          arg_print(i, 1);
          printf(", ");
          arg_print(i, 0);
          printf("\n");
		}
		else {
		  print_eff_addr(i);
		  printf("\n");
		}
      break;
    }

  case FLODPIPE:
    if ( (i->code == CLOAD_INDEXED) || (i->code ==CSLOAD_INDEXED) ) {
		if (PROC->id == KNC || PROC->id == KNC_SINGLE) {
		  print_eff_addr(i);
		}
		else {
		  printf("\t%s ",i->mnemonic);
		  arg_print(i,0);
		  printf(",");
		  arg_print(i,1);
		  printf(",");
		  arg_print(i,2);
		}
    } else {
      print_eff_addr(i);
    }
    printf("\n");
    break;
  case FSTRPIPE: /*Generate code for effective addresses*/
    if ( ( i->code == CSAVE_INDEXED ) || ( i->code == CSSAVE_INDEXED ) ){
	  if (PROC->id == KNC || PROC->id == KNC_SINGLE) {
		print_eff_addr(i);
	  }
	  else {
		printf("\t%s ",i->mnemonic);
		arg_print(i,0);
		printf(",");
		arg_print(i,1);
		printf(",");
		arg_print(i,2);
	  }
    } else {
      print_eff_addr(i);
    }
    printf("\n");
    break;

  case STORPIPE:
    print_eff_addr(i);
    printf("\n");
    break;

  case BRCHPIPE:
    printf("\t%s ",i->mnemonic);
    if ( PROC->switchargs[i->pipe] ){
      for(j=1;j<i->narg;j++) {
       arg_print(i,j);
       printf(" , ");
      }
      arg_print(i,0);
    }
    else {
      arg_print(i,0);
      if ( i->narg > 1 )
        printf(" , ");
      for(j=1;j<i->narg;j++) {
        arg_print(i,j);
        if ( j<i->narg-1 ) printf(" , ");
      }
    }
    if ( PROC->delayslot == 1 ) /*Insert a nop in branch delay slot */
      printf("\n\t%s",PROC->iset->instructions[DIRECTIVE][NOP].mnemonic);
    printf("\n");
    break;

  case SIMDHUMMERPIPE:

    printf("\t%s ",i->mnemonic);
    fprintf(stderr,"\t\t%s %d %d\n ",i->mnemonic,i->narg, PROC->switchargs[i->pipe]);

    // Complex madd
    // make_inst(SIMDHUMMERPIPE,VMADD_MII_IR,e,a,c,e);
    // make_inst(SIMDHUMMERPIPE,VMADD_RR_RI,e,a,c,e);
    if ( ((i->code == VMADD_MII_IR) || (i->code == VMADD_II_MIR))
         && PROC->switchargs[i->pipe] ) {
	// A+iB
        // 0 == dest
	// 1 == i
        // 2 == B
	// 3 == A
	// qvfxxnpmadd RT,RA, RC, RB <= -RAi RCi + RBr / RAr RCi + RBi
	//  => qvfxxnpmadd 0,2,1,3
	//
        // qvfxmadd RT,RA,RC,RB <= RAr*RCr + Rbr / RAr RCi + RBi
	//  => qvfxmadd e,a,c,e
	arg_print(i,0); printf(" , ");
	arg_print(i,2); printf(" , "); //
#define CORRECT_QPX
#ifdef CORRECT_QPX
	arg_print(i,1); printf(" , ");
	arg_print(i,3);
#else
	arg_print(i,3); printf(" , ");
	arg_print(i,1);
#endif
    } else if ( //((i->code == VMADD_RR_RI)
	       (i->narg == 4)
              && PROC->switchargs[i->pipe] ) {

        // qvfxmadd RT,RA,RC,RB <= RAr*RCr + Rbr / RAr RCi + RBi
	//  => qvfxmadd e,a,c,e
        arg_print(i,0); printf(" , ");
	arg_print(i,1); printf(" , "); //
#ifdef CORRECT_QPX
	arg_print(i,2); printf(" , ");
	arg_print(i,3);
#else
	arg_print(i,3); printf(" , ");
	arg_print(i,2);
#endif

    } else {

      arg_print(i,0);
      if ( i->narg > 1 )
        printf(" , ");
      for(j=1;j<i->narg;j++) {
        arg_print(i,j);
        if ( j<i->narg-1 ) printf(" , ");
      }
    }
    printf("\n");
    break;


  case SIMDMICPIPE:
	//fprintf(stderr, "\tDebug: %s %d %d\n ", i->mnemonic, i->narg, PROC->switchargs[i->pipe]);

	printf("\t%s ", i->mnemonic);

	if (PROC->switchargs[i->pipe] == 0)
	{
		puts("Error: unexpected PROC->switchargs for SIMDMICPIPE/SIMDMICSWIZPIPE");
		abort();
	}

	if ( i->code == MIC_AND_MEMOP ){
	  arg_print(i, 2);
	  printf("(");
	  arg_print(i, 3);
	  printf("){1to16} ");
	  printf(", ");
	  arg_print(i, 1);
	  printf(", ");
	  arg_print(i, 0);
	  printf("\n");
	} else 	if ( (i->code >= MIC_ADD_MEMOP) && (i->code <= MIC_FMADD231_MEMOP) ){
	  arg_print(i, 2);
	  printf("(");
	  arg_print(i, 3);
	  printf(") ");
	  printf(", ");
	  arg_print(i, 1);
	  printf(", ");
	  arg_print(i, 0);
	  printf("\n");
	}else{
	  for (j = (i->narg - 1); j >= 0; j--) {
	    arg_print(i, j);
	    if (j > 0) printf(", ");
	  }
	  printf("\n");
	}
	break;
  case SIMDMICSWIZPIPE:
	printf("\t%s ", i->mnemonic);
	if (PROC->switchargs[i->pipe] == 0) {
		puts("Error: unexpected PROC->switchargs for SIMDMICPIPE/SIMDMICSWIZPIPE"); abort();
	}

	if (i->code == MIC_VPERMF32X4 ){
	  arg_print(i,0);
	  arg_print(i,1);
	  printf(", ");
	  arg_print(i,2);
	  printf(", ");
	  arg_print(i,3);
	  printf("\n");
	} else if ( (i->code == CLOAD_SPLAT)||(i->code == CSLOAD_SPLAT) ){
	  arg_print(i,2);
	  printf("(");
	  arg_print(i,3);
	  printf("), ");
	  arg_print(i,0);
	  arg_print(i,1);
	  printf("\n");
	} else if (0 && (i->code == MIC_FMADD231_SWIZ) ) {
 // PAB ... I don't understand the operand ordering here. Result is first in manual????
 // Appears that VSUB_SWIZ, ADD, MUL place destination  last and that is counter to manual
		arg_print(i, 0);  // result vector reg
		printf(", ");
		arg_print(i, 1);
		printf(", ");
		arg_print(i, 2);  // swizzled vector reg
		arg_print(i, 3);  // swizzle pattern
		printf("\n");
	} else if (i->narg == 4) {  // no mask
		arg_print(i, 2);  // swizzled vector reg
		arg_print(i, 3);  // swizzle pattern
		printf(", ");
		arg_print(i, 1);
		printf(", ");
		arg_print(i, 0);  // result vector reg
		printf("\n");
	} else if (i->narg == 5) {  // mask
		arg_print(i, 3);  // swizzled vector reg
		arg_print(i, 4);  // swizzle pattern
		printf(", ");
		arg_print(i, 2);
		printf(", ");
		arg_print(i, 0);  // result vector reg
		arg_print(i, 1);  // mask
		printf("\n");
	} else {
		fprintf(stderr, "Error: Unexpected no. of operands!\n"); exit(-1);
	}
	break;

  case FMACPIPE:
  case FADDPIPE:
  case FMULPIPE:
  case FMOVPIPE:
    if (human_readable) {
      for (j=0;j<PROC->npipes;j++)
        if ( PROC->pipes[j].issuegroups & ( 1<<i->pipe ) ) p = j ;
      if ( p == 0 ) printf("\n");
    }

	if (PROC->id == KNC || PROC->id == KNC_SINGLE) {
		printf("\t%s ",i->mnemonic);
		for (j = (i->narg - 1); j > 1; j--) {
			arg_print(i, j);
			printf(", ");
		}
		arg_print(i, 0);
		arg_print(i, 1);
	}
	else {
		printf("\t%s ",i->mnemonic);

		if ( PROC->switchargs[i->pipe] ){
		  for(j=1;j<i->narg;j++) {
			arg_print(i,j);
			printf(", ");
		  }
		  arg_print(i,0);
		}
		else {
		  arg_print(i,0);
		  if ( i->narg > 1 )
			printf(", ");
		  for(j=1;j<i->narg;j++) {
			arg_print(i,j);
			if ( j<i->narg-1 ) printf(", ");
		  }
		}
	}
    printf("\n");
    break;


  default:
    //fprintf(stderr,"\tdefault %s %d %d %d\n ",i->mnemonic,i->narg,i->pipe, PROC->switchargs[i->pipe]);
    if ( human_readable ){
      for(j=0;j<PROC->npipes;j++)
        if ( PROC->pipes[j].issuegroups & ( 1<<i->pipe ) ) p = j ;
      if ( p == 0 ) printf("\n");
    }

    printf("\t%s ",i->mnemonic);


    if ( PROC->switchargs[i->pipe] ){
      for(j=1;j<i->narg;j++) {
        arg_print(i,j);
        printf(", ");
      }
      arg_print(i,0);
    }
    else {
      arg_print(i,0);
      if ( i->narg > 1 )
        printf(", ");
      for(j=1;j<i->narg;j++) {
        arg_print(i,j);
        if ( j<i->narg-1 ) printf(", ");
      }
    }

    printf("\n");
    break;
  }
  return;
}




void stall_inst_print(struct instruction *i)
{
  int j;

  /*
   * We've taken the knowledge of effective addressing out of the arg types to
   * make dependency checking logic simpler.
   * Consequently we need to check here whether the opcode demands construction
   * of n(reg) or n,reg arguments
   */
  if ( i->pipe == DIRECTIVE || i->pipe == CACHPIPE ) return;

  printf("\t\t\t/* Last Dependency:");
  switch ( i-> pipe ){
  case DIRECTIVE:
    break;
  case CACHPIPE:
    cache_print ( i);
    break;
  case LOADPIPE:
    if ( i->code == LOAD_IMM ) {
      printf("\t%s ",i->mnemonic);
      if ( PROC->switchargs[i->pipe] ){
         arg_print(i,1);
         printf(",");
         arg_print(i,0);
      }
      else{
         arg_print(i,0);
         printf(",");
         arg_print(i,1);
      }
      break;
    } else {
      print_eff_addr(i);
      break;
    }
  case STORPIPE:
  case FLODPIPE:
  case FSTRPIPE: /*Generate code for effective addresses*/
    print_eff_addr(i);
    break;

  case BRCHPIPE:
    printf("\t%s ",i->mnemonic);
    if ( PROC->switchargs[i->pipe] ){
      for(j=1;j<i->narg;j++) {
       arg_print(i,j);
       printf(" , ");
      }
      arg_print(i,0);
    }
    else {
      arg_print(i,0);
      if ( i->narg > 1 )
        printf(" , ");
      for(j=1;j<i->narg;j++) {
        arg_print(i,j);
        if ( j<i->narg-1 ) printf(" , ");
      }
    }
    if ( PROC->delayslot == 1 ) /*Insert a nop in branch delay slot */
      printf("\n\t%s",PROC->iset->instructions[DIRECTIVE][NOP].mnemonic);
    break;

  default:

    printf("\t");
    printf("%s ",i->mnemonic);

    if ( PROC->switchargs[i->pipe] ){
      for(j=1;j<i->narg;j++) {
        arg_print(i,j);
        printf(" , ");
      }
      arg_print(i,0);
    }
    else {
      arg_print(i,0);
      if ( i->narg > 1 )
        printf(" , ");
      for(j=1;j<i->narg;j++) {
        arg_print(i,j);
        if ( j<i->narg-1 ) printf(" , ");
      }
    }
    break;
  }
  printf("*/\n");
  return;
}



void inst_print_unknown(struct instruction *i)
{
  char arg0[80];
  char arg1[80];
  char arg2[80];
  char arg3[80];
  /*
   * We've taken the knowledge of effective addressing out of the arg types to
   * make dependency checking logic simpler.
   * Consequently we need to check here whether the opcode demands construction
   * of n(reg) or n,reg arguments
   */
  if ( i->code == Enter_Routine && i->pipe == DIRECTIVE) {
    set_label_prefix(i->mnemonic);
  }

  switch ( i-> pipe ){
  case DIRECTIVE:
    directive_print_unknown(i);
    break;
  case CACHPIPE:
    arg_sprint(i,0,arg0);
    arg_sprint(i,1,arg1);
    printf("/*preload %s %s */\n",arg0,arg1);
    break;
  case LOADPIPE:
    if ( i->code == ILD_EA_SHORT ) {
      fprintf(stderr,"Dumping a short load mnemonic: %s\n",i->mnemonic);
    }
    if ( i->code == LOAD_IMM ) {
      arg_print(i,0);
      printf("=");
      arg_print(i,1);
      printf(";\n");
    }
    else{
      arg_sprint(i,0,arg0);
      arg_sprint(i,1,arg1);
      arg_sprint(i,2,arg2);
      printf(i->mnemonic,arg0,arg1,arg2);
      printf("\n");
    }
    break;
  case PRAGMA:
    printf("/*%s %d*/\n",i->mnemonic,i->args[0]);
    break;
  case FLODPIPE:
      arg_sprint(i,0,arg0);
      arg_sprint(i,1,arg1);
      arg_sprint(i,2,arg2);
      printf(i->mnemonic,arg0,arg1,arg2);
      printf("\n");
      break;
  case STORPIPE:
  case FSTRPIPE: /*Generate code for effective addresses*/
      arg_sprint(i,0,arg0);
      arg_sprint(i,1,arg1);
      arg_sprint(i,2,arg2);
      printf(i->mnemonic,arg1,arg2,arg0);
      printf("\n");
    break;

  case BRCHPIPE:
    switch( i->code ){
    case BRANCH_GT:
      printf("  if ( (long) ");
      arg_print(i,0);
      printf(" >  ");
      printf("0 ) goto ");
      arg_print(i,1);
      printf(";\n");
      break;
    case BRANCH_EQ:
      printf("  if ( ");
      arg_print(i,0);
      printf(" ==  ");
      printf("0 ) goto ");
      arg_print(i,1);
      printf(";\n");
      break;
    case BRANCH_LE:
      printf("  if ( (long) ");
      arg_print(i,0);
      printf(" <= ");
      printf("0 ) goto ");
      arg_print(i,1);
      printf(";\n");
      break;
    case BRANCH_LT:
      printf("  if (  (long) ");
      arg_print(i,0);
      printf(" < ");
      printf("0 ) goto ");
      arg_print(i,1);
      printf(";\n");
      break;
    case BRANCH_GE:
      printf("  if (  (long) ");
      arg_print(i,0);
      printf(" >= ");
      printf("0 ) goto ");
      arg_print(i,1);
      printf(";\n");
      break;
    case BRANCH:
      printf(" goto ");
      arg_print(i,0);
      printf(";\n");
      break;
    case BRANCH_RET:
      printf("\nreturn;\n");
      break;
    }
    break;

  default:

    switch(i->narg){
    case 1:
      arg_sprint(i,0,arg0);
      printf(i->mnemonic,arg0);
      printf("\n");
      break;
    case 2:
      arg_sprint(i,0,arg0);
      arg_sprint(i,1,arg1);
      printf(i->mnemonic,arg0,arg1);
      printf("\n");
      break;
    case 3:
      arg_sprint(i,0,arg0);
      arg_sprint(i,1,arg1);
      arg_sprint(i,2,arg2);
      printf(i->mnemonic,arg0,arg1,arg2);
      printf("\n");
      break;
    case 4:
      arg_sprint(i,0,arg0);
      arg_sprint(i,1,arg1);
      arg_sprint(i,2,arg2);
      arg_sprint(i,3,arg3);
      printf(i->mnemonic,arg0,arg1,arg2,arg3);
      printf("\n");
      break;
    case 0:
      printf("%s", i->mnemonic);
      printf("\n");
      break;
    default:
      printf("Unkown inst format!!!");
      break;
    }
    break;
  }
  return;
}

void dump_regs( int kind)
{
  struct regs *set;
  int num;

  test_proc();
  switch ( kind ){
    case Fregs:
      set = PROC->fregs;
      break;
    case Iregs:
      set = PROC->iregs;
      break;
    default: puts("Bad Argument: Free reg"); exit(0);
  }
  for ( num = 0 ; num < set->regfile;num++){

    if ( set->allocated[num] ){
      printf( "#define %s %s\n",set->mnemonic[num],make_reg_name(kind,num));
    }
  }
  return;
}


void print_eff_addr(struct instruction *i)
{
  int j,p;

  if ( human_readable) {
    for(j=0;j<PROC->npipes;j++)
      if ( PROC->pipes[j].issuegroups & ( 1<<i->pipe ) ) p = j ;
    if ( p == 0 ) printf("/*New cycle*/\n");
    printf("\t");
    for(j=0;j<p;j++)
      printf("   ");
  } else printf("\t");



  switch ( PROC->id ){
  case USPARCII:
  case USPARCII_SINGLE:
    printf("%s ",i->mnemonic);
    if ( PROC->switchargs[i->pipe] ){
       printf("[");
       arg_print(i,2);
       printf("+");
       arg_print(i,1);
       printf("],");
       arg_print(i,0);
    }
    else{
       arg_print(i,0);
       printf(",[");
       arg_print(i,2);
       printf("+");
       arg_print(i,1);
       printf("]");
    }
    break;

  case KNC:
  case KNC_SINGLE:
	printf("%s ", i->mnemonic);

	if (PROC->switchargs[i->pipe] == 0) {
		fprintf(stderr, "Error: unexpected switchargs.\n");
		abort();
	}

	if (i->pipe == FSTRPIPE || i->pipe == STORPIPE)
	{
		arg_print(i, 0);
		printf(", ");
	}

	if (i->narg > 1)
	{
		arg_print(i, 1);
		printf("(");
		arg_print(i, 2);
		printf(")");
		if (i->narg > 3)
		{
			arg_print(i, 3);  // the MIC mask
		}
	}

	if (i->pipe == FLODPIPE || i->pipe == LOADPIPE)
	{
		printf(", ");
		arg_print(i, 0);
	}

	break;

  default:
    printf("%s ",i->mnemonic);
    if ( PROC->switchargs[i->pipe] ){
       arg_print(i,1);
       printf("(");
       arg_print(i,2);
       printf("),");
       arg_print(i,0);
    }
    else{
       arg_print(i,0);
       printf(",\t");
       arg_print(i,1);
       printf("(");
       arg_print(i,2);
       printf(")");
    }
    break;
  }
}

void arg_sprint(struct instruction *i, int j,char *buf)
{
  int atype=i->argtypes[j];

  if ( human_readable ) {
    switch( atype ) {
    case Fregs:
      if ( PROC->fregs->allocated[i->args[j]] ){
        sprintf(buf,"%s",PROC->fregs->mnemonic[i->args[j]]);
      }
      else{
        sprintf(buf,"%s",make_reg_name(Fregs,i->args[j]));
      }
      break;

    case Iregs:
      if ( PROC->iregs->allocated[i->args[j]] ){
        sprintf(buf,"%s",PROC->iregs->mnemonic[i->args[j]]);
      }
      else{
        sprintf(buf,"%s",make_reg_name(Iregs,i->args[j]));
      }
      break;
    case SIMM:
      sprintf(buf,"%s",get_offset_name(i->args[j]));break;
    case UIMM:
      sprintf(buf,"%s",get_offset_name(i->args[j]));break;
    case Label:
      sprintf(buf, "%s", make_label(i->args[j])); break;
    default:
      fprintf(stderr, "Error: Bad instruction\n");break;
    }
  }
  else {
      switch( atype ) {
      case Fregs:
        sprintf(buf,"%s",make_reg_name(Fregs,i->args[j]));break;
      case Iregs:
        sprintf(buf,"%s",make_reg_name(Iregs,i->args[j]));break;
      case SIMM:
        sprintf(buf,"%d",get_offset(i->args[j]));break;
      case UIMM:
        sprintf(buf,"0x%x",get_offset(i->args[j]));break;
      case Label:
        sprintf(buf, "%s", make_label(i->args[j])); break;
      default:
        fprintf(stderr, "Error: Bad instruction\n");break;
      }
  }

}



/* ...END OUTPUT routines}*/
