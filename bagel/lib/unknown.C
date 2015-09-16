/*
 *
 *  Copyright UKQCD Collaboration, November 2000.
 *  Written by Peter Boyle
 *
 *  It is provided under the GNU pubic License V2
 *  It is provided as is and is not guaranteed fit for any purpose.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "processor.h"
#include "unknown.h"
#include "bagel_int.h"

/*
 * File to produce C output from assembler!
 * Bit perverse, but aids debugging... and runs on my laptop
 */
extern struct processor *PROC;
extern int frame_membase;
extern int argcount;
enum {
  USE_FLOAT,
  USE_DOUBLE
};
int Ftype ;

void arg_sprint(struct instruction *i, int j,char *buf);

void set_unknown_type_float(void );
void set_unknown_type_double(void);
void set_unknown_type_float(void)
{
  Ftype = USE_FLOAT;
}
void set_unknown_type_double(void)
{
  Ftype = USE_DOUBLE;
}


int complex_simd_vector_length  =  1;
int integer_size = BAGEL_ISIZE;
int short_integer_size = BAGEL_ISSIZE;
const char *integer_name = BAGEL_ITYPE;
const char *short_integer_name = BAGEL_ISTYPE;
void set_unknown_bgq(void )
{
  complex_simd_vector_length = 2;
}
void set_unknown_vmx(void)
{
  complex_simd_vector_length = 4;
}


struct processor *create_unknown ()
{
  struct regs *iregs=new  regs;
  struct regs *fregs=new  regs;
  struct processor *proc = new  processor;
  struct pipeline *pipes = new  pipeline;
  /*Hardwire 1 pipe*/

  struct instructionset *iset = new  instructionset;

  int i;

  iregs->regfile=32;
  fregs->regfile=64*complex_simd_vector_length;
  iregs->fmt="R%d";
  fregs->fmt="F%d";

  for(i=0;i<32;i++){
    iregs->allocated[i] = 0;
    fregs->allocated[i] = 0;
    iregs->mnemonic[i] = NULL;
    fregs->mnemonic[i] = NULL;
  }

  proc->iregs = iregs;
  proc->fregs = fregs;
  proc->pipes = pipes;
  proc->nsimd = complex_simd_vector_length;
  proc->iset = iset;
  if ( Ftype == USE_DOUBLE ) {
    proc->FP_size= sizeof(double);
  } else {
    proc->FP_size= sizeof(float);
  }
  proc->FSP_size= sizeof(float);
  proc->I_size = integer_size;
  proc->IS_size = short_integer_size;
  proc->usenops = 0;
  proc->delayslot = 0;
  proc->RegWindows = 0;
  proc->Bias = 0;

  /*Stack area - this is a hack between us and start_routine*/
  proc->StackPointer=30;
  iregs->mnemonic[30] = "StackPointer";
  iregs->allocated[30] = 1;

   /*{... ABI definition ...*/
  proc->ABI_argbase = 0;
  proc->Isave = 0;
  proc->Fsave = 0;
  proc->CacheLine=32;
  /*Actually issue width is 2, but just pretend*/
  proc->npipes=1;
  for(i=0;i<N_PIPE_TYPES;i++){
    proc->switchargs[i] = 0;
  }

  pipes[0].issuegroups = FADDMASK | FMULMASK | FMACMASK | FMOVMASK | IALUMASK
                       | LOADMASK | STORMASK | FLODMASK | CACHMASK | FSTRMASK;

  proc->latency[FMACPIPE] = 5; /*Use ppc440 model*/
  proc->latency[FADDPIPE] = 5;
  proc->latency[FMULPIPE] = 5;
  proc->latency[FMOVPIPE] = 5;
  proc->latency[FLODPIPE] = 5;
  proc->latency[FSTRPIPE] = 5;
  proc->latency[IALUPIPE] = 1;
  proc->latency[LOADPIPE] = 3;
  proc->latency[STORPIPE] = 3;

  create_unknown_instructions(proc);
  return(proc);
}

void create_unknown_instructions( struct processor *proc)
{
  int i,j;

  struct instruction *fmacs = new  instruction [nfmacs];
  struct instruction *fadds = new  instruction [nfadds];
  struct instruction *fmuls = new  instruction [nfmuls];
  struct instruction *fmovs = new  instruction [nfmovs];
  struct instruction *fsaves= new  instruction [nfsaves];
  struct instruction *loads = new  instruction [nloads];
  struct instruction *floads= new  instruction [nfloads];
  struct instruction *cache = new  instruction [ncache];
  struct instruction *adds  = new  instruction [nadds];
  struct instruction *branch = new  instruction [nbrch];
  struct instruction *directive= new  instruction [ndirective];
  struct instruction *stores = new  instruction [nstore];

  /*Define the FMACS*/
  fmacs[FMADD].mnemonic  = (char *)"  %s = %s * %s + %s ;";
  fmacs[FMSUB].mnemonic  = (char *)"  %s = %s * %s - %s ;";
  fmacs[FNMADD].mnemonic = (char *)"  %s = - ( %s * %s + %s ) ;";
  fmacs[FNMSUB].mnemonic = (char *)"  %s = - ( %s * %s - %s );";
  for (i=0;i<nfmacs;i++){
    fmacs[i].narg  = 4;
    for(j=0;j<fmacs[i].narg;j++){
      fmacs[i].argtypes[j] = Fregs;
      fmacs[i].inout[j] = inarg;
    }
    fmacs[i].inout[0] = outarg;
    fmacs[i].pipe  = FMACPIPE;
    fmacs[i].code = i;
  }

  /*Define the FADDS*/
  fadds[FADD].mnemonic = (char *)"  %s = %s + %s ; ";
  fadds[FSUB].mnemonic = (char *)"  %s = %s - %s ; ";
  for (i=0;i<nfadds;i++){
    fadds[i].narg = 3;
    for(j=0;j<fadds[i].narg;j++){
      fadds[i].argtypes[j] = Fregs;
      fadds[i].inout[j] = inarg;
    }
    fadds[i].inout[0] = outarg;
    fadds[i].pipe  = FADDPIPE;
    fadds[i].code = i;
  }

  /*Define the FMOV*/
  fmovs[FMOV].mnemonic = (char *)"  %s = %s ;  ";
  fmovs[FNEG].mnemonic = (char *)"  %s = -%s ; ";
  for (i=0;i<nfmovs;i++){
    fmovs[i].narg = 2;
    for(j=0;j<fmovs[i].narg;j++){
      fmovs[i].argtypes[j] = Fregs;
      fmovs[i].inout[j] = inarg;
    }
    fmovs[i].inout[0] = outarg;
    fmovs[i].pipe  = FMOVPIPE;
    fmovs[i].code = i;
  }

  /*Define FMUL*/
  fmuls[FMUL].mnemonic = (char *)"  %s = %s * %s  ;";
  fmuls[FMUL].narg  = 3 ;
  for(j=0;j<fmuls[FMUL].narg;j++){
    fmuls[FMUL].argtypes[j] = Fregs;
    fmuls[FMUL].inout[j] = inarg;
  }
  fmuls[FMUL].inout[0] = outarg;
  fmuls[FMUL].pipe  = FMULPIPE;
  fmuls[FMUL].code = FMUL;


  /*Define loads*/
  loads[ILD_EA].mnemonic       = (char *)"  %s = * ( (integer *) ( (%s)+(%s) ) )  ;";
  loads[ILD_EA_SHORT].mnemonic = (char *)"  %s = * ( (short_integer *) ( (%s)+(%s) ) )  ;";
  loads[LOAD_REG_EA].mnemonic  = (char *)"  %s = * ( (integer *) ( (%s)+(%s) ) ) ;";
  loads[LOAD_ADDR].mnemonic    = (char *)"  %s =   ( %s ) + ( %s )  ;";
  loads[LOAD_IMM].mnemonic     = (char *)"  ( %s ) = ( %s ) ;  ";
  for (i=0;i<nloads;i++){
    loads[i].narg = 3;
    loads[i].argtypes[0] = Iregs;
    loads[i].argtypes[1] = SIMM;
    loads[i].argtypes[2] = Iregs;
    loads[i].inout[0] = outarg;
    loads[i].inout[1] = inarg;
    loads[i].inout[2] = inarg;
    loads[i].pipe  = LOADPIPE;
    loads[i].code  = i;
  }
  loads[LOAD_IMM].narg = 2;

  /*Define STORE - really the register save and restore*/
  /*Since these must be full register width we abstract them*/
  /* differently from the floating and integer which can be a lower precision*/
  stores[SAVE_REG_EA].mnemonic = (char *)"  *((integer *)( (%s) + (%s) )) = %s ;";
  stores[IST_EA].mnemonic      = (char *)"  *((integer *)( (%s) + (%s) )) = %s ;";
  for (i=0;i<nstore;i++){
    stores[i].narg = 3;
    stores[i].argtypes[0] = Iregs;
    stores[i].argtypes[1] = SIMM;
    stores[i].argtypes[2] = Iregs;
    stores[i].inout[0] = inarg;
    stores[i].inout[1] = inarg;
    stores[i].inout[2] = inarg;
    stores[i].pipe  = STORPIPE;
    stores[i].code  = i;
  }

  /*Define cache control codes*/
  cache[PREF_IMM].mnemonic = (char *)"BOLLOCKS "; /*Pre-fetch*/
  cache[WRITE_HINT].mnemonic = (char *)"BOLLOCKS "; /*Write hint - destructive*/
  cache[COUNTER_HINT].mnemonic = (char *)"BOLLOCKS ";
  for (i=0;i<ncache;i++){
    cache[i].narg = 1;
    cache[i].argtypes[0] = Iregs;
    cache[i].inout[0]  = inarg;
    cache[i].pipe  = CACHPIPE;
    cache[i].code  = i;
  }
  cache[PREF_IMM].narg = 2;
  cache[PREF_IMM].argtypes[0] = Iregs;
  cache[PREF_IMM].argtypes[1] = Iregs;
  cache[PREF_IMM].inout[1] = inarg;

  /*Floating point loads - does the APU interrupt the FPU arithmetic stream*/
  floads[FLOAD_EA].mnemonic     = (char *)"  %s = *( (Float *) ( (%s) + (%s) ) ) ; "  ;
  floads[FSLOAD_EA].mnemonic   = (char *)"  %s = *( (float *) ( (%s) + (%s) ) ) ; "  ;
  floads[FLOAD_REG_EA].mnemonic = (char *)"  %s = *( (Float *) ( (%s) + (%s) ) ) ; "  ;
  for (i=0;i<nfloads;i++){
    floads[i].narg = 3;
    floads[i].argtypes[0] = Fregs;
    floads[i].argtypes[1] = SIMM;
    floads[i].argtypes[2] = Iregs;

    floads[i].inout[0]  = outarg;
    floads[i].inout[1]  = inarg;
    floads[i].inout[2]  = inarg;

    floads[i].pipe  = FLODPIPE;
    floads[i].code  = i;
  }

  fsaves[FSAVE_EA].mnemonic     = (char *)"  *( (Float *) ( (%s) + (%s) ) ) = %s; "  ;
  fsaves[FSSAVE_EA].mnemonic   = (char *)"  *( (float *) ( (%s) + (%s) ) ) = %s; "  ;
  fsaves[FSAVE_REG_EA].mnemonic = (char *)"  *( (Float *) ( (%s) + (%s) ) ) = %s; "  ;
  for (i=0;i<nfsaves;i++){
    fsaves[i].narg = 3;
    fsaves[i].argtypes[0] = Fregs;
    fsaves[i].argtypes[1] = SIMM;
    fsaves[i].argtypes[2] = Iregs;

    fsaves[i].inout[0]  = inarg;
    fsaves[i].inout[1]  = inarg;
    fsaves[i].inout[2]  = inarg;

    fsaves[i].pipe  = FSTRPIPE;
    fsaves[i].code  = i;
  }


  /*Integer arithmetic pipe - Should only need add and subtract */
  /* Both with and without recording */
  adds[IADD].mnemonic = (char *)"  %s = %s + %s ; "  ;
  adds[ISUB].mnemonic = (char *)"  %s = %s - %s ; "  ;
  adds[IOR].mnemonic  = (char *)"  %s = %s | %s ; "  ;
  adds[CMOVGE].mnemonic = (char *)NULL;
  adds[CMOVGT].mnemonic = (char *)NULL;
  adds[CMOVLE].mnemonic = (char *)NULL;
  adds[CMOVLT].mnemonic = (char *)NULL;
  adds[ISUB_PREDICATE].mnemonic = NULL;
  adds[IOR_PREDICATE].mnemonic = NULL;
  adds[IMUL].mnemonic     = (char *)"  %s = %s * %s ; "  ;
  adds[IMUL_IMM].mnemonic = (char *)"  %s = %s * %s ; "  ;
  adds[IADD_IMM].mnemonic = (char *)"  %s = %s + %s ; "  ;
  adds[IAND_IMM].mnemonic = (char *)"  %s = %s & %s ; "  ;
  for (i=0;i<nadds;i++){
    adds[i].narg = 3;
    for(j=0;j<adds[i].narg;j++){
      adds[i].argtypes[j]= Iregs;
      adds[i].inout[j]=inarg;
    }
    adds[i].inout[0] = outarg;
    adds[i].pipe  = IALUPIPE;
    adds[i].code  = i;
  }

  /*Multiply immediate is a special case of IALU OPS*/
  adds[IMUL_IMM].argtypes[2] = SIMM;
  adds[IADD_IMM].argtypes[2] = SIMM;
  adds[IAND_IMM].argtypes[2] = SIMM;

  /*Branch instructions*/
  branch[BRANCH_GT].mnemonic = NULL;
  branch[BRANCH_EQ].mnemonic = NULL;
  branch[BRANCH_LT].mnemonic = NULL;
  branch[BRANCH].mnemonic    = NULL;
  branch[BRANCH_GE].mnemonic = NULL;
  branch[BRANCH_LE].mnemonic = NULL;
  branch[BRANCH_CTR].mnemonic= NULL;
  branch[BRANCH_RET].mnemonic= NULL;
  for (i=0;i<nbrch;i++){
    branch[i].narg = 2;
    branch[i].argtypes[0] = Iregs;
    branch[i].argtypes[1] = Label;
    branch[i].inout[0]  = outarg;
    branch[i].inout[1]  = inarg;
    branch[i].pipe  = BRCHPIPE;
    branch[i].code  = i;
  }
  branch[BRANCH_CTR].narg = 1;
  branch[BRANCH_CTR].argtypes[0] = Label;
  branch[BRANCH_CTR].inout[0]  = inarg;
  branch[BRANCH].narg = 1;
  branch[BRANCH].argtypes[0] = Label;
  branch[BRANCH].inout[0]  = inarg;

  /*Assembler directives*/
  directive[Target].mnemonic = (char *)"";
  directive[CmovTarget].mnemonic = (char *)"";
  directive[Cache_Align].mnemonic = (char *)"";
  directive[Enter_Routine].mnemonic = (char *)"";
  directive[Exit_Routine].mnemonic = (char *)"";
  directive[Pipe_Flush].mnemonic = (char *)"";
  directive[NOP].mnemonic = (char *)"nop  ";
  directive[FNOP].mnemonic = (char *)"fnop ";
  directive[DEBUG_MSG_INT].mnemonic = (char *)"oops";
  directive[DEBUG_MSG_FLOAT].mnemonic = (char *)"oops ";
  directive[DEBUG_MSG_COMPLEX].mnemonic = (char *)"oops ";

  for (i=0;i<ndirective;i++){
    directive[i].narg = 1;
    directive[i].argtypes[0] = Label;
    directive[i].pipe  = DIRECTIVE;
    directive[i].code  = i;
  }
  directive[DEBUG_MSG_INT].argtypes[0] = Iregs;
  directive[DEBUG_MSG_FLOAT].argtypes[0] = Fregs;
  directive[DEBUG_MSG_COMPLEX].argtypes[0] = Fregs;

  /*Now set these up as the instruction set*/
  proc->iset->instructions[FMACPIPE] = fmacs;
  proc->iset->instructions[FADDPIPE] = fadds;
  proc->iset->instructions[FMULPIPE] = fmuls;
  proc->iset->instructions[CACHPIPE] = cache;
  proc->iset->instructions[IALUPIPE] = adds;
  proc->iset->instructions[LOADPIPE] = loads;
  proc->iset->instructions[FLODPIPE] = floads;
  proc->iset->instructions[FSTRPIPE] = fsaves;
  proc->iset->instructions[BRCHPIPE] = branch;
  proc->iset->instructions[DIRECTIVE] = directive;
  proc->iset->instructions[STORPIPE] = stores;
  proc->iset->instructions[FMOVPIPE] = fmovs;

}

/*THE FOLLOWING THREE ROUTINES ARE POWER PC SPECIFIC */
static int Hdecrement;
int start_loop_unknown(int counter)
{
  (void)counter;  // suppresses "unused parameter" compiler warning

  int lab = get_target_label();
  Hdecrement = def_offset(-1,Byte,"minus1");

  make_inst(DIRECTIVE,Target,lab);
  return(lab);
}
void stop_loop_unknown(int branchno,int counter)
{
  /* Use the auto decrementing branch if possible */
  queue_load_addr(counter,Hdecrement,counter);
  make_inst(BRCHPIPE,BRANCH_GT,counter,branchno);

  return;
}
void check_iterations_unknown(int cntreg,int retno)
{
  make_inst(BRCHPIPE,BRANCH_LE,cntreg,retno);                    /*Branch if cntreg <= 0*/
}
/* ...} */

void directive_print_unknown(struct instruction *i)
{
  int k;
  char s[80];
    switch ( i->code ) {
    case CmovTarget:
    case Target:
      printf("%s:\n",make_label(i->args[0]));
      break;
    case Cache_Align:
      break;
    case Enter_Routine:
      if ( Ftype == USE_DOUBLE ) {
        printf("typedef double Float;\n");
      } else {
        printf("typedef float Float;\n");
      }
      printf("#include <stdio.h>\n");
      printf("#include <stdint.h>\n");
      printf("typedef %s integer;\n",integer_name);
      printf("typedef %s short_integer;\n",short_integer_name);
      printf("\n#ifdef __cplusplus\nextern \"C\" {\n#endif\n void %s (",
	     i->mnemonic);
      for(k=0;k<argcount;k++){
        printf("integer %s",PROC->iregs->mnemonic[k]);
        if (k<argcount-1) printf(",");
      }
      printf(")\n{\n\n");

      printf ("  /*Ensure natural alignment*/");
      printf ("  int iStackArea[%d+1];\n\n",frame_membase);
      printf ("  char *StackArea=(char *)iStackArea;\n\n");
      for ( k = argcount; k< PROC->iregs->regfile ; k++ ){
        if ( k != PROC->StackPointer ) {
          if ( PROC->iregs->allocated[k] ){
            printf ("  integer %s;\n",PROC->iregs->mnemonic[k]);
	  }
          else {
            printf ("  integer %s;\n",make_reg_name(Iregs,k));
	  }
	}
      }
      printf("  integer StackPointer = (integer) StackArea ;\n");
      for ( k = 0; k<PROC->fregs->regfile ; k++ ){
        if ( PROC->fregs->allocated[k] )
          printf ("  Float %s;\n",PROC->fregs->mnemonic[k]);
        else
          printf ("  Float %s;\n",make_reg_name(Fregs,k));
      }
      /*Temporary buffer now*/
      break;
    case Exit_Routine:
      printf("\n  return;\n}\n#ifdef __cplusplus\n}\n#endif\n");
      break;
    case Pipe_Flush:
    case LS_BARRIER:
      break;
    case NOP:
      break;
    case FNOP:
      break;
    case DEBUG_MSG_INT:
      arg_sprint(i, 0,s);
      printf("  printf(\"%s %%" BAGEL_IFMT "\\n\",%s);\n",s,s);
      break;
    case DEBUG_MSG_FLOAT:
      arg_sprint(i, 0,s);
      printf("  printf(\"%s %%le\\n\",(double)%s);\n",s,s);
      break;
    case DEBUG_MSG_COMPLEX:
      arg_sprint(i, 0,s);
      printf("  printf(\"%s \");\n",s);

      for (int j=0;j<complex_simd_vector_length*2;j++){

        if ( PROC->fregs->allocated[i->args[0]+j] ){
          sprintf(s,"%s",PROC->fregs->mnemonic[i->args[0]+j]);
        } else {
	  sprintf(s,"%s",make_reg_name(Fregs,i->args[0]+j));
	}

	if ( j>0 ) {
	  printf("  printf(\"_%%16.16le\",*((double *) &%s));\n",s);
	} else {
	  printf("  printf(\"%%16.16le\",*((double *) &%s));\n",s);
	}
      }
      printf("  printf(\"\\n\");\n");
      break;
    default:
      printf("Error: Unknown assembler directive requested\n");
      exit(0);
      break;
    }
}
void queue_cmovge_unknown(int cond,int src,int dst)
{
  int skipto = get_target_label();
  make_inst(BRCHPIPE,BRANCH_LT,cond,skipto);        /*Branch if cntreg < 0*/
  make_inst(IALUPIPE,IOR,dst,src,src);              /*Copy if cond >= 0 */
  make_inst(DIRECTIVE,CmovTarget,skipto);
}

void queue_cmovle_unknown(int cond,int src,int dst)
{
  int skipto = get_target_label();
  make_inst(BRCHPIPE,BRANCH_GT,cond,skipto);        /*Branch if cntreg > 0*/
  make_inst(IALUPIPE,IOR,dst,src,src);              /*Copy if cond <= 0 */
  make_inst(DIRECTIVE,CmovTarget,skipto);
}

void queue_cmovgt_unknown(int cond,int src,int dst)
{
  int skipto = get_target_label();
  make_inst(BRCHPIPE,BRANCH_LE,cond,skipto);        /*Branch if cntreg <= 0*/
  make_inst(IALUPIPE,IOR,dst,src,src);              /*Copy if cond > 0 */
  make_inst(DIRECTIVE,CmovTarget,skipto);

}

void queue_cmovlt_unknown(int cond,int src,int dst)
{
  int skipto = get_target_label();
  make_inst(BRCHPIPE,BRANCH_GE,cond,skipto);        /*Branch if cntreg >= 0*/
  make_inst(IALUPIPE,IOR,dst,src,src);              /*Copy if cond < 0 */
  make_inst(DIRECTIVE,CmovTarget,skipto);

}


void cache_print_unknown(struct instruction *i)
{

  switch(i->code){
  case WRITE_HINT:
      break;
  case PREF_IMM:
      break;
  case COUNTER_HINT:
      break;
  default:
      printf("Error: Unknown cache hint requested\n");
      exit(0);
      break;
  }
  return;
}


void debugF(int reg)
{
  make_inst(DIRECTIVE,DEBUG_MSG_FLOAT,reg);
}
void debugC(int reg)
{
  make_inst(DIRECTIVE,DEBUG_MSG_COMPLEX,reg);
}
void debugI(int reg)
{
  make_inst(DIRECTIVE,DEBUG_MSG_INT,reg);
}

