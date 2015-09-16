/*
 *
 *  Copyright Peter Boyle and Glasgow University 2000.
 *
 *  It is provided under the GNU pubic License V2
 *  It is provided as is and is not guaranteed fit for any purpose.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "processor.h"
#include "alpha.h"

struct processor *create_alpha164 (void)
{
  struct regs *iregs=new struct regs;
  struct regs *fregs=new struct regs;
  struct processor *proc = new struct processor;
  struct pipeline *pipes = new struct pipeline[4]; /*Hardwire 4 pipes for now!!*/
  struct instructionset *iset = new struct instructionset;

  int i;

  iregs->regfile=32;
  fregs->regfile=32;
  iregs->fmt="$%d";
  fregs->fmt="$f%d";
  for(i=0;i<32;i++){
    iregs->allocated[i] = 0;
    fregs->allocated[i] = 0;
    iregs->mnemonic[i] = NULL;
    fregs->mnemonic[i] = NULL;
  }


  proc->iregs = iregs;
  proc->fregs = fregs;
  proc->pipes = pipes;
  proc->iset = iset;
  proc->FP_size = 8;
  proc->I_size = 8;
  proc->usenops = 0;
  proc->delayslot = 0;
  proc->RegWindows = 0;
  proc->Bias = 0;

  /*{... ABI definition ...*/
  proc->ABI_argbase = 16;
    /*Allocate dedicated registers*/
  iregs->allocated[26] = 1;
  iregs->allocated[28] = 1;
  iregs->allocated[29] = 1;
  iregs->allocated[30] = 1;
  iregs->allocated[31] = 1;
  fregs->allocated[31] = 1;
  iregs->mnemonic[26] = "$ra";
  iregs->mnemonic[28] = "$at";
  iregs->mnemonic[29] = "$gp";
  iregs->mnemonic[30] = "$sp";
  iregs->mnemonic[31] = "zero";
  fregs->mnemonic[31] = "fzero";
  proc->StackPointer = 30;
  proc->ReturnPointer = 26;

    /*ABI defines a number of callee save registers*/
  proc->Isave=0;
  for ( i=9; i<=15 ; i ++){
    proc->Isave= proc->Isave | (1<<i);
  }
  proc->Isave= proc->Isave | (1<<26);
  proc->Isave= proc->Isave | (1<<30);
  proc->Fsave=0;
  for ( i=2; i<=9 ; i ++){
    proc->Fsave= proc->Fsave | (1<<i);
  }
  /*...END ABI}*/

  proc->CacheLine = 32; /*32 byte cache line*/
  proc->npipes=4;
  for(i=0;i<N_PIPE_TYPES;i++){
    proc->switchargs[i] = 0; /*Wrote this with PPC in mind so arg order is not always right*/
  }
  proc->switchargs[FMULPIPE] = 1;
  proc->switchargs[FADDPIPE] = 1;
  proc->switchargs[FMOVPIPE] = 1;
  proc->switchargs[IALUPIPE] = 1;

  pipes[1].issuegroups = FADDMASK|FMOVMASK;
  pipes[0].issuegroups = FMULMASK|FMOVMASK ;
  pipes[2].issuegroups = IALUMASK | BRCHMASK | LOADMASK
                       | FLODMASK | CACHMASK; /*Only one store per clock*/
  pipes[3].issuegroups = IALUMASK | BRCHMASK | LOADMASK  | FSTRMASK
                       | STORMASK | FLODMASK | CACHMASK;


  proc->latency[FMACPIPE] = 1000; /*No fmac on alpha*/
  proc->latency[FADDPIPE] = 4;
  proc->latency[FMULPIPE] = 4;
  proc->latency[FMOVPIPE] = 4;
  proc->latency[FLODPIPE] = 5; /*Guess a reasonable L1/2 latency*/
  proc->latency[FSTRPIPE] = 5;
  proc->latency[IALUPIPE] = 1;
  proc->latency[LOADPIPE] = 5;
  proc->latency[STORPIPE] = 5;

  create_alpha164_instructions(proc);
  return(proc);
}

void create_alpha164_instructions( struct processor *proc)
{
  int i,j;

  struct instruction *fmacs = new struct instruction [nfmacs];
  struct instruction *fadds = new struct instruction [nfadds];
  struct instruction *fmovs = new struct instruction [nfmovs];
  struct instruction *fmuls = new struct instruction [nfmuls];
  struct instruction *loads = new struct instruction [nloads];
  struct instruction *floads= new struct instruction [nfloads];
  struct instruction *fsaves= new struct instruction [nfsaves];
  struct instruction *cache = new struct instruction [ncache];
  struct instruction *adds  = new struct instruction [nadds];
  struct instruction *branch = new struct instruction [nbrch];
 struct instruction *directive = new struct instruction [ndirective];
  struct instruction *stores = new struct instruction [nstore];

  /*Define the FMACS*/
  fmacs[FMADD].mnemonic = NULL;
  fmacs[FMSUB].mnemonic = NULL;
  fmacs[FNMADD].mnemonic =NULL;
  fmacs[FNMSUB].mnemonic =NULL;
  for (i=0;i<nfmacs;i++){
    fmacs[i].narg  = 4;
    for(j=0;j<fmacs[i].narg;j++){
      fmacs[i].argtypes[j] = Fregs;
      fmacs[i].inout[j]   = inarg;
    }
    fmacs[i].inout[0]   = outarg;
    fmacs[i].pipe  = FMACPIPE;
    fmacs[i].code = i;
  }

  /*Define the FADDS*/
  fadds[FADD].mnemonic = "addt  ";
  fadds[FSUB].mnemonic = "subt  ";
  for (i=0;i<nfadds;i++){
    fadds[i].narg = 3;
    for(j=0;j<fadds[i].narg;j++){
      fadds[i].argtypes[j] = Fregs;
      fadds[i].inout[j]   = inarg;
    }
    fadds[i].inout[0]   = outarg;
    fadds[i].pipe  = FADDPIPE;
    fadds[i].code = i;
  }

  /*Define the FMOVS*/
  fmovs[FMOV].mnemonic = "fmov  ";
  fmovs[FNEG].mnemonic = "fneg  ";
  for (i=0;i<nfmovs;i++){
    fmovs[i].narg = 2;
    for(j=0;j<fmovs[i].narg;j++){
      fmovs[i].argtypes[j] = Fregs;
    }
    fmovs[i].inout[0]   = outarg;
    fmovs[i].inout[1]   = inarg;
    fmovs[i].pipe  = FMOVPIPE;
    fmovs[i].code = i;
  }

  /*Define FMUL*/
  fmuls[FMUL].mnemonic = "mult  ";
  fmuls[FMUL].narg  = 3 ;
  for(j=0;j<fmuls[FMUL].narg;j++){
    fmuls[FMUL].argtypes[j] = Fregs;
    fmuls[FMUL].inout[j]   = inarg;
  }
  fmuls[FMUL].inout[0]   = outarg;
  fmuls[FMUL].pipe  = FMULPIPE;
  fmuls[FMUL].code = FMUL;


  /*Define loads*/
  loads[ILD_EA].mnemonic = "ldq  ";
  loads[LOAD_ADDR].mnemonic = "lda  ";
  loads[LOAD_REG_EA].mnemonic = "ldq  ";
  for (i=0;i<nloads;i++){
    loads[i].narg = 3;
    loads[i].argtypes[0] = Iregs;
    loads[i].argtypes[1] = SIMM;
    loads[i].argtypes[2] = Iregs;
    loads[i].inout[0]  = outarg;
    loads[i].inout[1]  = inarg;
    loads[i].inout[2]  = inarg;
    loads[i].pipe  = LOADPIPE;
    loads[i].code  = i;
  }
  loads[LOAD_IMM].mnemonic = "ldq  ";
  loads[LOAD_IMM].narg = 2;

  /*Define STORE - really the register save and restore*/
  /*Since these must be full register width we abstract them*/
  /* differently from the floating and integer which can be a lower precision*/
  stores[SAVE_REG_EA].mnemonic = "stq  ";
  stores[IST_EA].mnemonic = "stq  ";
  for (i=0;i<nstore;i++){
    stores[i].narg = 3;
    stores[i].argtypes[0] = Iregs;
    stores[i].argtypes[1] = SIMM;
    stores[i].argtypes[2] = Iregs;
    stores[i].inout[0]  = inarg;
    stores[i].inout[1]  = inarg;
    stores[i].inout[2]  = inarg;
    stores[i].pipe  = STORPIPE;
    stores[i].code  = i;
  }

  /*Define cache control codes*/
  cache[PREF_IMM].mnemonic = "ldl $31, "; /*Pre-fetch*/
  cache[WRITE_HINT].mnemonic = "ldl $31, "; /*Write hint - destructive*/
  cache[COUNTER_HINT].mnemonic = NULL;
  for (i=0;i<ncache;i++){
    cache[i].narg = 1;
    cache[i].argtypes[0] = Iregs;
    cache[i].inout[0]  = inarg;
    cache[i].pipe  = CACHPIPE;
    cache[i].code  = i;
  }
  cache[PREF_IMM].narg = 2;
  cache[PREF_IMM].argtypes[0] = SIMM;
  cache[PREF_IMM].argtypes[1] = Iregs;

  /*Floating point loads - does the APU interrupt the FPU arithmetic stream*/
  floads[FLOAD_EA].mnemonic = "ldt "  ;
  floads[FLOAD_REG_EA].mnemonic = "ldt "  ;
  for (i=0;i<nfloads;i++){
    floads[i].narg = 3;
    floads[i].argtypes[0] = Fregs;
    floads[i].argtypes[1] = SIMM;
    floads[i].argtypes[2] = Iregs;
    floads[i].inout[0] = outarg;
    floads[i].inout[1] = inarg;
    floads[i].inout[2] = inarg;
    floads[i].pipe  = FLODPIPE;
    floads[i].code  = i;
  }

  fsaves[FSAVE_EA].mnemonic = "stt " ;
  fsaves[FSAVE_REG_EA].mnemonic = "stt " ;
  for (i=0;i<nfsaves;i++){
    fsaves[i].narg = 3;
    fsaves[i].argtypes[0] = Fregs;
    fsaves[i].argtypes[1] = SIMM;
    fsaves[i].argtypes[2] = Iregs;
    fsaves[i].inout[0] = inarg;
    fsaves[i].inout[1] = inarg;
    fsaves[i].inout[2] = inarg;
    fsaves[i].pipe  = FSTRPIPE;
    fsaves[i].code  = i;
  }


  /*Integer arithmetic pipe - Should only need add and subtract */
  /* Both with and without recording */
  adds[IADD].mnemonic = "addq  "  ;
  adds[ISUB].mnemonic = "subq  "  ;
  adds[IMUL].mnemonic = "mulq  "  ;
  adds[IMUL_IMM].mnemonic = "mulq  "  ;
  adds[IOR].mnemonic  = "or   "  ;
  adds[CMOVGE].mnemonic ="cmovge ";
  adds[CMOVGT].mnemonic ="cmovgt ";
  adds[CMOVLE].mnemonic ="cmovle ";
  adds[CMOVLT].mnemonic ="cmovlt ";
  adds[ISUB_PREDICATE].mnemonic = NULL;
  adds[IOR_PREDICATE].mnemonic = NULL;
  adds[IADD_IMM].mnemonic = "addq ";
  for (i=0;i<nadds;i++){
    adds[i].narg = 3;
    for(j=0;j<adds[i].narg;j++){
      adds[i].argtypes[j] = Iregs;
      adds[i].inout[j]  = inarg;
    }
    adds[i].inout[0]  = outarg;
    adds[i].pipe  = IALUPIPE;
    adds[i].code  = i;
  }
  adds[IMUL_IMM].argtypes[2] = SIMM;
  adds[IADD_IMM].argtypes[2] = SIMM;

  /*Branch instructions*/
  branch[BRANCH_GT].mnemonic = "bgt ";
  branch[BRANCH_EQ].mnemonic = "beq ";
  branch[BRANCH_LT].mnemonic = "blt ";
  branch[BRANCH].mnemonic    = NULL;
  branch[BRANCH_GE].mnemonic = "bge ";
  branch[BRANCH_LE].mnemonic = "ble ";
  branch[BRANCH_CTR].mnemonic= "bdnz   ";
  branch[BRANCH_RET].mnemonic= "blr    ";
  for (i=0;i<nbrch;i++){
    branch[i].narg = 2;
    branch[i].argtypes[0] = Iregs;
    branch[i].argtypes[1] = Label;
    branch[i].inout[0] = outarg;
    branch[i].inout[1] = inarg;
    branch[i].pipe  = BRCHPIPE;
    branch[i].code  = i;
  }

  /*Assembler directives*/
  directive[Target].mnemonic = "";
  directive[Cache_Align].mnemonic = "";
  directive[Enter_Routine].mnemonic = "";
  directive[Exit_Routine].mnemonic = "";
  directive[Pipe_Flush].mnemonic = "";
  directive[NOP].mnemonic = "nop  ";
  directive[FNOP].mnemonic = "fnop ";

  for (i=0;i<ndirective;i++){
    directive[i].narg = 1;
    directive[i].argtypes[0] = Label;
    directive[i].inout[0]  = inarg;
    directive[i].pipe  = DIRECTIVE;
    directive[i].code  = i;
  }

  /*Now set these up as the instruction set*/
  proc->iset->instructions[FMACPIPE] = fmacs;
  proc->iset->instructions[FADDPIPE] = fadds;
  proc->iset->instructions[FMULPIPE] = fmuls;
  proc->iset->instructions[CACHPIPE] = cache;
  proc->iset->instructions[IALUPIPE] = adds;
  proc->iset->instructions[LOADPIPE] = loads;
  proc->iset->instructions[FLODPIPE] = floads;
  proc->iset->instructions[BRCHPIPE] = branch;
  proc->iset->instructions[DIRECTIVE] = directive;
  proc->iset->instructions[STORPIPE] = stores;
  proc->iset->instructions[FMOVPIPE] = fmovs;
  proc->iset->instructions[FSTRPIPE] = fsaves;
}



/*THE FOLLOWING THREE ROUTINES ARE POWER PC SPECIFIC */
static int Hdecrement;
int start_loop_alpha164(int counter)
{
  (void)counter;  // suppresses "unused parameter" compiler warning
  int lab = get_target_label();
  Hdecrement = def_offset(-1,Byte,"minus1");
  make_inst(DIRECTIVE,Cache_Align,lab);
  make_inst(DIRECTIVE,Target,lab);
  return(lab);
}

void stop_loop_alpha164(int branchno,int counter)
{
  queue_load_addr(counter,Hdecrement,counter);
  make_inst(BRCHPIPE,BRANCH_GT,counter,branchno);
}
void check_iterations_alpha164(int cntreg,int retno)
{
  make_inst(BRCHPIPE,BRANCH_LE,cntreg,retno);      /*Branch if cntreg <= 0*/
}
/* ...} */

void directive_print_alpha164(struct instruction *i)
{
    switch ( i->code ) {
    case Target:
      printf("%s:\n",make_label(i->args[0]));
      break;
    case Cache_Align:
      printf(".align 4\n");
      break;
    case Enter_Routine:
      printf(".globl %s\n",i->mnemonic);
      printf(".ent %s\n",i->mnemonic);
      printf("%s:\n\n",i->mnemonic);
      printf(".set noreorder\n\n");
      printf(".set noat\n\n");
      break;
    case Exit_Routine:
      printf("\tret $31, ($26), 1\n");
      printf(".end %s\n",i->mnemonic);
      break;
    case Pipe_Flush:
      printf("/* Instruction order barrier inserted */\n");
      break;
    case NOP:
      printf("   nop\n");
      break;
    case FNOP:
      printf("   fnop \n");
      break;
    default:
      printf("Error: Unknown assembler directive requested\n");
      exit(0);
      break;
    }
}

void queue_cmovge_alpha164(int cond,int src,int dst)
{
  make_inst(IALUPIPE,CMOVGE,dst,cond,src);
}

void cache_print_alpha164(struct instruction *i)
{

  switch(i->code){
  case WRITE_HINT:
      printf("\t%s\t0(",i->mnemonic);
      arg_print(i,0);
      printf(")\n");
      break;
  case PREF_IMM:
      printf("\t%s\t",i->mnemonic);
      arg_print(i,0);
      printf("(");
      arg_print(i,1);
      printf(")\n");
      break;
  case COUNTER_HINT:
      printf("Counter SPR not present on alpha\n");
      exit(0);
      break;
  default:
      printf("Error: Unknown cache hint requested\n");
      exit(0);
      break;
  }
  return;
}







