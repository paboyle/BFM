
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
#include "usparc.h"

struct processor *create_usparcII (void)
{
  struct regs *iregs=new regs;
  struct regs *fregs=new regs;
  struct processor *proc = new processor;
  struct pipeline *pipes = new pipeline [4]; /*Hardwire 4 pipes for now!!*/
  struct instructionset *iset = new instructionset;
  char *tmp;
  int i;

  iregs->regfile=32;
  fregs->regfile=32;
  iregs->fmt="%%r%d";
  fregs->fmt="%%f%d";
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
  proc->FSP_size = 4;
  proc->I_size = 4;
  proc->usenops = 0;
  proc->delayslot = 1;
  proc->RegWindows = 1;
  proc->Bias = 0;

  /*{... ABI definition ...*/
  proc->ABI_argbase = 24;
    /*Mark global registers as allocated -> don't know compiler policy*/
  for ( i=0;i<=7;i++){
    if ( i!=1  && i != 4 && i!= 5 ){
      iregs->allocated[i] = 1;
      tmp = (char *)malloc(8);
      sprintf(tmp,"%%g%d",i);
      iregs->mnemonic[i] = tmp;
    }
  }

  iregs->allocated[14] = 1; /*Stack pointer is allocated*/
  iregs->allocated[30] = 1; /*Frame pointer is allocated*/
  iregs->allocated[31] = 1; /*Return pointer*/
  iregs->mnemonic[0] = "izero";
  iregs->mnemonic[14] = "%sp";
  iregs->mnemonic[30] = "%fp";
  iregs->mnemonic[31] = "retaddr";
  proc->StackPointer = 14;
  proc->ReturnPointer = 31;

    /*ABI defines a number of callee save registers*/
    /*
     * If I choose to not use the register windows and
     * save registers manually we omit globals (bottom 8)
     */
  proc->Isave=0xffffff00;
  proc->Fsave=0;   /*On sparc all Fregs are volatile (caller save)*/
  /*...END ABI}*/

  proc->CacheLine = 64; /*64 byte cache line*/
  proc->npipes=4;
  for(i=0;i<N_PIPE_TYPES;i++){
    proc->switchargs[i] = 0; /*Wrote this with PPC in mind so arg order is not always right*/
  }
  proc->switchargs[FMULPIPE] = 1;
  proc->switchargs[FADDPIPE] = 1;
  proc->switchargs[FMOVPIPE] = 1;
  proc->switchargs[IALUPIPE] = 1;
  proc->switchargs[LOADPIPE] = 1;
  proc->switchargs[FLODPIPE] = 1;

  pipes[1].issuegroups = FADDMASK|FMOVMASK;
  pipes[0].issuegroups = FMULMASK|FMOVMASK ;
  pipes[2].issuegroups = IALUMASK | BRCHMASK ;
  pipes[3].issuegroups = IALUMASK | BRCHMASK | LOADMASK | FSTRMASK
                       | STORMASK | FLODMASK | CACHMASK;


  proc->latency[FMACPIPE] = 1000; /*No fmac on alpha*/
  proc->latency[FADDPIPE] = 3;
  proc->latency[FMULPIPE] = 3;
  proc->latency[FMOVPIPE] = 3;
  proc->latency[FLODPIPE] = 5; /*Guess a reasonable L1/2 latency*/
  proc->latency[FSTRPIPE] = 5;
  proc->latency[IALUPIPE] = 1;
  proc->latency[LOADPIPE] = 5;
  proc->latency[STORPIPE] = 5;

  create_usparcII_instructions(proc);
  return(proc);
}

void create_usparcII_instructions( struct processor *proc)
{
  int i,j;

  struct instruction *fmacs = new instruction [nfmacs];
  struct instruction *fadds = new instruction [nfadds];
  struct instruction *fmovs = new instruction [nfmovs];
  struct instruction *fmuls = new instruction [nfmuls];
  struct instruction *loads = new instruction [nloads];
  struct instruction *floads= new instruction [nfloads];
  struct instruction *fsaves= new instruction [nfsaves];
  struct instruction *cache = new instruction [ncache];
  struct instruction *adds  = new instruction [nadds];
  struct instruction *branch = new instruction [nbrch];
  struct instruction *stores = new instruction [nstore];
  struct instruction *directive=new instruction [ndirective];

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
  fadds[FADD].mnemonic = "faddd  ";
  fadds[FSUB].mnemonic = "fsubd ";
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
  fmovs[FMOV].mnemonic = "fmovd  ";
  fmovs[FNEG].mnemonic = "fnegd  ";
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
  fmuls[FMUL].mnemonic = "fmuld  ";
  fmuls[FMUL].narg  = 3 ;
  for(j=0;j<fmuls[FMUL].narg;j++){
    fmuls[FMUL].argtypes[j] = Fregs;
    fmuls[FMUL].inout[j]   = inarg;
  }
  fmuls[FMUL].inout[0]   = outarg;
  fmuls[FMUL].pipe  = FMULPIPE;
  fmuls[FMUL].code = FMUL;


  /*Define loads*/
  loads[ILD_EA].mnemonic = "ld  ";
  loads[LOAD_REG_EA].mnemonic = "ld  ";
  loads[LOAD_ADDR].mnemonic = NULL; /*No short mnemonic on usparc bum!*/
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
  loads[LOAD_IMM].mnemonic = "mov  ";
  loads[LOAD_IMM].narg = 2;

  /*Define STORE - really the register save and restore*/
  /*Since these must be full register width we abstract them*/
  /* differently from the floating and integer which can be a lower precision*/
  stores[SAVE_REG_EA].mnemonic = "st  ";
  stores[IST_EA].mnemonic = "st  ";
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
  cache[PREF_IMM].mnemonic = "prefetch "; /*Pre-fetch*/
  cache[WRITE_HINT].mnemonic = "prefetch "; /*Write hint - different fcn
					     * codes insterted by cache_print
                                             */
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
  floads[FLOAD_EA].mnemonic = "ldd "  ;
  floads[FLOAD_REG_EA].mnemonic = "ldd "  ;
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

  fsaves[FSAVE_EA].mnemonic = "std " ;
  fsaves[FSAVE_REG_EA].mnemonic = "std " ;
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
  adds[IADD].mnemonic = "add  "  ;
  adds[ISUB].mnemonic = "sub  "  ;
  adds[IMUL].mnemonic = "umul  "  ;
  adds[IMUL_IMM].mnemonic = "umul  "  ;
  adds[IOR].mnemonic  = "or   "  ;
  adds[CMOVGE].mnemonic ="movrgez ";
  adds[CMOVGT].mnemonic ="movrgz ";
  adds[CMOVLE].mnemonic ="movrlez ";
  adds[CMOVLT].mnemonic ="movrlz ";
  adds[ISUB_PREDICATE].mnemonic ="subcc ";
  adds[IOR_PREDICATE].mnemonic = "orcc ";
  adds[IADD_IMM].mnemonic = "add  ";
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
  branch[BRANCH_GT].mnemonic = "bg ";
  branch[BRANCH_EQ].mnemonic = "be ";
  branch[BRANCH_LT].mnemonic = "ble ";
  branch[BRANCH].mnemonic    = "ba ";
  branch[BRANCH_GE].mnemonic = "bge ";
  branch[BRANCH_LE].mnemonic = "ble ";
  branch[BRANCH_CTR].mnemonic= NULL;
  branch[BRANCH_RET].mnemonic= NULL;
  for (i=0;i<nbrch;i++){
    branch[i].narg = 1;
    branch[i].argtypes[0] = Label;
    branch[i].inout[0] = inarg;
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
  directive[POP_WINDOW].mnemonic = "";
  directive[PUSH_WINDOW].mnemonic = "";
  directive[POP_WINDOW].argtypes[0] = SIMM;  /*Pass framesize as immediate offset*/
  directive[PUSH_WINDOW].argtypes[0]  = SIMM;


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



/*THE FOLLOWING THREE ROUTINES ARE PROCESSOR SPECIFIC */
static int Hdecrement;
int start_loop_usparcII(int counter)
{
  (void)counter;  // suppresses "unused parameter" compiler warning
  int lab = get_target_label();
  Hdecrement = def_offset(-1,Byte,"minus1");
  make_inst(BRCHPIPE,BRANCH,lab);
  make_inst(DIRECTIVE,Cache_Align,lab);
  make_inst(DIRECTIVE,Target,lab);
  return(lab);
}

void stop_loop_usparcII(int branchno,int counter)
{
  queue_load_addr(counter,Hdecrement,counter);
  queue_barrier();
  make_inst(IALUPIPE,IOR_PREDICATE,counter,counter,counter);
  make_inst(BRCHPIPE,BRANCH_GT,branchno);
}
void check_iterations_usparcII(int cntreg,int retno)
{
  make_inst(IALUPIPE,IOR_PREDICATE,cntreg,cntreg,cntreg);
  make_inst(BRCHPIPE,BRANCH_LE,retno);      /*Branch if cntreg <= 0*/
}
/* ...} */

void directive_print_usparcII(struct instruction *i)
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
      printf("%s:\n\n",i->mnemonic);

      break;
    case Exit_Routine:
      printf("\treturn %%i7+8\n\tnop\n");
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
    case POP_WINDOW:
      break;
    case PUSH_WINDOW:
      printf("\tsave %%sp,-176,%%sp\n",get_offset(i->args[0]));
      break;
    default:
      printf("Error: Unknown assembler directive requested\n");
      exit(0);
      break;
    }
}

void queue_cmovge_usparcII(int cond,int src,int dst)
{
  make_inst(IALUPIPE,CMOVGE,dst,cond,src);
}

void queue_cmovgt_usparcII(int cond,int src,int dst)
{
  make_inst(IALUPIPE,CMOVGT,dst,cond,src);
}

void queue_cmovle_usparcII(int cond,int src,int dst)
{
  make_inst(IALUPIPE,CMOVLE,dst,cond,src);
}

void queue_cmovlt_usparcII(int cond,int src,int dst)
{
  make_inst(IALUPIPE,CMOVLT,dst,cond,src);
}

void cache_print_usparcII(struct instruction *i)
{

  switch(i->code){
  case WRITE_HINT:
      printf("\t%s\t[",i->mnemonic);
      arg_print(i,0);
      printf("], 3\n");
      break;
  case PREF_IMM:
      printf("\t%s\t[",i->mnemonic);
      arg_print(i,1);
      printf(" + ");
      arg_print(i,0);
      printf("],1\n");
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



