/*
 *
 *  Copyright UKQCD Collaboration, March 2007.
 *  Written by Peter Boyle.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include "processor.h"
#include "registers.h"
extern struct processor *PROC;

void cayley_dperp( char *);

/*Options flags*/
int human = 0;
char name[80];
int dag=0;

Datum FourSpinType=Double;

char procname[80]="UNSPECIFIED";

void touch(int addr, int line);
void touch(int addr, int line)
{
      do_flush(addr,line);
      l2_touch(addr,line);
}

int main ( int argc, char *argv[])
{
  struct processor *myproc;
  int arg;
  char *c;

  name[0] = '\0';
  /*Process command line*/
  while ( ( arg = getopt(argc,argv,"Rn:P:sd")) != EOF){
    switch(arg){
    case 'R': human = 1; break;
    case 'n': if (strlen(optarg)<30) strcpy(name,optarg); break;
    case 'P': if (strlen(optarg)<30) strcpy(procname,optarg); break;
    case 'd': dag=1; break;
    case 's': FourSpinType= Single; break;
    default: fprintf(stderr,"Usage: %s -[Rd] -n routine_name\n",argv[0]); 
             fprintf(stderr,"\tR -> human readable .S format \n"); 
             fprintf(stderr,"\td -> dagger decompose\n"); 
             exit (1); break;
    }
  }

  /*Control output according to user set up args*/
  set_processor_optarg(procname);
  set_human_readable(human);
  setup_cmadds();

  /*Write the naive asm code*/
  cayley_dperp(name);

  /*Filter through an virtual out of order processor */
  schedule_for_proc();

  /*Dump the resulting code*/
  dump_instruction_queue();

  return(0);
}

void cayley_dperp( char *name)
{
  /*
   * This marks the argument registers as defined by ABI as off limits
   * to us until they are freed by "getarg()";
   */
  int dum = defargcount(1);

  /*Handles for the labels point*/
  int branchsite;
  int branchmu;

  /*-------------------------------------------------------------------------------
   * registers used 
   *-------------------------------------------------------------------------------
   */
  reg_array_2d(ChiA  ,Cregs,2,3); /*Neighbouring 2 spinor*/   // 6 
  reg_array_2d(ChiB  ,Cregs,2,3); /*Neighbouring 2 spinor*/   // 6
  reg_array_2d(Psi ,Cregs,4,3); /*Output 4-spinor - 12 regs*/ // 12 => 30
  alreg(ee,Cregs);
  alreg(eem,Cregs);
  alreg(dee,Cregs);

  offset_3d(CHIIMM,FourSpinType,4,3,2*nsimd());
  offset_2d(COEFF,Double,4,1);

  int t;

  /*
   * Integer registers
   */
  alreg(Psi_p,Iregs);   /*Pointer to current cpt in/out PSI field      */
  alreg(Chi_p,Iregs);   /*Pointer to current cpt in/out PSI field      */
  alreg(length,Iregs);  /*number of sites*/
  alreg(consts,Iregs);  /*Point to c[s],b[s],c[s] (Nsimd)*/
  alreg(consts_cur,Iregs);  /*Point to c[s],b[s],c[s] (Nsimd)*/
  alreg(Ls,Iregs);        

  alreg(sloop,Iregs);        
  alreg(args,Iregs);        
  alreg(s,Iregs);        
  alreg(s_o,Iregs);
  alreg(pref,Iregs);
  /*Useful integer immediate constants, in units of Fsize*/
  def_off( ZERO,Byte,0);
  def_off( MINUS_ONE,Byte,-1);
  def_off( CHI_IMM, FourSpinType,12*nsimd());
  def_off( PSI_IMM, FourSpinType,24*nsimd());
  def_off( NEG_PSI_IMM, FourSpinType,-24*nsimd());

  int Isize      = def_offset(PROC->I_size,Byte,"Isize");
  int Isize2      = def_offset(PROC->I_size*2,Byte,"Isize2");

  /*
   * Local variables for C++ code
   */
  int i,sp,co,j,k,nxt,ri;

  /*********************************************************************
   * Start of the "pseudo assembler proper.
   *********************************************************************/
  make_inst(DIRECTIVE,Enter_Routine,name);
  grab_stack(0);
  save_regs();

  /*
   * Define our arguments 
   */
  getarg(args); /*Pointer to arg list*/

  queue_iload(Chi_p, ZERO,args);  queue_load_addr(args,Isize,args);   // Input
  queue_iload(Psi_p, ZERO,args);  queue_load_addr(args,Isize,args);   // Output
  queue_iload(length,ZERO,args);  queue_load_addr(args,Isize,args);   
  queue_iload(Ls,    ZERO,args);  queue_load_addr(args,Isize,args);  
  queue_iload(consts,   ZERO,args); queue_load_addr(args,Isize,args);  
  
  /*Bugger 12 int regs for offsets as indexed loads only*/
  for ( i =0; i<12; i++ ) { 
    need_constant(i*2*SizeofDatum(FourSpinType)*nsimd());
  }

  /*
   * Site loop for(s=0;s<Ls;s++)
   */
  pragma(STORE_LIM,1);
  branchsite = start_loop(length);

  queue_iadd_imm(s,Chi_p,ZERO); 
  queue_iadd_imm(s_o,Psi_p,ZERO); 
  queue_iadd_imm(consts_cur,consts,ZERO);

  complex_rsplat(ee ,COEFF[0][0],consts_cur);
  complex_rsplat(eem,COEFF[1][0],consts_cur);

  // S=0 first to start up. This has coeff identity
  for(int sp=0;sp<2;sp++){
    int co=0;
      complex_load(Psi[sp][co],  CHIIMM[sp][co][0],s,FourSpinType);
      complex_load(Psi[sp+2][co],CHIIMM[sp+2][co][0],s,FourSpinType);
      co=2;
      complex_load(Psi[sp][co],  CHIIMM[sp][co][0],s,FourSpinType);
      complex_load(Psi[sp+2][co],CHIIMM[sp+2][co][0],s,FourSpinType);
      co=1;
      complex_load(Psi[sp][co],  CHIIMM[sp][co][0],s,FourSpinType);
      complex_load(Psi[sp+2][co],CHIIMM[sp+2][co][0],s,FourSpinType);
  }
  for(int sp=0;sp<2;sp++){ 
    for(int co=0;co<3;co++){
      complex_store(Psi[sp][co],  CHIIMM[sp][co][0],s_o,FourSpinType);
    }
  }
  for(int sp=0;sp<2;sp++){ 
    for(int co=0;co<3;co++){
      complex_store(Psi[sp+2][co],CHIIMM[sp+2][co][0],s_o,FourSpinType);
    }
  }
  for(int sp=0;sp<2;sp++){
    for(int co=0;co<3;co++){
      if ( dag )      complex_real_mul(ChiB[sp][co],eem,Psi[sp][co]);
      else            complex_real_mul(ChiB[sp][co],eem,Psi[sp+2][co]);
    }
  }
  
  queue_iadd_imm(sloop,Ls,MINUS_ONE);
  int branchls1 = start_loop(sloop);

  /////////////////////
  // Data pointers
  /////////////////////
  // Advance output pointer
  // Advance input pointer
  queue_iadd_imm(s,s,PSI_IMM);
  queue_iadd_imm(s_o,s_o,PSI_IMM);

  // Advance the constant pointer
  queue_iadd_imm(consts_cur,consts_cur,COEFF[2][0]);
  // Read the constants
  complex_rsplat(ee ,COEFF[0][0],consts_cur);
  complex_rsplat(eem,COEFF[1][0],consts_cur);
  
  if ( dag ) { 

    for(int co=0;co<3;co++){
      for(int sp=0;sp<2;sp++){
	complex_load(Psi[sp][co] ,  CHIIMM[sp  ][co][0],s,FourSpinType);    // P_+ Chi[s] -> Psi[s]
	complex_load(ChiA[sp][co],  CHIIMM[sp+2][co][0],s,FourSpinType); // P_- Chi[s] -> ChiA
      }
    }

    // Must load uee[s-1] and recurse for an upper inverse matrix
    for(int co=0;co<3;co++){
      for(int sp=0;sp<2;sp++){
	complex_real_madd(ChiB[sp][co],eem,Psi[sp][co]);
	complex_real_madd  (Psi[sp+2][co],ee,Psi[sp+2][co],ChiA[sp][co]);// P- recursion
      }
    }    

    for(int sp=0;sp<2;sp++){
      for(int co=0;co<3;co++){
	complex_store(Psi[sp][co],  CHIIMM[sp][co][0]  ,s_o,FourSpinType);
      }
    }  
    for(int sp=0;sp<2;sp++){
      for(int co=0;co<3;co++){
	complex_store(Psi[sp+2][co],CHIIMM[sp+2][co][0],s_o,FourSpinType);
      }
    }
  } else { 

    // Psi_-[s] = Chi_-[s]
    for(int co=0;co<3;co++){
      for(int sp=0;sp<2;sp++){
	complex_load(Psi[sp+2][co],  CHIIMM[sp+2][co][0],s,FourSpinType); // P_- Chi[s]
	complex_load(ChiA[sp][co],  CHIIMM[sp][co][0],s,FourSpinType);    // P_+ Chi[s]
      }
    }

    for(int co=0;co<3;co++){
      for(int sp=0;sp<2;sp++){
	complex_real_madd(ChiB[sp][co],eem,Psi[sp+2][co]); 
	complex_real_madd(Psi[sp][co],ee,Psi[sp][co],ChiA[sp][co]);  // P- recursion
      }
    }    
    // Accumulate the Um term for psi[Ls-1]. Last one [Ls-1] acquires the sum of all smaller s values
    // For Ls-1, we add with coefficient unity, and then ChiB contains the corrected term.
    for(int sp=0;sp<2;sp++){
      for(int co=0;co<3;co++){
	complex_store(Psi[sp][co],  CHIIMM[sp][co][0],s_o,FourSpinType);
      }
    }  
    for(int sp=0;sp<2;sp++){
      for(int co=0;co<3;co++){
	complex_store(Psi[sp+2][co],CHIIMM[sp+2][co][0],s_o,FourSpinType);
      }
    }

  }
  stop_loop(branchls1,sloop);

  // S=ls-1 first
  // We retain ChiB and Psi[Ls-1] in registers

  queue_iadd_imm(consts_cur,consts_cur,COEFF[2][0]);
  complex_rsplat(dee ,COEFF[0][0],consts_cur);
  queue_iadd_imm(consts_cur,consts_cur,COEFF[1][0]);

  for(int sp=0;sp<2;sp++){
    for(int co=0;co<3;co++){
      complex_real_mul(ChiB[sp][co],dee,ChiB[sp][co]);
    }
  }
  
  if(dag){
    for(int sp=0;sp<2;sp++){
      for(int co=0;co<3;co++){
	complex_real_mul(ChiA[sp][co],dee,Psi[sp+2][co]);
      }
    }
    for(int sp=0;sp<2;sp++){
      for(int co=0;co<3;co++){
	complex_store(ChiA[sp][co],CHIIMM[sp+2][co][0],s_o,FourSpinType);// ChiA is P-
      }
    }
    for(int sp=0;sp<2;sp++){
      for(int co=0;co<3;co++){
	complex_store(ChiB[sp][co],CHIIMM[sp+0][co][0],s_o,FourSpinType);// ChiB is P+
      }
    }
  } else { 
    for(int sp=0;sp<2;sp++){
      for(int co=0;co<3;co++){
	complex_real_mul(ChiA[sp][co],dee,Psi[sp][co]);
      }
    }
    for(int sp=0;sp<2;sp++){
      for(int co=0;co<3;co++){
	complex_store(ChiA[sp][co],CHIIMM[sp+0][co][0],s_o,FourSpinType);
      }
    }
    for(int sp=0;sp<2;sp++){
      for(int co=0;co<3;co++){
	complex_store(ChiB[sp][co],CHIIMM[sp+2][co][0],s_o,FourSpinType);
      
      }
    }
  }

  // Retain ChiB as this will enter every subsequent term, multiplied by eem.


  queue_iadd_imm(sloop,Ls,MINUS_ONE);
  int branchls2 = start_loop(sloop);

  complex_rsplat(dee ,COEFF[0][0],consts_cur);
  complex_rsplat(ee  ,COEFF[1][0],consts_cur);
  complex_rsplat(eem ,COEFF[2][0],consts_cur);
  queue_iadd_imm(consts_cur,consts_cur,COEFF[3][0]);
  
  queue_iadd_imm(s_o,s_o,NEG_PSI_IMM);
  for(int co=0;co<3;co++){
    for(int sp=0;sp<2;sp++){
      complex_load(Psi[sp][co],  CHIIMM[sp][co][0],s_o,FourSpinType);
      complex_load(Psi[sp+2][co],CHIIMM[sp+2][co][0],s_o,FourSpinType);
    }
  }

  for(int co=0;co<3;co++){
    for(int sp=0;sp<2;sp++){
      complex_real_mul(Psi[sp][co],dee,Psi[sp][co]);
      complex_real_mul(Psi[sp+2][co],dee,Psi[sp+2][co]);
    }
  }
  
  // eem terms pick up ChiA from last column always
  for(int co=0;co<3;co++){
    for(int sp=0;sp<2;sp++){
      if(dag) { 
	complex_real_madd(Psi[sp+2][co],eem,ChiA[sp][co],Psi[sp+2][co]); //P-
      } else { 
	complex_real_madd(Psi[sp][co],eem,ChiA[sp][co],Psi[sp][co]);
      }
    }
  }
  
  // ee terms use a recursion, via ChiB
  for(int co=0;co<3;co++){
    for(int sp=0;sp<2;sp++){
      if(dag){
	complex_real_madd(ChiB[sp][co],ee,ChiB[sp][co],Psi[sp][co]); //P+
      } else {
	complex_real_madd(ChiB[sp][co],ee,ChiB[sp][co],Psi[sp+2][co]);
      }
    }
  }  

  // Store
  if(dag){
    for(int sp=0;sp<2;sp++){
      for(int co=0;co<3;co++){
	complex_store(ChiB[sp][co],CHIIMM[sp+0][co][0],s_o,FourSpinType); // P+
      }
    }
    for(int sp=0;sp<2;sp++){
      for(int co=0;co<3;co++){
	complex_store(Psi[sp+2][co],CHIIMM[sp+2][co][0],s_o,FourSpinType);// P-
      }
    }
  } else {
    for(int sp=0;sp<2;sp++){
      for(int co=0;co<3;co++){
	complex_store(ChiB[sp][co],CHIIMM[sp+2][co][0],s_o,FourSpinType);
      }
    }
    for(int sp=0;sp<2;sp++){
      for(int co=0;co<3;co++){
	complex_store(Psi[sp+0][co],CHIIMM[sp+0][co][0],s_o,FourSpinType);
      }
    }
  }  
  
  stop_loop(branchls2,sloop);

  queue_imul_imm(s,Ls,PSI_IMM);
  queue_iadd(Chi_p,s,Chi_p); 
  queue_iadd(Psi_p,s,Psi_p); 

  stop_loop(branchsite,length);

  /*
   * EPILOGUE
   */
  restore_regs();
  free_stack();
  make_inst(DIRECTIVE,Exit_Routine,name);

  return;
}
  
