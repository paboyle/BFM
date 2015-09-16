/*
 *
 *  Copyright UKQCD Collaboration, March 2007.
 *  Written by Peter Boyle.
 *  This software is provided for NON-COMMERCIAL use only,
 *  and may not be redistributed without permission.
 *  It is provided as is and is not guaranteed fit for any purpose.
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


Datum FourSpinType=Double;

char procname[80]="UNSPECIFIED";

int main ( int argc, char *argv[])
{
  struct processor *myproc;
  int arg;
  char *c;

  name[0] = '\0';
  /*Process command line*/
  while ( ( arg = getopt(argc,argv,"Rn:P:s")) != EOF){
    switch(arg){
    case 'R': human = 1; break;
    case 'n': if (strlen(optarg)<30) strcpy(name,optarg); break;
    case 'P': if (strlen(optarg)<30) strcpy(procname,optarg); break;
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
  reg_array_2d(Chi ,Cregs,4,3); /*Neighbouring spinor*/
  reg_array_2d(Psi ,Cregs,4,3); /*Output 4-spinor - 12 regs*/
  alreg(c_m,Cregs);
  alreg(c_p,Cregs);
  alreg(b,Cregs);

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
  alreg(neigh ,Iregs);  /*Point to sm,sp for index s*/
  alreg(neigh_cur ,Iregs);  /*Point to sm,sp for index s*/

  alreg(Ls,Iregs);        

  alreg(sloop,Iregs);        
  alreg(args,Iregs);        
  alreg(s,Iregs);        
  alreg(s_m,Iregs);        

  /*Useful integer immediate constants, in units of Fsize*/
  def_off( ZERO_IMM,Byte,0);
  def_off( CHI_IMM, FourSpinType,12*nsimd());
  def_off( PSI_IMM, FourSpinType,24*nsimd());

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

  queue_iload(Chi_p, ZERO_IMM,args);  queue_load_addr(args,Isize,args);   // Input
  queue_iload(Psi_p, ZERO_IMM,args);  queue_load_addr(args,Isize,args);   // Output
  queue_iload(length,ZERO_IMM,args);  queue_load_addr(args,Isize,args);   
  queue_iload(Ls,    ZERO_IMM,args);  queue_load_addr(args,Isize,args);  
  queue_iload(consts,   ZERO_IMM,args); queue_load_addr(args,Isize,args);  
  queue_iload(neigh,   ZERO_IMM,args); queue_load_addr(args,Isize,args);  

  alreg(s_p,Iregs);        
  
  /*Bugger 12 int regs for offsets as indexed loads only*/
  for ( i =0; i<12; i++ ) { 
    need_constant(i*2*SizeofDatum(FourSpinType)*nsimd());
  }

  /*
   * Site loop for(s=0;s<Ls;s++)
   */
  branchsite = start_loop(length);

  /*
   * Fifth dimension loop
   */
  queue_iadd_imm(s,Chi_p,ZERO_IMM); 
  queue_iadd_imm(neigh_cur,neigh,ZERO_IMM);
  queue_iadd_imm(consts_cur,consts,ZERO_IMM);
  queue_iadd_imm(sloop,Ls,ZERO_IMM);
  int branchls = start_loop(sloop);

    /////////////////////
    // Data pointers
    /////////////////////
    // Address of Ident is "s"
    // Neighbour address for P+, P-
    queue_iload(s_m,ZERO_IMM,neigh_cur);
    queue_iload(s_p,Isize,neigh_cur);

    queue_imul_imm(s_m,s_m,PSI_IMM);
    queue_imul_imm(s_p,s_p,PSI_IMM);
  
    queue_iadd(s_m,s_m,Chi_p);
    queue_iadd(s_p,s_p,Chi_p);

    //////////////////////
    // Coefficients
    //////////////////////
    complex_rsplat(c_m,COEFF[0][0],consts_cur);// c_m P_ chi[s_m]
    complex_rsplat(b  ,COEFF[1][0],consts_cur);// b chi[s]
    complex_rsplat(c_p,COEFF[2][0],consts_cur);// c_p P_+ chi[s_p]

    for(int co=0;co<3;co++){
      for(int sp=0;sp<2;sp++){
	 complex_load(Chi[sp][co],  CHIIMM[sp][co][0],s_p,FourSpinType);
	 complex_load(Chi[sp+2][co],CHIIMM[sp+2][co][0],s_m,FourSpinType);
	 complex_load(Psi[sp][co],  CHIIMM[sp][co][0],s,FourSpinType);
	 complex_load(Psi[sp+2][co],CHIIMM[sp+2][co][0],s,FourSpinType);
      }
    }
    
    for(int co=0;co<3;co++){
      for(int sp=0;sp<4;sp++){
	complex_real_mul  (Psi[sp][co],b,Psi[sp][co]);
      }
    }    
    for(int co=0;co<3;co++){
      for(int sp=0;sp<2;sp++){
	complex_real_madd  (Psi[sp][co],c_p,Chi[sp][co]);  // P+
	complex_real_madd  (Psi[sp+2][co],c_m,Chi[sp+2][co]); // P-
      }
    }    

    pragma(STORE_LIM,1);
    for(int sp=0;sp<4;sp++){
      for(int co=0;co<3;co++){
	complex_store(Psi[sp][co],CHIIMM[sp][co][0],Psi_p,FourSpinType);
      }
    }

    queue_iadd_imm(s,s,PSI_IMM);         // Identity pointer
    queue_iadd_imm(Psi_p,Psi_p,PSI_IMM); // output pointer
    queue_iadd_imm(neigh_cur,neigh_cur,Isize2); 
    queue_iadd_imm(consts_cur,consts_cur,COEFF[3][0]);

  /*
   * TERMINATION point of the loop
   */

  stop_loop(branchls,sloop);
  queue_imul_imm(s,Ls,PSI_IMM);
  queue_iadd(Chi_p,s,Chi_p); 
  stop_loop(branchsite,length);

  /*
   * EPILOGUE
   */
  restore_regs();
  free_stack();
  make_inst(DIRECTIVE,Exit_Routine,name);

  return;
}
  
