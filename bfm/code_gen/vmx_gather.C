/*
 *
 *  Copyright UKQCD Collaboration, Feb 2002.
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

int swprefetch = 0;

void vmx_gather( char *);

void need_cache_line(int line);
void do_prefetch ( int ptr,int line );

int permute = 0;
char name[80]="";
char procname[80] = "UNSPECIFIED"; /*Default processor type*/
Datum FourSpinType = Double;

int main ( int argc, char *argv[])
{
  int arg;
  char *c;

  /*Process command line*/
  while ( ( arg = getopt(argc,argv,"n:P:hps")) != EOF){
    switch(arg){
    case 'n': if (strlen(optarg)<30) strcpy(name,optarg); break;
    case 'P': if (strlen(optarg)<30) strcpy(procname,optarg); break;
    case 'p': permute = 1; break;
    case 's': FourSpinType = Single ; break;
    default: fprintf(stderr,"Usage: %s -[hs] -n routine_name -P proc\t",argv[0]); 
             exit (1); break;
    }
  }

  /*Control output according to user set up args*/
  set_processor_optarg(procname);

  /*queue the naive asm code*/
  vmx_gather(name);
  /*Filter through an virtual out of order processor */
  //  schedule_for_proc();
  /*Dump the resulting code*/
  dump_instruction_queue();
  return(0);
}
/*
void vmx_gather(integer *sites,
                 Float *psi
                 integer comm_off, 
                 Float *comm_buf,
                 integer nface,
                );
*/

void vmx_gather( char *name)
{
  int dum = defargcount(6);

  /*Integer register usage*/
  alreg(sites,Iregs);
  alreg(comm_offset,Iregs);
  alreg(nface,Iregs);
  alreg(psi_p,Iregs);
  alreg(comm_p,Iregs);
  alreg(Ls,Iregs);

  /*Floating register usage*/
  reg_array_1d(T,Cregs,12);
  def_off(ZERO_IMM  ,FourSpinType,0);
  def_off(PSI_IMM  ,FourSpinType,24*nsimd());
  offset_2d(PSIIMM,FourSpinType,12,2*nsimd());

  int Isize      = def_offset(PROC->I_size,Byte,"Isize");

  int brchf,brchs,retno;       /*Branch target handles*/
  int i;

  make_inst(DIRECTIVE,Enter_Routine,name);
  grab_stack(0);
  save_regs();

  getarg(sites);           /*Get args*/
  getarg(psi_p);           /*Get args*/
  getarg(comm_offset);     /*Get args*/
  getarg(comm_p);           /*Get args*/
  getarg(nface);           /*Get args*/
  getarg(Ls);           /*Get args*/

  alreg(in_p,Iregs);
  alreg(out_p,Iregs);
  alreg(s,Iregs);
  alreg(permreg,Cregs);
  alreg(pref,Iregs);
  for( i = 0;i<12;i++) { 
    need_constant(i*2*SizeofDatum(FourSpinType)*nsimd());
  }
  if (nsimd() == 1) permute = 0;

  /*Branch to stack restore if length <1*/
  pragma(DCBT_SPACE,0);
  retno = get_target_label();
  check_iterations(nface,retno); 

  if ( permute ) 
    complex_simd_init(permreg);

  queue_imul_imm(out_p,comm_offset,PSI_IMM);
  queue_imul(out_p,out_p,Ls);
  queue_iadd(out_p,out_p,comm_p);

  pragma(LOAD_LIM,1);
  /* for i=0, i < nface */
  brchf = start_loop(nface);
  queue_iadd_imm(s,Ls,ZERO_IMM);
  // Gather address calculation
  queue_iload(in_p,ZERO_IMM,sites);
  queue_imul(in_p,in_p,Ls);
  queue_imul_imm(in_p,in_p,PSI_IMM);
  queue_iadd(in_p,in_p,psi_p);

  brchs = start_loop(s);

  int sc_idx[12]={0,2,4,6,8,10,1,3,5,7,9,11};
  for(int sc=0;sc<12;sc++) { 
    int spinco=sc;
    complex_load(T[spinco],PSIIMM[spinco][0],in_p,FourSpinType); 
  }
  for(int sc=0;sc<12;sc++) { 
    int spinco=sc;
    if ( permute ) { 
      complex_simd_permute(3,T[spinco],T[spinco]);
    }
  }  
  for(int sc=0;sc<12;sc++) { 
    int spinco=sc;
    complex_store(T[spinco],PSIIMM[spinco][0],out_p,FourSpinType);
  }

  queue_iadd_imm(out_p,out_p,PSI_IMM);
  queue_iadd_imm(in_p,in_p,PSI_IMM);

  // end loop
  stop_loop(brchs,s);
  queue_load_addr(sites,Isize,sites);
  stop_loop(brchf,nface);

  // where to jump once loop is over
  make_inst(DIRECTIVE,Target,retno);

  // Cleanup
  restore_regs();
  free_stack();

  // finish
  make_inst(DIRECTIVE,Exit_Routine,name);

  return;
}








