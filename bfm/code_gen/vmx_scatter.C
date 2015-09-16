
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


void vmx_scatter( char *);

void need_cache_line(int line);
void do_prefetch ( int ptr,int line );
void do_writehint ( int ptr,int line );

int permute = 0;
char name[80]="";
char procname[80] = "UNSPECIFIED"; /*Default processor type*/
Datum FourSpinType = Double;
Datum TwoSpinType = Double;

int main ( int argc, char *argv[])
{
  int arg;
  char *c;

  /*Process command line*/
  while ( ( arg = getopt(argc,argv,"n:P:psS")) != EOF){
    switch(arg){
    case 'n': if (strlen(optarg)<30) strcpy(name,optarg); break;
    case 'P': if (strlen(optarg)<30) strcpy(procname,optarg); break;
    case 'p': permute = 1; break;
    case 'S': TwoSpinType = Single ; break;
    case 's': FourSpinType = TwoSpinType = Single ; break;
    default: fprintf(stderr,"Usage: %s -[hs] -n routine_name -P proc\t",argv[0]); 
             exit (1); break;
    }
  }

  /*Control output according to user set up args*/
  set_processor_optarg(procname);

  /*queue the naive asm code*/
  vmx_scatter(name);
  /*Filter through an virtual out of order processor */
  schedule_for_proc();
  /*Dump the resulting code*/
  dump_instruction_queue();
  return(0);
}
/*
void vmx_scatter(
                 Float *comm_buf,
                 Float *psi,
                 integer nface,
		 int Ls
                );
*/

void vmx_scatter( char *name)
{
  int dum = defargcount(4);

  /*Integer register usage*/
  alreg(in_p,Iregs);
  alreg(out_p,Iregs);
  alreg(Ls,Iregs);
  alreg(nface,Iregs);
  alreg(s,Iregs);

  /*Floating register usage*/
  reg_array_1d(T,Cregs,12);

  def_off(ZERO_IMM ,FourSpinType,0);
  def_off(PSI_IMM  ,FourSpinType,24*nsimd());
  offset_2d(PSIIMM,FourSpinType,6,2*nsimd());
  def_off(CHI_IMM  ,TwoSpinType,12*nsimd());
  offset_2d(CHIIMM,TwoSpinType,6,2*nsimd());

  int brchf,brchs,retno;       /*Branch target handles*/
  int i;

  make_inst(DIRECTIVE,Enter_Routine,name);
  grab_stack(0);
  save_regs();

  getarg(in_p);           /*Get args*/
  getarg(out_p);           /*Get args*/
  getarg(nface);           /*Get args*/
  getarg(Ls);           /*Get args*/


  for( i = 0;i<6;i++) { 
    need_constant(i*2*SizeofDatum(FourSpinType)*nsimd());
    need_constant(i*2*SizeofDatum(TwoSpinType)*nsimd());
  }

  /*Branch to stack restore if length <1*/
  retno = get_target_label();
  check_iterations(nface,retno); 
  //  pragma(LOAD_LIM,0);
  pragma(DCBT_SPACE,1);
  pragma(DCBT_POST,0);
  /* for i=0, i < nface */
  brchf = start_loop(nface);
  queue_iadd_imm(s,Ls,ZERO_IMM);
  brchs = start_loop(s);

    for(int spinco=0;spinco<6;spinco++) {
      complex_load(T[spinco],CHIIMM[spinco][0],in_p,TwoSpinType); 
    }  
    queue_iadd_imm(in_p,in_p,CHI_IMM);

    for(int spinco=0;spinco<6;spinco++) {
      complex_store(T[spinco],PSIIMM[spinco][0],out_p,FourSpinType);
    }

    queue_iadd_imm(out_p,out_p,PSI_IMM);
  // end loop
  stop_loop(brchs,s);
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








