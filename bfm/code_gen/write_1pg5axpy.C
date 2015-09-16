
/*
 *
 *  Copyright UKQCD Collaboration, November 2000.
 *  Written by Peter Boyle.
 *  This software is provided for NON-COMMERCIAL use only,
 *  and may not be redistributed without permission.
 *  It is provided as is and is not guaranteed fit for any purpose.
 *
 */

extern "C" {
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
}

#include "processor.h"
#include "registers.h"

void qcdoc_1pG5axpy( char *);

char name[80]="";

char procname[80] = "UNSPECIFIED"; /*Default processor type*/

Datum FourSpinType = Double;

int main ( int argc, char *argv[])
{
  int arg;
  char *c;

  /*Process command line*/
  while ( ( arg = getopt(argc,argv,"n:P:s")) != EOF){
    switch(arg){
    case 'n': if (strlen(optarg)<70) strcpy(name,optarg); break;
    case 'P': if (strlen(optarg)<70) strcpy(procname,optarg); break;
    case 's': FourSpinType = Single; break;
    default: fprintf(stderr,"Usage: %s -[h] -n routine_name -P proc\t",argv[0]); 
             exit (1); break;
    }
  }

  /*Control output according to user set up args*/
  set_processor_optarg(procname);
  setup_cmadds();

  /*queue the naive asm code*/
  qcdoc_1pG5axpy(name);

  /*Filter through an virtual out of order processor */
  schedule_for_proc();
  /*Dump the resulting code*/
  dump_instruction_queue();
  return(0);

}

  /*
   * 1pG5vaxpy ( vec3, a, vec1, vec2 , nvec )
   * fpoint a
   * fpoint vec1[nvec][4][Ncol][rei] 
   * fpoint vec2[nvec][4][Ncol][rei] 
   * fpoint vec2[nvec][4][Ncol][rei] 
   * asmint nvec
   *
   * for(i=0;i<nvec;i++)
   * for(c=0;c<Ncol;c++)
   * for(r=0;r<rei;r++)
   *   vec3[i][0][c][r] = A * vec1[i][0][c][r] + vec2[i][0][c][r];
   *   vec3[i][1][c][r] = A * vec1[i][1][c][r] + vec2[i][1][c][r];
   *
   */

void qcdoc_1pG5axpy( char *name)
{
  int dum = defargcount(5);

  /*Integer register usage*/
  alreg(vec1ptr,Iregs);
  alreg(vec2ptr,Iregs);
  alreg(vec1pre,Iregs);
  alreg(vec2pre,Iregs);
  alreg(counter,Iregs);
  alreg(Aptr,Iregs);
  alreg(outptr,Iregs);

  /*Floating register usage*/
  alreg(A,Fregs);
  reg_array_3d(vec1,Fregs,2,3,2);
  reg_array_3d(vec2,Fregs,2,3,2);
  reg_array_2d(vec3,Fregs,3,2);

  def_off(ZERO,FourSpinType,0);
  def_off(VEC_ATOM,FourSpinType,12);
  def_off(VEC_SKIP,FourSpinType,24);
  def_off(VEC_SKIP_TWO,FourSpinType,48);
  offset_3d(VEC_IMM,FourSpinType,4,3,2);

  struct stream *PreVec1;
  struct stream *PreVec2;

  int brchno,retno; /*Branch target handles*/
  int co,rei;

  make_inst(DIRECTIVE,Enter_Routine,name);
  grab_stack(0);
  save_regs();

  getarg(vec2ptr);           /*Get args*/
  getarg(Aptr);
  getarg(vec1ptr);           /*Get args*/
  getarg(counter);

  queue_fload(A,ZERO,Aptr); /*Find the complex scale factor*/

  PreVec1= create_stream(VEC_ATOM,vec1ptr ,counter,STREAM_IN ,STRIDED,VEC_SKIP);
  PreVec2= create_stream(VEC_ATOM,vec2ptr ,counter,STREAM_IN ,STRIDED,VEC_SKIP);

  /*Branch to stack restore if length <1*/
  retno = get_target_label();
  check_iterations(counter,retno); 


  for(co=0;co<3;co++) {
  for(rei=0;rei<2;rei++){
    queue_fload(vec1[0][co][rei],VEC_IMM[0][co][rei],vec1ptr,FourSpinType);
    queue_fload(vec2[0][co][rei],VEC_IMM[0][co][rei],vec2ptr,FourSpinType);
  }
  }

  brchno = start_loop(counter);

  pragma(DCBT_SPACE,6);
  pragma(STORE_LIM,2);

  //  queue_iadd_imm(vec1pre,vec1ptr,VEC_SKIP_TWO);
  //  queue_iadd_imm(vec2pre,vec2ptr,VEC_SKIP_TWO);

  //  do_prefetch(vec1pre,0);
  //  do_prefetch(vec2pre,0);
  //  do_prefetch(vec1pre,1);
  //  do_prefetch(vec1pre,2);
  //  do_prefetch(vec2pre,1);
  //  do_prefetch(vec2pre,2);
  
  for(co=0;co<3;co++) {
  for(rei=0;rei<2;rei++){
    queue_fload(vec1[1][co][rei],VEC_IMM[1][co][rei],vec1ptr,FourSpinType);
    queue_fload(vec2[1][co][rei],VEC_IMM[1][co][rei],vec2ptr,FourSpinType);
    queue_fmadd(vec3[co][rei],A,vec1[0][co][rei],vec2[0][co][rei]);
  }
  }
  for(co=0;co<3;co++) {
  for(rei=0;rei<2;rei++){
    queue_fstore(vec3[co][rei],VEC_IMM[0][co][rei],vec2ptr,FourSpinType);
  }
  }

  for(co=0;co<3;co++) {
  for(rei=0;rei<2;rei++){
    queue_fmadd(vec3[co][rei],A,vec1[1][co][rei],vec2[1][co][rei]);
  }
  }
  for(co=0;co<3;co++) {
  for(rei=0;rei<2;rei++){
    queue_fstore(vec3[co][rei],VEC_IMM[1][co][rei],vec2ptr,FourSpinType);
  }
  }

  iterate_stream(PreVec1);
  iterate_stream(PreVec2);
  //  queue_prefetch(PreVec1);

  for(co=0;co<3;co++) {
  for(rei=0;rei<2;rei++){
    queue_fload(vec1[0][co][rei],VEC_IMM[0][co][rei],vec1ptr,FourSpinType);
    queue_fload(vec2[0][co][rei],VEC_IMM[0][co][rei],vec2ptr,FourSpinType);
  }
  }

  stop_loop(brchno,counter);

  make_inst(DIRECTIVE,Target,retno);

  restore_regs();
  free_stack();
  make_inst(DIRECTIVE,Exit_Routine,name);

  return;
}








