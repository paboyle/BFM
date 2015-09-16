
/*
 *
 *  Copyright UKQCD Collaboration, November 2000.
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
extern struct processor *PROC;

#include "registers.h"


void qcdoc_vaxpy_norm( char *);
int do_vaxmy = 0;
int do_vaxby = 0;
int overwrite1 = 0;
int overwrite2 = 0;
int swprefetch = 0;
Datum SpinorType = Double;
char name[80]="";

char procname[80] = "UNSPECIFIED"; /*Default processor type*/

int main ( int argc, char *argv[])
{
  int arg;
  char *c;

  /*Process command line*/
  while ( ( arg = getopt(argc,argv,"lrbn:P:ms")) != EOF){
    switch(arg){
    case 'n': if (strlen(optarg)<40) strcpy(name,optarg); break;
    case 'P': if (strlen(optarg)<40) strcpy(procname,optarg); break;
    case 'b': do_vaxby = 1; break;
    case 'r': overwrite2 = 1; break;
    case 'l': overwrite1 = 1; break;
    case 'm': do_vaxmy = 1; break;
    case 's': SpinorType = Single; break;
    default: fprintf(stderr,"Usage: %s -[h] -n routine_name -P proc\t",argv[0]); 
             exit (1); break;
    }
  }

  if ( overwrite1 && overwrite2 ) {
    printf("Cannot overwrite both x and y\n");
    exit(-1);
  }

  /*Control output according to user set up args*/
  set_processor_optarg(procname);
  setup_cmadds();

  /*queue the naive asm code*/
  qcdoc_vaxpy_norm(name);

  /*Filter through an virtual out of order processor */
  schedule_for_proc();
  /*Dump the resulting code*/
  dump_instruction_queue();
  return(0);

}

  /*
   * vaxpy ( a, vec1, vec2 , nvec )
   * fpoint a
   * fpoint vec1[nvec][Ncol][rei] 
   * fpoint vec2[nvec][Ncol][rei] 
   * asmint nvec
   *
   * for(i=0;i<nvec;i++)
   *   vec1[i] = A * vec1[i] + vec2 [i];
   *
   */

void qcdoc_vaxpy_norm( char *name)
{
  int dum = defargcount(7);

  /*Integer register usage*/
  alreg(vec1ptr,Iregs);
  alreg(vec2ptr,Iregs);
  alreg(vec3ptr,Iregs);
  alreg(counter,Iregs);
  alreg(Aptr,Iregs);
  alreg(Bptr,Iregs);
  alreg(outptr,Iregs);
  alreg(normp,Iregs);

  /*Floating register usage*/
  alreg(A,Fregs);
  alreg(B,Fregs);
  reg_array_2d(vec1,Fregs,3,2);
  reg_array_2d(vec2,Fregs,3,2);
  reg_array_2d(vec3,Fregs,3,2);
  reg_array_2d(dot ,Fregs,3,2);
  //  alreg(dot,Fregs);
  def_off(ZERO,SpinorType,0);
  def_off(VEC_ATOM,SpinorType,6);
  int Isize = PROC->I_size;
  int word_size = def_offset(Isize,Byte,"word_size");
  offset_2d(VEC_IMM,SpinorType,3,2);

  struct stream *PreVec1;
  struct stream *PreVec2;
  struct stream *PreVec3;

  int brchno,retno; /*Branch target handles*/
  int co,rei;

  make_inst(DIRECTIVE,Enter_Routine,name);
  grab_stack(0);
  save_regs();

  if ( (!overwrite1) && (!overwrite2) ) {
    getarg(vec3ptr);
  }
  getarg(Aptr);
  if ( do_vaxby) getarg(Bptr);
  getarg(vec1ptr);           /*Get args*/
  getarg(vec2ptr);           /*Get args*/
  getarg(counter);
  getarg(normp);

  alreg(tmp,Iregs);

  for (int i =0; i<12; i++ ) { 
    need_constant(i*2*SizeofDatum(SpinorType)*nsimd());
  }
  
  queue_fload(A,ZERO,Aptr,Double); /*Find the complex scale factor*/
  if ( do_vaxby ) queue_fload(B,ZERO,Bptr,Double); /*Find the complex scale factor*/
  queue_iload_imm(tmp,ZERO);

  /*
   * Ugly hack to get floating zero
   */
  queue_istore(tmp,ZERO,normp);
  if ( PROC->I_size < PROC->FP_size ) {
    queue_istore(tmp,word_size,normp);
  }
  /*
   * Insert a label to prevent reordering
   */
  make_inst(DIRECTIVE,Target,get_target_label());

  for ( co = 0 ; co < 3; co ++ ) {
    for ( rei=0;rei<2;rei++){
      queue_fload(dot[co][rei],ZERO,normp,Double);
    }
  }

  //  queue_fload(dot,ZERO,normp);


  PreVec1= create_stream(VEC_ATOM,vec1ptr ,counter,STREAM_IN ,LINEAR);
  PreVec2= create_stream(VEC_ATOM,vec2ptr ,counter,STREAM_IN ,LINEAR);
  if ( (!overwrite1) && (!overwrite2) ) {
    PreVec3 = create_stream(VEC_ATOM,vec3ptr,counter,STREAM_OUT,LINEAR);
  }

  /*Branch to stack restore if length <1*/
  retno = get_target_label();
  check_iterations(counter,retno); 

  /*
   * Start software pipeline
   */
  
  co = 0;
  for (rei=0;rei<2;rei++){
    queue_fload(vec1[co][rei],VEC_IMM[co][rei],vec1ptr,SpinorType);
    queue_fload(vec2[co][rei],VEC_IMM[co][rei],vec2ptr,SpinorType);
  }
  co = 1;
  for (rei=0;rei<2;rei++){
    queue_fload(vec1[co][rei],VEC_IMM[co][rei],vec1ptr,SpinorType);
    queue_fload(vec2[co][rei],VEC_IMM[co][rei],vec2ptr,SpinorType);
  }

  brchno = start_loop(counter);

  if ( swprefetch ) {
    pragma(DCBT_SPACE,5);
    do_prefetch(vec1ptr,0);
    do_prefetch(vec2ptr,0);
    do_prefetch(vec1ptr,1);
    do_prefetch(vec2ptr,1);
    do_prefetch(vec1ptr,2);
    do_prefetch(vec2ptr,2);
    if ( SpinorType == Double ){
    do_prefetch(vec1ptr,3);
    do_prefetch(vec2ptr,3);
    do_prefetch(vec1ptr,4);
    do_prefetch(vec2ptr,4);
    do_prefetch(vec1ptr,5);
    do_prefetch(vec2ptr,5);
    }
  }

  if ( overwrite1 ) {
    queue_iadd_imm(outptr,vec1ptr,ZERO);
  } else if ( overwrite2 ) {
    queue_iadd_imm(outptr,vec2ptr,ZERO);
  } else {
    queue_iadd_imm(outptr,vec3ptr,ZERO);
  }

  co=0;
  for(rei=0;rei<2;rei++){
    if ( do_vaxmy ) { 
      make_inst(FMACPIPE,FMSUB,vec3[co][rei],A,vec1[co][rei],vec2[co][rei]);
    } else if ( do_vaxby ) { 
      queue_fmul(vec2[co][rei],B,vec2[co][rei]);
      queue_fmadd(vec3[co][rei],A,vec1[co][rei],vec2[co][rei]);
    } else {
      queue_fmadd(vec3[co][rei],A,vec1[co][rei],vec2[co][rei]);
    }    
    queue_fmadd(dot[co][rei],vec3[co][rei],vec3[co][rei],dot[co][rei]);
    //    queue_fmadd(dot,vec3[co][rei],vec3[co][rei],dot);
    queue_fstore(vec3[co][rei],VEC_IMM[co][rei],outptr,SpinorType);
  }
  co = 1;
  int s=0;
  for(rei=0;rei<2;rei++){
    if ( do_vaxmy ) { 
      make_inst(FMACPIPE,FMSUB,vec3[co][rei],A,vec1[co][rei],vec2[co][rei]);
    } else if ( do_vaxby ) { 
      queue_fmul(vec2[co][rei],B,vec2[co][rei]);
      queue_fmadd(vec3[co][rei],A,vec1[co][rei],vec2[co][rei]);
    } else {
      queue_fmadd(vec3[co][rei],A,vec1[co][rei],vec2[co][rei]);
    }    
    queue_fmadd(dot[co][rei],vec3[co][rei],vec3[co][rei],dot[co][rei]);
    //    queue_fmadd(dot,vec3[co][rei],vec3[co][rei],dot);
    queue_fstore(vec3[co][rei],VEC_IMM[co][rei],outptr,SpinorType);
  }
  co = 2;
  for(rei=0;rei<2;rei++){
    queue_fload(vec1[co][rei],VEC_IMM[co][rei+s*2],vec1ptr,SpinorType);
    queue_fload(vec2[co][rei],VEC_IMM[co][rei+s*2],vec2ptr,SpinorType);
    if ( do_vaxmy ) { 
      make_inst(FMACPIPE,FMSUB,vec3[co][rei],A,vec1[co][rei],vec2[co][rei]);
    } else if ( do_vaxby ) { 
      queue_fmul (vec2[co][rei],B,vec2[co][rei]);
      queue_fmadd(vec3[co][rei],A,vec1[co][rei],vec2[co][rei]);
    } else {
      queue_fmadd(vec3[co][rei],A,vec1[co][rei],vec2[co][rei]);
    }    
    queue_fmadd(dot[co][rei],vec3[co][rei],vec3[co][rei],dot[co][rei]);
    //    queue_fmadd(dot,vec3[co][rei],vec3[co][rei],dot);
    queue_fstore(vec3[co][rei],VEC_IMM[co][rei],outptr,SpinorType);
  }


  /*
  if ( (!overwrite1) && (!overwrite2) ) {
    queue_prefetch(PreVec3);
  }
  */
  
    iterate_stream(PreVec1);
    iterate_stream(PreVec2);

    if ( (!overwrite1) && (!overwrite2) ) {
      iterate_stream(PreVec3);
    }
    queue_prefetch(PreVec1);
    queue_prefetch(PreVec2);

  co = 0;
  for (rei=0;rei<2;rei++){
    queue_fload(vec1[co][rei],VEC_IMM[co][rei],vec1ptr,SpinorType);
    queue_fload(vec2[co][rei],VEC_IMM[co][rei],vec2ptr,SpinorType);
  }
  co = 1;
  for (rei=0;rei<2;rei++){
    queue_fload(vec1[co][rei],VEC_IMM[co][rei],vec1ptr,SpinorType);
    queue_fload(vec2[co][rei],VEC_IMM[co][rei],vec2ptr,SpinorType);
  }


  stop_loop(brchno,counter);

  queue_fadd(dot[0][0],dot[0][0],dot[0][1]);
  queue_fadd(dot[1][0],dot[1][0],dot[1][1]);
  queue_fadd(dot[2][0],dot[2][0],dot[2][1]);
  queue_fadd(dot[0][0],dot[0][0],dot[1][0]);
  queue_fadd(dot[0][0],dot[0][0],dot[2][0]);
  queue_fstore(dot[0][0],ZERO,normp,Double);

  //  queue_fstore(dot,ZERO,normp);

  make_inst(DIRECTIVE,Target,retno);

  restore_regs();
  free_stack();
  make_inst(DIRECTIVE,Exit_Routine,name);

  return;
}








