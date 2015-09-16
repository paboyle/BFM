
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


void qcdoc_merge( char *);
Datum SpinorType = Double;
int half_precision=0;
char name[80]="";

char procname[80] = "UNSPECIFIED"; /*Default processor type*/

int main ( int argc, char *argv[])
{
  int arg;
  char *c;

  /*Process command line*/
  while ( ( arg = getopt(argc,argv,"hn:P:ms")) != EOF){
    switch(arg){
    case 'h': half_precision=1; break;
    case 'n': if (strlen(optarg)<40) strcpy(name,optarg); break;
    case 'P': if (strlen(optarg)<40) strcpy(procname,optarg); break;
    case 's': SpinorType = Single; break;
    default: fprintf(stderr,"Usage: %s -[h] -n routine_name -P proc\t",argv[0]); 
             exit (1); break;
    }
  }

  if ( half_precision )  SpinorType=Half;

  /*Control output according to user set up args*/
  set_processor_optarg(procname);
  setup_cmadds();

  /*queue the naive asm code*/
  qcdoc_merge(name);

  schedule_for_proc();

  /*Dump the resulting code*/
  dump_instruction_queue();
  return(0);

}

  /*
   * vmx_merge ( out, vec1, vec2 , nvec )
   */

void qcdoc_merge( char *name)
{
  int dum = defargcount(5);

  /*Integer register usage*/
  alreg(outptr,Iregs);
  alreg(vec1ptr,Iregs);
  alreg(vec2ptr,Iregs);
  alreg(counter,Iregs);

  /*Floating register usage*/
  reg_array_1d(vec1,Cregs,3);
  reg_array_1d(vec2,Cregs,3);
  reg_array_1d(oreg,Cregs,6);
  alreg(permreg,Cregs);

  def_off(ZERO,SpinorType,0);
  def_off (IN_ATOM,SpinorType,6*nsimd()); // 2spins worth, 3 colors x complex 
  def_off (OUT_ATOM,SpinorType,12*nsimd());// 2spins worth, 3 colors x complex x simd  
  def_off(bits16,Byte,0xFFFF);
  def_off(thirtytwo,Byte,32);
  def_off(sixteen,Byte,16);


  offset_2d(CHI_IMM,SpinorType,6,2*nsimd());

  int Isize = PROC->I_size;
  int word_size = def_offset(Isize,Byte,"word_size");

  struct stream *PreOut;
  struct stream *PreVec1;
  struct stream *PreVec2;

  int brchno,retno; /*Branch target handles*/
  int co;

  make_inst(DIRECTIVE,Enter_Routine,name);
  int bias = grab_stack(64);
  save_regs();
  queue_iadd_imm(PROC->StackPointer,PROC->StackPointer,bias);

  getarg(outptr);           /*Get args*/
  getarg(vec1ptr);  
  getarg(vec2ptr);
  getarg(counter);

  alreg(Mask,Iregs);
  alreg(Convert1,Iregs);
  alreg(Convert2,Iregs);
  int memory = PROC->StackPointer;

  for (int i =0; i<6; i++ ) { 
    need_constant(i*2*SizeofDatum(SpinorType)*nsimd());
  }
  need_constant(64);
  complex_simd_init(permreg);

  if ( half_precision ) {
    queue_iload_imm(Mask,ZERO);
    queue_ori(Mask,Mask,bits16);
    queue_lshift(Mask,Mask,thirtytwo);
    queue_ori(Mask,Mask,bits16);
    queue_lshift(Mask,Mask,sixteen);
  }

  /*
   * Insert a label to prevent reordering
   */
  make_inst(DIRECTIVE,Target,get_target_label());

  PreVec1= create_stream(IN_ATOM,vec1ptr ,counter,STREAM_IN ,LINEAR);
  PreVec2= create_stream(IN_ATOM,vec2ptr ,counter,STREAM_IN ,LINEAR);
  PreOut = create_stream(OUT_ATOM,outptr  ,counter,STREAM_OUT ,LINEAR);

  /*Branch to stack restore if length <1*/
  retno = get_target_label();
  check_iterations(counter,retno); 

  /*
   * Start software pipeline
   */

  brchno = start_loop(counter);

  int indco[3]={0,1,2};
  int permute_mu=3;

  for(int ico=0;ico<3;ico++){
    co = indco[ico];
    // Could do entirely in integer unit for half precision to accelerate this
    if ( half_precision ) { 
      complex_load_half(vec1[co],CHI_IMM[co][0],vec1ptr,memory,Convert1,Convert2,Mask);
      complex_load_half(vec2[co],CHI_IMM[co][0],vec2ptr,memory,Convert1,Convert2,Mask);
    } else { 
      complex_load(vec1[co],CHI_IMM[co][0],vec1ptr,SpinorType);
      complex_load(vec2[co],CHI_IMM[co][0],vec2ptr,SpinorType);
    }
  }

  {
    // Merge the vectors
    for(co=0;co<3;co++) complex_simd_merge (0,permute_mu,oreg[co*2]  ,vec1[co],vec2[co]);
    for(co=0;co<3;co++) complex_simd_merge (1,permute_mu,oreg[co*2+1],vec1[co],vec2[co]);
  }
  //  make_inst(DIRECTIVE,LS_BARRIER);
  for(int i=0;i<6;i++){ // 2 SIMD sites, 3 colors, 2 spins 2 complex == 24 floats
    if ( half_precision ) { 
      complex_store_half(oreg[i],CHI_IMM[i][0],outptr,memory,Convert1,Convert2,Mask);
    } else { 
      complex_store(oreg[i],CHI_IMM[i][0],outptr,SpinorType);
    }
  }

  iterate_stream(PreVec1);
  iterate_stream(PreVec2);

  do_prefetch(vec1ptr,0);
  do_prefetch(vec2ptr,0);
  do_prefetch(vec1ptr,1);
  do_prefetch(vec2ptr,1);


  iterate_stream(PreOut);

  stop_loop(brchno,counter);
  
  make_inst(DIRECTIVE,Target,retno);

  queue_isub_imm(PROC->StackPointer,PROC->StackPointer,bias);
  restore_regs();
  free_stack();
  make_inst(DIRECTIVE,Exit_Routine,name);

  return;
}








