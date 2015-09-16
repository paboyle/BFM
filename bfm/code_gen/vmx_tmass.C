
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


void qcdoc_tmass( char *);
Datum SpinorType = Double;
char name[80]="";

char procname[80] = "UNSPECIFIED"; /*Default processor type*/

int main ( int argc, char *argv[])
{
  int arg;
  char *c;

  /*Process command line*/
  while ( ( arg = getopt(argc,argv,"n:P:ms")) != EOF){
    switch(arg){
    case 'n': if (strlen(optarg)<40) strcpy(name,optarg); break;
    case 'P': if (strlen(optarg)<40) strcpy(procname,optarg); break;
    case 's': SpinorType = Single; break;
    default: fprintf(stderr,"Usage: %s -[h] -n routine_name -P proc\t",argv[0]); 
             exit (1); break;
    }
  }

  /*Control output according to user set up args*/
  set_processor_optarg(procname);
  setup_cmadds();

  /*queue the naive asm code*/
  qcdoc_tmass(name);

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

void qcdoc_tmass( char *name)
{
  int dum = defargcount(5);

  /*Integer register usage*/
  alreg(vec1ptr,Iregs);
  alreg(vec2ptr,Iregs);
  alreg(counter,Iregs);
  alreg(Aptr,Iregs);
  alreg(Bptr,Iregs);
  alreg(outptr,Iregs);
  alreg(normp,Iregs);

  /*Floating register usage*/
  alreg(A,Cregs);
  alreg(B,Cregs);
  reg_array_1d(vec1,Cregs,3);
  reg_array_1d(vec2,Cregs,3);
  def_off(ZERO,SpinorType,0);
  def_off(VEC_ATOM,SpinorType,6*nsimd());
  offset_2d(VEC_IMM,SpinorType,3,2*nsimd());
  offset_2d(VEC_IMMD,Double,3,2*nsimd());

  int Isize = PROC->I_size;
  int word_size = def_offset(Isize,Byte,"word_size");

  struct stream *PreVec1;
  struct stream *PreVec2;

  int brchno,retno; /*Branch target handles*/
  int co;

  make_inst(DIRECTIVE,Enter_Routine,name);
  grab_stack(0);
  save_regs();

  getarg(vec2ptr);
  getarg(Aptr);
  getarg(Bptr);
  getarg(vec1ptr);           /*Get args*/
  getarg(counter);

  alreg(tmp,Iregs);

  for (int i =0; i<12; i++ ) { 
    need_constant(i*2*SizeofDatum(SpinorType)*nsimd());
  }
  
  complex_rsplat(A,ZERO,Aptr,Double); /*Find the complex scale factor*/
  complex_rsplat(B,ZERO,Bptr,Double); /*Find the complex scale factor*/

  /*
   * Insert a label to prevent reordering
   */
  make_inst(DIRECTIVE,Target,get_target_label());

  PreVec1= create_stream(VEC_ATOM,vec1ptr ,counter,STREAM_IN ,LINEAR);
  PreVec2= create_stream(VEC_ATOM,vec2ptr ,counter,STREAM_OUT ,LINEAR);

  /*Branch to stack restore if length <1*/
  retno = get_target_label();
  check_iterations(counter,retno); 

  /*
   * Start software pipeline
   */

  brchno = start_loop(counter);


  int indco[3]={0,2,1};

  for(int spin=0;spin<4;spin++){
  for(int ico=0;ico<3;ico++){
    co = indco[ico];
    complex_load(vec1[co],VEC_IMM[co][0],vec1ptr,SpinorType);
  }
  for(int ico=0;ico<3;ico++){
    co = indco[ico];

    // v2 = a V1 + i g5 V1
    simd_mul(vec2[co],A,vec1[co]); // Real mass term (coeff_ident)
    simd_mul(vec1[co],B,vec1[co]); // imag mass term (coeff_ident)
    if ( (spin <  2) ) {
      complex_ApiB(vec2[co],vec2[co],vec1[co]);
    } else {
      complex_AmiB(vec2[co],vec2[co],vec1[co]);
    }
  }
  for(int ico=0;ico<3;ico++){
    co = indco[ico];
    complex_store(vec2[co],VEC_IMM[co][0],vec2ptr,SpinorType);
  }

  iterate_stream(PreVec1);
  iterate_stream(PreVec2);
  }
  stop_loop(brchno,counter);
  
  make_inst(DIRECTIVE,Target,retno);

  restore_regs();
  free_stack();
  make_inst(DIRECTIVE,Exit_Routine,name);

  return;
}








