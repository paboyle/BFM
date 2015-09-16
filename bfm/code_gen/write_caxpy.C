

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

#include "registers.h"


void qcdoc_caxpy( char *);

int overwrite1 = 0;
int overwrite2 = 0;
int do_vaxby = 0;

char name[80]="";

char procname[80] = "UNSPECIFIED"; /*Default processor type*/
Datum FourSpinType = Double;

int main ( int argc, char *argv[])
{
  int arg;
  char *c;

  int human = 0;

  /*Process command line*/
  while ( ( arg = getopt(argc,argv,"lrbRn:P:ms")) != EOF){
    switch(arg){
    case 'n': if (strlen(optarg)<20) strcpy(name,optarg); break;
    case 'P': if (strlen(optarg)<20) strcpy(procname,optarg); break;
    case 'r': overwrite2 = 1; break;
    case 'l': overwrite1 = 1; break;
    case 'R': human = 1; break;
    case 'b': do_vaxby = 1; break;
    case 's': FourSpinType= Single; break;
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
  qcdoc_caxpy(name);

  /*Filter through an virtual out of order processor */
  //schedule_for_proc();
  set_human_readable(human);
  /*Dump the resulting code*/
  dump_instruction_queue();
  return(0);

}

  /*
   * caxpy ( a, vec1, vec2 , nvec )
   * fpoint a
   * fpoint vec1[nvec][Ncol][rei] 
   * fpoint vec2[nvec][Ncol][rei] 
   * asmint nvec
   *
   * for(i=0;i<nvec;i++)
   *   vec1[i] = A * vec1[i] + vec2 [i];
   *
   */

void qcdoc_caxpy( char *name)
{
  int dum = defargcount(6);

  /*Integer register usage*/
  alreg(vec1ptr,Iregs);
  alreg(vec2ptr,Iregs);
  alreg(vec3ptr,Iregs);
  alreg(counter,Iregs);
  alreg(Aptr,Iregs);
  alreg(Bptr,Iregs);
  alreg(outptr,Iregs);

  /*Floating register usage*/
  alreg(A,Cregs);
  alreg(B,Cregs);
  reg_array_1d(vec1,Cregs,6);
  reg_array_1d(vec2,Cregs,6);
  reg_array_1d(vec3,Cregs,6);

  def_off(ZERO,FourSpinType,0);
  def_off(VEC_ATOM,FourSpinType,12*nsimd());
  offset_2d(VEC_IMM,FourSpinType,6,2*nsimd());

  int brchno,retno; /*Branch target handles*/
  int co,rei;

  make_inst(DIRECTIVE,Enter_Routine,name);
  grab_stack(0);
  save_regs();

  if ( (!overwrite1) && (!overwrite2) ) {
    getarg(vec3ptr);
  }
  getarg(Aptr);
  if ( do_vaxby )  getarg(Bptr);
  getarg(vec1ptr);           /*Get args*/
  getarg(vec2ptr);           /*Get args*/
  getarg(counter);

  for ( int i =0; i<6; i++ ) { 
    need_constant(i*2*SizeofDatum(FourSpinType)*nsimd());
  }
  for ( int i =0; i<3; i++ ) { 
    need_constant(i*32);
  }


  complex_csplat(A,ZERO,Aptr,Double);                 /*Find the complex scale factor - always use double*/
  if ( do_vaxby ) complex_csplat(B,ZERO,Bptr,Double); /*Find the complex scale factor*/

  /*Branch to stack restore if length <1*/
  retno = get_target_label();
  check_iterations(counter,retno); 

  /*
   * Start software pipeline
   */


  brchno = start_loop(counter);

  if ( overwrite1 ) {
    queue_iadd_imm(outptr,vec1ptr,ZERO);
  } else if ( overwrite2 ) {
    queue_iadd_imm(outptr,vec2ptr,ZERO);
  } else {
    queue_iadd_imm(outptr,vec3ptr,ZERO);
  }

  for(int i=0;i<3;i++){
    complex_load(vec1[i*2],VEC_IMM[i*2][0],vec1ptr,FourSpinType);
    complex_load(vec2[i*2],VEC_IMM[i*2][0],vec2ptr,FourSpinType);
  }

  for(int i=0;i<3;i++){
    complex_load(vec1[i*2+1],VEC_IMM[i*2+1][0],vec1ptr,FourSpinType);
    complex_load(vec2[i*2+1],VEC_IMM[i*2+1][0],vec2ptr,FourSpinType);
  }

  for(co=0;co<6;co++){
    if ( do_vaxby ) { 
      complex_mul(vec3[co] ,B,vec2[co]);
      complex_maddto(vec3[co],vec3[co],A,vec1[co]); 
    } else {
      complex_maddto(vec3[co],vec2[co],A,vec1[co]); // v3 = a v1 + v2
    }
  }

  queue_iadd_imm(vec1ptr,vec1ptr,VEC_ATOM);
  queue_iadd_imm(vec2ptr,vec2ptr,VEC_ATOM);
  do_prefetch(vec1ptr,0);
  do_prefetch(vec2ptr,0);

  for(co=0;co<6;co++){
    complex_store(vec3[co],VEC_IMM[co][0],outptr,FourSpinType);
  }

  queue_iadd_imm(vec3ptr,vec3ptr,VEC_ATOM);

  stop_loop(brchno,counter);


  make_inst(DIRECTIVE,Target,retno);

  restore_regs();
  free_stack();
  make_inst(DIRECTIVE,Exit_Routine,name);

  return;
}








