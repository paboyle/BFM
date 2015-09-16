

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


extern struct processor *PROC;
void qcdoc_vaxpy( char *);

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
    case 'R': human = 1; break;
    case 's': FourSpinType= Single; break;
    default: fprintf(stderr,"Usage: %s -[h] -n routine_name -P proc\t",argv[0]); 
             exit (1); break;
    }
  }

  /*Control output according to user set up args*/
  set_processor_optarg(procname);
  setup_cmadds();

  /*queue the naive asm code*/
  qcdoc_vaxpy(name);

  /*Filter through an virtual out of order processor */
  //schedule_for_proc();
  set_human_readable(human);
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

void qcdoc_vaxpy( char *name)
{
  int dum = defargcount(1);

  /*Integer register usage*/
  alreg(args,Iregs);
  alreg(vec1ptr,Iregs);
  alreg(vec2ptr,Iregs);
  alreg(counter,Iregs);
  alreg(skip,Iregs);
  alreg(Aptr,Iregs); // Coeffs
  alreg(outptr,Iregs);
  
  /*Floating register usage*/
  alreg(Al,Cregs);
  alreg(Bl,Cregs);
  alreg(Ah,Cregs);
  alreg(Bh,Cregs);

  reg_array_1d(vec1,Cregs,6);
  reg_array_1d(vec2,Cregs,6);
  reg_array_1d(vec3,Cregs,6);

  def_off(ZERO_IMM,Byte,0);
  def_off(ZERO,FourSpinType,0);
  def_off(DOUBLE_IMM,Double,1);
  def_off(VEC_ATOM,FourSpinType,12*nsimd());
  offset_2d(VEC_IMM,FourSpinType,12,2*nsimd());
  offset_2d(CONST_IMM,Double,4,2*nsimd());

  int Isize      = def_offset(PROC->I_size,Byte,"Isize");

  int brchno,retno; /*Branch target handles*/
  int co,rei;

  make_inst(DIRECTIVE,Enter_Routine,name);
  grab_stack(0);
  save_regs();

  getarg(args);
  queue_iload(outptr, ZERO_IMM,args);  queue_load_addr(args,Isize,args);  
  queue_iload(vec1ptr, ZERO_IMM,args);  queue_load_addr(args,Isize,args);  
  queue_iload(vec2ptr, ZERO_IMM,args);  queue_load_addr(args,Isize,args);  
  queue_iload(counter, ZERO_IMM,args);  queue_load_addr(args,Isize,args);  
  queue_iload(skip, ZERO_IMM,args);  queue_load_addr(args,Isize,args);  
  queue_iload(Aptr, ZERO_IMM,args); 

  for ( int i =0; i<12; i++ ) { 
    need_constant(i*2*SizeofDatum(FourSpinType)*nsimd());
  }
  for ( int i =0; i<3; i++ ) { 
    need_constant(i*32);
  }

  complex_rsplat(Al,ZERO_IMM,Aptr,Double);queue_load_addr(Aptr,DOUBLE_IMM,Aptr);  
  complex_rsplat(Ah,ZERO_IMM,Aptr,Double);queue_load_addr(Aptr,DOUBLE_IMM,Aptr);  
  complex_rsplat(Bl,ZERO_IMM,Aptr,Double);queue_load_addr(Aptr,DOUBLE_IMM,Aptr);  
  complex_rsplat(Bh,ZERO_IMM,Aptr,Double);queue_load_addr(Aptr,DOUBLE_IMM,Aptr);  

  /*Branch to stack restore if length <1*/
  retno = get_target_label();
  check_iterations(counter,retno); 

  /*
   * Start software pipeline
   */

  brchno = start_loop(counter);

  // Low halves
  for(int i=0;i<3;i++){
    complex_load(vec1[i*2],VEC_IMM[i*2][0],vec1ptr,FourSpinType);
    complex_load(vec2[i*2],VEC_IMM[i*2][0],vec2ptr,FourSpinType);
  }
  for(int i=0;i<3;i++){
    complex_load(vec1[i*2+1],VEC_IMM[i*2+1][0],vec1ptr,FourSpinType);
    complex_load(vec2[i*2+1],VEC_IMM[i*2+1][0],vec2ptr,FourSpinType);
  }
  for(co=0;co<6;co++){
    simd_mul(vec2[co] ,Bl,vec2[co]);
    simd_madd(vec3[co],Al,vec1[co],vec2[co]);
    complex_store(vec3[co],VEC_IMM[co][0],outptr,FourSpinType);
  }

  // High halves
  for(int i=0;i<3;i++){
    complex_load(vec1[i*2],VEC_IMM[i*2+6][0],vec1ptr,FourSpinType);
    complex_load(vec2[i*2],VEC_IMM[i*2+6][0],vec2ptr,FourSpinType);
  }
  for(int i=0;i<3;i++){
    complex_load(vec1[i*2+1],VEC_IMM[i*2+7][0],vec1ptr,FourSpinType);
    complex_load(vec2[i*2+1],VEC_IMM[i*2+7][0],vec2ptr,FourSpinType);
  }
  for(co=0;co<6;co++){
    simd_mul( vec2[co],Bh,vec2[co]);
    simd_madd(vec3[co],Ah,vec1[co],vec2[co]);
    complex_store(vec3[co],VEC_IMM[co+6][0],outptr,FourSpinType);
  }

  queue_iadd(vec1ptr,vec1ptr,skip);
  queue_iadd(vec2ptr,vec2ptr,skip);
  queue_iadd(outptr,outptr,skip);
  do_prefetch(vec1ptr,0);
  do_prefetch(vec2ptr,0);


  stop_loop(brchno,counter);


  make_inst(DIRECTIVE,Target,retno);

  restore_regs();
  free_stack();
  make_inst(DIRECTIVE,Exit_Routine,name);

  return;
}








