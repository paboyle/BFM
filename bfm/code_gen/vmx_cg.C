

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

char name[80]="";

char procname[80] = "UNSPECIFIED"; /*Default processor type*/
Datum FourSpinType = Double;

int main ( int argc, char *argv[])
{
  int arg;
  char *c;

  int human = 0;

  /*Process command line*/
  while ( ( arg = getopt(argc,argv,"Rn:P:s")) != EOF){
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
  //  schedule_for_proc();
  set_human_readable(human);
  /*Dump the resulting code*/
  dump_instruction_queue();
  return(0);

}

/* Optimised CG update
    // Can fuse these two together. Further, producing
    // "p" in a cache-blocked fashion could allow no reload
    // of p from main memory (~8 timeslices at a time??)
    // Psi update
    cp = axpy_norm(r,mmp,r,-a);
    axpy(psi,p,psi,a);
    axpy(p,p,r,b);
   * for(i=0;i<nvec;i++)
   *   vec3[i] =-A * vec4[i] + vec3 [i] 
   *   vec1[i] = A * vec2[i] + vec1 [i];
   *   vec2[i] = B * vec2[i] + vec3 [i];
   */

void qcdoc_vaxpy( char *name)
{
  int dum = defargcount(1);

  /*Integer register usage*/
  alreg(vec1ptr,Iregs);
  alreg(vec2ptr,Iregs);
  alreg(vec3ptr,Iregs);
  alreg(vec4ptr,Iregs);
  alreg(counter,Iregs);
  alreg(Aptr,Iregs);
  alreg(Bptr,Iregs);
  alreg(outptr1,Iregs);
  alreg(outptr2,Iregs);
  alreg(outptr3,Iregs);
  alreg(normp,Iregs);
  alreg(args,Iregs);

  /*Floating register usage*/
  alreg(A,Cregs);
  alreg(B,Cregs);
  reg_array_1d(vec1,Cregs,6);
  reg_array_1d(vec2,Cregs,6);
  reg_array_1d(vec3,Cregs,6);
  reg_array_1d(vec4,Cregs,6);
  alreg(nrm,Cregs);
  alreg(stmp,Fregs);
  alreg(sdot,Fregs);

  def_off( ZERO_IMM,Byte,0);
  def_off(ZERO,FourSpinType,0);


  def_off(VEC_ATOM,FourSpinType,12*nsimd());
  offset_2d(VEC_IMM,FourSpinType,6,2*nsimd());
  offset_2d(VEC_IMMD,Double,6,2*nsimd());
  int Isize      = def_offset(PROC->I_size,Byte,"Isize");

  int brchno,retno; /*Branch target handles*/
  int co,rei;

  make_inst(DIRECTIVE,Enter_Routine,name);
  int handle=grab_stack(SizeofDatum(Double)*2*nsimd());
  save_regs();

  getarg(args);   // alpha
  alreg(buf,Iregs);
  queue_iadd_imm(buf,PROC->StackPointer,handle);

  queue_iload(Aptr, ZERO_IMM,args);   queue_load_addr(args,Isize,args);   //alpha
  queue_iload(Bptr, ZERO_IMM,args);   queue_load_addr(args,Isize,args);   //beta
  queue_iload(vec1ptr, ZERO_IMM,args);   queue_load_addr(args,Isize,args);   //psi
  queue_iload(vec2ptr, ZERO_IMM,args);   queue_load_addr(args,Isize,args);   //p
  queue_iload(vec3ptr, ZERO_IMM,args);   queue_load_addr(args,Isize,args);   //r
  queue_iload(vec4ptr, ZERO_IMM,args);   queue_load_addr(args,Isize,args);   //mmp
  queue_iload(counter, ZERO_IMM,args);   queue_load_addr(args,Isize,args);   //length
  queue_iload(normp, ZERO_IMM,args);   queue_load_addr(args,Isize,args);   //nrm

  for ( int i =0; i<6; i++ ) { 
    need_constant(i*2*SizeofDatum(FourSpinType)*nsimd());
  }
  for ( int i =0; i<3; i++ ) { 
    need_constant(i*32);
    need_constant(i*64);
  }


  complex_rsplat(nrm,ZERO,normp,Double);
  complex_rsplat(A,ZERO,Aptr,Double); /*Find the complex scale factor - always use double*/
  complex_rsplat(B,ZERO,Bptr,Double); /*Find the complex scale factor*/

  /*Branch to stack restore if length <1*/
  retno = get_target_label();
  check_iterations(counter,retno); 

  /*
   * Start software pipeline
   */

  brchno = start_loop(counter);

  for(int i=0;i<3;i++){
    complex_load(vec1[i*2],VEC_IMM[i*2][0],vec1ptr,FourSpinType);
    complex_load(vec2[i*2],VEC_IMM[i*2][0],vec2ptr,FourSpinType);
    complex_load(vec3[i*2],VEC_IMM[i*2][0],vec3ptr,FourSpinType);
    complex_load(vec4[i*2],VEC_IMM[i*2][0],vec4ptr,FourSpinType);
  }
  queue_iadd_imm(outptr1,vec1ptr,ZERO);
  queue_iadd_imm(outptr2,vec2ptr,ZERO);
  queue_iadd_imm(outptr3,vec3ptr,ZERO);

  queue_iadd_imm(vec1ptr,vec1ptr,VEC_ATOM);
  queue_iadd_imm(vec2ptr,vec2ptr,VEC_ATOM);
  queue_iadd_imm(vec3ptr,vec3ptr,VEC_ATOM);

  do_prefetch(vec1ptr,0);
  do_prefetch(vec2ptr,0);
  do_prefetch(vec3ptr,0);

  for(int i=0;i<3;i++){
    complex_load(vec1[i*2+1],VEC_IMM[i*2+1][0],outptr1,FourSpinType);
    complex_load(vec2[i*2+1],VEC_IMM[i*2+1][0],outptr2,FourSpinType);
    complex_load(vec3[i*2+1],VEC_IMM[i*2+1][0],outptr3,FourSpinType);
    complex_load(vec4[i*2+1],VEC_IMM[i*2+1][0],vec4ptr,FourSpinType);
  }


  do_prefetch(vec1ptr,1);
  do_prefetch(vec2ptr,1);
  do_prefetch(vec3ptr,1);

  
  /*
   *   vec3[i] =-A * vec4[i] + Vec3 [i];
   *   vec1[i] = A * vec2[i] + vec1 [i];
   *   vec2[i] = B * vec2[i] + vec3 [i];
   */
  for(co=0;co<6;co++)  simd_nmadd(vec3[co],A,vec4[co],vec3[co]);
  for(co=0;co<6;co++){
    simd_madd (vec1[co],A,vec2[co],vec1[co]);
    complex_store(vec3[co],VEC_IMM[co][0],outptr3,FourSpinType);
  }

  queue_iadd_imm(vec4ptr,vec4ptr,VEC_ATOM);
  do_prefetch(vec4ptr,0);
  do_prefetch(vec4ptr,1);

  if ( FourSpinType == Double ) do_prefetch(vec1ptr,2);
  for(co=0;co<6;co++){
   simd_madd(vec2[co],B,vec2[co],vec3[co]);
   complex_store(vec1[co],VEC_IMM[co][0],outptr1,FourSpinType);
  }
  if ( FourSpinType == Double ) do_prefetch(vec2ptr,2);
  if ( FourSpinType == Double ) do_prefetch(vec3ptr,2);
  for(co=0;co<6;co++){
   simd_madd (nrm,vec3[co],vec3[co],nrm);
   complex_store(vec2[co],VEC_IMM[co][0],outptr2,FourSpinType);
  }
  if ( FourSpinType == Double ) do_prefetch(vec4ptr,2);


  stop_loop(brchno,counter);

  complex_store(nrm,VEC_IMMD[0][0],buf,Double);
  for(int i=0;i<2*nsimd();i++) {
    if ( i==0){
      queue_fload(sdot,VEC_IMMD[0][i],buf,Double);
    }else{
      queue_fload(stmp,VEC_IMMD[0][i],buf,Double);
      queue_fadd(sdot,sdot,stmp);
    }
  }
  queue_fstore(sdot,ZERO,normp,Double);

  make_inst(DIRECTIVE,Target,retno);

  restore_regs();
  free_stack();
  make_inst(DIRECTIVE,Exit_Routine,name);

  return;
}








