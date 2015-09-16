
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


void qcdoc_inner( char *);
int do_vaxmy = 0;
int do_vaxby = 0;
int overwrite1 = 0;
int overwrite2 = 0;
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
  qcdoc_inner(name);

  /*Filter through an virtual out of order processor */
  //  schedule_for_proc();

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

void qcdoc_inner( char *name)
{
  int dum = defargcount(5);

  /*Integer register usage*/
  alreg(vec1ptr,Iregs);
  alreg(vec2ptr,Iregs);
  //alreg(vec3ptr,Iregs);
  int outptr;
  alreg(counter,Iregs);
  //alreg(Aptr,Iregs);
  //alreg(Bptr,Iregs);
  alreg(normp,Iregs);
  alreg(tmp,Iregs);

  /*Floating register usage*/
  //alreg(A,Cregs);
  //alreg(B,Cregs);
  reg_array_1d(vec1,Cregs,6);
  reg_array_1d(vec2,Cregs,6);
  //reg_array_1d(vec3,Cregs,6);
  reg_array_1d(dot ,Cregs,6);
  alreg(sdot,Cregs);
  alreg(stmp,Cregs);
  reg_array_1d(fdot ,Fregs,2);
  reg_array_1d(ftmp ,Fregs,2);


  def_off(ZERO,SpinorType,0);
  def_off(VEC_ATOM,SpinorType,12*nsimd());
  offset_2d(VEC_IMM,SpinorType,6,2*nsimd());
  offset_2d(VEC_IMMD,Double,6,2*nsimd());

  int Isize = PROC->I_size;
  int word_size = def_offset(Isize,Byte,"word_size");

  struct stream *PreVec1;
  struct stream *PreVec2;
  //struct stream *PreVec3;

  int brchno,retno; /*Branch target handles*/
  int co;

  make_inst(DIRECTIVE,Enter_Routine,name);
  grab_stack(0);
  save_regs();

  //if ( (!overwrite1) && (!overwrite2) ) {
  //  getarg(vec3ptr);
  //}
  //getarg(Aptr);
  //if ( do_vaxby) getarg(Bptr);
  getarg(vec1ptr);           /*Get args*/
  getarg(vec2ptr);           /*Get args*/
  getarg(counter);
  getarg(normp);

  for (int i =0; i<12; i++ ) { 
    need_constant(i*2*SizeofDatum(SpinorType)*nsimd());
  }
  
  //complex_load(A,ZERO,Aptr,Double); /*Find the complex scale factor*/
  //if ( do_vaxby ) complex_load(B,ZERO,Bptr,Double); /*Find the complex scale factor*/

  /*
   * Insert a label to prevent reordering
   */
  make_inst(DIRECTIVE,Target,get_target_label());

  for ( co = 0 ; co < 6; co ++ ) {
    complex_load(dot[co],ZERO,normp,Double);
  }
  complex_load(sdot,ZERO,normp,Double);
  queue_fload(fdot[0],ZERO,normp,Double);
  queue_fload(fdot[1],ZERO,normp,Double);

  PreVec1= create_stream(VEC_ATOM,vec1ptr ,counter,STREAM_IN ,LINEAR);
  PreVec2= create_stream(VEC_ATOM,vec2ptr ,counter,STREAM_IN ,LINEAR);
  //if ( (!overwrite1) && (!overwrite2) ) {
  //  PreVec3 = create_stream(VEC_ATOM,vec3ptr,counter,STREAM_OUT,LINEAR);
  //}

  /*Branch to stack restore if length <1*/
  retno = get_target_label();
  check_iterations(counter,retno); 

  /*
   * Start software pipeline
   */

  brchno = start_loop(counter);

  if ( overwrite1 ) {
    outptr=vec1ptr;
  } else if ( overwrite2 ) {
    outptr=vec2ptr;
  } else {
    //outptr=vec3ptr;
  }

  int indco[6]={0,2,4,1,3,5};

  for(int ico=0;ico<6;ico++){
    co = indco[ico];
    complex_load(vec1[co],VEC_IMM[co][0],vec1ptr,SpinorType);
    complex_load(vec2[co],VEC_IMM[co][0],vec2ptr,SpinorType);
  }

  for(int ico=0;ico<6;ico++){
    co = ico;
    if ( do_vaxmy ) { 
      exit(-1);
    } else if ( do_vaxby ) { 
      //simd_mul(vec2[co],B,vec2[co]);
      //simd_madd(vec3[co],A,vec1[co],vec2[co]);
    } else {
      //simd_madd(vec3[co],A,vec1[co],vec2[co]);
	//complex_add(vec3[co],vec1[co],vec2[co]);
    }    
  }
  // In case we are overwriting on arg,
  // an extra prefetch here prevents refetch from L2 by L1p
  // Due to store colliding with data ahead in L1p
  iterate_stream(PreVec1);
  iterate_stream(PreVec2);
  do_prefetch(vec1ptr,0);
  do_prefetch(vec2ptr,0);

  for(int ico=0;ico<6;ico++){
    co = ico;
    //simd_madd(dot[co],vec3[co],vec3[co],dot[co]);
    complex_conjmaddto(dot[co],dot[co],vec1[co],vec2[co]);
  }
  //for(int ico=0;ico<6;ico++){
    //co = ico;
    //complex_store(vec3[co],VEC_IMM[co][0],outptr,SpinorType);
  //}


  //if ( (!overwrite1) && (!overwrite2) ) {
  //  iterate_stream(PreVec3);
  //}

  stop_loop(brchno,counter);
  
  for(co=1;co<6;co++) complex_add(dot[0],dot[0],dot[co]);
  complex_store(dot[0],ZERO,normp,Double);
  for(int i=0;i<nsimd();i++) {
    queue_fload(ftmp[0],VEC_IMMD[0][2*i],normp,Double);
    queue_fadd(fdot[0],fdot[0],ftmp[0]);
    queue_fload(ftmp[1],VEC_IMMD[0][2*i+1],normp,Double);
    queue_fadd(fdot[1],fdot[1],ftmp[1]);
    //complex_load(stmp,VEC_IMMD[i][0],normp,Double);
    //complex_add(sdot,sdot,stmp);
  }
  queue_fstore(fdot[0],ZERO,normp,Double);
  queue_fstore(fdot[1],VEC_IMMD[0][1],normp,Double);
  //complex_store(sdot,ZERO,normp,Double);
  make_inst(DIRECTIVE,Target,retno);

  restore_regs();
  free_stack();
  make_inst(DIRECTIVE,Exit_Routine,name);

  return;
}








