
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


void queue_face_scatter( char *);

void need_cache_line(int line);
void do_prefetch ( int ptr,int line );


char name[80]="";
char procname[80] = "UNSPECIFIED"; /*Default processor type*/
Datum TwoSpinType = Double;
int bgl=0;
int pad=16;
int caller_save  = 0;

int main ( int argc, char *argv[])
{
  int arg;
  char *c;

  /*Process command line*/
  while ( ( arg = getopt(argc,argv,"bVsn:P:")) != EOF){
    switch(arg){
    case 'n': if (strlen(optarg)<30) strcpy(name,optarg); break;
    case 'P': if (strlen(optarg)<30) strcpy(procname,optarg); break;
    case 'V': caller_save = 1; break;
    case 's': TwoSpinType = Single; break;
    case 'b': bgl = 1; pad=12; break;
    default: fprintf(stderr,"Usage: %s -[hs] -n routine_name -P proc\t",argv[0]); 
             exit (1); break;
    }
  }

  /*Control output according to user set up args*/
  set_processor_optarg(procname);
  if ( caller_save ) set_abi_caller_save_all();

  /*queue the naive asm code*/
  queue_face_scatter(name);
  /*Filter through an virtual out of order processor */
  schedule_for_proc();
  /*Dump the resulting code*/
  dump_instruction_queue();
  return(0);
}
/*
 * void face_scatter(Float *Chip,Float *rbp,Float **face_table, unsigned long nbound)
 * {
 *  for ( i = 0; i < nbound; i++) {
 *    face_p = (Float *)(Chip + PAD_HALF_SPINOR_SIZE * face_table[i] );
 *    for( atom = 0 ; atom < HALF_SPINOR_SIZE; atom++){
 *      face_p[atom] = rbp[atom];
 *    }
 *    rbp += PAD_HALF_SPINOR_SIZE;
 *  }
 * }
 */

void queue_face_scatter( char *name)
{
  int dum = defargcount(6);

  /*Integer register usage*/
  alreg(Chip,Iregs);
  alreg(rbp,Iregs);
  alreg(facetable,Iregs);
  alreg(nbound,Iregs);

  /*Floating register usage*/
  reg_array_1d(T,Cregs,6);

  def_off(CHI_ATOM  ,TwoSpinType,12);
  def_off(CHI_PADDED,TwoSpinType,pad); /*Skip a bit more for pec alignment  */
  offset_1d(CHI_IMM ,TwoSpinType,12);  /*Dont care about internal layout*/

  struct stream *ReceiveBuf;
  struct stream *Face;

  int brchno,retno;       /*Branch target handles*/
  int t1,t2,t3,t4,t5,t6;  /*integer vars for routine proper*/
  int i;

  make_inst(DIRECTIVE,Enter_Routine,name);
  grab_stack(0);
  save_regs();

  getarg(Chip);           /*Get args*/
  getarg(rbp);           /*Get args*/

  need_cache_line(0);
  need_cache_line(1);
  need_cache_line(2);

  need_cache_line(4);
  need_cache_line(5);
  need_cache_line(6);

  do_prefetch(rbp,0);
  do_prefetch(rbp,1);
  do_prefetch(rbp,2);

  getarg(facetable);     /*Get args*/
  getarg(nbound);


  ReceiveBuf = create_stream(CHI_PADDED,rbp ,nbound,STREAM_IN ,STRIDED,CHI_PADDED);
  Face       = create_stream(CHI_PADDED,Chip,nbound,STREAM_OUT,LOOKUP,facetable);

  for ( i =0; i<6; i++ ) { 
    need_constant(i*2*SizeofDatum(TwoSpinType));
  }

  /*Branch to stack restore if length <1*/
  retno = get_target_label();
  check_iterations(nbound,retno); 
  pragma(STORE_LIM,2);
  pragma(DCBT_SPACE,7);


  /* for i=0, i < nbound */
  brchno = start_loop(nbound);

  do_prefetch(rbp,0);
  do_prefetch(rbp,1);
  do_prefetch(rbp,2);
/*
  pragma(LOAD_LIM,6);
  pragma(STORE_LIM,2);

*/
  
  // for atom=0; atom < HALF_SPINOR_SIZE; atom++
  for ( i = 0 ; i < 6 ; i ++ ) {
    //    printf("i=%d\n",i);fflush(stdout);
    //freg T[i] = rpb+offset for ith element of spinor
    complex_load(T[i], CHI_IMM[i*2], rbp,TwoSpinType);
    complex_store(T[i], CHI_IMM[i*2], Chip,TwoSpinType);

  }
  
  iterate_stream(ReceiveBuf);
  iterate_stream(Face);


  // end for
  stop_loop(brchno,nbound);

  // where to jump once loop is over
  make_inst(DIRECTIVE,Target,retno);

  // Cleanup
  restore_regs();
  free_stack();

  // finish
  make_inst(DIRECTIVE,Exit_Routine,name);

  return;
}








