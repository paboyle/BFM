/*
 *
 *  Copyright UKQCD Collaboration, October 2000.
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
void do_prefetch ( int ptr,int line );
void do_writehint ( int ptr,int line );
#include "registers.h"


void need_cache_line(int line);


Datum FourSpinType = Double;
Datum TwoSpinType = Double;
Datum GaugeType = Double;

void qcdoc_su3_recon( char *);

int human = 0;
int dagger=0;

char name[80]="";
char procname[80]="UNSPECIFIED";

int main ( int argc, char *argv[])
{
  int arg;
  char *c;

  /*Process command line*/
  while ( ( arg = getopt(argc,argv,"sRn:P:d")) != EOF){
    switch(arg){
    case 'R': human = 1; break;
    case 'n': if (strlen(optarg)<30) strcpy(name,optarg); break;
    case 'P': if (strlen(optarg)<30) strcpy(procname,optarg);break;
    case 'd': dagger = 1; break;
    case 's': TwoSpinType = FourSpinType = GaugeType = Single; break;

    default: fprintf(stderr,"Usage: %s -[Rd] -n routine_name\n",
                     argv[0]); 
             fprintf(stderr,"\tR -> human readable .S format \n"); 
             fprintf(stderr,"\td -> dagger reconstruct\n"); 
             fprintf(stderr,"\ts -> Sloppy precision on 2-spinors\n"); 
             exit (1); break;
    }
  }

  set_processor_optarg(procname);
  set_human_readable(human);
  setup_cmadds(); /*Loses 6 regs*/

  /*Write the naive asm code*/
  qcdoc_su3_recon(name);

  /*Filter through an virtual out of order processor */
  schedule_for_proc();

  /*Dump the resulting code*/
  dump_instruction_queue();

  return(0);
}

void qcdoc_su3_recon( char *name)
{
  /****  This section defines all the registers and offsets I need ****/

  /*
   * This marks the argument registers as defined by ABI as off limits
   * to us until they are freed by "getarg()";
   */
  int dum = defargcount(4);

  /*Handle for the loop entry point*/
  int branchsite;
  int branchmu;
  int retno ;

  /*------------------------------------------------------------------
   * Floating point registers
   *------------------------------------------------------------------
   */

  // Reconstruct 8 registers for 4 spinor
  //  reg_array_2d(PSI,Fregs,4,2); 
  reg_array_3d(PSI,Fregs,3,4,2); 
  offset_3d(PSI_IMM,FourSpinType,4,3,2);    /*Offsets within 4 spinor*/

  // Reconstruct 2 spinor registers
#define  NEO 2
  reg_array_3d(Atmp,Fregs,1,2,2); /*CHIplus  regs */
  reg_array_3d(Btmp,Fregs,1,2,2); /*CHIminus regs */
  int A[NEO][2][2] = {
    Atmp[0][0][0],    Atmp[0][0][1],
    Atmp[0][1][0],    Atmp[0][1][1],
    -1,-1,-1,-1
  };
  int B[NEO][2][2] = {
    Btmp[0][0][0],    Btmp[0][0][1],
    Btmp[0][1][0],    Btmp[0][1][1],
    -1,-1,-1,-1
  };

  /*Regs for SU3 two spinor multiply ... overlap with the reconstruct*/
  /*                                                      registers  */
  int CHIR[3][2][2] = {
                 A[0][0][0],A[0][0][1],
                 A[0][1][0],A[0][1][1],
                 B[0][0][0],B[0][0][1],
                 B[0][1][0],B[0][1][1],
                 PSI[0][0][0],PSI[0][0][1],
                 PSI[0][1][0],PSI[0][1][1]
                };
  offset_3d(CHI_IMM,TwoSpinType,3,2,2);

  /*Registers for the gauge link (2 rows)*/
  int UA[3][2] = { 
                 {PSI[0][2][0],PSI[0][2][1]},
		 {PSI[2][1][0],PSI[2][1][1]},
                 {PSI[1][0][0],PSI[1][0][1]}
                };
  int UB[3][2] = { 
                 {PSI[1][1][0],PSI[1][1][1]},
		 {PSI[2][0][0],PSI[2][0][1]},
                 {PSI[1][2][0],PSI[1][2][1]},
                };
  offset_3d(GIMM    , GaugeType, 3, 3 ,2 );

  // Other 8 registers used for reduction variables in SU3.
  // Could use these in reconstruct??
  int E[2] = { PSI[2][2][0],PSI[2][2][1]};

  /*
   * FCD used for drain of Chi 
   * Overlap with PSI[*][3][*]
   */
  int F[2] = {PSI[0][3][0],PSI[0][3][1]};
  int C[2] = {PSI[1][3][0],PSI[1][3][1]};
  int D[2] = {PSI[2][3][0],PSI[2][3][1]};

  /*
   * Integer registers
   */
  alreg(psi,Iregs); 
  alreg(Umu,Iregs);
  alreg(Ufetch,Iregs);

  alreg(Chiin,Iregs);
  alreg(Chiout,Iregs);
  
  alreg(Chifetch,Iregs);

  reg_array_1d(Chiplus,Iregs,4);/*Pointers to the 8 2-spinors for recombination*/
  reg_array_1d(Chiminus,Iregs,4);

  alreg(mu,Iregs);
  alreg(Chidrain,Iregs);
  alreg(pref,Iregs);

  alreg(mem,Iregs);
  alreg(length,Iregs);

  int Isize = PROC->I_size;
  int Fsize = PROC->FP_size;

  def_off( ZERO_IMM, Byte,0);
  def_off( PSI_ATOM, FourSpinType, 24);
  def_off( CHI_ATOM, TwoSpinType, 12);
  def_off( PAD_CHI_ATOM, TwoSpinType, 16);
  def_off( MAT_IMM, GaugeType, 18);

  int Ndim   = def_offset(4,Byte,"Ndim");
  int Ndimm1 = def_offset(3,Byte,"Ndimm1");
  int hbias,bias;

  /*Offsets handles to stack*/
  int hbitbucket = def_offset(16*Isize,Byte,"hbitbucket");
  int Tsize;
  if ( TwoSpinType == Double ) Tsize = PROC->FP_size;
  else Tsize = PROC->FSP_size;
  int hstk0   = def_offset(16*Isize+12*Tsize  ,Byte,"hstk0");
  int hstk1   = def_offset(16*Isize+2*12*Tsize,Byte,"hstk1");
  int hstk2   = def_offset(16*Isize+3*12*Tsize,Byte,"hstk2");
  int hstk3   = def_offset(16*Isize+4*12*Tsize,Byte,"hstk3");

  int hIsize  = def_offset(Isize,Byte,"Isize");

  int i,co,j,k,nxt,ri,sp,nxtco,eop,eo_a,eo_b;

  /***********************************************************************/

                 /*
		  * PROLOGUE
		  */

  make_inst(DIRECTIVE,Enter_Routine,name);
  
  /*Allocate stack save any callee save registers we need etc...*/
  int stack_buf_size;
  stack_buf_size = 16*Isize + 
                   12*Fsize * 5 ;

  hbias = grab_stack(stack_buf_size);
  bias = get_offset(hbias);
  save_regs();
  queue_iadd_imm(mem,PROC->StackPointer,hbias); /*Pointer to buf on stack*/

  /*Define our arguments - all pointers ala fortran*/
  getarg(psi);
  getarg(Umu);
  getarg(Chiin);
  getarg(length);
  /*{... Process arguments ...*/

  queue_iload(length,ZERO_IMM,length);      /*Load in sx counter*/

  retno = get_target_label(); /*Branch to exit if yzt <1*/
  check_iterations(length,retno); 

  need_cache_line(0);
  need_cache_line(1);
  need_cache_line(2);
  need_cache_line(3);
  need_cache_line(4);

  pragma(DCBT_SPACE,5);
  pragma(DCBT_POST,1);

#define LOAD_U(comin,comax)\
  /*Load two link rows*/\
  for( i = comin;i<=comax;i++ ){\
    for( ri=0;ri<2;ri++){  \
      queue_fload(UA[i][ri],GIMM[i][0][ri],Umu,GaugeType);\
      queue_fload(UB[i][ri],GIMM[i][1][ri],Umu,GaugeType);\
    } \
  }

#define PRELOAD_U  LOAD_U(0,1)
#define POSTLOAD_U  LOAD_U(2,2)

  PRELOAD_U

#define LOAD_CHI(comin,comax) \
    /*Load Chi column*/\
    for( i = comin;i<=comax;i++ ){\
      for( ri=0;ri<2;ri++){\
        queue_fload(CHIR[i][0][ri],CHI_IMM[i][0][ri],Chiin,TwoSpinType);\
      } \
      for( ri=0;ri<2;ri++){\
        queue_fload(CHIR[i][1][ri],CHI_IMM[i][1][ri],Chiin,TwoSpinType);\
      } \
    }

#define PRELOAD_CHI  LOAD_CHI(0,1)
#define POSTLOAD_CHI  LOAD_CHI(2,2)

#define POSTLOAD \
    POSTLOAD_CHI \
    POSTLOAD_U 

  do_prefetch(Chiin,0);
  do_prefetch(Chiin,1);
  if ( SizeofDatum(TwoSpinType) == 8 ) do_prefetch(Chiin,2);

  PRELOAD_CHI
    
  /*
   * Start site loop
   */

  queue_iadd_imm(Chidrain,mem,hbitbucket);

  branchsite = start_loop(length);

  queue_iadd_imm(Chiout,mem,hstk0);

  /*
   * Loop over mu in asm
   */
  queue_iload_imm(mu,Ndimm1);

#define CHIDRAIN \
      queue_fstore(F[0],CHI_IMM[1][1][0],Chidrain,TwoSpinType);\
      queue_fstore(F[1],CHI_IMM[1][1][1],Chidrain,TwoSpinType);\
      queue_fstore(C[0],CHI_IMM[2][0][0],Chidrain,TwoSpinType);\
      queue_fstore(C[1],CHI_IMM[2][0][1],Chidrain,TwoSpinType);\
      queue_fstore(D[0],CHI_IMM[2][1][0],Chidrain,TwoSpinType);\
      queue_fstore(D[1],CHI_IMM[2][1][1],Chidrain,TwoSpinType);


#define PREFETCH_CHI \
  queue_iadd_imm(Chifetch,Chiin,PAD_CHI_ATOM);\
  do_prefetch(Chifetch,0);\
  do_prefetch(Chifetch,1);\
  if ( SizeofDatum(TwoSpinType) == 8 ) do_prefetch(Chifetch,2);

#define PREFETCH_CHIF \
  queue_iadd_imm(Chifetch,Chifetch,PAD_CHI_ATOM);\
  do_prefetch(Chifetch,0);\
  do_prefetch(Chifetch,1);\
  if ( SizeofDatum(TwoSpinType) == 8 ) do_prefetch(Chifetch,2);



  for ( int unroll=0;unroll<2;unroll++ ) { 

  if ( unroll==0 ) { 
  branchmu = start_loop(mu);
  pragma(DCBT_SPACE,5);
  pragma(STORE_LIM,1);
  pragma(LOAD_LIM,2);
  } else { 
  pragma(STORE_LIM,2);
  pragma(DCBT_SPACE,5);
  pragma(DCBT_POST,1);
  pragma(DCBT_PRE,0);
  pragma(LOAD_LIM,2);
  }

  CHIDRAIN
  POSTLOAD

  if ( unroll == 0 ) { 
    PREFETCH_CHI
    queue_iadd_imm(Ufetch,Umu,MAT_IMM);
    do_prefetch(Ufetch,0);
    do_prefetch(Ufetch,1);
    do_prefetch(Ufetch,2);
    if ( GaugeType == Double ) { 
      do_prefetch(Ufetch,3);
      do_prefetch(Ufetch,4);
    }
  } else {
    pragma(DCBT_SPACE,3);
    PREFETCH_CHI
    PREFETCH_CHIF
    PREFETCH_CHIF
    PREFETCH_CHIF
  }




    j=0;
    queue_three_cmuls(C[0],C[1],UA[j][0],UA[j][1],CHIR[j][0][0],CHIR[j][0][1],
		      D[0],D[1],UA[j][0],UA[j][1],CHIR[j][1][0],CHIR[j][1][1],
		      E[0],E[1],UB[j][0],UB[j][1],CHIR[j][0][0],CHIR[j][0][1]);
    j=1;
    queue_three_cmadds(C[0],C[1],UA[j][0],UA[j][1],CHIR[j][0][0],CHIR[j][0][1],
		       D[0],D[1],UA[j][0],UA[j][1],CHIR[j][1][0],CHIR[j][1][1],
		       E[0],E[1],UB[j][0],UB[j][1],CHIR[j][0][0],CHIR[j][0][1]);
    j=2;
                                                                                                                                                     
    queue_three_cmadds(C[0],C[1],UA[j][0],UA[j][1],CHIR[j][0][0],CHIR[j][0][1],
		       D[0],D[1],UA[j][0],UA[j][1],CHIR[j][1][0],CHIR[j][1][1],
		       E[0],E[1],UB[j][0],UB[j][1],CHIR[j][0][0],CHIR[j][0][1]);
                                                                                                                                                     
    /*Store the first three results*/
    queue_fstore(C[0],CHI_IMM[0][0][0],Chiout,TwoSpinType);
    queue_fstore(C[1],CHI_IMM[0][0][1],Chiout,TwoSpinType);
    queue_fstore(D[0],CHI_IMM[0][1][0],Chiout,TwoSpinType);
    queue_fstore(D[1],CHI_IMM[0][1][1],Chiout,TwoSpinType);
    queue_fstore(E[0],CHI_IMM[1][0][0],Chiout,TwoSpinType);
    queue_fstore(E[1],CHI_IMM[1][0][1],Chiout,TwoSpinType);

    /*Load the third row*/
    for(j=0;j<3;j++){
      for(ri=0;ri<2;ri++){
	queue_fload(UA[j][ri],GIMM[j][2][ri],Umu,GaugeType);
      }
    }
    /*Gauge layout is linear, mu faster than site*/
    queue_iadd_imm(Umu,Umu,MAT_IMM);


    /*Now the second set of three cdots*/
                                                                                                                                                     
    j=0;
    queue_three_cmuls(F[0],F[1],UB[j][0],UB[j][1],CHIR[j][1][0],CHIR[j][1][1],
		      C[0],C[1],UA[j][0],UA[j][1],CHIR[j][0][0],CHIR[j][0][1],
		      D[0],D[1],UA[j][0],UA[j][1],CHIR[j][1][0],CHIR[j][1][1]);
    j=1;
    queue_three_cmadds(F[0],F[1],UB[j][0],UB[j][1],CHIR[j][1][0],CHIR[j][1][1],
		       C[0],C[1],UA[j][0],UA[j][1],CHIR[j][0][0],CHIR[j][0][1],
		       D[0],D[1],UA[j][0],UA[j][1],CHIR[j][1][0],CHIR[j][1][1]);
    j=2;
    queue_three_cmadds(F[0],F[1],UB[j][0],UB[j][1],CHIR[j][1][0],CHIR[j][1][1],
		       C[0],C[1],UA[j][0],UA[j][1],CHIR[j][0][0],CHIR[j][0][1],
		       D[0],D[1],UA[j][0],UA[j][1],CHIR[j][1][0],CHIR[j][1][1]);

    /**************END SU3 CODE *************/
                                                                                                                                                     
    queue_iadd_imm(Chiin,Chiin,PAD_CHI_ATOM);
    queue_iadd_imm(Chidrain,Chiout,ZERO_IMM);
    queue_iadd_imm(Chiout,Chiout,CHI_ATOM);

    if ( unroll == 0 ) {

      PRELOAD_U
      PRELOAD_CHI

    }

  /*********************************************************/
  /****************** END OF SU3 MULTIPLY ******************/
  /*********************************************************/

  if ( unroll== 0 ) { 
  stop_loop(branchmu,mu); /* End loop over mu*/
  make_inst(DIRECTIVE,Target,get_target_label() ); /*delineate the sections*/
  }
  }

  
  /*********************************************************/
  /****************** START OF RECONSTRUCT *****************/
  /*********************************************************/

  //Address calculation...
  // Chiminus -> Stack  and  ChiPlus -> Chiin

  pragma(STORE_INORDER,1);
  queue_iadd_imm(Chiminus[0],mem,hstk0);

  /*For register use reasons loop over colour outermost*/

#define LOAD_CHI_MU0(eo,co) \
    for( sp = 0; sp<2;sp++ ){\
      for( ri = 0; ri<2;ri++ ){\
	queue_fload(A[eo][sp][ri],CHI_IMM[co][sp][ri],Chiminus[0],TwoSpinType);\
	if ( co == 0 ) {\
	  queue_fload(B[eo][sp][ri],CHI_IMM[co][sp][ri],Chiin,TwoSpinType);\
	  queue_iadd_imm(Chiplus[0],Chiin,ZERO_IMM);\
	} else {\
	  queue_fload(B[eo][sp][ri],CHI_IMM[co][sp][ri],Chiplus [0],TwoSpinType);\
	}\
      }}


  pragma(LOAD_LIM,2);
  LOAD_CHI_MU0(0,0)
  pragma(DCBT_POST,1);

  CHIDRAIN

  int neo_a = NEO;
  int neo_b = NEO;
  eo_a = 0;
  eo_b = 0;

  for ( co = 0; co <3 ; co ++ ){

    pragma(LOAD_LIM,1);
    if ( co == 0 ) {
      // Use the third colour for unrolling the loads
      A[1][0][0] = PSI[2][0][0];
      A[1][0][1] = PSI[2][0][1];
      A[1][1][0] = PSI[2][1][0];
      A[1][1][1] = PSI[2][1][1];
      B[1][0][0] = PSI[2][2][0];
      B[1][0][1] = PSI[2][2][1];
      B[1][1][0] = PSI[2][3][0];
      B[1][1][1] = PSI[2][3][1];
      queue_iadd_imm(Chiminus[1],mem,hstk1); // This is invariant of loop
                                             // Take out
      queue_iadd_imm(Chiplus[1],Chiin     ,PAD_CHI_ATOM);
    }

  /***************************************************************
  * MU = 0 reconstruct                                           *
  ****************************************************************/



    if ( co == 2 ) { 
      // Flip to not unrolled due to register pressure
      neo_b = 1;
      neo_a = 2;

      A[1][0][0] = PSI[0][0][0];
      A[1][0][1] = PSI[0][0][1];
      A[1][1][0] = PSI[1][0][0];
      A[1][1][1] = PSI[1][0][1];

      pragma(DCBT_POST,0);
      pragma(DCBT_SPACE,1);
      queue_iadd_imm(Ufetch,Umu,ZERO_IMM);
      //      do_prefetch(Ufetch,0);
      do_prefetch(Ufetch,1);
      do_prefetch(Ufetch,2);
      if ( GaugeType == Double ) {
	do_prefetch(Ufetch,3);
	do_prefetch(Ufetch,4);
      }
    }
    /* psi_0 =  Chiplus[0] + Chiminus[0] */
    /* psi_1 =  Chiplus[1] + Chiminus[1] */

    queue_fadd(PSI[co][0][0],B[eo_b][0][0],A[eo_a][0][0]);
    queue_fadd(PSI[co][0][1],B[eo_b][0][1],A[eo_a][0][1]);
    queue_fadd(PSI[co][1][0],B[eo_b][1][0],A[eo_a][1][0]);
    queue_fadd(PSI[co][1][1],B[eo_b][1][1],A[eo_a][1][1]);

    // Dagger = 0:
    /* psi_2 =-iChiplus[1] +iChiminus[1] */
    /* psi_3 =-iChiplus[0] +iChiminus[0] */
    // Dagger = 1:
    /* psi_2 = iChiplus[1] -iChiminus[1] */
    /* psi_3 = iChiplus[0] -iChiminus[0] */
    if ( dagger == 0 ) {
      queue_fsub(PSI[co][2][0],B[eo_b][1][1],A[eo_a][1][1]);
      queue_fsub(PSI[co][2][1],A[eo_a][1][0],B[eo_b][1][0]);
      queue_fsub(PSI[co][3][0],B[eo_b][0][1],A[eo_a][0][1]);
      queue_fsub(PSI[co][3][1],A[eo_a][0][0],B[eo_b][0][0]);
    } else { 
      queue_fsub(PSI[co][2][0],A[eo_a][1][1],B[eo_b][1][1]);
      queue_fsub(PSI[co][2][1],B[eo_b][1][0],A[eo_a][1][0]);
      queue_fsub(PSI[co][3][0],A[eo_a][0][1],B[eo_b][0][1]);
      queue_fsub(PSI[co][3][1],B[eo_b][0][0],A[eo_a][0][0]);
    }

  /***************************************************************
  * MU = 1 reconstruct                                           *
  ****************************************************************/
 
    eo_a = (eo_a+1)%neo_a;
    eo_b = (eo_b+1)%neo_b;
    for( sp = 0; sp<2;sp++ ){
      for( ri = 0; ri<2;ri++ ){
 
	queue_fload(A[eo_a][sp][ri],CHI_IMM[co][sp][ri],Chiminus[1],TwoSpinType);
	queue_fload(B[eo_b][sp][ri],CHI_IMM[co][sp][ri],Chiplus [1],TwoSpinType);

      }
    }

    if ( co == 0 ) {
      queue_iadd_imm(Chiminus[2],mem,hstk2);
      queue_iadd_imm(Chiminus[3],mem,hstk3);
      queue_iadd_imm(Chiplus[2],Chiplus[1],PAD_CHI_ATOM);
      queue_iadd_imm(Chiplus[3],Chiplus[2],PAD_CHI_ATOM);
    }

    /* psi_0 +=  Chiplus[0] + Chiminus[0] */
    /* psi_1 +=  Chiplus[1] + Chiminus[1] */

    queue_fadd(PSI[co][0][0],PSI[co][0][0],B[eo_b][0][0]);
    queue_fadd(PSI[co][0][1],PSI[co][0][1],B[eo_b][0][1]);
    queue_fadd(PSI[co][1][0],PSI[co][1][0],B[eo_b][1][0]);
    queue_fadd(PSI[co][1][1],PSI[co][1][1],B[eo_b][1][1]);

    queue_fadd(PSI[co][0][0],PSI[co][0][0],A[eo_a][0][0]);
    queue_fadd(PSI[co][0][1],PSI[co][0][1],A[eo_a][0][1]);
    queue_fadd(PSI[co][1][0],PSI[co][1][0],A[eo_a][1][0]);
    queue_fadd(PSI[co][1][1],PSI[co][1][1],A[eo_a][1][1]);

    //Dagger == 0
    /* psi_2 +=  Chiplus[1] - Chiminus[1] */
    /* psi_3 += -Chiplus[0] + Chiminus[0] */
    //Dagger == 1
    /* psi_2 -=  Chiplus[1] - Chiminus[1] */
    /* psi_3 -= -Chiplus[0] + Chiminus[0] */
    if ( dagger == 0 ) {
      queue_fadd(PSI[co][2][0],PSI[co][2][0],B[eo_b][1][0]);
      queue_fadd(PSI[co][2][1],PSI[co][2][1],B[eo_b][1][1]);
      queue_fsub(PSI[co][2][0],PSI[co][2][0],A[eo_a][1][0]);
      queue_fsub(PSI[co][2][1],PSI[co][2][1],A[eo_a][1][1]);

      queue_fsub(PSI[co][3][0],PSI[co][3][0],B[eo_b][0][0]);
      queue_fsub(PSI[co][3][1],PSI[co][3][1],B[eo_b][0][1]);
      queue_fadd(PSI[co][3][0],PSI[co][3][0],A[eo_a][0][0]);
      queue_fadd(PSI[co][3][1],PSI[co][3][1],A[eo_a][0][1]);
    } else { 
      queue_fsub(PSI[co][2][0],PSI[co][2][0],B[eo_b][1][0]);
      queue_fsub(PSI[co][2][1],PSI[co][2][1],B[eo_b][1][1]);
      queue_fadd(PSI[co][2][0],PSI[co][2][0],A[eo_a][1][0]);
      queue_fadd(PSI[co][2][1],PSI[co][2][1],A[eo_a][1][1]);

      queue_fadd(PSI[co][3][0],PSI[co][3][0],B[eo_b][0][0]);
      queue_fadd(PSI[co][3][1],PSI[co][3][1],B[eo_b][0][1]);
      queue_fsub(PSI[co][3][0],PSI[co][3][0],A[eo_a][0][0]);
      queue_fsub(PSI[co][3][1],PSI[co][3][1],A[eo_a][0][1]);
    }

  /***************************************************************
  * MU = 2 reconstruct                                           *
  ****************************************************************/
    eo_a = (eo_a+1)%neo_a;
    eo_b = (eo_b+1)%neo_b;
    for( sp = 0; sp<2;sp++ ){
      for( ri = 0; ri<2;ri++ ){
 
	queue_fload(A[eo_a][sp][ri],CHI_IMM[co][sp][ri],Chiminus[2],TwoSpinType);
	queue_fload(B[eo_b][sp][ri],CHI_IMM[co][sp][ri],Chiplus [2],TwoSpinType);

      }
    }
 
    /* psi_0 +=  Chiplus[0] + Chiminus[0] */
    /* psi_1 +=  Chiplus[1] + Chiminus[1] */

    queue_fadd(PSI[co][0][0],PSI[co][0][0],B[eo_b][0][0]);
    queue_fadd(PSI[co][0][1],PSI[co][0][1],B[eo_b][0][1]);
    queue_fadd(PSI[co][1][0],PSI[co][1][0],B[eo_b][1][0]);
    queue_fadd(PSI[co][1][1],PSI[co][1][1],B[eo_b][1][1]);

    queue_fadd(PSI[co][0][0],PSI[co][0][0],A[eo_a][0][0]);
    queue_fadd(PSI[co][0][1],PSI[co][0][1],A[eo_a][0][1]);
    queue_fadd(PSI[co][1][0],PSI[co][1][0],A[eo_a][1][0]);
    queue_fadd(PSI[co][1][1],PSI[co][1][1],A[eo_a][1][1]);

    //Dagger == 0
    /* psi_2 +=-iChiplus[0] +iChiminus[0] */
    /* psi_3 += iChiplus[1] -iChiminus[1] */
    //Dagger == 1

    /* psi_2 -=-iChiplus[0] +iChiminus[0] */
    /* psi_3 -= iChiplus[1] -iChiminus[1] */
    if ( dagger == 0 ) { 
      queue_fadd(PSI[co][2][0],PSI[co][2][0],B[eo_b][0][1]);
      queue_fsub(PSI[co][2][1],PSI[co][2][1],B[eo_b][0][0]);
      queue_fsub(PSI[co][2][0],PSI[co][2][0],A[eo_a][0][1]);
      queue_fadd(PSI[co][2][1],PSI[co][2][1],A[eo_a][0][0]);

      queue_fsub(PSI[co][3][0],PSI[co][3][0],B[eo_b][1][1]);
      queue_fadd(PSI[co][3][1],PSI[co][3][1],B[eo_b][1][0]);
      queue_fadd(PSI[co][3][0],PSI[co][3][0],A[eo_a][1][1]);
      queue_fsub(PSI[co][3][1],PSI[co][3][1],A[eo_a][1][0]);
    } else { 
      queue_fsub(PSI[co][2][0],PSI[co][2][0],B[eo_b][0][1]);
      queue_fadd(PSI[co][2][1],PSI[co][2][1],B[eo_b][0][0]);
      queue_fadd(PSI[co][2][0],PSI[co][2][0],A[eo_a][0][1]);
      queue_fsub(PSI[co][2][1],PSI[co][2][1],A[eo_a][0][0]);

      queue_fadd(PSI[co][3][0],PSI[co][3][0],B[eo_b][1][1]);
      queue_fsub(PSI[co][3][1],PSI[co][3][1],B[eo_b][1][0]);
      queue_fsub(PSI[co][3][0],PSI[co][3][0],A[eo_a][1][1]);
      queue_fadd(PSI[co][3][1],PSI[co][3][1],A[eo_a][1][0]);
    }


  /***************************************************************
  * MU = 3 reconstruct                                           *
  ****************************************************************/
    pragma(LOAD_LIM,2);
 
    eo_a = (eo_a+1)%neo_a;
    eo_b = (eo_b+1)%neo_b;
    for( sp = 0; sp<2;sp++ ){
      for( ri = 0; ri<2;ri++ ){
 	queue_fload(A[eo_a][sp][ri],CHI_IMM[co][sp][ri],Chiminus[3],TwoSpinType);
	queue_fload(B[eo_b][sp][ri],CHI_IMM[co][sp][ri],Chiplus [3],TwoSpinType );
      }
    }

    /* psi_0 +=  Chiplus[0] + Chiminus[0] */
    /* psi_1 +=  Chiplus[1] + Chiminus[1] */

    queue_fadd(PSI[co][0][0],PSI[co][0][0],B[eo_b][0][0]);
    queue_fadd(PSI[co][0][1],PSI[co][0][1],B[eo_b][0][1]);
    queue_fadd(PSI[co][1][0],PSI[co][1][0],B[eo_b][1][0]);
    queue_fadd(PSI[co][1][1],PSI[co][1][1],B[eo_b][1][1]);


    //Dagger == 0
    /* psi_2 +=  Chiplus[0] - Chiminus[0] */
    /* psi_3 +=  Chiplus[1] - Chiminus[1] */
    //Dagger == 1
    /* psi_2 -=  Chiplus[0] - Chiminus[0] */
    /* psi_3 -=  Chiplus[1] - Chiminus[1] */
    if ( dagger == 0 ) { 
      queue_fadd(PSI[co][2][0],PSI[co][2][0],B[eo_b][0][0]);
      queue_fadd(PSI[co][2][1],PSI[co][2][1],B[eo_b][0][1]);
      queue_fadd(PSI[co][3][0],PSI[co][3][0],B[eo_b][1][0]);
      queue_fadd(PSI[co][3][1],PSI[co][3][1],B[eo_b][1][1]);
    } else { 
      queue_fsub(PSI[co][2][0],PSI[co][2][0],B[eo_b][0][0]);
      queue_fsub(PSI[co][2][1],PSI[co][2][1],B[eo_b][0][1]);
      queue_fsub(PSI[co][3][0],PSI[co][3][0],B[eo_b][1][0]);
      queue_fsub(PSI[co][3][1],PSI[co][3][1],B[eo_b][1][1]);
    }

    queue_fadd(PSI[co][0][0],PSI[co][0][0],A[eo_a][0][0]);
    queue_fadd(PSI[co][0][1],PSI[co][0][1],A[eo_a][0][1]);
    queue_fadd(PSI[co][1][0],PSI[co][1][0],A[eo_a][1][0]);
    queue_fadd(PSI[co][1][1],PSI[co][1][1],A[eo_a][1][1]);

    if ( dagger == 0 ) { 
      queue_fsub(PSI[co][2][0],PSI[co][2][0],A[eo_a][0][0]);
      queue_fsub(PSI[co][2][1],PSI[co][2][1],A[eo_a][0][1]);
      queue_fsub(PSI[co][3][0],PSI[co][3][0],A[eo_a][1][0]);
      queue_fsub(PSI[co][3][1],PSI[co][3][1],A[eo_a][1][1]);
    } else { 
      queue_fadd(PSI[co][2][0],PSI[co][2][0],A[eo_a][0][0]);
      queue_fadd(PSI[co][2][1],PSI[co][2][1],A[eo_a][0][1]);
      queue_fadd(PSI[co][3][0],PSI[co][3][0],A[eo_a][1][0]);
      queue_fadd(PSI[co][3][1],PSI[co][3][1],A[eo_a][1][1]);
    }
    /*
     * Store the spinors. If this is problematic
     * in terms of PEC WriteBuf misses, I could
     * store to the stack and copy out later.
     */ 

    if ( co != 2 ) {
      LOAD_CHI_MU0(0,co+1)
	eo_a=0;
	eo_b=0;
    }

    queue_fstore(PSI[co][0][0],PSI_IMM[0][co][0],psi,FourSpinType);
    queue_fstore(PSI[co][0][1],PSI_IMM[0][co][1],psi,FourSpinType);

  }

  /*
   * Store out in linear order now
   */
  pragma(STORE_LIM,2);
  pragma(DCBT_SPACE,8);

  for ( co=0;co<3;co ++ ) {
    queue_fstore(PSI[co][1][0],PSI_IMM[1][co][0],psi,FourSpinType);
    queue_fstore(PSI[co][1][1],PSI_IMM[1][co][1],psi,FourSpinType);
  }
  for ( co=0;co<3;co ++ ) {
    queue_fstore(PSI[co][2][0],PSI_IMM[2][co][0],psi,FourSpinType);
    queue_fstore(PSI[co][2][1],PSI_IMM[2][co][1],psi,FourSpinType);
  }
  if ( TwoSpinType == FourSpinType ) {
       queue_iadd_imm(Chidrain,psi,CHI_ATOM);
  } else { 
    queue_iadd_imm(Chidrain,mem,hbitbucket);
    for ( co=0;co<3;co ++ ) {
      queue_fstore(PSI[co][3][0],PSI_IMM[3][co][0],psi,FourSpinType);
      queue_fstore(PSI[co][3][1],PSI_IMM[3][co][1],psi,FourSpinType);
    }
  }

  queue_iadd_imm(psi,psi,PSI_ATOM);
    /*
     * Put in an artificial dependency here
     * to try to stop the preloads getting above the last load of 
     * reconstruct.
     */
  queue_iadd_imm(Chiplus[3],Chiplus[3],ZERO_IMM);
  queue_iadd_imm(Chiin     ,Chiplus[3],PAD_CHI_ATOM);
  pragma(DCBT_SPACE,0);
  do_prefetch(Chiin,0);
  do_prefetch(Chiin,1);
  if ( SizeofDatum(TwoSpinType) == 8 )do_prefetch(Chiin,2);
  PRELOAD_U
  PRELOAD_CHI

  /* TERMINATION point of the loop*/
  stop_loop(branchsite,length);

  CHIDRAIN

  make_inst(DIRECTIVE,Target,retno);

                  /*
		   *
		   * EPILOGUE
		   *
		   */

  restore_regs();
  free_stack();
  make_inst(DIRECTIVE,Exit_Routine,name);

  return;

}









