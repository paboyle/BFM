
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

#include "registers.h"

/*
 * This version assumes that we have 32 registers and will overflow 
 * the register pool on alpha (31 wired to zero)
 * Therefore Highly tweaked for QCDOC.
 * The reason is that the 5 cycle fpu latency forces us into a 
 * 3 interleaved CMADD rather than 2 interleaved CMADD mode to 
 * get the pipelining *absolutely* optimal.
 *
 * Now transpose the 2-spinor structure so that the su3-2spinor
 * writes come out in linear order
 */
extern struct rotating_reg *CMADregs;

void qcdoc_decom_su3( char *);

/*Options flags*/
int human = 0;
char name[80];
int dagger = 0;

Datum FourSpinType = Double;
Datum TwoSpinType  = Double;
Datum GaugeType    = Double;


void need_cache_line(int line);
void do_prefetch ( int ptr,int line );
void do_writehint ( int ptr,int line );
#define do_WH do_prefetch

char procname[80]="UNSPECIFIED";

int main ( int argc, char *argv[])
{
  struct processor *myproc;
  int arg;
  char *c;

  name[0] = '\0';
  /*Process command line*/
  while ( ( arg = getopt(argc,argv,"sRL:n:P:UdV")) != EOF){
    switch(arg){
    case 'R': human = 1; break;
    case 'n': if (strlen(optarg)<30) strcpy(name,optarg); break;
    case 'P': if (strlen(optarg)<30) strcpy(procname,optarg); break;
    case 'd': dagger = 1; break;
    case 's': TwoSpinType = FourSpinType = GaugeType = Single; break;
    
    default: fprintf(stderr,"Usage: %s -[Rd] -n routine_name\n",argv[0]); 
             fprintf(stderr,"\tR -> human readable .S format \n"); 
             fprintf(stderr,"\td -> dagger decompose\n"); 
             fprintf(stderr,"\ts -> sloppy precision on 2-spinors\n"); 
             exit (1); break;
    }
  }

  /*Control output according to user set up args*/
  set_processor_optarg(procname);
  set_human_readable(human);
  setup_cmadds();

  /*Write the naive asm code*/
  qcdoc_decom_su3(name);

  /*Filter through an virtual out of order processor */
  schedule_for_proc();

  /*Dump the resulting code*/
  dump_instruction_queue();

  return(0);
}

// Operation is:
//
// dagger == 0 
// Chi0..3     -> mem(i)
// Chi4..7     -> stack
// U Chi4..7   -> mem(i)
//
// dagger == 1
// Chi0..3     -> stack
// Chi4..7     -> mem(i-4)
// U Chi0..3   -> mem(i+4)
//
// Where:
//
// Chi0 => 1-gamma_0
// Chi1 => 1-gamma_1
// Chi2 => 1-gamma_2
// Chi3 => 1-gamma_3
// Chi4 => 1+gamma_0
// Chi5 => 1+gamma_1
// Chi6 => 1+gamma_2
// Chi7 => 1+gamma_3
//
// mem(i) is derived from pointer table tab in the following way:
// 
// dagger == 0:
//
//      (1-gamma_i) -> table(i)
//   U_i(1+gamma_i) -> table(i+4)
//
//   (1+gamma_i)    -> table(i)
//   U_i(1-gamma_i) -> table(i+4)
//
void qcdoc_decom_su3( char *name)
{
  /****  This section defines all the registers and offsets I need ****/

  /*
   * This marks the argument registers as defined by ABI as off limits
   * to us until they are freed by "getarg()";
   */
  int dum = defargcount(4);

  /*Handles for the labels point*/
  int branchsite;
  int branchmu;
  int retno;

  /*-------------------------------------------------------------------------------
   * Floating point registers used in decompose phase.
   *-------------------------------------------------------------------------------
   */


  /* 24 Registers used for the PSI spinor  (8 remaining are used below)*/
  reg_array_2d(U,Fregs,3,2); /*Input 4-spinor*/
  reg_array_2d(V,Fregs,3,2);
  reg_array_2d(W,Fregs,3,2);
  reg_array_2d(X,Fregs,3,2);
  /*Offsets within 4 spinor used for loading 4-spinor*/
  offset_2d(UIMM,FourSpinType,3,2);
  offset_2d_shft(VIMM,FourSpinType,3,2,6);
  offset_2d_shft(WIMM,FourSpinType,3,2,12);
  offset_2d_shft(XIMM,FourSpinType,3,2,18);

  int t;


  /*-------------------------------------------------------------------------------
   * Floating point registers used in SU3 phase.
   *-------------------------------------------------------------------------------
   */


  /*First row of SU3 matrix -> reuse "W" -> 6 regs*/
  int UA[3][2] = {
                 W[0][0],W[0][1],
                 W[1][0],W[1][1],
                 W[2][0],W[2][1],
                 };
  /*Second row of SU3 matrix -> reuse "X" -> 6 regs*/
  int UB[3][2] = {
                 X[0][0],X[0][1],
                 X[1][0],X[1][1],
                 X[2][0],X[2][1],
                 };
  offset_3d(GIMM    , GaugeType, 3, 3 ,2 );


  /*
   * Registers used for the 2 spinor. Layout is reim, spin, colour.
   * This is transposed from usual so that we produce the 2 spinors
   * in a linearly moving fashion with the natural loop ordering.
   */
  int Chir[3][2][2] = { 
                 U[0][0],U[0][1],
                 U[1][0],U[1][1],
                 U[2][0],U[2][1],
                 V[0][0],V[0][1],
                 V[1][0],V[1][1],
                 V[2][0],V[2][1]
                 };
  offset_3d(CHIIMM,TwoSpinType,3,2,2);


  /*
   * Reduction variables C,D,E,F - 8 of these for SU3 phase,.
   * Register pressure forces need to overlap with the tmp regs in decompose phase
   * 
   */
  reg_array_1d(C,Fregs,2);
  reg_array_1d(D,Fregs,2);
  reg_array_1d(E,Fregs,2);
  reg_array_1d(F,Fregs,2);

  /*
   * Integer registers
   */
  alreg(psi,Iregs);    /*Pointer to current cpt of PSI field      */
  alreg(Umu,Iregs);    /*Pointer to the current cpt of gauge field*/

  alreg(Chiin,Iregs);  /*Pointer to input 2spinor during SU3      */
  alreg(Chiout,Iregs); /*Pointer to the output to chi field       */

  alreg(mem,Iregs);    /*Pointer to our temporary storage area on stack*/
  alreg(length,Iregs); /*number of sites*/

  alreg(Chiplus,Iregs);     /*Pointers to projected Chi's */
  alreg(Chitmp1,Iregs);   
  alreg(Chitmp2,Iregs);   
  alreg(Chiminus,Iregs);    /*One is to stack, one is to mem
			     *which gets assigned where depends on sign of decompose
			     */
  
  alreg(pref,Iregs);     /*Pointers to projected Chi's */
  alreg(mu,Iregs);

  alreg(tab,Iregs);       /*Pointer to current entry in offset table*/

  int Isize = PROC->I_size;
  int Fsize = PROC->FP_size;
  if ( TwoSpinType == Single ) { 
   Fsize = PROC->FSP_size;
  }
 
  /*Useful integer immediate constants, in units of Fsize*/
  def_off( ZERO_IMM, Byte, 0);
  def_off( PSI_IMM, FourSpinType, 24);
  def_off( CHI_IMM, TwoSpinType, 12);
  def_off( PAD_CHI_IMM, TwoSpinType, 16);
  def_off( MAT_IMM, GaugeType, 18);

  /*
   *Useful immediate constants
   */
  int Ndim = def_offset(4,Byte,"Ndim");
  int Ndimm1 = def_offset(3,Byte,"Ndimm1");
  int hIsize  = def_offset(Isize,Byte,"hIsize"); 
  int hbias,bias;

  int hbitbucket = def_offset(16*Isize,Byte,"hbitbucket");
  int hstk0   = def_offset(16*Isize+1*12*Fsize,Byte,"hstk0");
  int hstk1   = def_offset(16*Isize+2*12*Fsize,Byte,"hstk1");
  int hstk2   = def_offset(16*Isize+3*12*Fsize,Byte,"hstk2");
  int hstk3   = def_offset(16*Isize+4*12*Fsize,Byte,"hstk3");
  int mult;

  /*
   * Lookup tables
   */
  struct stream *PreChi;

  /*
   * Local variables for C++ code
   */
  int i,co,j,k,nxt,ri;

  /*********************************************************************/


  /*--------------------------------------------------------------------
   * Start of the "pseudo assembler proper.
   *--------------------------------------------------------------------
   */

  make_inst(DIRECTIVE,Enter_Routine,name);

  /*Allocate stack, save any callee save registers we need etc...*/
  /* Layout of our region is as follows                   */
  /* 48*Fsize <- temporary spinor area for each mu        */
  /* 12*Fsize <- Bit bucket for writing on 2spin of trash */
  /* 16*Isize  <- Pointer storage                          */

  int stack_buf_size;
  stack_buf_size = Isize * 16
    + Fsize * 12 
    + Fsize * 12 * 4;

  hbias = grab_stack(stack_buf_size);
  bias  = get_offset(hbias);

  save_regs();

  queue_iadd_imm(mem,PROC->StackPointer,hbias);  /*Pointer to buf on stack*/

  /*
   * Define our arguments 
   * QCDOC_ChDecom_[p,m]_[mhct]su3 (
   *                                fpoint PSI[xyzt][sx][4][3][2],
   *                                fpoint GAUGE[xyzt][sx][4][3][3][2],
   *                                asmint  *length,
   *                                void *tab
   *                               );
   */
  need_cache_line(0);
  need_cache_line(1);
  need_cache_line(2);
  need_cache_line(3);
  need_cache_line(4);

  getarg(psi);
  getarg(Umu);
  getarg(length);
  getarg(tab);


  /*
   *  Layout of data section on Stack
   *
   *  mem + 4*Isize  <- &Psi ...We arrangefor the decompose to go to the (cached)
   *  mem + 3*Isize  <- &stk3...stack and the su3 to read back from the stack 
   *  mem + 2*Isize  <- &stk2...so that we can store gather the writes after the su3
   *  mem + 1*Isize  <- &stk1...and not have them ever in the cache.
   *  mem + 0*Isize  <- &stk0
   */

  /*{... Process arguments ...*/

  queue_iload(length,ZERO_IMM,length);  /*Load in counter*/

  retno = get_target_label(); /*Branch to exit if length <1*/
  check_iterations(length,retno); 


  PreChi = create_stream(CHI_IMM,Chiout,length,STREAM_OUT,POINTER,tab);


  /*
   * Start site loop
   */
  queue_iadd_imm(mem,PROC->StackPointer,hbias);    /*Locations on stack*/

  /* 
   * Pipeline the loading of the PSI vector
   */
#define LOAD_PSI_UPPER(comin,comax) \
  for( co = comin; co<=comax;co++ ){ \
    for( ri = 0; ri<2;ri++ ){\
      queue_fload(U[co][ri],UIMM[co][ri],psi,FourSpinType);\
    }\
    for( ri = 0; ri<2;ri++ ){\
      queue_fload(V[co][ri],VIMM[co][ri],psi,FourSpinType);\
    }\
  } 

#define LOAD_PSI_LOWER(comin,comax) \
  for( co = comin; co<=comax;co++ ){ \
    for( ri = 0; ri<2;ri++ ){\
      queue_fload(W[co][ri],WIMM[co][ri],psi,FourSpinType);\
    }\
    for( ri = 0; ri<2;ri++ ){\
      queue_fload(X[co][ri],XIMM[co][ri],psi,FourSpinType);\
    }\
  } 


#define PRELOAD_PSI \
  LOAD_PSI_UPPER(0,0)\
  LOAD_PSI_LOWER(0,0)\
  LOAD_PSI_UPPER(1,1)\
  LOAD_PSI_LOWER(1,1)

#define POSTLOAD_PSI\
  LOAD_PSI_UPPER(2,2)\
  LOAD_PSI_LOWER(2,2)

  PRELOAD_PSI

  int Utmp = Chitmp2;

  branchsite = start_loop(length);

  
  /***************************************************************
  * Load 4 spinor *
  ****************************************************************/
  pragma(LOAD_LIM,1);
  pragma(DCBT_SPACE,8);

  POSTLOAD_PSI


  /***************************************************************
  * MU = 0 decompose                                       *
  ****************************************************************/

  /*
   * Address computation:
   *
   * if ( dagger ) 
   *    Chiminus -> stack           & Chiplus -> tab(0) = Chiout
   * else
   *    Chiminus -> tab(0) = chiout & Chiplus -> stack
   */ 

  if ( dagger ) { 
    queue_iadd_imm(Chiminus,mem,hstk0);
    queue_iadd_imm(Chiplus,Chiout,ZERO_IMM); 
  } else { 
    queue_iadd_imm(Chiplus,mem,hstk0);  
    queue_iadd_imm(Chiminus,Chiout,ZERO_IMM); 
  }
  iterate_stream(PreChi);
  if ( dagger ) { /*tmp1==minus[mu=1]*/
    queue_iadd_imm(Chitmp1,mem,hstk1);
    queue_iadd_imm(Chitmp2,Chiout,ZERO_IMM); 
  } else { 
    queue_iadd_imm(Chitmp2,mem,hstk1);  
    queue_iadd_imm(Chitmp1,Chiout,ZERO_IMM); 
  }
  iterate_stream(PreChi);

  for( co=0; co < 3 ; co ++){

    // Minus decompose (2 stores = 1 gather line) for this color
     /* chia_lower =  psi(0) -i psi(3) */
    t = C[0];
    queue_fadd(t,U[co][0],X[co][1]);
    queue_fstore(t,CHIIMM[co][0][0],Chiminus,TwoSpinType);

    t = C[1];
    queue_fsub(t,U[co][1],X[co][0]);
    queue_fstore(t,CHIIMM[co][0][1],Chiminus,TwoSpinType);


    // Plus decompose (2 stores = 1 gather line) for this color
    /* chib_upper = psi(0) + i psi(3) */
    t = D[0];
    queue_fsub(t,U[co][0],X[co][1]);
    queue_fstore(t,CHIIMM[co][0][0],Chiplus,TwoSpinType);

    t = D[1];
    queue_fadd(t,U[co][1],X[co][0]);
    queue_fstore(t,CHIIMM[co][0][1],Chiplus,TwoSpinType);


    // Minus decompose (2 stores = 1 gather line) for this color
     /* chib_lower = psi(1) - i psi(2) */
    t = E[0];
    queue_fadd(t,V[co][0],W[co][1]);
    queue_fstore(t,CHIIMM[co][1][0],Chiminus,TwoSpinType);

    t = E[1];
    queue_fsub(t,V[co][1],W[co][0]);
    queue_fstore(t,CHIIMM[co][1][1],Chiminus,TwoSpinType);

    // Plus decompose (2 stores = 1 gather line) for this color
     /* chib_lower = psi(1) + i psi(2) */
    t = F[0];
    queue_fsub(t,V[co][0],W[co][1]);
    queue_fstore(t,CHIIMM[co][1][0],Chiplus,TwoSpinType);

    t = F[1];
    queue_fadd(t,V[co][1],W[co][0]);
    queue_fstore(t,CHIIMM[co][1][1],Chiplus,TwoSpinType);

  }


  /***************************************************************
  * MU = 1 decompose                                             *
  ****************************************************************/

  /*
   * Address computation:
   *
   * if ( dagger ) 
   *    Chiminus -> stack           & Chiplus -> tab(1) = Chiout
   * else
   *    Chiminus -> tab(1) = chiout & Chiplus -> stack
   */ 

  for( co=0; co < 3 ; co ++){

    /*Write on gather line of the ChiMinus decompose*/
    /* chib_upper = psi(0) + psi(3) */
    t = C[0];
    queue_fadd(t,U[co][0],X[co][0]);
    queue_fstore(t,CHIIMM[co][0][0],Chitmp1,TwoSpinType);

    t = C[1];
    queue_fadd(t,U[co][1],X[co][1]);
    queue_fstore(t,CHIIMM[co][0][1],Chitmp1,TwoSpinType);


    /*Write one gather line of the ChiPlus decompose*/
    /* chib_upper = psi(0) - psi(3) */
    t = D[0];
    queue_fsub(t,U[co][0],X[co][0]);
    queue_fstore(t,CHIIMM[co][0][0],Chitmp2,TwoSpinType);

    t = D[1];
    queue_fsub(t,U[co][1],X[co][1]);
    queue_fstore(t,CHIIMM[co][0][1],Chitmp2,TwoSpinType);


    /*Write one gather line of the ChiMinus decompose*/
     /* chib_upper = psi(1) - psi(2) */
    t = E[0];
    queue_fsub(t,V[co][0],W[co][0]);
    queue_fstore(t,CHIIMM[co][1][0],Chitmp1,TwoSpinType);

    t = E[1];
    queue_fsub(t,V[co][1],W[co][1]);
    queue_fstore(t,CHIIMM[co][1][1],Chitmp1,TwoSpinType);


    /*Write one gather line of the ChiPlus decompose*/
     /* chib_upper = psi(1) + psi(2) */
    t = F[0];
    queue_fadd(t,V[co][0],W[co][0]);
    queue_fstore(t,CHIIMM[co][1][0],Chitmp2,TwoSpinType);

    t = F[1];
    queue_fadd(t,V[co][1],W[co][1]);
    queue_fstore(t,CHIIMM[co][1][1],Chitmp2,TwoSpinType);

  }
  //Load in the Next link 
  pragma(DCBT_SPACE,4);
  do_prefetch(Umu,0);
  do_prefetch(Umu,1);
  do_prefetch(Umu,2);
  if ( GaugeType == Double ) {
    do_prefetch(Umu,3);
    do_prefetch(Umu,4);
  }

  /***************************************************************
  * MU = 2 decompose                                             *
  ****************************************************************/

  /*
   * Address computation:
   *
   * if ( dagger ) 
   *    Chiminus -> stack           & Chiplus -> tab(2) = Chiout
   * else
   *    Chiminus -> tab(2) = chiout & Chiplus -> stack
   */ 

  if ( dagger ) { 
    queue_iadd_imm(Chiminus,mem,hstk2);
    queue_iadd_imm(Chiplus,Chiout,ZERO_IMM); 
  } else { 
    queue_iadd_imm(Chiplus,mem,hstk2);  
    queue_iadd_imm(Chiminus,Chiout,ZERO_IMM); 
  }
  iterate_stream(PreChi);

  for( co=0; co < 3 ; co ++){
    

    /* chib_upper = psi(0) - i psi(2) */
    t = C[0];
    queue_fadd(t,U[co][0],W[co][1]);
    queue_fstore(t,CHIIMM[co][0][0],Chiminus,TwoSpinType);

    t = C[1];
    queue_fsub(t,U[co][1],W[co][0]);
    queue_fstore(t,CHIIMM[co][0][1],Chiminus,TwoSpinType);

    /* chib_upper = psi(0) + i psi(2) */
    t = D[0];
    queue_fsub(t,U[co][0],W[co][1]);
    queue_fstore(t,CHIIMM[co][0][0],Chiplus,TwoSpinType);

    t = D[1];
    queue_fadd(t,U[co][1],W[co][0]);
    queue_fstore(t,CHIIMM[co][0][1],Chiplus,TwoSpinType);

       /* chib_lower = psi(1) + i psi(3) */
    t = E[0];
    queue_fsub(t,V[co][0],X[co][1]);
    queue_fstore(t,CHIIMM[co][1][0],Chiminus,TwoSpinType);

    t = E[1];
    queue_fadd(t,V[co][1],X[co][0]);
    queue_fstore(t,CHIIMM[co][1][1],Chiminus,TwoSpinType);


       /* chib_lower = psi(1) - i psi(3) */
    t = F[0];
    queue_fadd(t,V[co][0],X[co][1]);
    queue_fstore(t,CHIIMM[co][1][0],Chiplus,TwoSpinType);

    t = F[1];
    queue_fsub(t,V[co][1],X[co][0]);
    queue_fstore(t,CHIIMM[co][1][1],Chiplus,TwoSpinType);



  }


  /***************************************************************
  * MU = 3 decompose                                             *
  ****************************************************************/

  /*
   * Address computation:
   *
   * if ( dagger ) 
   *    Chiminus -> stack           & Chiplus -> tab(3) = Chiout
   * else
   *    Chiminus -> tab(3) = chiout & Chiplus -> stack
   */ 

  if ( dagger ) { 
    queue_iadd_imm(Chiminus,mem,hstk3);
    queue_iadd_imm(Chiplus,Chiout,ZERO_IMM); 
  } else { 
    queue_iadd_imm(Chiplus,mem,hstk3);  
    queue_iadd_imm(Chiminus,Chiout,ZERO_IMM); 
  }
  iterate_stream(PreChi);

#define DO_M_0(C,co,store)\
    /*Chib_upper = psi(0) - psi(2)*/\
    for ( ri = 0; ri < 2; ri ++){\
      t = C[ri];\
      queue_fsub(t,U[co][ri],W[co][ri]);\
      if(store) queue_fstore(t,CHIIMM[co][0][ri],Chiminus,TwoSpinType);\
    }

#define DO_P_0(D,co)\
    /*Chib_upper = psi(0) + psi(2)*/\
    for ( ri = 0; ri < 2; ri ++){\
      t = D[ri];\
      queue_fadd(t,U[co][ri],W[co][ri]);\
      queue_fstore(t,CHIIMM[co][0][ri],Chiplus,TwoSpinType);\
    }

#define DO_M_1(E,co,store) \
    /*Chib_lower = psi(1) - psi(3)*/\
    for ( ri = 0; ri < 2; ri ++){\
      t = E[ri];\
      queue_fsub(t,V[co][ri],X[co][ri]);\
      if (store) queue_fstore(t,CHIIMM[co][1][ri],Chiminus,TwoSpinType);\
    }

#define DO_P_1(F,co)\
    /*Chib_lower = psi(1) + psi(3)*/\
    for ( ri = 0; ri < 2; ri ++){\
      t = F[ri];\
      queue_fadd(t,V[co][ri],X[co][ri]);\
      queue_fstore(t,CHIIMM[co][1][ri],Chiplus,TwoSpinType);\
    }
  DO_M_0(C,0,1)
  DO_P_0(D,0)
  DO_M_1(E,0,1)
  DO_P_1(F,0)

  pragma(STORE_LIM,1);
  DO_M_0(C,1,1)
  DO_P_0(D,1)

  DO_P_1(F,1)
  DO_P_0(C,2)
  DO_P_1(D,2)
    /*
     * Want to finish with M_1(F1)
     * Want to finish with M_0(C2)
     * Want to finish with M_1(D2)
     */
  DO_M_1(F,1,0)
  DO_M_0(C,2,0)
  DO_M_1(D,2,0)

  /*
   * Now reload the chi's from cache,
   * multiply by gauge links and store to pointers tab(4-7)
   */

  /*Set the pointers to the stack arrays*/
  queue_iload_imm(mu,Ndimm1);
  queue_iadd_imm(Chiin,mem,hstk0);
  queue_iadd_imm(Chitmp1,Chiminus,ZERO_IMM);

#define SAVE_CHI_DRAIN\
    queue_fstore(F[0],CHIIMM[1][1][0],Chitmp1,TwoSpinType);\
    queue_fstore(F[1],CHIIMM[1][1][1],Chitmp1,TwoSpinType);\
    queue_fstore(C[0],CHIIMM[2][0][0],Chitmp1,TwoSpinType);\
    queue_fstore(C[1],CHIIMM[2][0][1],Chitmp1,TwoSpinType);\
    queue_fstore(D[0],CHIIMM[2][1][0],Chitmp1,TwoSpinType);\
    queue_fstore(D[1],CHIIMM[2][1][1],Chitmp1,TwoSpinType);

#define LOAD_U(comin,comax) \
    for( i = comin;i<=comax;i++ ){\
      for( ri=0;ri<2;ri++){  \
	queue_fload(UA[i][ri],GIMM[0][i][ri],Umu,GaugeType);\
	queue_fload(UB[i][ri],GIMM[1][i][ri],Umu,GaugeType);\
      } \
    }

#define LOAD_CHI(comin,comax) \
    for( i = comin;i<=comax;i++ ){\
      for( ri=0;ri<2;ri++){\
        queue_fload(Chir[i][0][ri],CHIIMM[i][0][ri],Chiin,TwoSpinType);\
        queue_fload(Chir[i][1][ri],CHIIMM[i][1][ri],Chiin,TwoSpinType);\
      } \
    }

  pragma(LOAD_LIM,10);

  LOAD_CHI(0,0)
  LOAD_U(0,0)
  LOAD_CHI(1,1)

  /*
   * Loop over mu in asm
   */
  for ( int unroll=0;unroll<2;unroll++) { 

  pragma(DCBT_SPACE,4);
  pragma(DCBT_POST,1);
  if ( unroll == 0 ) {
  branchmu = start_loop(mu);
  }

  SAVE_CHI_DRAIN
  LOAD_U(1,1)

  if ( unroll == 1 ) {
    pragma(DCBT_SPACE,4);
    pragma(DCBT_POST,0);
    queue_iadd_imm(psi,psi,PSI_IMM); /*Increment psi pointer*/
    do_prefetch(psi,0);
    do_prefetch(psi,1);
    do_prefetch(psi,2);
    if ( FourSpinType == Double ) {
      queue_iadd_imm(pref,psi,CHI_IMM); /*Increment psi pointer*/
      do_prefetch(pref,0);
      do_prefetch(pref,1);
      do_prefetch(pref,2);
    }

  }


  pragma(LOAD_LIM,1);
  LOAD_CHI(2,2)
  LOAD_U(2,2)
    

/*******************SU3dag CODE PROPER *****************/
    int fr1,fr2,fr3,fr4;
    
    j=0;
    queue_three_conjmuls(C[0],C[1],UA[j][0],UA[j][1],Chir[j][0][0],Chir[j][0][1],
			 D[0],D[1],UA[j][0],UA[j][1],Chir[j][1][0],Chir[j][1][1],
			 E[0],E[1],UB[j][0],UB[j][1],Chir[j][0][0],Chir[j][0][1]);
    j=1;
    queue_three_conjmadds(C[0],C[1],UA[j][0],UA[j][1],Chir[j][0][0],Chir[j][0][1],
			  D[0],D[1],UA[j][0],UA[j][1],Chir[j][1][0],Chir[j][1][1],
			  E[0],E[1],UB[j][0],UB[j][1],Chir[j][0][0],Chir[j][0][1]);
    j=2;
    queue_three_conjmadds(C[0],C[1],UA[j][0],UA[j][1],Chir[j][0][0],Chir[j][0][1],
			  D[0],D[1],UA[j][0],UA[j][1],Chir[j][1][0],Chir[j][1][1],
			  E[0],E[1],UB[j][0],UB[j][1],Chir[j][0][0],Chir[j][0][1]);

    
    /*Store the first three results*/
    queue_fstore(C[0],CHIIMM[0][0][0],Chiout,TwoSpinType);
    queue_fstore(C[1],CHIIMM[0][0][1],Chiout,TwoSpinType);
    queue_fstore(D[0],CHIIMM[0][1][0],Chiout,TwoSpinType);
    queue_fstore(D[1],CHIIMM[0][1][1],Chiout,TwoSpinType);
    queue_fstore(E[0],CHIIMM[1][0][0],Chiout,TwoSpinType);
    queue_fstore(E[1],CHIIMM[1][0][1],Chiout,TwoSpinType);
    pragma(DCBT_SPACE,4);
    /*Load the third row*/
    for(j=0;j<3;j++){
      for(ri=0;ri<2;ri++){
	queue_fload(UA[j][ri],GIMM[2][j][ri],Umu,GaugeType);
      }
    }
    if ( unroll == 0 ) {
    queue_iadd_imm(pref,Umu,MAT_IMM);
    do_prefetch(pref,0);
    do_prefetch(pref,1);
    do_prefetch(pref,2);
    if ( GaugeType == Double ) {
      queue_iadd_imm(pref,pref,CHI_IMM);
      do_prefetch(pref,0);
      do_prefetch(pref,1);
    }
    }
    /*Gauge layout is linear, mu faster than site*/
    queue_iadd_imm(Umu,Umu,MAT_IMM);

    /*Now the second set of three cdots*/

    j=0;
    queue_three_conjmuls(F[0],F[1],UB[j][0],UB[j][1],Chir[j][1][0],Chir[j][1][1],
			 C[0],C[1],UA[j][0],UA[j][1],Chir[j][0][0],Chir[j][0][1],
			 D[0],D[1],UA[j][0],UA[j][1],Chir[j][1][0],Chir[j][1][1]);
    j=1;
    queue_three_conjmadds(F[0],F[1],UB[j][0],UB[j][1],Chir[j][1][0],Chir[j][1][1], 
			  C[0],C[1],UA[j][0],UA[j][1],Chir[j][0][0],Chir[j][0][1],
			  D[0],D[1],UA[j][0],UA[j][1],Chir[j][1][0],Chir[j][1][1]);
    j=2;
    queue_three_conjmadds(F[0],F[1],UB[j][0],UB[j][1],Chir[j][1][0],Chir[j][1][1],
			  C[0],C[1],UA[j][0],UA[j][1],Chir[j][0][0],Chir[j][0][1],
			  D[0],D[1],UA[j][0],UA[j][1],Chir[j][1][0],Chir[j][1][1]);


/**************END SU3 CODE *************/
    
    queue_iadd_imm(Chitmp1,Chiout,ZERO_IMM);
    queue_iadd_imm(Chiin,Chiin,CHI_IMM);

    iterate_stream(PreChi);

    if ( unroll == 0 ) { 

    LOAD_CHI(0,0)
    LOAD_U(0,0)
    LOAD_CHI(1,1)

    stop_loop(branchmu,mu); /* End loop over mu*/
    make_inst(DIRECTIVE,Target,get_target_label() ); /*delineate the sections*/
    }
  }
  pragma(STORE_LIM,1);

    /*
     * Update address values and fetch 6 cache lines. This is far too
     * early - the 16 lines for gauge fields probably evict this.
     * Want to only fetch this on the fourth link down below.
     * This would also free up 6 load store slots and save around 10 cycles
     */
  /*
   * NO reorder to keep the prefetches ahead of the loads
   */
  make_inst(DIRECTIVE,Target,get_target_label());
  pragma(LOAD_LIM,10);
  SAVE_CHI_DRAIN
  PRELOAD_PSI

  /* TERMINATION point of the loop*/
  stop_loop(branchsite,length);
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









