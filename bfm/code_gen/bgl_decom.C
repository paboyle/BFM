
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
 * This version assumes that we have 32r + 32i registers and will overflow 
 * the register pool on everything except BG/L (and itanium)
 */
extern struct rotating_reg *CMADregs;

void qcdoc_decom_su3( char *);

/*Options flags*/
int human = 0;
char name[80];
int dagger = 0;


Datum FourSpinType=Double;
Datum TwoSpinType=Double;
Datum GaugeType=Double;

int StoreWithoutAllocate = 0;

void need_half_line(int line);
void need_cache_line(int line);
void do_prefetch ( int ptr,int line );
void do_writehint ( int ptr,int line );

void do_WH(int A,int B)
{
  if (StoreWithoutAllocate ) return;
  if (SizeofDatum(TwoSpinType) == 8 )  do_writehint(A,B);
  else do_prefetch(A,B); /*In single precision case, half cacheline prevents write hinting*/
}
char procname[80]="UNSPECIFIED";

int main ( int argc, char *argv[])
{
  struct processor *myproc;
  int arg;
  char *c;

  name[0] = '\0';
  /*Process command line*/
  while ( ( arg = getopt(argc,argv,"xsRL:n:P:UdV")) != EOF){
    switch(arg){
    case 'R': human = 1; break;
    case 'n': if (strlen(optarg)<30) strcpy(name,optarg); break;
    case 'P': if (strlen(optarg)<30) strcpy(procname,optarg); break;
    case 'd': dagger = 1; break;
    case 's': TwoSpinType = FourSpinType = GaugeType = Single; break;

    case 'x': StoreWithoutAllocate = 1; break;
    
    default: fprintf(stderr,"Usage: %s -[Rd] -n routine_name\n",argv[0]); 
             fprintf(stderr,"\tR -> human readable .S format \n"); 
             fprintf(stderr,"\td -> dagger decompose\n"); 
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
// Chi0123     -> mem(0,2,4,6)
// U Chi4567   -> mem(1,3,5,7)
//
// dagger == 1
// Chi4567     -> mem(0,2,4,6)
// U Chi0..3   -> mem(1,3,5,7)
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

/*
  Upointer, 
  Uregs, pointer register for chiout, Chi temp regs, Chi input regs 
 */
void hsu3_mult (int Up, 
               struct rotating_reg *Ureg,
               int chip,
               register_array_2d &ChiResult,
               register_array_2d &ChiIn);

void qcdoc_decom_su3( char *name)
{
  /****  This section defines all the registers and offsets I need ****/

  /*
   * This marks the argument registers as defined by ABI as off limits
   * to us until they are freed by "getarg()";
   */
  int dum = defargcount(5);

  /*Handles for the labels point*/
  int branchsite;
  int branchmu;
  int retno;

  /*-------------------------------------------------------------------------------
   * registers used in decompose phase.
   *-------------------------------------------------------------------------------
   */
  /* 12 cmplx Registers used for the PSI spinor  */
  reg_array_2d(Psi,Cregs,4,3); /*Input 4-spinor*/
  offset_3d(PSIIMM,FourSpinType,4,3,2);

  int t;

  /*-------------------------------------------------------------------------------
   * cmplx rotating registers used in SU3 phase.
   *-------------------------------------------------------------------------------
   */
  struct rotating_reg *Ureg = create_rotating_reg(Cregs,8,"Ureg"); /*Input 4-spinor*/

  /*
   * 6+6 cmplx Registers used for 2 spinors. Layout is reim, spin, colour.
   * This is transposed from usual so that we produce the 2 spinors
   * in a linearly moving fashion with the natural loop ordering.
   */
  reg_array_2d(ChiM ,Cregs,2,3); 
  reg_array_2d(ChiP ,Cregs,2,3); 
  offset_3d(CHIIMM,TwoSpinType,2,3,2);

  /*
   * Integer registers
   */
  alreg(psi_p,Iregs);    /*Pointer to current cpt of PSI field      */
  alreg(U_p,Iregs);    /*Pointer to the current cpt of gauge field*/
  alreg(Chi_p,Iregs); /*Pointer to the output to chi field       */
  alreg(length,Iregs); /*number of sites*/
  alreg(tab,Iregs);       /*Pointer to current entry in offset table*/
  alreg(Complex_i,Iregs);       /*Pointer to current entry in offset table*/
  alreg(pref,Iregs);       /*Pointer to current entry in offset table*/
  alreg(pref1,Iregs);       /*Pointer to current entry in offset table*/

  /*Useful integer immediate constants, in units of Fsize*/
  def_off( ZERO_IMM, Byte,0);
  def_off( PSI_IMM, FourSpinType,24);
  def_off( CHI_IMM, TwoSpinType,12);
  def_off( PAD_CHI_IMM, TwoSpinType,12);
  def_off( MAT_IMM, GaugeType,18);

  /*
   * Lookup tables
   */
  struct stream *PreChi;

  /*
   * Local variables for C++ code
   */
  int i,sp,co,j,k,nxt,ri;

  /*********************************************************************/


  /*--------------------------------------------------------------------
   * Start of the "pseudo assembler proper.
   *--------------------------------------------------------------------
   */
  pragma(STORE_LIM,1);
  pragma(DCBT_SPACE,3);
  pragma(DCBT_POST,1);
  make_inst(DIRECTIVE,Enter_Routine,name);
  grab_stack(0);
  save_regs();

  /*
   * Define our arguments 
   * QCDOC_ChDecom_[p,m]_[mhct]su3 (
   *                                fpoint PSI[xyzt][sx][4][3][2],
   *                                fpoint GAUGE[xyzt][sx][4][3][3][2],
   *                                asmint  *length,
   *                                void *tab
   *                               );
   */
  getarg(psi_p);
  getarg(U_p);
  getarg(length);
  getarg(tab);
  getarg(Complex_i);


  /*
   * Load common constants into Iregs
   */

  /*Bugger 12 int regs for offsets*/
  for ( i =0; i<12; i++ ) { 
    need_constant(i*2*SizeofDatum(FourSpinType));
  }
  for ( i =0; i<6; i++ ) { 
    need_constant(i*2*SizeofDatum(TwoSpinType));
  }
  for ( i =0; i<9; i++ ) { 
    need_constant(i*2*SizeofDatum(GaugeType));
  }

  /*{... Process arguments ...*/
  queue_iload(length,ZERO_IMM,length);  /*Load in counter*/


  retno = get_target_label(); /*Branch to exit if length <1*/
  check_iterations(length,retno); 

  PreChi = create_stream(CHI_IMM,Chi_p,length,STREAM_OUT,POINTER,tab);


  /*
   * Start site loop
   */

  for( co = 0; co<3;co++ ){ 
    for(  sp = 0; sp<4;sp++ ){ 
      complex_load(Psi[sp][co],PSIIMM[sp][co][0],psi_p,FourSpinType);
    }
  }
  int creg;/*Cycle through the whole set*/

  creg = get_rotating_register(Ureg);
  complex_constants_prepare(creg,Complex_i);

  branchsite = start_loop(length);


  /***************************************************************
  * Load 4 spinor *
  ****************************************************************/
  queue_iadd_imm(pref,U_p,ZERO_IMM);


  for ( int mu=0;mu<4;mu++ ) {

    if ( mu > 0 ) {
      creg = get_rotating_register(Ureg);
      complex_constants_prepare(creg,Complex_i);
    }

    switch(mu) { 
    case 0: 
      /* Gx
       *  0 0  0  i    [0]+-i[3]
       *  0 0  i  0    [1]+-i[2]
       *  0 -i 0  0
       * -i 0  0  0
       *
       */
      for( co=0; co < 3 ; co ++) complex_ApiB(ChiP[0][co],Psi[0][co],Psi[3][co]);
      for( co=0; co < 3 ; co ++) complex_AmiB(ChiM[0][co],Psi[0][co],Psi[3][co]);
      for( co=0; co < 3 ; co ++) complex_ApiB(ChiP[1][co],Psi[1][co],Psi[2][co]);
      for( co=0; co < 3 ; co ++) complex_AmiB(ChiM[1][co],Psi[1][co],Psi[2][co]);
      break;
    case 1:
      /*Gy
       *  0 0  0  -1  [0] -+ [3]
       *  0 0  1  0   [1] +- [2]
       *  0 1  0  0
       * -1 0  0  0
       */
      for( co=0; co < 3 ; co ++) complex_sub(ChiP[0][co],Psi[0][co],Psi[3][co]);
      for( co=0; co < 3 ; co ++) complex_add(ChiM[0][co],Psi[0][co],Psi[3][co]);
      for( co=0; co < 3 ; co ++) complex_add(ChiP[1][co],Psi[1][co],Psi[2][co]);
      for( co=0; co < 3 ; co ++) complex_sub(ChiM[1][co],Psi[1][co],Psi[2][co]);
      break;
    case 2:
      /*Gz
       *  0 0  i  0   [0]+-i[2]
       *  0 0  0 -i   [1]-+i[3]
       * -i 0  0  0
       *  0 i  0  0
       */
      for( co=0; co < 3 ; co ++) complex_ApiB(ChiP[0][co],Psi[0][co],Psi[2][co]);
      for( co=0; co < 3 ; co ++) complex_AmiB(ChiM[0][co],Psi[0][co],Psi[2][co]);
      for( co=0; co < 3 ; co ++) complex_AmiB(ChiP[1][co],Psi[1][co],Psi[3][co]);
      for( co=0; co < 3 ; co ++) complex_ApiB(ChiM[1][co],Psi[1][co],Psi[3][co]);
      break;
    case 3:
      /*Gt
       *  0 0  1  0 [0]+-[2]
       *  0 0  0  1 [1]+-[3]
       *  1 0  0  0
       *  0 1  0  0
       */
      for( co=0; co < 3 ; co ++) complex_add(ChiP[0][co],Psi[0][co],Psi[2][co]);
      for( co=0; co < 3 ; co ++) complex_sub(ChiM[0][co],Psi[0][co],Psi[2][co]);
      for( co=0; co < 3 ; co ++) complex_add(ChiP[1][co],Psi[1][co],Psi[3][co]);
      for( co=0; co < 3 ; co ++) complex_sub(ChiM[1][co],Psi[1][co],Psi[3][co]);
      break;
    default:
      exit(-1);
    }

    if ( dagger ) {
      for( sp=0; sp < 2 ; sp ++)
	for( co=0; co < 3 ; co ++)
	  complex_store(ChiP[sp][co],CHIIMM[sp][co][0],Chi_p,TwoSpinType);
      iterate_stream(PreChi);
      hsu3_mult(U_p,Ureg,Chi_p,ChiP,ChiM);
    } else { 
      for( sp=0; sp < 2 ; sp ++)
	for( co=0; co < 3 ; co ++)
	  complex_store(ChiM[sp][co],CHIIMM[sp][co][0],Chi_p,TwoSpinType);
      iterate_stream(PreChi);
      do_WH(Chi_p,0);
      do_WH(Chi_p,1);
      if ( SizeofDatum(TwoSpinType) == 8 ) do_WH(Chi_p,2);

      queue_iadd_imm(pref1,Chi_p,ZERO_IMM);
      hsu3_mult(U_p,Ureg,pref1,ChiM,ChiP);
    }

    iterate_stream(PreChi); 
    do_WH(Chi_p,0);
    do_WH(Chi_p,1);
    if ( SizeofDatum(TwoSpinType) == 8 ) do_WH(Chi_p,2);

    queue_iadd_imm(U_p,U_p,MAT_IMM);
    queue_iadd_imm(pref,pref,MAT_IMM);
    if ( mu & 0x0 )
      do_prefetch(pref,0);
    do_prefetch(pref,1);
    do_prefetch(pref,2);
    if( SizeofDatum(FourSpinType) == 8 ) {
      do_prefetch(pref,3);
      do_prefetch(pref,4);
    }

  }
  queue_iadd_imm(psi_p,psi_p,PSI_IMM);
/*
  queue_iadd_imm(pref,psi_p,PSI_IMM);
  do_prefetch(pref,0);
  do_prefetch(pref,1);
  do_prefetch(pref,2);
  do_prefetch(pref,3);
  do_prefetch(pref,4);
  do_prefetch(pref,5);
*/  
  for( co = 0; co<3;co++ ){ 
    for(  sp = 0; sp<4;sp++ ){ 
      complex_load(Psi[sp][co],PSIIMM[sp][co][0],psi_p,FourSpinType);
    }
  }
  creg = get_rotating_register(Ureg);
  complex_constants_prepare(creg,Complex_i);
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

void hsu3_mult (int Up, struct rotating_reg *Ureg,int chip,	
		register_array_2d &ChiResult,
		register_array_2d &ChiIn)
{

  offset_3d(GIMM    , GaugeType, 3, 3 ,2 );
  offset_3d(XIMM,TwoSpinType,2,3,2);

  /*Up is the pointer to the gfield*/
  int u[3][3];

  for(int i=0;i<3;i++) {
  for(int j=0;j<3;j++) {
    u[i][j] = get_rotating_register(Ureg);
  }
  }
  
  for(int j=0;j<3;j++) {
    complex_load(u[0][j],GIMM[0][j][0],Up,GaugeType);
    complex_load(u[1][j],GIMM[1][j][0],Up,GaugeType);
  }

  {
    int j=0; /*Chi[c][s] = U[c][j] Chi [j][s]*/
    complex_three_conjmuls (ChiResult[0][0],u[0][j],ChiIn[0][j],
			    ChiResult[1][0],u[0][j],ChiIn[1][j],
			    ChiResult[0][1],u[1][j],ChiIn[0][j]);
    j=1;
    complex_three_conjmadds(ChiResult[0][0],u[0][j],ChiIn[0][j],
			    ChiResult[1][0],u[0][j],ChiIn[1][j],
			    ChiResult[0][1],u[1][j],ChiIn[0][j]);
    j=2;
    complex_three_conjmadds(ChiResult[0][0],u[0][j],ChiIn[0][j],
  			    ChiResult[1][0],u[0][j],ChiIn[1][j],
			    ChiResult[0][1],u[1][j],ChiIn[0][j]);
  }

  for(int j=0;j<3;j++) {
    complex_load(u[2][j],GIMM[2][j][0],Up,GaugeType);
  }

  {
    int j=0; /*Chi[c][s] = U[c][j] Chi [j][s]*/
    complex_three_conjmuls (ChiResult[1][1],u[1][j],ChiIn[1][j],
			    ChiResult[0][2],u[2][j],ChiIn[0][j],
			    ChiResult[1][2],u[2][j],ChiIn[1][j]);
    j=1;
    complex_three_conjmadds(ChiResult[1][1],u[1][j],ChiIn[1][j],
			    ChiResult[0][2],u[2][j],ChiIn[0][j],
			    ChiResult[1][2],u[2][j],ChiIn[1][j]);
    j=2;
    complex_three_conjmadds(ChiResult[1][1],u[1][j],ChiIn[1][j],
			    ChiResult[0][2],u[2][j],ChiIn[0][j],
			    ChiResult[1][2],u[2][j],ChiIn[1][j]);
  }

  /*Save the output*/
  for(int sp=0;sp<2;sp++) {
   for(int j=0;j<3;j++) {
    complex_store(ChiResult[sp][j],XIMM[sp][j][0],chip,TwoSpinType);
   }
  }
}







