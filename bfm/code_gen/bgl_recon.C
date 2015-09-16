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
#define do_WH do_prefetch
#include "registers.h"


void need_cache_line(int line);
void need_half_line(int line);

Datum FourSpinType = Double;
Datum TwoSpinType  = Double;
Datum GaugeType    = Double;

void su3_mult(int Up, 
               struct rotating_reg *Ureg,
               int chip,
	       register_array_2d &ChiResult,
	       register_array_2d &ChiIn);


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
    case 's': TwoSpinType=FourSpinType=GaugeType=Single; break;

    default: fprintf(stderr,"Usage: %s -[Rd] -n routine_name\n",
                     argv[0]); 
             fprintf(stderr,"\tR -> human readable .S format \n"); 
             fprintf(stderr,"\td -> dagger reconstruct\n"); 
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
  int dum = defargcount(5);

  /*Handle for the loop entry point*/
  int branchsite;
  int retno ;

  /*------------------------------------------------------------------
   * Floating point registers
   *------------------------------------------------------------------
   */

  // Reconstruct 12 registers for 4 spinor
  reg_array_2d(PSI,Cregs,4,3); 
  offset_3d(PSI_IMM,FourSpinType,4,3,2);    /*Offsets within 4 spinor*/

  // Reconstruct 2 spinor registers
  reg_array_2d(ChiPlus ,Cregs,2,3); /*CHIplus  6-regs */
  reg_array_2d(ChiMinus,Cregs,2,3); /*CHIminus 6-regs */
  offset_3d(CHI_IMM,TwoSpinType,2,3,2);

  /*Registers for the gauge link (2 rows)*/
  struct rotating_reg *Ureg = create_rotating_reg(Cregs,8,"Ureg"); /*Input 4-spinor*/

  /*
   * Integer registers
   */
  alreg(psi_p,Iregs); 
  alreg(U_p,Iregs);
  alreg(Chi_p,Iregs);
  alreg(length,Iregs);
  alreg(Complex_i,Iregs);

  def_off( ZERO_IMM, Byte,0);
  def_off( PSI_ATOM, FourSpinType,24);
  def_off( CHI_ATOM, TwoSpinType,12);
  def_off( PAD_CHI_ATOM, TwoSpinType,12);
  def_off( MAT_IMM, GaugeType,18);

  int i,co,j,k,ri,sp;

  /***********************************************************************/

  
  /*Allocate stack save any callee save registers we need etc...*/
  make_inst(DIRECTIVE,Enter_Routine,name);
  grab_stack(0);
  save_regs();

  /*Define our arguments - all pointers ala fortran*/
  getarg(psi_p);
  getarg(U_p);
  getarg(Chi_p);
  getarg(length);
  getarg(Complex_i);

  /*
   * Load common constants into Iregs for indexed loads
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
  queue_iload(length,ZERO_IMM,length);      /*Load in sx counter*/

  retno = get_target_label(); /*Branch to exit if yzt <1*/
  check_iterations(length,retno); 

  /*
   * Start site loop
   */
  branchsite = start_loop(length);

  for( int mu=0;mu<4;mu++) {

    register_array_2d *UChi;
    register_array_2d *Chi;

    if( dagger == 1 ) { 
      UChi = & ChiPlus;
      Chi  = & ChiMinus;
    } else { 
      Chi  = & ChiPlus;
      UChi = & ChiMinus;
    }

    su3_mult(U_p, Ureg, Chi_p, *UChi, *Chi);

    queue_iadd_imm(Chi_p,Chi_p,PAD_CHI_ATOM);
    do_prefetch(Chi_p,1);
    do_prefetch(Chi_p,2);
    if ( SizeofDatum(TwoSpinType) == 8 ) do_prefetch(Chi_p,3);

    queue_iadd_imm(U_p,U_p,MAT_IMM);
    do_prefetch(U_p,1);
    do_prefetch(U_p,2);
    if(SizeofDatum(GaugeType)==8){
      do_prefetch(U_p,3);
      do_prefetch(U_p,4);
      do_prefetch(U_p,5);
    }

    for(int sp=0;sp<2;sp++) {
     for(int j=0;j<3;j++) {
      complex_load((*Chi)[sp][j],CHI_IMM[sp][j][0],Chi_p,TwoSpinType);
     }
    }
    queue_iadd_imm(Chi_p,Chi_p,PAD_CHI_ATOM);
    do_prefetch(Chi_p,1);
    do_prefetch(Chi_p,2);
    if ( SizeofDatum(TwoSpinType) == 8 ) do_prefetch(Chi_p,3);

    int tmp_reg  = get_rotating_register(Ureg);
    complex_constants_prepare(tmp_reg,Complex_i);

    /*Upper two components are simple, by virtue of our choice of projection*/      
    if ( mu == 0 ) {
      for(co=0;co<3;co++) complex_add(PSI[0][co],ChiPlus[0][co],ChiMinus[0][co]);
      for(co=0;co<3;co++) complex_add(PSI[1][co],ChiPlus[1][co],ChiMinus[1][co]);
    } else { 
      for(co=0;co<3;co++) complex_add(PSI[0][co],PSI[0][co],ChiPlus[0][co]);
      for(co=0;co<3;co++) complex_add(PSI[1][co],PSI[1][co],ChiPlus[1][co]);
      for(co=0;co<3;co++) complex_add(PSI[0][co],PSI[0][co],ChiMinus[0][co]);
      for(co=0;co<3;co++) complex_add(PSI[1][co],PSI[1][co],ChiMinus[1][co]);
    }


    /*Lower two components are complex, by virtue of our choice of projection*/      
    switch(mu) { 
      case 0:
        for(co=0;co<3;co++) complex_sub(PSI[2][co],ChiPlus[1][co],ChiMinus[1][co]);
        for(co=0;co<3;co++) complex_sub(PSI[3][co],ChiPlus[0][co],ChiMinus[0][co]);
	break;
      case 1:
        for(co=0;co<3;co++) {
          int tmp = get_rotating_register(Ureg);
	  if ( tmp == tmp_reg ) tmp = get_rotating_register(Ureg);

	  complex_AmiB(tmp, ChiPlus[1][co],PSI[2][co]);
	  complex_sub (PSI[2][co],tmp,ChiMinus[1][co]);
	}
        for(co=0;co<3;co++) {
          int tmp = get_rotating_register(Ureg);
	  if ( tmp == tmp_reg ) tmp = get_rotating_register(Ureg);

	  complex_AmiB(tmp,ChiMinus[0][co],PSI[3][co]);
	  complex_sub (PSI[3][co],tmp, ChiPlus[0][co]);
	}
	break;
      case 2:
    /* psi_2 +=-iChiplus[0] +iChiminus[0] */
    /* psi_3 += iChiplus[1] -iChiminus[1] */
        for(co=0;co<3;co++) complex_AmiB(PSI[2][co],PSI[2][co],ChiPlus[0][co]);
        for(co=0;co<3;co++) complex_ApiB(PSI[3][co],PSI[3][co],ChiPlus[1][co]);
        for(co=0;co<3;co++) complex_ApiB(PSI[2][co],PSI[2][co],ChiMinus[0][co]);
        for(co=0;co<3;co++) complex_AmiB(PSI[3][co],PSI[3][co],ChiMinus[1][co]);
	break;
      case 3:
    /* psi_2 +=  Chiplus[0] - Chiminus[0] */
    /* psi_3 +=  Chiplus[1] - Chiminus[1] */
        for(co=0;co<3;co++) complex_add(PSI[2][co],PSI[2][co],ChiPlus[0][co]);
        for(co=0;co<3;co++) complex_add(PSI[3][co],PSI[3][co],ChiPlus[1][co]);
        for(co=0;co<3;co++) complex_sub(PSI[2][co],PSI[2][co],ChiMinus[0][co]);
        for(co=0;co<3;co++) complex_sub(PSI[3][co],PSI[3][co],ChiMinus[1][co]);
	break;
    }     
  }

  for(sp=0;sp<4;sp++){
    for(co=0;co<3;co++){
      complex_store(PSI[sp][co],PSI_IMM[sp][co][0],psi_p,FourSpinType);
    }
  }

  queue_iadd_imm(psi_p,psi_p,PSI_ATOM);

  do_WH(psi_p,1); 
  do_WH(psi_p,2); 
  if ( SizeofDatum(FourSpinType) == 8 ) {
    do_WH(psi_p,3); 
    do_WH(psi_p,4); 
    do_WH(psi_p,5); 
  }

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

void su3_mult(int Up, 
               struct rotating_reg *Ureg,
               int chip,
	       register_array_2d &ChiResult,
	       register_array_2d &ChiIn)
{
  offset_3d(GIMM    , GaugeType, 3, 3 ,2);
  offset_3d(XIMM,TwoSpinType,2,3,2);

  /*Up is the pointer to the gfield*/
  int u[3][3];

  for(int i=0;i<3;i++) {
  for(int j=0;j<3;j++) {
    u[i][j] = get_rotating_register(Ureg);
  }
  }
  
  for(int j=0;j<3;j++) {
    complex_load(u[0][j],GIMM[j][0][0],Up,GaugeType);
    /*Load the input 2spinor*/
    for(int sp=0;sp<2;sp++) {
     complex_load(ChiIn[sp][j],XIMM[sp][j][0],chip,TwoSpinType);
    }
    complex_load(u[1][j],GIMM[j][1][0],Up,GaugeType);
  }


  {
    int j=0; /*Chi[c][s] = U[c][j] Chi [j][s]*/
    complex_three_cmuls (ChiResult[0][0],u[0][j],ChiIn[0][j],
			 ChiResult[1][0],u[0][j],ChiIn[1][j],
			 ChiResult[0][1],u[1][j],ChiIn[0][j]);
    j=1;
    complex_three_cmadds(ChiResult[0][0],u[0][j],ChiIn[0][j],
			 ChiResult[1][0],u[0][j],ChiIn[1][j],
			 ChiResult[0][1],u[1][j],ChiIn[0][j]);
    j=2;
    complex_three_cmadds(ChiResult[0][0],u[0][j],ChiIn[0][j],
			 ChiResult[1][0],u[0][j],ChiIn[1][j],
			 ChiResult[0][1],u[1][j],ChiIn[0][j]);
  }

  for(int j=0;j<3;j++) {
    complex_load(u[2][j],GIMM[j][2][0],Up,GaugeType);
  }

  {
    int j=0; /*Chi[c][s] = U[c][j] Chi [j][s]*/
    complex_three_cmuls (ChiResult[1][1],u[1][j],ChiIn[1][j],
			 ChiResult[0][2],u[2][j],ChiIn[0][j],
			 ChiResult[1][2],u[2][j],ChiIn[1][j]);
    j=1;
    complex_three_cmadds(ChiResult[1][1],u[1][j],ChiIn[1][j],
			 ChiResult[0][2],u[2][j],ChiIn[0][j],
			 ChiResult[1][2],u[2][j],ChiIn[1][j]);
    j=2;
    complex_three_cmadds(ChiResult[1][1],u[1][j],ChiIn[1][j],
			 ChiResult[0][2],u[2][j],ChiIn[0][j],
			 ChiResult[1][2],u[2][j],ChiIn[1][j]);
  }


}
