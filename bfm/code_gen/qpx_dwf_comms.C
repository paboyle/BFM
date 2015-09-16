/*
 *
 *  Copyright UKQCD Collaboration, March 2007.
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

void dwf_dslash( char *);
void set_label_fmt(char *);
/*Options flags*/
int human = 0;
char name[80];
int dagger = 0;
int dd2 = 0;

void do_writehint(int,int);
#undef LOCK_GAUGE
#define do_WH(A,B) ;
#define   lock l1_lock
#define unlock do_flush

//#define   lock do_prefetch
//#define unlock l1_unlock

void jump(int where);
void touch(int addr, int line);

void jump(int where)
{
  make_inst(BRCHPIPE,BRANCH,where); 
}
void touch(int addr, int line)
{
      do_flush(addr,line);
      if ( dd2 ) l2_touch(addr,line);
      else do_prefetch(addr,line);
}


Datum FourSpinType=Double;
Datum TableType=ShortInteger;
Datum GaugeType=Double;
int addto = 0;
char procname[80]="UNSPECIFIED";

bool do_interior=true;
bool do_exterior=true;
bool do_dperp   =false;
bool do_axpy    =true;
bool do_addto   =true;

int main ( int argc, char *argv[])
{
  struct processor *myproc;
  int arg;
  char *c;

  name[0] = '\0';
  /*Process command line*/
  while ( ( arg = getopt(argc,argv,"aRn:P:dsfiepx")) != EOF){
    switch(arg){
    case 'R': human = 1; break;
    case 'n': if (strlen(optarg)<30) strcpy(name,optarg); break;
    case 'P': if (strlen(optarg)<30) strcpy(procname,optarg); break;
    case 'd': dagger = 1; break;
    case 'a': addto = 1; break;
    case 'i': do_interior=1;do_exterior=0;break;
    case 'e': do_exterior=1;do_interior=0;break;
    case 'p': do_dperp=1;break;
    case 'x': do_axpy=1;break;
    case 'f': dd2 = 1; break;
    case 's': FourSpinType=Single; GaugeType=Single; break;
    default: fprintf(stderr,"Usage: %s -[Rd] -n routine_name\n",argv[0]); 
             fprintf(stderr,"\tR -> human readable .S format \n"); 
             fprintf(stderr,"\td -> dagger decompose\n"); 
             exit (1); break;
    }
  }

  if ( do_interior ) { 
    printf("/*Interior links*/\n");
  }
  if ( do_exterior ) { 
    printf("/*Exterior links*/\n");
  }
  if ( addto ) { 
    printf("/*Add to result field*/\n");
  }
  if ( do_dperp ) { 
    printf("/*5d hopping term*/\n");
  }

  /*Control output according to user set up args*/
  set_processor_optarg(procname);
  set_human_readable(human);
  setup_cmadds();

  /*Write the naive asm code*/
  dwf_dslash(name);

  /*Filter through an virtual out of order processor */
  schedule_for_proc();

  /*Dump the resulting code*/
  dump_instruction_queue();

  return(0);
}

void complex_six_cmaddtos(
		       int h1, int e1,int a1,int c1,
                       int h2, int e2,int a2,int c2,
                       int h3, int e3,int a3,int c3,
		       int h4, int e4,int a4,int c4,
                       int h5, int e5,int a5,int c5,
                       int h6, int e6,int a6,int c6);

void dwf_dslash( char *name)
{
  /*
   * This marks the argument registers as defined by ABI as off limits
   * to us until they are freed by "getarg()";
   */
  int dum = defargcount(1);
  int retno;

  /*------------------------------------------------------------------
   * registers used 
   *------------------------------------------------------------------
   */
  reg_array_2d(Psi,Cregs,4,3);  /* Accumulated 4-spinor - 12 regs*/
  reg_array_2d(Chi ,Cregs,2,3); // 2 spinor      - 6 regs 
  reg_array_2d(UChi,Cregs,2,3); // 2 spinors     - 6 regs
  reg_array_1d(tmp,Cregs,8); 
  int creg = tmp[7];

  int Chimu[4][3] = {
    {Chi[0][0],Chi[0][1],Chi[0][2]},
    {Chi[1][0],Chi[1][1],Chi[1][2]},
    {UChi[0][0],UChi[0][1],UChi[0][2]},
    {UChi[1][0],UChi[1][1],UChi[1][2]}
  };

  int U[3][3] = {
    {tmp[6],tmp[2],tmp[5]}, // Note reuse of tmp[6,0]. Care taken 
    {tmp[0],tmp[3],tmp[6]}, // U[1][2], U[2][2] are the danger elements
    {tmp[1],tmp[4],tmp[0]}
  };

  int UChiR[2][3] = { // Hold result for UChi 
    {tmp[1],tmp[2],tmp[3]},
    {tmp[4],tmp[5],tmp[6]}
  };
  offset_3d(CHIIMM,FourSpinType,2,3,2*nsimd());
  offset_3d(PSIIMM,FourSpinType,4,3,2*nsimd());

  int t;

  /*
   * Integer registers
   */
  alreg(U_p,Iregs);    /*Pointer to the current cpt of gauge field    */
  alreg(U_p_s,Iregs);  /*Pointer to the current cpt of gauge field    */
  alreg(Chi_p,Iregs);  /*Pointer to the input four spinor             */
  alreg(Psi_p,Iregs);  /*Pointer to current cpt output PSI field      */
  alreg(Chi_mu,Iregs); /*Pointer to the input four spinor            */
  alreg(length,Iregs);  /*number of sites*/
  alreg(tab,Iregs);     /*Pointer to current entry in offset table*/
  alreg(Complex_i,Iregs);  /*Point to (0,1)x Nsimd*/
  alreg(Ls,Iregs);        
  alreg(s,Iregs);        

  alreg(dperp_ss,Iregs);        
  alreg(dperp_s,Iregs);        
  alreg(dperp_XX,Iregs);

  alreg(s_offset,Iregs);        
  alreg(args,Iregs);        
  int ptr;

  /*Useful integer immediate constants, in units of Fsize*/
  def_off( ZERO_IMM,Byte,0);
  def_off( minusone,Byte,-1);
  def_off( one,Byte,1);
  def_off( m64,Byte,-64);
  def_off( PSI_IMM, FourSpinType,24*nsimd());
  def_off( CHI_IMM, FourSpinType,12*nsimd());
  def_off( MAT_IMM, GaugeType,18*nsimd());
  offset_3d(GIMM    , GaugeType, 3, 3 ,2*nsimd() );
  offset_1d(TAB_IMM,TableType,17);

  // Mask bits for predicating directions
  def_off( mask_0,Byte,1);
  def_off( mask_1,Byte,2);
  def_off( mask_2,Byte,4);
  def_off( mask_3,Byte,8);
  def_off( mask_4,Byte,16);
  def_off( mask_5,Byte,32);
  def_off( mask_6,Byte,64);
  def_off( mask_7,Byte,128);
  int mask_imm[8] = { 
    mask_0,
    mask_1,
    mask_2,
    mask_3,
    mask_4,
    mask_5,
    mask_6,
    mask_7
  };

  // Integer sizes
  int Isize      = def_offset(PROC->I_size,Byte,"Isize");
  int ISsize     = def_offset(PROC->IS_size,Byte,"ISsize");
  int i,j,co,sp;

  /*********************************************************************/

  pragma(DCBT_SPACE,0);
  pragma(LOAD_LIM,1);
  pragma(LOAD_SPACE,1);
  make_inst(DIRECTIVE,Enter_Routine,name);

  grab_stack(0);
  save_regs();

  /*********************************************
   * our arguments 
   *********************************************
   */
  getarg(args); /*Pointer to arg list*/

  queue_iload(Psi_p, ZERO_IMM,args);  queue_load_addr(args,Isize,args);  
  queue_iload(Chi_p, ZERO_IMM,args);  queue_load_addr(args,Isize,args);  
  queue_iload(U_p,   ZERO_IMM,args);  queue_load_addr(args,Isize,args);  
  queue_iload(length,ZERO_IMM,args);  queue_load_addr(args,Isize,args);  
  queue_iload(Ls,    ZERO_IMM,args);  queue_load_addr(args,Isize,args);  
  queue_iload(tab,   ZERO_IMM,args);  queue_load_addr(args,Isize,args);  
  queue_iload(Complex_i,ZERO_IMM,args);   

  /**************************************************
   * Load common constants into Iregs
   **************************************************
   */

  for (int i =0; i<12; i++ ) { 
    need_constant(i*2*SizeofDatum(FourSpinType)*nsimd());
  }
  for (int i =0; i<9; i++ ) { 
    need_constant(i*2*SizeofDatum(GaugeType)*nsimd());
  }
  complex_constants_prepare(creg,Complex_i);

  /**********************************
   * Site loop 
   **********************************
   */

  retno = get_target_label(); 
  check_iterations(length,retno); 
  int branchsite = start_loop(length);
#ifdef LOCK_GAUGE
  if ( do_interior ) {
  lock(tab,0);
  if ( GaugeType == Double) {
  queue_load_addr(s_offset,ZERO_IMM,U_p); // lock from new U_p
  for(int i=0;i<5;i++) { 
    lock(s_offset,0); // lock next set of links (32 lines, 2304 bytes)
    lock(s_offset,1); // Exactly 18*2*8*8 = 2304
    lock(s_offset,2);
    lock(s_offset,3);
    lock(s_offset,4);
    lock(s_offset,5);
    queue_load_addr(s_offset,PSI_IMM,s_offset); 
  }
  lock(s_offset,0); 
  lock(s_offset,1); 
  }
  }
#endif

  //////////////////////////////////////////////////
  // lock next set of links (32 lines, 2304 bytes)
  // lock from new U_p
  // Exactly 18*2*8*8 = 2304
  // Simply suppress if doing face only. Could predicate it.
  ///////////////////////////////////////////////////
  queue_iadd_imm(U_p_s,U_p,ZERO_IMM);
  
  /**********************************
   *  s-loop
   **********************************
   */
  queue_iadd_imm (s,Ls,ZERO_IMM);
  queue_iload_imm(dperp_s, ZERO_IMM);
  queue_iload_imm(s_offset,ZERO_IMM);        

  // Mask holding the directions on halo
  int mask = args;
  queue_iload_short(mask,TAB_IMM[10],tab);

  int branchls   = start_loop(s); 

  queue_iadd_imm(U_p_s,U_p,ZERO_IMM);

  /////////////////////////////////////////////////////////
  // Setup the Psi field for accumulation
  // Could scale Psi (axpy).
  /////////////////////////////////////////////////////////
  if ( addto ) {
    complex_load(Psi[0][0],PSIIMM[0][0][0],Psi_p,FourSpinType);
    complex_load(Psi[0][2],PSIIMM[0][2][0],Psi_p,FourSpinType);
    complex_load(Psi[1][1],PSIIMM[1][1][0],Psi_p,FourSpinType);
    complex_load(Psi[2][0],PSIIMM[2][0][0],Psi_p,FourSpinType);
    complex_load(Psi[2][2],PSIIMM[2][2][0],Psi_p,FourSpinType);
    complex_load(Psi[3][1],PSIIMM[3][1][0],Psi_p,FourSpinType);

    complex_load(Psi[0][1],PSIIMM[0][1][0],Psi_p,FourSpinType);
    complex_load(Psi[1][0],PSIIMM[1][0][0],Psi_p,FourSpinType);
    complex_load(Psi[1][2],PSIIMM[1][2][0],Psi_p,FourSpinType);
    complex_load(Psi[2][1],PSIIMM[2][1][0],Psi_p,FourSpinType);
    complex_load(Psi[3][0],PSIIMM[3][0][0],Psi_p,FourSpinType);
    complex_load(Psi[3][2],PSIIMM[3][2][0],Psi_p,FourSpinType);
  } else if ( do_dperp ) {

    /*checkerboarded site 5th dim hopping term
     *result_cb is parity par of result
     *If result_cb == cb4d then neighbours on source_cb are s-1,s 
     *If result_cb != cb4d then neighbours on source_cb are s,s+1
     *
     * Need negtwomass and two as constants.
     *
     * Apply negtwomass when straddle boundary
     */

    int coeffp=U[2][0];
    int coeffm=U[2][1];

    int two_offset        = PSIIMM[1][0][0];
    int negtwomass_offset = PSIIMM[1][1][0];

    // Base pointer
    queue_iload_short(Chi_mu,TAB_IMM[9],tab);  
    queue_iadd       (Chi_mu,Chi_mu,Chi_p);    

    int cb_zero_branch  = get_target_label();
    int cb_done_branch = get_target_label();

    queue_iload_short(dperp_XX,TAB_IMM[11],tab); 
    check_iterations (dperp_XX,cb_zero_branch);

    ///////////
    // CB==1
    ///////////

    queue_iadd_imm(dperp_ss,dperp_s,one);
    queue_isub(dperp_ss,dperp_ss,Ls);

    complex_load(coeffm,two_offset,Complex_i);
    complex_load(coeffp,negtwomass_offset,Complex_i);

    int wrap = get_target_label();
    conditional_branch_cmpzero ( BRANCH_EQ, dperp_ss, wrap );
      complex_load(coeffp,two_offset,Complex_i);
      queue_iadd_imm(dperp_ss,dperp_s,one);
    make_inst(DIRECTIVE,Target,wrap);

    queue_imul_imm(dperp_XX,dperp_ss,PSI_IMM);  // XX<- s+1 offset
    queue_imul_imm(dperp_ss,dperp_s,PSI_IMM);   // ss<- s-1 offset

    queue_iadd    (dperp_ss,dperp_ss,Chi_mu);  
    queue_iadd    (Chi_mu,dperp_XX,Chi_mu);    
// ss     -> s-1  (checkerboarded-> s) 
// Chi_mu -> s+1  (checkerboarded-> s+1)

    jump(cb_done_branch);
    
    make_inst(DIRECTIVE,Target,cb_zero_branch);

    ///////////
    // CB==0
    ///////////
    int nowrap = get_target_label();

    complex_load(coeffp,two_offset,Complex_i);
    complex_load(coeffm,two_offset,Complex_i);

    queue_iadd_imm(dperp_ss,dperp_s,minusone);
    conditional_branch_cmpzero ( BRANCH_GE, dperp_ss, nowrap );
      complex_load(coeffm,negtwomass_offset,Complex_i);
      queue_iadd    (dperp_ss,dperp_ss,Ls);
    make_inst(DIRECTIVE,Target,nowrap);

    queue_imul_imm(dperp_ss,dperp_ss,PSI_IMM);
    queue_iadd(dperp_ss,dperp_ss,Chi_mu); 

    queue_imul_imm(dperp_XX,dperp_s,PSI_IMM);
    queue_iadd(Chi_mu,dperp_XX,Chi_mu);   

// ss     -> s-1  (checkerboarded-> s-1)
// Chi_mu -> s+1  (checkerboarded-> s)


    make_inst(DIRECTIVE,Target,cb_done_branch);

    {
	int op=2;
	int om=0;
	if ( dagger ) {
	  op=0;
	  om=2;
	}

	complex_load(Psi[0+op][0],PSIIMM[0+op][0][0],Chi_mu,FourSpinType);
	complex_load(Psi[0+om][0],PSIIMM[0+om][0][0],dperp_ss,FourSpinType);

	complex_load(Psi[0+op][2],PSIIMM[0+op][2][0],Chi_mu,FourSpinType);
	complex_load(Psi[0+om][2],PSIIMM[0+om][2][0],dperp_ss,FourSpinType);

	complex_load(Psi[1+op][1],PSIIMM[1+op][1][0],Chi_mu,FourSpinType);
	complex_load(Psi[1+om][1],PSIIMM[1+om][1][0],dperp_ss,FourSpinType);

	complex_load(Psi[0+op][1],PSIIMM[0+op][1][0],Chi_mu,FourSpinType);
	complex_load(Psi[0+om][1],PSIIMM[0+om][1][0],dperp_ss,FourSpinType);

	complex_load(Psi[1+op][0],PSIIMM[1+op][0][0],Chi_mu,FourSpinType);
	complex_load(Psi[1+om][0],PSIIMM[1+om][0][0],dperp_ss,FourSpinType);

	complex_load(Psi[1+op][2],PSIIMM[1+op][2][0],Chi_mu,FourSpinType);
	complex_load(Psi[1+om][2],PSIIMM[1+om][2][0],dperp_ss,FourSpinType);

	for(int sp=0;sp<2;sp++){
	for(int co=0;co<3;co++){
	  complex_real_mul  (Psi[sp+op][co],coeffp,Psi[sp+op][co]);
	  complex_real_mul  (Psi[sp+om][co],coeffm,Psi[sp+om][co]);
	}}
	queue_iadd_imm(dperp_s,dperp_s,one);
    }

  } else { // painfully loading zero... 
    for(co=0;co<3;co++) complex_load(Psi[0][co],PSIIMM[0][2][0],Complex_i);
    for(co=0;co<3;co++) complex_load(Psi[1][co],PSIIMM[0][2][0],Complex_i);
    for(co=0;co<3;co++) complex_load(Psi[2][co],PSIIMM[0][2][0],Complex_i);
    for(co=0;co<3;co++) complex_load(Psi[3][co],PSIIMM[0][2][0],Complex_i);
  }

  /////////////////////////////////////////////////////////
  // Loop over directions
  /////////////////////////////////////////////////////////

  for ( int mu=0;mu<4;mu++ ) {
   for ( int pmdir=0;pmdir<2;pmdir++ ) {

     /* spinor pointer  */
     int pm;
     if ( dagger ) pm = 1-pmdir;
     else          pm = pmdir;

     int dir = pmdir+mu*2;

     queue_iload_short(Chi_mu,TAB_IMM[dir],tab);  
     queue_iadd       (Chi_mu,Chi_mu,Chi_p);    
     queue_iadd       (Chi_mu,Chi_mu,s_offset); // offset for this "s"

     touch(Chi_mu,0);

     // A. Prefetching... can compute chi_nxt also 
     // B. Think about 5d hopping term here.

     ////////////////////////////////////////////////////////
     // Branch structure for interior/exterior neighbours
     ////////////////////////////////////////////////////////
     int lab_skip_mu     = get_target_label(); 
     int lab_proj_mu     = get_target_label(); 
     int lab_su3_mu      = get_target_label(); 

     queue_iand_imm(dperp_ss,mask,mask_imm[pmdir+mu*2]); // non-zero if exterior
     check_iterations(dperp_ss,lab_proj_mu);

     ////////////////////////////////////////////////////////
     // Exterior points are already projected. Just load.
     ////////////////////////////////////////////////////////

     if (do_exterior){
       complex_load(Chi[0][0],PSIIMM[0][0][0],Chi_mu,FourSpinType);
       complex_load(Chi[0][2],PSIIMM[0][2][0],Chi_mu,FourSpinType);
       complex_load(Chi[1][1],PSIIMM[1][1][0],Chi_mu,FourSpinType);
       complex_load(Chi[0][1],PSIIMM[0][1][0],Chi_mu,FourSpinType);
       complex_load(Chi[1][0],PSIIMM[1][0][0],Chi_mu,FourSpinType);
       complex_load(Chi[1][2],PSIIMM[1][2][0],Chi_mu,FourSpinType);
       jump(lab_su3_mu);
     } else {
       jump(lab_skip_mu);
     }

     ////////////////////////////////////////////////////////
     // Interior points are not already projected.
     ////////////////////////////////////////////////////////
     make_inst(DIRECTIVE,Target,lab_proj_mu);
     if ( do_interior ) { 
       
       queue_iadd       (Chi_mu,Chi_mu,s_offset); // offset for this "s"

       complex_load(Chimu[0][0],PSIIMM[0][0][0],Chi_mu,FourSpinType);
       complex_load(Chimu[0][2],PSIIMM[0][2][0],Chi_mu,FourSpinType);
       complex_load(Chimu[1][1],PSIIMM[1][1][0],Chi_mu,FourSpinType);
       complex_load(Chimu[2][0],PSIIMM[2][0][0],Chi_mu,FourSpinType);
       complex_load(Chimu[2][2],PSIIMM[2][2][0],Chi_mu,FourSpinType);
       complex_load(Chimu[3][1],PSIIMM[3][1][0],Chi_mu,FourSpinType);

       complex_load(Chimu[0][1],PSIIMM[0][1][0],Chi_mu,FourSpinType);
       complex_load(Chimu[1][0],PSIIMM[1][0][0],Chi_mu,FourSpinType);
       complex_load(Chimu[1][2],PSIIMM[1][2][0],Chi_mu,FourSpinType);
       complex_load(Chimu[2][1],PSIIMM[2][1][0],Chi_mu,FourSpinType);
       complex_load(Chimu[3][0],PSIIMM[3][0][0],Chi_mu,FourSpinType);
       complex_load(Chimu[3][2],PSIIMM[3][2][0],Chi_mu,FourSpinType);

       /****************************************************************
	* Spin project 4 spinor                                        *
	****************************************************************/
       if ( mu == 0 ) {
	 /* Gx
	  *  0 0  0  i    [0]+-i[3]
	  *  0 0  i  0    [1]+-i[2]
	  *  0 -i 0  0
	  * -i 0  0  0
	  */
	 if ( pm ==0 ) {
	   for(co=0;co<3;co++) complex_ApiB(Chi[0][co],Chimu[0][co],Chimu[3][co]);
	   for(co=0;co<3;co++) complex_ApiB(Chi[1][co],Chimu[1][co],Chimu[2][co]);
	 } else {
	   for(co=0;co<3;co++) complex_AmiB(Chi[0][co],Chimu[0][co],Chimu[3][co]);
	   for(co=0;co<3;co++) complex_AmiB(Chi[1][co],Chimu[1][co],Chimu[2][co]);
	 }
       } else if ( mu == 1 ) {
	 /*Gy
	  *  0 0  0  -1  [0] -+ [3]
	  *  0 0  1  0   [1] +- [2]
	  *  0 1  0  0
	  * -1 0  0  0
	  */
	 if ( pm ==0 ) {
	   for(co=0;co<3;co++) complex_sub(Chi[0][co],Chimu[0][co],Chimu[3][co]);
	   for(co=0;co<3;co++) complex_add(Chi[1][co],Chimu[1][co],Chimu[2][co]);
	 } else {
	   for(co=0;co<3;co++) complex_add(Chi[0][co],Chimu[0][co],Chimu[3][co]);
	   for(co=0;co<3;co++) complex_sub(Chi[1][co],Chimu[1][co],Chimu[2][co]);
	 }
       } else if ( mu == 2 ) {
	 /*Gz
	  *  0 0  i  0   [0]+-i[2]
	  *  0 0  0 -i   [1]-+i[3]
	  * -i 0  0  0
	  *  0 i  0  0
	  */
	 if ( pm ==0 ) {
	   for(co=0;co<3;co++) complex_ApiB(Chi[0][co],Chimu[0][co],Chimu[2][co]);
	   for(co=0;co<3;co++) complex_AmiB(Chi[1][co],Chimu[1][co],Chimu[3][co]);
	 } else {
	   for(co=0;co<3;co++) complex_AmiB(Chi[0][co],Chimu[0][co],Chimu[2][co]);
	   for(co=0;co<3;co++) complex_ApiB(Chi[1][co],Chimu[1][co],Chimu[3][co]);
	 }
       } else if ( mu == 3 ) {
	 /*Gt
	  *  0 0  1  0 [0]+-[2]
	  *  0 0  0  1 [1]+-[3]
	  *  1 0  0  0
	  *  0 1  0  0
	  */
	 if ( pm ==0 ) {
	   for(co=0;co<3;co++) complex_add(Chi[0][co],Chimu[0][co],Chimu[2][co]);
	   for(co=0;co<3;co++) complex_add(Chi[1][co],Chimu[1][co],Chimu[3][co]);
	 } else {
	   for(co=0;co<3;co++) complex_sub(Chi[0][co],Chimu[0][co],Chimu[2][co]); 
	   for(co=0;co<3;co++) complex_sub(Chi[1][co],Chimu[1][co],Chimu[3][co]);
	 }
       } 
     } else { 
       jump(lab_skip_mu);
     }

     make_inst(DIRECTIVE,Target,lab_su3_mu);

     /************************************************************
      * SU3 multiply.
      ************************************************************
      */
     complex_load(U[0][0],GIMM[0][0][0],U_p_s,GaugeType);
     complex_load(U[0][2],GIMM[0][2][0],U_p_s,GaugeType);
     complex_load(U[1][1],GIMM[1][1][0],U_p_s,GaugeType);
     complex_load(U[2][0],GIMM[2][0][0],U_p_s,GaugeType);

     complex_load(U[0][1],GIMM[0][1][0],U_p_s,GaugeType);
     complex_load(U[1][0],GIMM[1][0][0],U_p_s,GaugeType);
     complex_load(U[2][1],GIMM[2][1][0],U_p_s,GaugeType);

     j=0;
     complex_six_cmuls (UChi[0][0],U[0][j],Chi[0][j],
			UChi[0][1],U[1][j],Chi[0][j],
			UChi[0][2],U[2][j],Chi[0][j],
			UChi[1][0],U[0][j],Chi[1][j],
			UChi[1][1],U[1][j],Chi[1][j],
			UChi[1][2],U[2][j],Chi[1][j]);

     j=1;
     complex_load(U[1][2],GIMM[1][2][0],U_p_s,GaugeType);
     complex_load(U[2][2],GIMM[2][2][0],U_p_s,GaugeType);

     complex_six_cmadds (UChi[0][0],U[0][j],Chi[0][j],
			 UChi[0][1],U[1][j],Chi[0][j],
			 UChi[0][2],U[2][j],Chi[0][j],
			 UChi[1][0],U[0][j],Chi[1][j],
			 UChi[1][1],U[1][j],Chi[1][j],
			 UChi[1][2],U[2][j],Chi[1][j]);
     j=2;
     complex_six_cmaddtos (UChiR[0][0],UChi[0][0],U[0][j],Chi[0][j],
			   UChiR[0][1],UChi[0][1],U[1][j],Chi[0][j],
			   UChiR[0][2],UChi[0][2],U[2][j],Chi[0][j],
			   UChiR[1][0],UChi[1][0],U[0][j],Chi[1][j],
			   UChiR[1][1],UChi[1][1],U[1][j],Chi[1][j],
			   UChiR[1][2],UChi[1][2],U[2][j],Chi[1][j]);


     /*****************************************************************
      * Reconstruct / accumulate into Psi
      *****************************************************************
      */

     /*Upper two components are simple, by virtue of our choice of projection*/      
     for(co=0;co<3;co++) complex_add(Psi[0][co],Psi[0][co],UChiR[0][co]);
     for(co=0;co<3;co++) complex_add(Psi[1][co],Psi[1][co],UChiR[1][co]);
      
     /*Lower two components are complex, by virtue of our choice of projection*/      
     if ( (mu==0) && (pm==0) ) for(co=0;co<3;co++) complex_AmiB(Psi[2][co],Psi[2][co],UChiR[1][co]);
     if ( (mu==0) && (pm==1) ) for(co=0;co<3;co++) complex_ApiB(Psi[2][co],Psi[2][co],UChiR[1][co]);
     if ( (mu==1) && (pm==0) ) for(co=0;co<3;co++) complex_add (Psi[2][co],Psi[2][co],UChiR[1][co]);
     if ( (mu==1) && (pm==1) ) for(co=0;co<3;co++) complex_sub (Psi[2][co],Psi[2][co],UChiR[1][co]);
     if ( (mu==2) && (pm==0) ) for(co=0;co<3;co++) complex_AmiB(Psi[2][co],Psi[2][co],UChiR[0][co]);
     if ( (mu==2) && (pm==1) ) for(co=0;co<3;co++) complex_ApiB(Psi[2][co],Psi[2][co],UChiR[0][co]);
     if ( (mu==3) && (pm==0) ) for(co=0;co<3;co++) complex_add (Psi[2][co],Psi[2][co],UChiR[0][co]);
     if ( (mu==3) && (pm==1) ) for(co=0;co<3;co++) complex_sub (Psi[2][co],Psi[2][co],UChiR[0][co]);

     if ( (mu==0) && (pm==0) ) for(co=0;co<3;co++) complex_AmiB(Psi[3][co],Psi[3][co],UChiR[0][co]);
     if ( (mu==0) && (pm==1) ) for(co=0;co<3;co++) complex_ApiB(Psi[3][co],Psi[3][co],UChiR[0][co]);
     if ( (mu==1) && (pm==0) ) for(co=0;co<3;co++) complex_sub (Psi[3][co],Psi[3][co],UChiR[0][co]);
     if ( (mu==1) && (pm==1) ) for(co=0;co<3;co++) complex_add (Psi[3][co],Psi[3][co],UChiR[0][co]);
     if ( (mu==2) && (pm==0) ) for(co=0;co<3;co++) complex_ApiB(Psi[3][co],Psi[3][co],UChiR[1][co]);
     if ( (mu==2) && (pm==1) ) for(co=0;co<3;co++) complex_AmiB(Psi[3][co],Psi[3][co],UChiR[1][co]);
     if ( (mu==3) && (pm==0) ) for(co=0;co<3;co++) complex_add (Psi[3][co],Psi[3][co],UChiR[1][co]);
     if ( (mu==3) && (pm==1) ) for(co=0;co<3;co++) complex_sub (Psi[3][co],Psi[3][co],UChiR[1][co]);

     if ( (do_exterior==0) || (do_interior==0) ){
       make_inst(DIRECTIVE,Target,lab_skip_mu);
     }
     queue_load_addr(U_p_s,MAT_IMM,U_p_s); 

   }      /*END MU/PM loops*/
  }
     
  /* Store result Psi & Pointer update*/    
  for(sp=0;sp<4;sp++){
    for(co=0;co<3;co++){
      complex_store(Psi[sp][co],PSIIMM[sp][co][0],Psi_p,FourSpinType);
    }
  }
  queue_iadd_imm(Psi_p,Psi_p,PSI_IMM);
  queue_iadd_imm(s_offset,s_offset,CHI_IMM);

  stop_loop(branchls,s); 

#ifdef LOCK_GAUGE
  if ( do_interior ) {
  if ( GaugeType == Double ) {
  queue_load_addr(s_offset,ZERO_IMM,U_p); // unlock from s_offset
  for(int i=0;i<5;i++) { // 9*8 32 bytes => 36 clines => 6*6
    unlock(s_offset,0);  // Unlock our old U's.
    unlock(s_offset,1);  
    unlock(s_offset,2);  
      unlock(s_offset,3);  
      unlock(s_offset,4);  
      unlock(s_offset,5); 
    }
    queue_load_addr(s_offset,PSI_IMM,s_offset); 
  unlock(s_offset,0); 
  unlock(s_offset,1); 
  unlock(tab,0); 
  }
  }
#endif

  queue_load_addr(s_offset,ZERO_IMM,U_p);
  queue_iadd_imm(U_p,U_p_s,ZERO_IMM);
  queue_load_addr(tab,TAB_IMM[16],tab);
  touch(U_p,0);
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



