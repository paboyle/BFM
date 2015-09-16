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

void dwf_deriv( char *);
void set_label_fmt(char *);
/*Options flags*/
int human = 0;
char name[80];
int dagger = 0;
int dd2 = 0;

void do_writehint(int,int);

int LOCK_GAUGE=0;
#define do_WH(A,B) ;
#define   lock l1_lock
#define unlock do_flush

void jump(int where);
void jump(int where)
{
  make_inst(BRCHPIPE,BRANCH,where); 
}
void touch(int addr, int line);
void touch(int addr, int line)
{
      do_flush(addr,line);
      if ( dd2 ) l2_touch(addr,line);
      else do_prefetch(addr,line);
}


Datum ConstantType=Double;
Datum FourSpinType=Double;
Datum TableType=ShortInteger;
Datum GaugeType=Double;
char procname[80]="UNSPECIFIED";

int main ( int argc, char *argv[])
{
  struct processor *myproc;
  int arg;
  char *c;

  name[0] = '\0';
  /*Process command line*/
  while ( ( arg = getopt(argc,argv,"aRn:P:dskfiloerpx")) != EOF){
    switch(arg){
    case 'R': human = 1; break;
    case 'n': if (strlen(optarg)<30) strcpy(name,optarg); break;
    case 'P': if (strlen(optarg)<30) strcpy(procname,optarg); break;
    case 'd': dagger = 1; break;
    case 'f': dd2 = 1; break;
    case 's': FourSpinType=Single; GaugeType=Single;  break;
    default: fprintf(stderr,"Usage: %s -[Rd] -n routine_name\n",argv[0]); 
             fprintf(stderr,"\tR -> human readable .S format \n"); 
             fprintf(stderr,"\td -> dagger decompose\n"); 
             exit (1); break;
    }
  }
  /*Control output according to user set up args*/
  streams_reg_alloc(0);
  set_processor_optarg(procname);
  set_human_readable(human);
  setup_cmadds();

  /*Write the naive asm code*/
  dwf_deriv(name);

  /*Filter through an virtual out of order processor */
  schedule_for_proc();

  /*Dump the resulting code*/
  dump_instruction_queue();

  return(0);
}

void dwf_deriv( char *name)
{
  /*
   * This marks the argument registers as defined by ABI as off limits
   * to us until they are freed by "getarg()";
   */
  int dum = defargcount(1);
  int retno;

  /*
   * S=phi^dag (MdagM)^-1 phi
   * 
   * dS = phi^dag (MdagM)^-1 [ dMdag M + Mdag dM ] (MdagM)^-1 phi
   *
   * Let X = (MdagM)^-1 phi
   *     Y = M X = M^-dag phi
   *
   * Want terms:    Ydag dM X
   *                Xdag dMdag Y 
   * 
   * Take Xdag 1-gamma Y
   *
   * Still a bit confused about the 1+g 1-g terms; but this may be simply a factor of two as we add +h.c.
   * Will continue to follow Chroma's routine
   */
  reg_array_2d(Y,Cregs,4,3);  // 4 spinor - 24 regs 
  reg_array_2d(X,Cregs,4,3);  // 4 spinor - 12 regs 
  reg_array_1d(F,Cregs,3);    // Force
  alreg(Z,Cregs);             // Zero
  alreg(creg,Cregs);             // Zero

  offset_3d(CHIIMM,FourSpinType,2,3,2*nsimd());
  offset_3d(PSIIMM,FourSpinType,4,3,2*nsimd());
  offset_3d(GIMM  ,GaugeType, 3, 3 ,2*nsimd() );

  def_off( GAUGE_SITE_IMM, FourSpinType,4*18*nsimd());
  def_off( MAT_IMM  , GaugeType,18*nsimd());
  def_off( PSI_IMM  , FourSpinType,24*nsimd());
  def_off( CHI_IMM  , FourSpinType,12*nsimd());
  def_off( CONST_ZERO_OFFSET,Double,2*2*nsimd());

  /*
   * Integer registers
   */
  alreg(F_p,Iregs);  /*Pointer to the current cpt of force field    */
  alreg(F_p_s,Iregs);
  alreg(Y_mu,Iregs);  
  alreg(Y_p,Iregs);  
  alreg(X_p,Iregs);  

  alreg(length,Iregs);   /*number of sites*/
  alreg(tab,Iregs);      /*Pointer to current entry in offset table*/
  alreg(Complex_i,Iregs);/*Point to (0,1)x Nsimd*/
  alreg(Ls,Iregs);        
  alreg(s,Iregs);        

  alreg(recbuf_base,Iregs);
  alreg(args,Iregs);        
  alreg(s_offset,Iregs);        
  /*Useful integer immediate constants, in units of Fsize*/
  def_off( ZERO_IMM,Byte,0);
  def_off( minusone,Byte,-1);
  def_off( one,Byte,1);

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
  alreg(mask ,Iregs);

  offset_1d(TAB_IMM,TableType,17);

  // Integer sizes
  int Isize      = def_offset(PROC->I_size,Byte,"Isize");
  int ISsize     = def_offset(PROC->IS_size,Byte,"ISsize");
  int i,j,co,sp;

  /*********************************************************************/

  make_inst(DIRECTIVE,Enter_Routine,name);
  grab_stack(0);
  save_regs();

  /*********************************************
   * our arguments 
   *********************************************
   */
  getarg(args); /*Pointer to arg list*/

  queue_iload(X_p, ZERO_IMM,args);   queue_load_addr(args,Isize,args);   //0
  queue_iload(Y_p, ZERO_IMM,args);   queue_load_addr(args,Isize,args);   //1
  queue_iload(F_p,   ZERO_IMM,args);   queue_load_addr(args,Isize,args);   //2
  queue_iload(length,ZERO_IMM,args);   queue_load_addr(args,Isize,args);   //3
  queue_iload(Ls,    ZERO_IMM,args);   queue_load_addr(args,Isize,args);   //4
  queue_iload(tab,   ZERO_IMM,args);   queue_load_addr(args,Isize,args);   //5
  queue_iload(Complex_i,ZERO_IMM,args);queue_load_addr(args,Isize,args);   //6
                                       queue_load_addr(args,Isize,args);   //7  
  queue_iload(recbuf_base,ZERO_IMM,args);queue_load_addr(args,Isize,args); //8  

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
  complex_load(Z,CONST_ZERO_OFFSET,Complex_i,Double);

  // Site loop 
  retno = get_target_label(); 
  check_iterations(length,retno); 
  int branchsite = start_loop(length);
  
  // S loop
  queue_iload_short(mask,TAB_IMM[10],tab);
  queue_iadd_imm (s,Ls,ZERO_IMM);
  queue_iload_imm(s_offset,ZERO_IMM);        
  int branchls   = start_loop(s); 
  queue_iadd_imm(F_p_s,F_p,ZERO_IMM);
  //  debugI(s);

  // Loop over directions
  for ( int mu=0;mu<4;mu++ ) {

    int dir = mu*2+1;  // Always in forward dir

     // Complex branch structure for interior/exterior neighbours
     int lab_proj_mu     = get_target_label(); 
     int lab_continue    = get_target_label(); 
     queue_iand_imm  (Y_mu,mask,mask_imm[dir]); // non-zero if exterior
     check_iterations(Y_mu,lab_proj_mu);

     // Exterior points are already projected. Just load.
       queue_iload_short(Y_mu,TAB_IMM[dir],tab);  
       queue_iadd       (Y_mu,Y_mu,recbuf_base);  
       //       debugI(Y_mu);
       //debugI(recbuf_base);
       queue_iadd       (Y_mu,Y_mu,s_offset); 
       for(int sp=0;sp<2;sp++){
	 for(int co=0;co<3;co++){
	   complex_load(Y[sp][co],PSIIMM[sp][co][0],Y_mu,FourSpinType);
	 }
       }

  jump(lab_continue);
  make_inst(DIRECTIVE,Target,lab_proj_mu); 

      // Interior points are not already projected.
      // * Spin project 4 spinor                                        

     queue_iload_short(Y_mu,TAB_IMM[dir],tab);  
     //     debugI(tab);
     //     debugI(Y_mu);
     queue_iadd       (Y_mu,Y_mu,Y_p);    
     queue_iadd       (Y_mu,Y_mu,s_offset); // offset for this "s"
     queue_iadd       (Y_mu,Y_mu,s_offset); // offset for this "s"

     for(int sp=0;sp<4;sp++){
       for(int co=0;co<3;co++){
	 complex_load(X[sp][co],PSIIMM[sp][co][0],Y_mu,FourSpinType);
	 //	 debugC(X[sp][co]);
       }
     }

    int pm = 1;       // pm=0 == 1+gamma, pm=1 => 1-gamma
    if ( dagger ) pm = 0;

     if ( mu == 0 ) {
	 if ( pm ==0 ) {
	   for(co=0;co<3;co++) complex_ApiB(Y[0][co],X[0][co],X[3][co]);
	   for(co=0;co<3;co++) complex_ApiB(Y[1][co],X[1][co],X[2][co]);
	 } else {
	   for(co=0;co<3;co++) complex_AmiB(Y[0][co],X[0][co],X[3][co]);
	   for(co=0;co<3;co++) complex_AmiB(Y[1][co],X[1][co],X[2][co]);
	 }
       } else if ( mu == 1 ) {
	 if ( pm ==0 ) {
	   for(co=0;co<3;co++) complex_sub(Y[0][co],X[0][co],X[3][co]);
	   for(co=0;co<3;co++) complex_add(Y[1][co],X[1][co],X[2][co]);
	 } else {
	   for(co=0;co<3;co++) complex_add(Y[0][co],X[0][co],X[3][co]);
	   for(co=0;co<3;co++) complex_sub(Y[1][co],X[1][co],X[2][co]);
	 }
       } else if ( mu == 2 ) {
	 if ( pm ==0 ) {
	   for(co=0;co<3;co++) complex_ApiB(Y[0][co],X[0][co],X[2][co]);
	   for(co=0;co<3;co++) complex_AmiB(Y[1][co],X[1][co],X[3][co]);
	 } else {
	   for(co=0;co<3;co++) complex_AmiB(Y[0][co],X[0][co],X[2][co]);
	   for(co=0;co<3;co++) complex_ApiB(Y[1][co],X[1][co],X[3][co]);
	 }
       } else if ( mu == 3 ) {
	 if ( pm ==0 ) {
	   for(co=0;co<3;co++) complex_add(Y[0][co],X[0][co],X[2][co]);
	   for(co=0;co<3;co++) complex_add(Y[1][co],X[1][co],X[3][co]);
	 } else {
	   for(co=0;co<3;co++) complex_sub(Y[0][co],X[0][co],X[2][co]); 
	   for(co=0;co<3;co++) complex_sub(Y[1][co],X[1][co],X[3][co]);
	 }
     }

  make_inst(DIRECTIVE,Target,lab_continue);

       ///////////////////////////////////////////////////////////////
       // Y contains spin projection of forward neighbour in mu direction
       // Repromote to Y to 4 spinor
       ///////////////////////////////////////////////////////////////
       for(int co_y=0;co_y<3;co_y++){
	 
	 if ( (mu==0) && (pm==0) )  complex_AmiB(Y[2][co_y],Z,Y[1][co_y]);
	 if ( (mu==0) && (pm==1) )  complex_ApiB(Y[2][co_y],Z,Y[1][co_y]);
	 if ( (mu==1) && (pm==0) )  complex_add (Y[2][co_y],Z,Y[1][co_y]);
	 if ( (mu==1) && (pm==1) )  complex_sub (Y[2][co_y],Z,Y[1][co_y]);
	 if ( (mu==2) && (pm==0) )  complex_AmiB(Y[2][co_y],Z,Y[0][co_y]);
	 if ( (mu==2) && (pm==1) )  complex_ApiB(Y[2][co_y],Z,Y[0][co_y]);
	 if ( (mu==3) && (pm==0) )  complex_add (Y[2][co_y],Z,Y[0][co_y]);
	 if ( (mu==3) && (pm==1) )  complex_sub (Y[2][co_y],Z,Y[0][co_y]);

	 if ( (mu==0) && (pm==0) ) complex_AmiB(Y[3][co_y],Z,Y[0][co_y]);
	 if ( (mu==0) && (pm==1) ) complex_ApiB(Y[3][co_y],Z,Y[0][co_y]);
	 if ( (mu==1) && (pm==0) ) complex_sub (Y[3][co_y],Z,Y[0][co_y]);
	 if ( (mu==1) && (pm==1) ) complex_add (Y[3][co_y],Z,Y[0][co_y]);
	 if ( (mu==2) && (pm==0) ) complex_ApiB(Y[3][co_y],Z,Y[1][co_y]);
	 if ( (mu==2) && (pm==1) ) complex_AmiB(Y[3][co_y],Z,Y[1][co_y]);
	 if ( (mu==3) && (pm==0) ) complex_add (Y[3][co_y],Z,Y[1][co_y]);
	 if ( (mu==3) && (pm==1) ) complex_sub (Y[3][co_y],Z,Y[1][co_y]);

       }

       ///////////////////////////////////////////////////////////////
       // Load X
       ///////////////////////////////////////////////////////////////
       for(int co_x=0;co_x<3;co_x++){
	 for(int sp=0;sp<4;sp++) {
	   complex_load(X[sp][co_x],PSIIMM[sp][co_x][0],X_p,FourSpinType);
	 }
       }

       ///////////////////////////////////////////////////////////////
       // Spin trace tensor product
       ///////////////////////////////////////////////////////////////
       for(int co_x=0;co_x<3;co_x++){
	 // Spin trace outer product
	 for ( int co_y=0;co_y<3;co_y++) complex_load (F[co_y],GIMM[co_y][co_x][0],F_p_s);  
	 for ( int co_y=0;co_y<3;co_y++) complex_conjmadd(F[co_y],X[0][co_x],Y[0][co_y]);  
	 for ( int co_y=0;co_y<3;co_y++) complex_conjmadd(F[co_y],X[1][co_x],Y[1][co_y]);  
	 for ( int co_y=0;co_y<3;co_y++) complex_conjmadd(F[co_y],X[2][co_x],Y[2][co_y]);  
	 for ( int co_y=0;co_y<3;co_y++) complex_conjmadd(F[co_y],X[3][co_x],Y[3][co_y]);  
	 
	 for ( int co_y=0;co_y<3;co_y++) complex_store(F[co_y],GIMM[co_y][co_x][0],F_p_s);  
	 
       }

     queue_load_addr(F_p_s,MAT_IMM,F_p_s); 

  }
  queue_iadd_imm(X_p,X_p,PSI_IMM);
  queue_iadd_imm(s_offset,s_offset,CHI_IMM);
  stop_loop(branchls,s); 
  queue_iadd_imm(F_p,F_p_s,ZERO_IMM);
  queue_load_addr(tab,TAB_IMM[16],tab);
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



