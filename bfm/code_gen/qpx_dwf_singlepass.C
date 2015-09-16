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

void do_writehint(int,int);
#define do_WH(A,B) ;

Datum FourSpinType=Double;
Datum GaugeType=Double;
int addto = 0;
char procname[80]="UNSPECIFIED";

int main ( int argc, char *argv[])
{
  struct processor *myproc;
  int arg;
  char *c;

  name[0] = '\0';
  /*Process command line*/
  while ( ( arg = getopt(argc,argv,"aRn:P:ds")) != EOF){
    switch(arg){
    case 'R': human = 1; break;
    case 'n': if (strlen(optarg)<30) strcpy(name,optarg); break;
    case 'P': if (strlen(optarg)<30) strcpy(procname,optarg); break;
    case 'd': dagger = 1; break;
    case 'a': addto = 1; break;
    case 's': FourSpinType=Single; GaugeType=Single; break;
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
  dwf_dslash(name);

  /*Filter through an virtual out of order processor */
  //schedule_for_proc();

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
  alreg(Chi_nxt,Iregs);/*Pointer to the input four spinor            */
  alreg(length,Iregs);  /*number of sites*/
  alreg(tab,Iregs);     /*Pointer to current entry in offset table*/
  alreg(tab_s,Iregs);   /*Pointer to current entry in offset table*/
  alreg(Complex_i,Iregs);  /*Point to (0,1)x Nsimd*/
  alreg(Ls,Iregs);        
  alreg(s,Iregs);        
  alreg(ss,Iregs);        
  alreg(s_offset,Iregs);        
  alreg(args,Iregs);        
  int ptr;

  /*Useful integer immediate constants, in units of Fsize*/
  def_off( ZERO_IMM,Byte,0);
  def_off( m1,Byte,-1);
  def_off( m64,Byte,-64);
  def_off( PSI_IMM, FourSpinType,24*nsimd());
  def_off( CHI_IMM, FourSpinType,12*nsimd());
  def_off( MAT_IMM, GaugeType,18*nsimd());
  offset_3d(GIMM    , GaugeType, 3, 3 ,2*nsimd() );

  int Isize      = def_offset(PROC->I_size,Byte,"Isize");
  int i,j,co,sp;
  /*********************************************************************/

  /*--------------------------------------------------------------------
   * Start of the "pseudo assembler proper.
   *--------------------------------------------------------------------
   */
  make_inst(DIRECTIVE,Enter_Routine,name);

  grab_stack(0);
  save_regs();

  /*
   * Define our arguments 
   */
  getarg(args); /*Pointer to arg list*/

  queue_iload(Psi_p, ZERO_IMM,args);  queue_load_addr(args,Isize,args);  
  queue_iload(Chi_p, ZERO_IMM,args);  queue_load_addr(args,Isize,args);  
  queue_iload(U_p,   ZERO_IMM,args);  queue_load_addr(args,Isize,args);  
  queue_iload(length,ZERO_IMM,args);  queue_load_addr(args,Isize,args);  
  queue_iload(Ls,    ZERO_IMM,args);  queue_load_addr(args,Isize,args);  
  queue_iload(tab,   ZERO_IMM,args);  queue_load_addr(args,Isize,args);  
  queue_iload(Complex_i,ZERO_IMM,args);   
  int pref = args;


  /*
   * Load common constants into Iregs
   */
  /*Bugger 12 int regs for offsets as indexed loads only*/
  for (int i =0; i<12; i++ ) { 
    need_constant(i*2*SizeofDatum(FourSpinType)*nsimd());
  }
  for (int i =0; i<9; i++ ) { 
    need_constant(i*2*SizeofDatum(GaugeType)*nsimd());
  }
  complex_constants_prepare(creg,Complex_i);

  retno = get_target_label(); /*Branch to exit if length <1*/
  check_iterations(length,retno); 


  /*Initial locking phase : i, table, gauge*/
  //  l1_lock(Complex_i,0); // Unlocked in epilogue
  queue_iadd_imm(U_p_s,U_p,ZERO_IMM);

  /*Initial 4-spinor preload*/
  queue_iload    (Chi_mu,ZERO_IMM,tab);  
  queue_imul     (Chi_mu,Chi_mu,Ls);       // These can be eliminated
  queue_imul_imm (Chi_mu,Chi_mu,PSI_IMM);  // By premultiplying table
  queue_iadd     (Chi_mu,Chi_mu,Chi_p);    

  for(sp=0;sp<4;sp++ ) {
    for(co=0;co<3;co++){ 
      complex_load(Chimu[sp][co],PSIIMM[sp][co][0],Chi_mu,FourSpinType);
    }
  }

  /*
   * Site loop & s-loop
   */

  int branchsite = start_loop(length);


  queue_load_addr(s_offset,ZERO_IMM,U_p); // lock from new U_p
  for(int i=0;i<5;i++) { 
    l1_lock(s_offset,0); // lock next set of links (32 lines, 2304 bytes)
    l1_lock(s_offset,1); // Exactly 18*2*8*8 = 2304
    l1_lock(s_offset,2); 
    l1_lock(s_offset,3); 
    l1_lock(s_offset,4); 
    l1_lock(s_offset,5); 
    queue_load_addr(s_offset,PSI_IMM,s_offset); 
  }
  l1_lock(s_offset,0); 
  l1_lock(s_offset,1); 

  queue_iadd_imm(s,Ls,ZERO_IMM);
  queue_iload_imm(s_offset,ZERO_IMM);        

  l1_lock(tab,0); // 8 * 8 bytes per site; aligned => 1 line

  int branchls   = start_loop(s); 

  pragma(DCBT_SPACE,0);
  //  pragma(LOAD_LIM,1);
  //  complex_constants_assert(cmplx_i);
  queue_iadd_imm(U_p_s,U_p,ZERO_IMM);
  queue_iadd_imm(tab_s,tab,ZERO_IMM);
  
  for ( int mu=0;mu<4;mu++ ) {
   for ( int pmdir=0;pmdir<2;pmdir++ ) {

    // pm is the sign for the gamma matrix
    // pmdir is the sign for the spatial displacement (0=plus, 1=minus)
    // Matrix is 1-gammamu Psi(x+mu)
    // So in non-dagger case, first element in two spinor is ChiMinus
    int pm;
    if ( dagger ) pm = 1-pmdir;
    else          pm = pmdir;
  
    /***************************
     * Compute and prefetch chinxt
     * Raise even higher?
     ***************************
     */
    queue_load_addr(tab_s,Isize,tab_s); // increment table inside the mu/pm loop
    if ((mu==3) && (pmdir==1)){
      queue_iadd_imm(ss,s,m1); // ss=s-1 (ls-1...0)
                               // s < Ls-1 : revert back to same 4d site, bump s_offset
      queue_cmovgt(ss,tab,tab_s);
      queue_iadd_imm (s_offset,s_offset,PSI_IMM);
                                    // s == Ls-1 : use s_offset = 0 & let tab_s advance
      queue_cmovle(ss,ss,s_offset); // ss==0, move ss in.
    }
    queue_iload    (Chi_nxt,ZERO_IMM,tab_s);  
    queue_imul     (Chi_nxt,Chi_nxt,Ls);         
    queue_imul_imm (Chi_nxt,Chi_nxt,PSI_IMM);    
    queue_iadd     (Chi_nxt,Chi_nxt,Chi_p);    
    queue_iadd     (Chi_nxt,Chi_nxt,s_offset); // offset for this "s"

    /****************************************
     * First prefetch kicked off (psi or chi_mu)
     ****************************************
     */

    complex_load(Chimu[2][1],PSIIMM[2][1][0],Chi_mu,FourSpinType);
    complex_load(Chimu[3][0],PSIIMM[3][0][0],Chi_mu,FourSpinType);
    complex_load(Chimu[3][2],PSIIMM[3][2][0],Chi_mu,FourSpinType);

    do_flush(Chi_nxt,0);
    if ( (mu==0) && (pmdir==0) && addto ) do_prefetch(Psi_p,0);
    else do_prefetch(Chi_nxt,0);

    //Try to only have one outstanding long latency load per thread on BGQ
    //leave lmq slots for progress in hits in l1p

    /**************************
     * Begin load gauge fields
     * Even only as LMQ conflicts
     **************************
     */
    complex_load(U[0][0],GIMM[0][0][0],U_p_s,GaugeType);
    complex_load(U[0][2],GIMM[0][2][0],U_p_s,GaugeType);
    complex_load(U[1][1],GIMM[1][1][0],U_p_s,GaugeType);
    complex_load(U[2][0],GIMM[2][0][0],U_p_s,GaugeType);
    queue_iadd_imm(Chi_mu,Chi_nxt,ZERO_IMM);

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

      /************************************************************
       * SU3 multiply. Collect odd loads for U
       ************************************************************
       */
#undef DEBUG
#ifdef DEBUG
      debugI(Chi_mu);
      for(sp=0;sp<2;sp++)
      for(co=0;co<3;co++)
	debugC(Chi[sp][co]);
#endif

      complex_load(U[0][1],GIMM[0][1][0],U_p_s,GaugeType);
      complex_load(U[1][0],GIMM[1][0][0],U_p_s,GaugeType);
      complex_load(U[2][1],GIMM[2][1][0],U_p_s,GaugeType);

      for(j=0;j<3;j++){
	if ( j==1 ){ // Couldn't load earlier due to reg pressure
	  complex_load(U[1][2],GIMM[1][2][0],U_p_s,GaugeType);
	  complex_load(U[2][2],GIMM[2][2][0],U_p_s,GaugeType);
	}
	if ( j==2) {
	  if ( (mu == 0 ) && (pmdir == 0) && (addto) ) {
	    complex_load(Psi[0][2],PSIIMM[0][2][0],Psi_p,FourSpinType);
	    complex_load(Psi[1][1],PSIIMM[1][1][0],Psi_p,FourSpinType);
	    complex_load(Psi[2][0],PSIIMM[2][0][0],Psi_p,FourSpinType);
	    complex_load(Psi[2][2],PSIIMM[2][2][0],Psi_p,FourSpinType);
	    complex_load(Psi[3][1],PSIIMM[3][1][0],Psi_p,FourSpinType);
	    complex_load(Psi[0][0],PSIIMM[0][0][0],Psi_p,FourSpinType);
	  }
	}
	if ( j==0 ){
	  complex_six_cmuls (UChi[0][0],U[0][j],Chi[0][j],
			     UChi[0][1],U[1][j],Chi[0][j],
			     UChi[0][2],U[2][j],Chi[0][j],
			     UChi[1][0],U[0][j],Chi[1][j],
			     UChi[1][1],U[1][j],Chi[1][j],
			     UChi[1][2],U[2][j],Chi[1][j]);
	} else { 
	  if (j==2) {

	    if ( (mu != 0 ) || (pmdir != 0) ) {
	      /********************************************************
	       *Begin load of Chi_mu for next iter - could lift these
	       ********************************************************/
	      complex_load(Chimu[1][0],PSIIMM[1][0][0],Chi_mu,FourSpinType);
	      complex_load(Chimu[1][1],PSIIMM[1][1][0],Chi_mu,FourSpinType);
	    } else {
	      do_prefetch(Chi_mu,0);
	    }
	    
	    complex_six_cmaddtos (UChiR[0][0],UChi[0][0],U[0][j],Chi[0][j],
				UChiR[0][1],UChi[0][1],U[1][j],Chi[0][j],
				UChiR[0][2],UChi[0][2],U[2][j],Chi[0][j],
				UChiR[1][0],UChi[1][0],U[0][j],Chi[1][j],
				UChiR[1][1],UChi[1][1],U[1][j],Chi[1][j],
				UChiR[1][2],UChi[1][2],U[2][j],Chi[1][j]);
	  } else {
	    complex_six_cmadds (UChi[0][0],U[0][j],Chi[0][j],
				UChi[0][1],U[1][j],Chi[0][j],
				UChi[0][2],U[2][j],Chi[0][j],
				UChi[1][0],U[0][j],Chi[1][j],
				UChi[1][1],U[1][j],Chi[1][j],
				UChi[1][2],U[2][j],Chi[1][j]);
	  }
	}
      }

#ifdef DEBUG
      for(sp=0;sp<2;sp++)
      for(co=0;co<3;co++)
	debugC(UChiR[sp][co]);
#endif
      /*****************************
       *Psi loading
       *****************************/    
      if ( (mu == 0 ) && (pmdir == 0) ) {
	if ( addto ) {
	  complex_load(Psi[0][1],PSIIMM[0][1][0],Psi_p,FourSpinType);
	  complex_load(Psi[1][0],PSIIMM[1][0][0],Psi_p,FourSpinType);
	  complex_load(Psi[1][2],PSIIMM[1][2][0],Psi_p,FourSpinType);
	  complex_load(Psi[2][1],PSIIMM[2][1][0],Psi_p,FourSpinType);
	  complex_load(Psi[3][0],PSIIMM[3][0][0],Psi_p,FourSpinType);
	  complex_load(Psi[3][2],PSIIMM[3][2][0],Psi_p,FourSpinType);
	} else { // why loading zero...
	  for(co=0;co<3;co++) complex_load(Psi[0][co],PSIIMM[0][2][0],Complex_i);
	  for(co=0;co<3;co++) complex_load(Psi[1][co],PSIIMM[0][2][0],Complex_i);
	  for(co=0;co<3;co++) complex_load(Psi[2][co],PSIIMM[0][2][0],Complex_i);
	  for(co=0;co<3;co++) complex_load(Psi[3][co],PSIIMM[0][2][0],Complex_i);
	}
	/********************************************************
	 *Begin load of Chi_mu for next iter - could lift these
	 ********************************************************/
	complex_load(Chimu[1][0],PSIIMM[1][0][0],Chi_mu,FourSpinType);
	complex_load(Chimu[1][1],PSIIMM[1][1][0],Chi_mu,FourSpinType);
      }

      

      /***********************
       * Gauge pointer update
       ************************
       */
      queue_load_addr(U_p_s,MAT_IMM,U_p_s); 
      if ((mu==3) && (pmdir==1)){
        queue_cmovgt(ss,U_p,U_p_s); // s<Ls-1 : loop back to U_p
      }
      if ((mu==3) && (pmdir==0)){ // Chi-mu touch is in now
	do_prefetch(U_p_s,1);     // Send gauge link touch if needed
      }

  /*****************************************************************
   * Reconstruct / accumulate into Psi
   *****************************************************************
   */
    /*Upper two components are simple, by virtue of our choice of projection*/      
    /*If addto, load Psi, else must do complex move for mu==0*/
      for(co=0;co<3;co++) complex_add(Psi[0][co],Psi[0][co],UChiR[0][co]);
      for(co=0;co<3;co++) complex_add(Psi[1][co],Psi[1][co],UChiR[1][co]);
      
    /*Lower two components are complex, by virtue of our choice of projection*/      
      if ( (mu==0) && (pm==0) )	for(co=0;co<3;co++) complex_AmiB(Psi[2][co],Psi[2][co],UChiR[1][co]);
      if ( (mu==0) && (pm==1) )	for(co=0;co<3;co++) complex_ApiB(Psi[2][co],Psi[2][co],UChiR[1][co]);
      if ( (mu==1) && (pm==0) )	for(co=0;co<3;co++) complex_add (Psi[2][co],Psi[2][co],UChiR[1][co]);
      if ( (mu==1) && (pm==1) ) for(co=0;co<3;co++) complex_sub (Psi[2][co],Psi[2][co],UChiR[1][co]);
      if ( (mu==2) && (pm==0) ) for(co=0;co<3;co++) complex_AmiB(Psi[2][co],Psi[2][co],UChiR[0][co]);
      if ( (mu==2) && (pm==1) ) for(co=0;co<3;co++) complex_ApiB(Psi[2][co],Psi[2][co],UChiR[0][co]);
      if ( (mu==3) && (pm==0) ) for(co=0;co<3;co++) complex_add (Psi[2][co],Psi[2][co],UChiR[0][co]);
      if ( (mu==3) && (pm==1) ) for(co=0;co<3;co++) complex_sub (Psi[2][co],Psi[2][co],UChiR[0][co]);

      complex_load(Chimu[0][0],PSIIMM[0][0][0],Chi_mu,FourSpinType);
      complex_load(Chimu[0][1],PSIIMM[0][1][0],Chi_mu,FourSpinType);


      complex_load(Chimu[2][0],PSIIMM[2][0][0],Chi_mu,FourSpinType);
      complex_load(Chimu[2][2],PSIIMM[2][2][0],Chi_mu,FourSpinType);
      complex_load(Chimu[3][1],PSIIMM[3][1][0],Chi_mu,FourSpinType);

      complex_load(Chimu[0][2],PSIIMM[0][2][0],Chi_mu,FourSpinType);
      complex_load(Chimu[1][2],PSIIMM[1][2][0],Chi_mu,FourSpinType);

      if ( (mu==0) && (pm==0) )	for(co=0;co<3;co++) complex_AmiB(Psi[3][co],Psi[3][co],UChiR[0][co]);
      if ( (mu==0) && (pm==1) ) for(co=0;co<3;co++) complex_ApiB(Psi[3][co],Psi[3][co],UChiR[0][co]);
      if ( (mu==1) && (pm==0) )	for(co=0;co<3;co++) complex_sub(Psi[3][co],Psi[3][co],UChiR[0][co]);
      if ( (mu==1) && (pm==1) ) for(co=0;co<3;co++) complex_add(Psi[3][co],Psi[3][co],UChiR[0][co]);
      if ( (mu==2) && (pm==0) )	for(co=0;co<3;co++) complex_ApiB(Psi[3][co],Psi[3][co],UChiR[1][co]);
      if ( (mu==2) && (pm==1) ) for(co=0;co<3;co++) complex_AmiB(Psi[3][co],Psi[3][co],UChiR[1][co]);
      if ( (mu==3) && (pm==0) ) for(co=0;co<3;co++) complex_add(Psi[3][co],Psi[3][co],UChiR[1][co]);
      if ( (mu==3) && (pm==1) ) for(co=0;co<3;co++) complex_sub(Psi[3][co],Psi[3][co],UChiR[1][co]);


      /* Store result Psi & Pointer update*/    
      if ( (mu==3) && (pmdir==1) ) {
	for(sp=0;sp<4;sp++){
	  for(co=0;co<3;co++){
#ifdef DEBUG
	    debugC(Psi[sp][co]);
#endif
	    complex_store(Psi[sp][co],PSIIMM[sp][co][0],Psi_p,FourSpinType);
	  }
	}
      }

      /*END MU/PM loops*/

   }
  }


  queue_iadd_imm(Psi_p,Psi_p,PSI_IMM);

  stop_loop(branchls,s); /*Loop over Ls*/   /* TERMINATION point of the loop*/

  queue_load_addr(s_offset,ZERO_IMM,U_p); // unlock from s_offset
  for(int i=0;i<5;i++) { // 9*8 32 bytes => 36 clines => 6*6
    l1_unlock(s_offset,0);  // Unlock our old U's.
    l1_unlock(s_offset,1);  
    l1_unlock(s_offset,2);  
    l1_unlock(s_offset,3);  
    l1_unlock(s_offset,4);  
    l1_unlock(s_offset,5); 
    queue_load_addr(s_offset,PSI_IMM,s_offset); 
  }
  l1_unlock(s_offset,0); 
  l1_unlock(s_offset,1); 
  queue_iadd_imm(U_p,U_p_s,ZERO_IMM);

  l1_unlock(tab,0);
  queue_load_addr(tab,ZERO_IMM,tab_s); // Advance

  stop_loop(branchsite,length);

  //l1_unlock(Complex_i,0);  

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




