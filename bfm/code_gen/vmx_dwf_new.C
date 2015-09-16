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
  schedule_for_proc();

  /*Dump the resulting code*/
  dump_instruction_queue();

  return(0);
}

void su3_mult ( int U[3][3],
		register_array_2d &ChiResult,
		register_array_2d &ChiIn);

void dwf_dslash( char *name)
{
  /*
   * This marks the argument registers as defined by ABI as off limits
   * to us until they are freed by "getarg()";
   */
  int dum = defargcount(1);

  /*Handles for the labels point*/
  int branchsite;
  int branchmu;
  int retno;

  /*-------------------------------------------------------------------------------
   * registers used 
   *-------------------------------------------------------------------------------
   */
  
  reg_array_2d(Chi ,Cregs,2,3); // 2 spionr      - 6 regs 
  reg_array_2d(UChi,Cregs,2,3); // 2 spinors     - 6 regs
  reg_array_2d(Psi,Cregs,4,3); /*Output 4-spinor - 12 regs*/
  reg_array_1d(tmp,Cregs,6); 
  alreg(creg,Cregs);
  alreg(creg1,Cregs); /*This is spare for now 
                       *plan to use to ease thrashing of constant
		       */

  /*Ideally need to hold 12+9+6+6 regs = 33 for the 
   *input 4-spin, 2 spin * 2 + gauge field
   *Use 6 spare to ease the pressure. 
   */

  int U[3][3] = {  /*Overlap the gauge field regs with 9 of the result 4 spinor*/
    {Psi[0][0],Psi[0][1],Psi[0][2]},
    {Psi[1][0],Psi[1][1],Psi[1][2]},
    {Psi[2][0],Psi[2][1],Psi[2][2]}
  };

  offset_3d(CHIIMM,FourSpinType,2,3,2*nsimd());
  offset_3d(PSIIMM,FourSpinType,4,3,2*nsimd());

  int t;

  /*
   * Integer registers
   */
  alreg(U_p,Iregs);    /*Pointer to the current cpt of gauge field    */
  alreg(Chi_p,Iregs);  /*Pointer to the input four spinor             */
  alreg(Psi_p,Iregs);  /*Pointer to current cpt output PSI field      */
  alreg(Chi_mu,Iregs);  /*Pointer to the input four spinor            */
  alreg(twospinor_base,Iregs); /*Two spinor intermediate array               */
  alreg(length,Iregs);  /*number of sites*/
  alreg(tab,Iregs);     /*Pointer to current entry in offset table*/

  alreg(Complex_i,Iregs);  /*Point to (0,1)x Nsimd*/
  alreg(Ls,Iregs);        
  alreg(Ls_skip,Iregs);    /*Increment to hop to next s-address in two spinor*/
  alreg(s,Iregs);        
  alreg(args,Iregs);        

  /*Useful integer immediate constants, in units of Fsize*/
  def_off( ZERO_IMM,Byte,0);
  def_off( two_ND,Byte,8);
  def_off( m1,Byte,-1);
  def_off( m2,Byte,-2);
  def_off( PSI_IMM, FourSpinType,24*nsimd());
  def_off( TWO_PSI_IMM, FourSpinType,48*nsimd());
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
  queue_iload(twospinor_base,ZERO_IMM,args);    queue_load_addr(args,Isize,args);   
  queue_iload(Complex_i,ZERO_IMM,args);   


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

  alreg(twospinor,Iregs);        
  alreg(twospinor_dir,Iregs);
  alreg(pref,Iregs);

  /*Reuse the argument register*/
  int chinxt   = Complex_i; // Dont need it again (I hope)


  retno = get_target_label(); /*Branch to exit if length <1*/
  check_iterations(length,retno); 

    int Psimu[4][3];
    for(co = 0; co<3;co++ ){ 
      if (co==0) { 
	Psimu[0][co] = tmp[2];
	Psimu[1][co] = tmp[3];
	Psimu[2][co] = tmp[4];
	Psimu[3][co] = tmp[5];
      } else if ( co==1 ) { 
	Psimu[0][co] = tmp[0];
	Psimu[1][co] = tmp[1];
	Psimu[2][co] = Psi[3][0]; 
	Psimu[3][co] = Psi[3][1];
      }else if ( co == 2 ) { 
	Psimu[0][co] = Psi[3][2];
	Psimu[1][co] = creg1;
	Psimu[2][co] = UChi[0][0];
	Psimu[3][co] = UChi[0][1];
      }
    }


  /*
   * Site loop
   */
  branchsite = start_loop(length);



  /****************************************************************
   * Spin project 4 spinor                                        *
   * and su3 multiply it                                          *
   ****************************************************************/
  queue_iload(Chi_mu,ZERO_IMM,tab);  
  queue_imul(Chi_mu,Chi_mu,Ls);         

  queue_iadd_imm(twospinor_dir,twospinor_base,ZERO_IMM);

  /* 
   * Unroll loop over 8 directions
   */


  for ( int mu=0;mu<4;mu++ ) {
   for ( int pmdir=0;pmdir<2;pmdir++ ) {


    queue_load_addr(tab,Isize,tab);       // Always increment table inside the mu/pm loop
    queue_iload    (chinxt,ZERO_IMM,tab); // Chinxt may not be used if we break to the 
    queue_load_addr(pref,Isize,tab);      // Touch the next element
    do_prefetch    (pref,0);

    /*
     * Load the gauge link and iterate pointer
     * This is not prefetched so I order it to expose as many distinct sequential line fetches  
     * as possible to the memory system to fill the three LFB's as efficiently as possible.
     * A second pass pulls in to registers from L1 the "other halves" of the cache lines.
     */

    complex_load(U[0][0],GIMM[0][0][0],U_p,GaugeType);
    complex_load(U[0][2],GIMM[0][2][0],U_p,GaugeType);
    complex_load(U[1][1],GIMM[1][1][0],U_p,GaugeType);
    complex_load(U[2][0],GIMM[2][0][0],U_p,GaugeType);
    complex_load(U[2][2],GIMM[2][2][0],U_p,GaugeType);


    pragma(LOAD_SPACE,0);
    pragma(STORE_LIM,10);
    pragma(DCBT_SPACE,0);
    pragma(DCBT_POST,0);
    // pm is the sign for the gamma matrix
    // pmdir is the sign for the spatial displacement (0=plus, 1=minus)
    // Matrix is 1-gammamu Psi(x+mu)
    // So in non-dagger case, first element in two spinor is ChiMinus
    int pm;
    if ( dagger ) pm = 1-pmdir;
    else          pm = pmdir;




    if ( mu == 0 && pmdir == 0){
      /*address calculation for hopping term : Base + nbrsite*Ls*PSI_IMM */
      queue_imul_imm(Chi_mu,Chi_mu,PSI_IMM);    
    }

    if ( mu == 0 && pmdir == 0){
      queue_iadd(Chi_mu,Chi_mu,Chi_p);    
      for(sp=0;sp<4;sp++ ) {
	for(co=0;co<3;co++){ 
	  complex_load(Psimu[sp][co],PSIIMM[sp][co][0],Chi_mu,FourSpinType);
	}
      }
      // pre-increment Chi
      queue_iadd_imm(Chi_mu,Chi_mu,PSI_IMM);
    }

    queue_imul(chinxt,chinxt,Ls);    

    queue_iload_imm(Ls_skip,CHI_IMM);

    queue_imul_imm(chinxt,chinxt,PSI_IMM);    

    queue_iadd_imm(s,Ls,ZERO_IMM);
    queue_iadd_imm(twospinor,twospinor_dir,ZERO_IMM); 

    queue_iadd(chinxt,chinxt,Chi_p);    


    if ( (mu == 3) && (pmdir == 1) ) { 
      do_prefetch(Psi_p,0);
      do_prefetch(tab,1);            
    } else {
      do_prefetch(chinxt,0);
    }

    //    complex_load(U[2][2],GIMM[2][2][0],U_p,GaugeType);
    queue_iadd_imm(U_p,U_p,ZERO_IMM); // Artificial dep to prevent load
                                      // reorder for LSQ flush issues
                                      // Could use the result of a load ...
    queue_fmov(U[2][2],U[2][2]);
    complex_load(U[0][1],GIMM[0][1][0],U_p,GaugeType);
    complex_load(U[1][0],GIMM[1][0][0],U_p,GaugeType);
    complex_load(U[1][2],GIMM[1][2][0],U_p,GaugeType);
    complex_load(U[2][1],GIMM[2][1][0],U_p,GaugeType);

    queue_imul_imm(Ls_skip,Ls_skip,two_ND);


  /*
   * Ls loop
   */

    int branchls = start_loop(s);

    /*
     * On last s iteration pointer swizzle
     */
    int ss = pref;
    queue_iadd_imm(ss,s,m1); // ss = s -1;
    if ( (mu == 3) && (pmdir == 1) ) { 
      queue_cmovle(ss,twospinor_base,Chi_mu);
      if ( addto ) do_prefetch(Psi_p,0); // touch the Psi_p stream
                                         // if we are adding to it
    } else { 
      queue_cmovle(ss,chinxt,Chi_mu);
    }

    // Prefetch the next Psi. Better late than never.
    do_prefetch(Chi_mu,0);
    do_prefetch(Chi_mu,1);
    do_prefetch(Chi_mu,2);
    queue_iadd_imm(pref,Chi_mu,CHI_IMM);
    do_prefetch(pref,0);
    do_prefetch(pref,1);
    do_prefetch(pref,2);



      if ( mu == 0 ) {

      /* Gx
       *  0 0  0  i    [0]+-i[3]
       *  0 0  i  0    [1]+-i[2]
       *  0 -i 0  0
       * -i 0  0  0
       *
       */

        if ( pm ==0 ) {
          for(co=0;co<3;co++) complex_ApiB(Chi[0][co],Psimu[0][co],Psimu[3][co]);
          for(co=0;co<3;co++) complex_ApiB(Chi[1][co],Psimu[1][co],Psimu[2][co]);
        } else {
          for(co=0;co<3;co++) complex_AmiB(Chi[0][co],Psimu[0][co],Psimu[3][co]);
          for(co=0;co<3;co++) complex_AmiB(Chi[1][co],Psimu[1][co],Psimu[2][co]);
        }

      } else if ( mu == 1 ) {
      
      /*Gy
       *  0 0  0  -1  [0] -+ [3]
       *  0 0  1  0   [1] +- [2]
       *  0 1  0  0
       * -1 0  0  0
       */
       
        if ( pm ==0 ) {
          for(co=0;co<3;co++) complex_sub(Chi[0][co],Psimu[0][co],Psimu[3][co]);
          for(co=0;co<3;co++) complex_add(Chi[1][co],Psimu[1][co],Psimu[2][co]);
        } else {
          for(co=0;co<3;co++) complex_add(Chi[0][co],Psimu[0][co],Psimu[3][co]);
          for(co=0;co<3;co++) complex_sub(Chi[1][co],Psimu[1][co],Psimu[2][co]);
        }

      } else if ( mu == 2 ) {

      /*Gz
       *  0 0  i  0   [0]+-i[2]
       *  0 0  0 -i   [1]-+i[3]
       * -i 0  0  0
       *  0 i  0  0
       */
        if ( pm ==0 ) {
          for(co=0;co<3;co++) complex_ApiB(Chi[0][co],Psimu[0][co],Psimu[2][co]);
	  for(co=0;co<3;co++) complex_AmiB(Chi[1][co],Psimu[1][co],Psimu[3][co]);
        } else {
          for(co=0;co<3;co++) complex_AmiB(Chi[0][co],Psimu[0][co],Psimu[2][co]);
          for(co=0;co<3;co++) complex_ApiB(Chi[1][co],Psimu[1][co],Psimu[3][co]);
        }

      } else if ( mu == 3 ) {

      /*Gt
       *  0 0  1  0 [0]+-[2]
       *  0 0  0  1 [1]+-[3]
       *  1 0  0  0
       *  0 1  0  0
       */
        if ( pm ==0 ) {
          for(co=0;co<3;co++) complex_add(Chi[0][co],Psimu[0][co],Psimu[2][co]);
	  for(co=0;co<3;co++) complex_add(Chi[1][co],Psimu[1][co],Psimu[3][co]);
        } else {
          for(co=0;co<3;co++) complex_sub(Chi[0][co],Psimu[0][co],Psimu[2][co]); 
          for(co=0;co<3;co++) complex_sub(Chi[1][co],Psimu[1][co],Psimu[3][co]);
        }

      } 

    su3_mult(U,UChi,Chi);

    pragma(DCBT_SPACE,0);
    pragma(STORE_LIM,20);



// Address calculation for two spinor much more complicated as reversed loop ordering
//
// Address should be (s * 8 + mu*2+pmdir)*12*sizeof(double)*Nsimd 
// As loop Ls innermost need several pointers
// twospinor_base
// twospinor_dir
// twospinor
//   for (site){
//     twospinor_dir = twospinor_base
//     for ( dir ){
//        twospinor = twospinor_dir
//        for ( s ){
//          blah
//          twospinor = twospinor + Ls*8*12*sizeof(double)
//        }
//        twospinor_dir = twospinor_dir + 12*sizeof(double)
//      }
//    }
//
    /*
     * Drain the stores
     */
    for(co=0;co<3;co++) complex_store(UChi[0][co],CHIIMM[0][co][0],twospinor,FourSpinType);
    for(co=0;co<3;co++) complex_store(UChi[1][co],CHIIMM[1][co][0],twospinor,FourSpinType);
    queue_iadd(twospinor,twospinor,Ls_skip);

    queue_mbar();
    for(sp=0;sp<4;sp++){
    for(co=0;co<3;co++){
      complex_load(Psimu[sp][co],PSIIMM[sp][co][0],Chi_mu,FourSpinType);
    }}

    // pre-increment Chi
    queue_iadd_imm(Chi_mu,Chi_mu,PSI_IMM);
    
    
    complex_constants_assert(cmplx_i);

    stop_loop(branchls,s); /*Loop over Ls*/

    queue_iadd_imm(twospinor_dir,twospinor_dir,CHI_IMM);

   }
  }


  /*
   * Immediately get the fetching in place for next loop
   * In addto case, we haven't already pulled in the Psi field.
   * Note we pointer swizzled to force the two spinors to be pulled in in L1 above
   */
  if (addto ){
    queue_iadd_imm(pref,Psi_p,ZERO_IMM);
    do_prefetch(pref,0);
    do_prefetch(pref,1);
    do_prefetch(pref,2);
    queue_iadd_imm(pref,pref,CHI_IMM);    
    do_prefetch(pref,0);
    do_prefetch(pref,1);
    do_prefetch(pref,2);
  }


  pragma(LOAD_SPACE,0);

   /* Now we reconstruct this site from the two spinor array very much like the old bgl_recon.C */

  /*
   * Ls loop
   */
  pragma(LOAD_LIM,10);
  pragma(DCBT_SPACE,0);
  pragma(DCBT_POST,0);

  /*
   * Set up the two spinor pointer, and start the prefetch
   * This could be commenced with the mu=3 case above?
   * For each mu we have six cache lines to fetch
   */
  queue_iadd_imm(twospinor,twospinor_base,ZERO_IMM);
  queue_iadd_imm(s,Ls,ZERO_IMM);

    /*ChiPlus and ChiMinus refer to the dirac matrix sign*/
    /*Chi is from the negative direction, UChi from the positive direction*/
  int ChiPlus[2][3];
  int ChiMinus[2][3];
  for(sp=0;sp<2;sp++){
    for(co=0;co<3;co++){
        ChiPlus [sp][co] = tmp[(1-sp)*3+co];
        ChiMinus[sp][co] = Chi[sp][co];
    }
  }


  /*****************************************************************
   * Reconstruct main loop
   *****************************************************************
   */
  int branchlsrec = start_loop(s);

  pragma(DCBT_SPACE,0);
  pragma(DCBT_POST,0);

  /*Optional prefetch of Psi*/
  if ( addto ) {
   
      for(sp=0;sp<4;sp++) {
	for(co=0;co<3;co++) {
	  complex_load(Psi[sp][co],PSIIMM[sp][co][0],Psi_p,FourSpinType);
	}
      }
    queue_iadd_imm(pref,Psi_p,CHI_IMM);	
    do_prefetch(pref,0);
    do_prefetch(pref,1);
    do_prefetch(pref,2);
    queue_iadd_imm(pref,pref,CHI_IMM);    
    do_prefetch(pref,0);
    do_prefetch(pref,1);
    do_prefetch(pref,2);
    queue_iadd_imm(pref,pref,CHI_IMM);
  }

  /*
   * Loop over directions
   */
  for ( int mu=0;mu<4;mu++ ) {

    /*
     * Register selection
     */  
    for(sp=0;sp<2;sp++){
    for(co=0;co<3;co++){
      if ( mu == 0 ) {
        ChiPlus [sp][co] = tmp[(1-sp)*3+co];
        ChiMinus[sp][co] = Chi[sp][co];
      }
      if ( mu == 1 ) {
        ChiPlus [sp][co] = UChi[sp][co];
        ChiMinus[sp][co] = Chi[sp][co];
      } 
      if ( mu == 2 ) {
        ChiPlus [sp][co] = tmp[sp*3+co];
        ChiMinus[sp][co] = Chi[sp][co];
      } 
      if ( mu == 3 ) {
        ChiPlus [sp][co] = tmp[sp*3+co];
        ChiMinus[sp][co] = UChi[sp][co];
      } 
    }
    }

    /*
     * Load this direction's CHI's
     */
    for(sp=0;sp<2;sp++) {
      for(j=0;j<3;j++) {
	if ( dagger )  complex_load( ChiMinus[sp][j],CHIIMM[sp][j][0],twospinor,FourSpinType);
	else           complex_load( ChiPlus[sp][j] ,CHIIMM[sp][j][0],twospinor,FourSpinType);
      }
    }
    queue_iadd_imm(twospinor,twospinor,CHI_IMM);
    for(sp=0;sp<2;sp++) {
      for(j=0;j<3;j++) {
	if ( dagger )  complex_load( ChiPlus[sp][j] ,CHIIMM[sp][j][0],twospinor,FourSpinType);
	else           complex_load( ChiMinus[sp][j],CHIIMM[sp][j][0],twospinor,FourSpinType);
      }
    }
    queue_iadd_imm(twospinor,twospinor,CHI_IMM);

    /*
     * And prefetch the two spinor
     * NB. Should swizzle in the next site's fields for mu==3 and s=1.
     */
    queue_iadd_imm(pref,twospinor,CHI_IMM);
    do_prefetch(twospinor,0);
    do_prefetch(twospinor,1);
    do_prefetch(twospinor,2);
    do_prefetch(pref,0);
    do_prefetch(pref,1);
    do_prefetch(pref,2);

    /*Upper two components are simple, by virtue of our choice of projection*/      
   if ( (mu == 0) && (addto == 0) ) {
	for(co=0;co<3;co++) complex_add(Psi[0][co],ChiPlus[0][co],ChiMinus[0][co]);
	for(co=0;co<3;co++) complex_add(Psi[1][co],ChiPlus[1][co],ChiMinus[1][co]);
    } else {
      if ( dagger ) {
      for(co=0;co<3;co++) complex_add(Psi[0][co],Psi[0][co],ChiMinus[0][co]);
      for(co=0;co<3;co++) complex_add(Psi[1][co],Psi[1][co],ChiMinus[1][co]);
      for(co=0;co<3;co++) complex_add(Psi[0][co],Psi[0][co],ChiPlus[0][co]);
      for(co=0;co<3;co++) complex_add(Psi[1][co],Psi[1][co],ChiPlus[1][co]);
      } else {
      for(co=0;co<3;co++) complex_add(Psi[0][co],Psi[0][co],ChiPlus[0][co]);
      for(co=0;co<3;co++) complex_add(Psi[1][co],Psi[1][co],ChiPlus[1][co]);
      for(co=0;co<3;co++) complex_add(Psi[0][co],Psi[0][co],ChiMinus[0][co]);
      for(co=0;co<3;co++) complex_add(Psi[1][co],Psi[1][co],ChiMinus[1][co]);
      }
    }

    /*Lower two components are complex, by virtue of our choice of projection*/      
    switch(mu) { 
      case 0:
	if ( addto ) { 
	  for(co=0;co<3;co++) complex_AmiB(Psi[2][co],Psi[2][co],ChiPlus[1][co]);
	  for(co=0;co<3;co++) complex_AmiB(Psi[3][co],Psi[3][co],ChiPlus[0][co]);
	  for(co=0;co<3;co++) complex_ApiB(Psi[2][co],Psi[2][co],ChiMinus[1][co]);
	  for(co=0;co<3;co++) complex_ApiB(Psi[3][co],Psi[3][co],ChiMinus[0][co]);
	} else { 
	  for(co=0;co<3;co++) complex_sub(tmp[co]  ,ChiPlus[1][co],ChiMinus[1][co]);//Need to mult by -i -> Psi[2]
	  for(co=0;co<3;co++) complex_sub(tmp[3+co],ChiPlus[0][co],ChiMinus[0][co]);// -> Psi[3]
	  //debugF(tmp[3]);
	}
	break;
      case 1:
	if ( addto ) {
	    for(co=0;co<3;co++) complex_add(Psi[3][co],Psi[3][co],ChiMinus[0][co]);
	    for(co=0;co<3;co++) complex_add(Psi[2][co],Psi[2][co],ChiPlus [1][co]);
	    for(co=0;co<3;co++) complex_sub(Psi[3][co],Psi[3][co],ChiPlus [0][co]);
	    for(co=0;co<3;co++) complex_sub(Psi[2][co],Psi[2][co],ChiMinus[1][co]);
	} else { 
	    for(co=0;co<3;co++) complex_AmiB(Psi[3][co],ChiMinus[0][co],tmp[co+3]);
	    for(co=0;co<3;co++) complex_AmiB(Psi[2][co],ChiPlus [1][co],tmp[co]);
	    for(co=0;co<3;co++) complex_sub(Psi[3][co],Psi[3][co],ChiPlus [0][co]);
	    for(co=0;co<3;co++) complex_sub(Psi[2][co],Psi[2][co],ChiMinus[1][co]);
	}
	
	break;
      case 2:
    /* psi_2 +=-iChiplus[0] +iChiminus[0] */
    /* psi_3 += iChiplus[1] -iChiminus[1] */
	if ( dagger ) { 
	  for(co=0;co<3;co++) complex_ApiB(Psi[2][co],Psi[2][co],ChiMinus[0][co]);
	  for(co=0;co<3;co++) complex_AmiB(Psi[3][co],Psi[3][co],ChiMinus[1][co]);
	  for(co=0;co<3;co++) complex_AmiB(Psi[2][co],Psi[2][co],ChiPlus[0][co]);
	  for(co=0;co<3;co++) complex_ApiB(Psi[3][co],Psi[3][co],ChiPlus[1][co]);
	} else {
	  for(co=0;co<3;co++) complex_AmiB(Psi[2][co],Psi[2][co],ChiPlus[0][co]);
	  for(co=0;co<3;co++) complex_ApiB(Psi[3][co],Psi[3][co],ChiPlus[1][co]);
	  for(co=0;co<3;co++) complex_ApiB(Psi[2][co],Psi[2][co],ChiMinus[0][co]);
	  for(co=0;co<3;co++) complex_AmiB(Psi[3][co],Psi[3][co],ChiMinus[1][co]);
	}
	break;
      case 3:
    /* psi_2 +=  Chiplus[0] - Chiminus[0] */
    /* psi_3 +=  Chiplus[1] - Chiminus[1] */
	if ( dagger ) {
	  for(co=0;co<3;co++) complex_sub(Psi[2][co],Psi[2][co],ChiMinus[0][co]);
	  for(co=0;co<3;co++) complex_sub(Psi[3][co],Psi[3][co],ChiMinus[1][co]);
	  for(co=0;co<3;co++) complex_add(Psi[2][co],Psi[2][co],ChiPlus[0][co]);
	  for(co=0;co<3;co++) complex_add(Psi[3][co],Psi[3][co],ChiPlus[1][co]);
	}else{
	  for(co=0;co<3;co++) complex_add(Psi[2][co],Psi[2][co],ChiPlus[0][co]);
	  for(co=0;co<3;co++) complex_add(Psi[3][co],Psi[3][co],ChiPlus[1][co]);
	  for(co=0;co<3;co++) complex_sub(Psi[2][co],Psi[2][co],ChiMinus[0][co]);
	  for(co=0;co<3;co++) complex_sub(Psi[3][co],Psi[3][co],ChiMinus[1][co]);
	}
	break;
    }     
    //debugF(Psi[3][0]);

  } /*Loops over direction*/

  /* Store result Psi & Pointer update*/
  for(sp=0;sp<4;sp++){
    for(co=0;co<3;co++){
      complex_store(Psi[sp][co],PSIIMM[sp][co][0],Psi_p,FourSpinType);
    }
  }
  queue_iadd_imm(Psi_p,Psi_p,PSI_IMM);

  complex_constants_assert(cmplx_i);
  stop_loop(branchlsrec,s); /*Loop over Ls*/


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

void su3_mult (int u[3][3], 
		register_array_2d &ChiResult,
		register_array_2d &ChiIn)

{
  int j,sp;

  for(j=0;j<3;j++){
      if ( j==0 ){
        complex_six_cmuls (ChiResult[0][0],u[0][j],ChiIn[0][j],
			   ChiResult[0][1],u[1][j],ChiIn[0][j],
			   ChiResult[0][2],u[2][j],ChiIn[0][j],
			   ChiResult[1][0],u[0][j],ChiIn[1][j],
			   ChiResult[1][1],u[1][j],ChiIn[1][j],
			   ChiResult[1][2],u[2][j],ChiIn[1][j]);
      } else { 
        complex_six_cmadds (ChiResult[0][0],u[0][j],ChiIn[0][j],
			    ChiResult[0][1],u[1][j],ChiIn[0][j],
			    ChiResult[0][2],u[2][j],ChiIn[0][j],
			    ChiResult[1][0],u[0][j],ChiIn[1][j],
			    ChiResult[1][1],u[1][j],ChiIn[1][j],
			    ChiResult[1][2],u[2][j],ChiIn[1][j]);
      }
  }

}


