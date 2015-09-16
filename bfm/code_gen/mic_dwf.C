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
int half_precision=0;
void do_writehint(int,int);

int LOCK_GAUGE=0;
#define do_WH(A,B) ;
#define   lock l1_lock
#define unlock do_flush

//#define   lock do_prefetch
//#define unlock l1_unlock

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
int addto = 0;
int scaleoutput = 0;
int sgny=0;
char procname[80]="UNSPECIFIED";

bool do_interior=true;
bool do_exterior=true;
bool do_dperp   =false;
bool do_axpy    =true;
bool do_addto   =true;
bool do_norm    =false;
bool do_scale    =false;

int main ( int argc, char *argv[])
{
  struct processor *myproc;
  int arg;
  char *c;

  name[0] = '\0';
  /*Process command line*/
  while ( ( arg = getopt(argc,argv,"ahRn:P:dskfiloerpx")) != EOF){
    switch(arg){
    case 'R': human = 1; break;
    case 'n': if (strlen(optarg)<30) strcpy(name,optarg); break;
    case 'P': if (strlen(optarg)<30) strcpy(procname,optarg); break;
    case 'd': dagger = 1; break;
    case 'y': sgny = 1; break;
    case 'a': addto = 1; break;
    case 'h': half_precision = 1; break;
    case 'k': do_scale=1;break;
    case 'i': do_interior=1;do_exterior=0;break;
    case 'e': do_exterior=1;do_interior=0;break;
    case 'o': scaleoutput = 1; break;
    case 'l': LOCK_GAUGE  = 1; break;
    case 'p': do_dperp=1;break;
    case 'r': do_norm=1;break;
    case 'x': do_axpy=1;break;
    case 'f': dd2 = 1; break;
    case 's': FourSpinType=Single; GaugeType=Single;  break;
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
  streams_reg_alloc(0);
  set_processor_optarg(procname);
  set_human_readable(human);
  setup_cmadds();

  /*Write the naive asm code*/
  dwf_dslash(name);

  /*Filter through an virtual out of order processor */
  //  schedule_for_proc();

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

#if 0
  reg_array_1d(tmp,Cregs,8); 
  int creg = tmp[7];
  int nreg = tmp[6];
  int U[3][3] = {
    {tmp[0],tmp[3],tmp[0]}, // Note reuse of tmp[6,0]. Care taken 
    {tmp[1],tmp[4],tmp[1]}, // U[1][2], U[2][2] are the danger elements
    {tmp[2],tmp[5],tmp[2]}
  };
#else
  reg_array_1d(tmp,Cregs,4); 
  int creg = tmp[3];
  int nreg = tmp[0];
  int U[3][3] = {
    {tmp[0],tmp[0],tmp[0]}, 
    {tmp[1],tmp[1],tmp[1]}, 
    {tmp[2],tmp[2],tmp[2]}
  };
#endif
  int Chimu[4][3] = {
    {Chi[0][0],Chi[0][1],Chi[0][2]},
    {Chi[1][0],Chi[1][1],Chi[1][2]},
    {UChi[0][0],UChi[0][1],UChi[0][2]},
    {UChi[1][0],UChi[1][1],UChi[1][2]}
  };

  offset_3d(CHIIMM,FourSpinType,2,3,2*nsimd());
  offset_3d(PSIIMM,FourSpinType,4,3,2*nsimd());
  offset_3d(PSIIMMH,Half,4,3,2*nsimd());

  int t;

  /*
   * Integer registers -- 12 regs
   */
  // Following are persistent state across s-loop iterations
  alreg(U_p,Iregs);    /*Pointer to the current cpt of gauge field    */
  alreg(Chi_p,Iregs);  /*Pointer to the input four spinor             */
  alreg(tab,Iregs);     /*Pointer to current entry in offset table*/
  alreg(Complex_i,Iregs);  /*Point to (0,1)x Nsimd*/
  alreg(Ls,Iregs);        
  alreg(s,Iregs);        
  alreg(recbuf_base,Iregs);
  alreg(mask,Iregs);
  alreg(s_offset,Iregs);  // Compute as multiple of s ?

  // Temporaries
  // Have to spill length, Psi_p, outptr to stack. As these are outside direction
  // loop the overhead is still strongly suppressed.
  alreg(temp1,Iregs);
  alreg(temp2,Iregs);
  alreg(temp3,Iregs);
  alreg(temp4,Iregs);
  alreg(temp5,Iregs);
  int length = temp3;
  int Psi_p  = temp4;
  int outptr = temp5;

  // Alias args to mask
  int args = mask;

  /*Useful integer immediate constants, in units of Fsize*/
  def_off( ZERO_IMM,Byte,0);
  def_off( minusone,Byte,-1);
  def_off( one,Byte,1);
  def_off( m64,Byte,-64);
  def_off( GAUGE_SITE_IMM, FourSpinType,8*18*nsimd());
  def_off( PSI_IMM  , FourSpinType,24*nsimd());
  def_off( CHI_IMM, FourSpinType,12*nsimd());
  def_off( MAT_IMM, GaugeType,18*nsimd());
  offset_3d(GIMM    , GaugeType, 3, 3 ,2*nsimd() );
  offset_1d(TAB_IMM,TableType,17);
  offset_1d(SPILL_IMM,Integer,16);
  const int SpillLength=8;
  const int SpillPsiP  =9;
  const int SpillOutptr=10;
  const int SpillMask  =11;
  
  //  int cmplx_i_offset    = CONSTIMM[0][0];// i
  //  int cmplx_mi_offset   = CONSTIMM[1][0];// -i
  def_off(zero_offset,ConstantType,2*nsimd()*2);
  def_off(two_offset ,ConstantType,3*nsimd()*2);
  def_off(negtwomass_offset,ConstantType,4*nsimd()*2);
  def_off(A_offset         ,ConstantType,5*nsimd()*2);
  def_off(B_offset         ,ConstantType,6*nsimd()*2);
  def_off(norm_offset      ,ConstantType,7*nsimd()*2);

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

  def_off(bits16,Byte,0xFFFF);
  def_off(thirtytwo,Byte,32);
  def_off(sixteen,Byte,16);


  int i,j,co,sp;

  /*********************************************************************/

  pragma(DCBT_SPACE,0);
  pragma(LOAD_LIM,1);
  pragma(LOAD_SPACE,0);
  make_inst(DIRECTIVE,Enter_Routine,name);

  int hbias = grab_stack(128); // 128 bytes on stack frame
  save_regs();
  queue_iadd_imm(PROC->StackPointer,PROC->StackPointer,hbias);

  /*********************************************
   * our arguments 
   *********************************************
   */
  getarg(args); /*Pointer to arg list*/

  queue_iload(Psi_p, ZERO_IMM,args);   queue_load_addr(args,Isize,args);   //0
  queue_iload(Chi_p, ZERO_IMM,args);   queue_load_addr(args,Isize,args);   //1
  queue_iload(U_p,   ZERO_IMM,args);   queue_load_addr(args,Isize,args);   //2
  queue_iload(length,ZERO_IMM,args);   queue_load_addr(args,Isize,args);   //3
  queue_iload(Ls,    ZERO_IMM,args);   queue_load_addr(args,Isize,args);   //4
  queue_iload(tab,   ZERO_IMM,args);   queue_load_addr(args,Isize,args);   //5
  queue_iload(Complex_i,ZERO_IMM,args);queue_load_addr(args,Isize,args);   //6
                                       queue_load_addr(args,Isize,args);   //7  
  queue_iload(recbuf_base,ZERO_IMM,args);queue_load_addr(args,Isize,args); //8  
  queue_iload(outptr,ZERO_IMM,args);                                    //9

  queue_istore(Psi_p,SPILL_IMM[SpillPsiP],PROC->StackPointer);
  queue_istore(outptr,SPILL_IMM[SpillOutptr],PROC->StackPointer);


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

  complex_load(nreg,zero_offset,Complex_i,ConstantType); // Zero the norm reg

  if ( half_precision ) {
    int Mask = temp4;
    queue_iload_imm(Mask,ZERO_IMM);
    queue_ori(Mask,Mask,bits16);
    queue_lshift(Mask,Mask,thirtytwo);
    queue_ori(Mask,Mask,bits16);
    queue_lshift(Mask,Mask,sixteen);
    queue_istore(Mask,SPILL_IMM[SpillMask],PROC->StackPointer);
  }

  /**********************************
   * Site loop 
   **********************************
   */

  retno = get_target_label(); 
  check_iterations(length,retno); 
  int branchsite = start_loop(length);
  queue_istore(length,SPILL_IMM[SpillLength],PROC->StackPointer);
  if( LOCK_GAUGE) {
  lock(tab,0);
  if ( do_interior ) {
    queue_load_addr(s_offset,ZERO_IMM,U_p); // lock from new U_p
    if ( GaugeType == Double) {
      for(int i=0;i<6;i++) { 
	lock(s_offset,0); // lock next set of links (32 lines, 2304 bytes)
	lock(s_offset,1); // Exactly 18*2*8*8 = 2304
	lock(s_offset,2);
	lock(s_offset,3);
	lock(s_offset,4);
	lock(s_offset,5);
	queue_load_addr(s_offset,PSI_IMM,s_offset); 
      }
    } else { 
      for(int i=0;i<6;i++) { 
	lock(s_offset,0); // lock next set of links (18 lines)
	lock(s_offset,1); // Exactly 18*2*8*4 = 1152
	lock(s_offset,2);
	queue_load_addr(s_offset,PSI_IMM,s_offset); 
      }
    }
  }
  }

  //////////////////////////////////////////////////
  // lock next set of links (32 lines, 2304 bytes)
  // lock from new U_p
  // Exactly 18*2*8*8 = 2304
  // Simply suppress if doing face only. Could predicate it.
  ///////////////////////////////////////////////////
  
  /**********************************
   *  s-loop
   **********************************
   */
  queue_iadd_imm (s,Ls,ZERO_IMM);
  queue_iload_imm(s_offset,ZERO_IMM);        
  //  queue_iload_imm(dperp_s, ZERO_IMM);

  // Mask holding the directions on halo
  queue_iload_short(mask,TAB_IMM[10],tab);
  
  int lab_nowork;
  if ( ! do_interior ) { 
    //    lab_nowork = get_target_label(); 
    //    check_iterations(mask,lab_nowork);
  }
  int branchls   = start_loop(s); 


  if ( do_dperp ) {

    int dperp_chi_mu  = temp1;
    int dperp_chi_mu_p= temp2;
    int dperp_chi_mu_m= temp3;

    int dperp_s      = temp1;  // carefully alias with dperp_chi_mu
    int dperp_sp     = temp2;  // s+1-Ls or s+1 
    int dperp_sm     = temp3;  

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

    int cb_zero_branch  = get_target_label();
    int cb_done_branch = get_target_label();

    { 
      int dperp_cb = temp2;
      queue_iload_short(dperp_cb,TAB_IMM[11],tab);  // dperp_cb == cb value for this site
      check_iterations (dperp_cb,cb_zero_branch );  // conditional branch based on cb
    }

    ///////////
    // CB==1
    ///////////

    complex_load(coeffm,two_offset,Complex_i,ConstantType);
    complex_load(coeffp,negtwomass_offset,Complex_i,ConstantType);


    queue_isub( dperp_s, Ls,s);           // dperp_s = true s index

    // if ( s==Ls-1 ) ss = s+1-Ls
    // else           ss = s+1
    queue_iadd_imm(dperp_sp,dperp_s,one); 
    queue_isub(dperp_sp,dperp_sp,Ls); 
    int wrap = get_target_label();
    conditional_branch_cmpzero ( BRANCH_EQ, dperp_sp, wrap );
      complex_load(coeffp,two_offset,Complex_i,ConstantType);
      queue_iadd_imm(dperp_sp,dperp_s,one);
    make_inst(DIRECTIVE,Target,wrap);

    queue_imul_imm(dperp_sp,dperp_sp,PSI_IMM);  // XX<- s+1 offset
    queue_imul_imm(dperp_sm,dperp_s,PSI_IMM);   // ss<- s-1 offset

      // sm     -> s-1  (checkerboarded-> s) 
      // sp     -> s+1  (checkerboarded-> s+1)
    // Base pointer
    queue_iload_short(dperp_chi_mu_p,TAB_IMM[9],tab);  
    queue_iadd       (dperp_chi_mu_p,dperp_chi_mu_p,Chi_p);    // Aliased with dperp_s which is no longer live
    queue_iadd    (dperp_chi_mu_m,dperp_sm,dperp_chi_mu);  
    queue_iadd    (dperp_chi_mu_p,dperp_sp,dperp_chi_mu);    

    jump(cb_done_branch);

    make_inst(DIRECTIVE,Target,cb_zero_branch);

    ///////////
    // CB==0
    ///////////
    int nowrap = get_target_label();

    complex_load(coeffp,two_offset,Complex_i,ConstantType);
    complex_load(coeffm,two_offset,Complex_i,ConstantType);

    queue_iadd_imm(dperp_sm,dperp_s,minusone);
    conditional_branch_cmpzero ( BRANCH_GE, dperp_sm, nowrap );
      complex_load(coeffm,negtwomass_offset,Complex_i,ConstantType);
      queue_iadd    (dperp_sm,dperp_sm,Ls);
    make_inst(DIRECTIVE,Target,nowrap);

    // dperp_sm     -> s-1  (checkerboarded-> s-1)
    // dperp_sp     -> s+1  (checkerboarded-> s)
    queue_imul_imm(dperp_sp,dperp_s,PSI_IMM);
    queue_imul_imm(dperp_sm,dperp_sm,PSI_IMM);
    
    queue_iload_short(dperp_chi_mu,TAB_IMM[9],tab);  
    queue_iadd       (dperp_chi_mu,dperp_chi_mu,Chi_p);    
    
    queue_iadd(dperp_chi_mu_m,dperp_sm,dperp_chi_mu); 
    queue_iadd(dperp_chi_mu_p,dperp_sp,dperp_chi_mu);   
    
    make_inst(DIRECTIVE,Target,cb_done_branch);

    {
	int op=2;
	int om=0;
	if ( dagger ) {
	  op=0;
	  om=2;
	}

	complex_load(Psi[0+op][0],PSIIMM[0+op][0][0],dperp_chi_mu_p,FourSpinType);
	complex_load(Psi[0+om][0],PSIIMM[0+om][0][0],dperp_chi_mu_m,FourSpinType);

	complex_load(Psi[0+op][2],PSIIMM[0+op][2][0],dperp_chi_mu_p,FourSpinType);
	complex_load(Psi[0+om][2],PSIIMM[0+om][2][0],dperp_chi_mu_m,FourSpinType);

	complex_load(Psi[1+op][1],PSIIMM[1+op][1][0],dperp_chi_mu_p,FourSpinType);
	complex_load(Psi[1+om][1],PSIIMM[1+om][1][0],dperp_chi_mu_m,FourSpinType);

	complex_load(Psi[0+op][1],PSIIMM[0+op][1][0],dperp_chi_mu_p,FourSpinType);
	complex_load(Psi[0+om][1],PSIIMM[0+om][1][0],dperp_chi_mu_m,FourSpinType);

	complex_load(Psi[1+op][0],PSIIMM[1+op][0][0],dperp_chi_mu_p,FourSpinType);
	complex_load(Psi[1+om][0],PSIIMM[1+om][0][0],dperp_chi_mu_m,FourSpinType);

	complex_load(Psi[1+op][2],PSIIMM[1+op][2][0],dperp_chi_mu_p,FourSpinType);
	complex_load(Psi[1+om][2],PSIIMM[1+om][2][0],dperp_chi_mu_m,FourSpinType);

	for(int sp=0;sp<2;sp++){
	for(int co=0;co<3;co++){
	  complex_real_mul  (Psi[sp+op][co],coeffp,Psi[sp+op][co]);
	  complex_real_mul  (Psi[sp+om][co],coeffm,Psi[sp+om][co]);
	}}
    }

  }

  /////////////////////////////////////////////////////////
  // Setup the Psi field for accumulation
  // Could scale Psi (axpy).
  /////////////////////////////////////////////////////////
  if ( addto && (do_dperp == 0) ) {

    int A=U[2][1];
    
    queue_iload(Psi_p,SPILL_IMM[SpillPsiP],PROC->StackPointer);
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

    if ( do_scale ) { 

      complex_load(A,A_offset,Complex_i,ConstantType);

	simd_mul(Psi[0][0],A,Psi[0][0]);
	simd_mul(Psi[0][2],A,Psi[0][2]);
	simd_mul(Psi[1][1],A,Psi[1][1]);
	simd_mul(Psi[2][0],A,Psi[2][0]);
	simd_mul(Psi[2][2],A,Psi[2][2]);
	simd_mul(Psi[3][1],A,Psi[3][1]);


	simd_mul(Psi[0][1],A,Psi[0][1]);
	simd_mul(Psi[1][0],A,Psi[1][0]);
	simd_mul(Psi[1][2],A,Psi[1][2]);
	simd_mul(Psi[2][1],A,Psi[2][1]);
	simd_mul(Psi[3][0],A,Psi[3][0]);
	simd_mul(Psi[3][2],A,Psi[3][2]);

    }

  } 

  if ( addto && (do_dperp == 1) ) {

    int A = U[2][1];
    queue_iload(Psi_p,SPILL_IMM[SpillPsiP],PROC->StackPointer);
    complex_load(A,A_offset,Complex_i,ConstantType);

      complex_load(Chimu[0][0],PSIIMM[0][0][0],Psi_p,FourSpinType);
      complex_load(Chimu[0][2],PSIIMM[0][2][0],Psi_p,FourSpinType);
      complex_load(Chimu[1][1],PSIIMM[1][1][0],Psi_p,FourSpinType);

      complex_load( UChi[0][0],PSIIMM[2][0][0],Psi_p,FourSpinType);
      complex_load( UChi[0][2],PSIIMM[2][2][0],Psi_p,FourSpinType);
      complex_load( UChi[1][1],PSIIMM[3][1][0],Psi_p,FourSpinType);

      complex_load(Chimu[0][1],PSIIMM[0][1][0],Psi_p,FourSpinType);
      complex_load(Chimu[1][0],PSIIMM[1][0][0],Psi_p,FourSpinType);
      complex_load(Chimu[1][2],PSIIMM[1][2][0],Psi_p,FourSpinType);
      complex_load(UChi[0][1] ,PSIIMM[2][1][0],Psi_p,FourSpinType);
      complex_load(UChi[1][0] ,PSIIMM[3][0][0],Psi_p,FourSpinType);
      complex_load(UChi[1][2] ,PSIIMM[3][2][0],Psi_p,FourSpinType);

      simd_madd(Psi[0][0],A,Chimu[0][0],Psi[0][0]);
      simd_madd(Psi[1][0],A,Chimu[1][0],Psi[1][0]);
      simd_madd(Psi[2][0],A,UChi[0][0],Psi[2][0]);
      simd_madd(Psi[3][0],A,UChi[1][0],Psi[3][0]);

      simd_madd(Psi[0][1],A,Chimu[0][1],Psi[0][1]);
      simd_madd(Psi[1][1],A,Chimu[1][1],Psi[1][1]);
      simd_madd(Psi[2][1],A, UChi[0][1],Psi[2][1]);
      simd_madd(Psi[3][1],A, UChi[1][1],Psi[3][1]);

      simd_madd(Psi[0][2],A,Chimu[0][2],Psi[0][2]);
      simd_madd(Psi[1][2],A,Chimu[1][2],Psi[1][2]);
      simd_madd(Psi[2][2],A,UChi[0][2],Psi[2][2]);
      simd_madd(Psi[3][2],A,UChi[1][2],Psi[3][2]);

  } else  if ( (addto == 0) && (do_dperp == 0 ) )  { // painfully loading zero... 

    for(co=0;co<3;co++) complex_load(Psi[0][co],zero_offset,Complex_i,ConstantType);
    for(co=0;co<3;co++) complex_load(Psi[1][co],zero_offset,Complex_i,ConstantType);
    for(co=0;co<3;co++) complex_load(Psi[2][co],zero_offset,Complex_i,ConstantType);
    for(co=0;co<3;co++) complex_load(Psi[3][co],zero_offset,Complex_i,ConstantType);

  }


  /////////////////////////////////////////////////////////
  // Loop over directions
  /////////////////////////////////////////////////////////

  { 
    // Scope the temporary regs U_p_s remains  live across the direction loop
    int U_p_s = temp1;
    int Chi_mu =temp2; // live only during direction loop
    int temp   =temp3;
    queue_iadd_imm(U_p_s,U_p,ZERO_IMM);

    for ( int mu=0;mu<4;mu++ ) {
      for ( int pmdir=0;pmdir<2;pmdir++ ) {
	
	/* spinor pointer  */
	int pm;
	if ( dagger ) pm = 1-pmdir;
	else          pm = pmdir;
	
	int dir = pmdir+mu*2;
	
	
	// A. Prefetching... can compute chi_nxt also 
	
	////////////////////////////////////////////////////////
	// Branch structure for interior/exterior neighbours
	////////////////////////////////////////////////////////
	int lab_skip_mu     = get_target_label(); 
	int lab_proj_mu     = get_target_label(); 
	int lab_su3_mu      = get_target_label(); 

	// If we decompress first ca 10% cycle overhead might do better
	// Address computation
#define COMPUTE_ADDRESS(mu,pmdir)			\
	queue_iand_imm(temp,mask,mask_imm[pmdir+mu*2]); \
	queue_iload_short(Chi_mu,TAB_IMM[dir],tab);  \
	queue_iadd(temp5,s_offset,s_offset);\
	if ( half_precision ) { \
	  queue_rshift(temp4,s_offset,one);\
	} else { \
	  queue_mov(s_offset,temp4);\
	}\
	queue_cmovgt(temp,temp4,temp5);\
	queue_iadd       (Chi_mu,Chi_mu,temp5);    \
	queue_mov(Chi_p,temp5);\
	queue_cmovgt(temp,recbuf_base,temp5);\
	queue_iadd(Chi_mu,Chi_mu,temp5);\

#define PREFETCH_SPINOR \
	do_prefetch(Chi_mu,0);\
	do_prefetch(Chi_mu,1);\
	do_prefetch(Chi_mu,2); \
	if ( SizeofDatum(FourSpinType)==8){\
	  do_prefetch(Chi_mu,3);\
	  do_prefetch(Chi_mu,4);\
	  do_prefetch(Chi_mu,5);\
	}

	COMPUTE_ADDRESS(mu,pmdir);

	PREFETCH_SPINOR;
	
	complex_load(Chi[0][0],PSIIMM[0][0][0],Chi_mu,FourSpinType);
	complex_load(Chi[1][0],PSIIMM[1][0][0],Chi_mu,FourSpinType);

	complex_load(Chi[0][2],PSIIMM[0][2][0],Chi_mu,FourSpinType);
	complex_load(Chi[0][1],PSIIMM[0][1][0],Chi_mu,FourSpinType);
	complex_load(Chi[1][2],PSIIMM[1][2][0],Chi_mu,FourSpinType);
	complex_load(Chi[1][1],PSIIMM[1][1][0],Chi_mu,FourSpinType);

	// branch leq so if interior goto proj
	check_iterations(temp,lab_proj_mu);

	// Exterior points are already projected. Just load.
	if (do_exterior){
	  
	  if ( half_precision ) { 
	    int memory=PROC->StackPointer;
	    int Convert1=temp4;
	    int Convert2=temp5;
	    int Mask    =temp;

	    queue_iload(Mask,SPILL_IMM[SpillMask],PROC->StackPointer);
	    // four words (32B double, 16B single) are stored as 8 bytes in mem
	    complex_load_half(Chi[0][0],PSIIMMH[0][0][0],Chi_mu,memory,Convert1,Convert2,Mask);
	    complex_load_half(Chi[0][1],PSIIMMH[0][1][0],Chi_mu,memory,Convert1,Convert2,Mask);
	    complex_load_half(Chi[0][2],PSIIMMH[0][2][0],Chi_mu,memory,Convert1,Convert2,Mask);
	    complex_load_half(Chi[1][0],PSIIMMH[1][0][0],Chi_mu,memory,Convert1,Convert2,Mask);
	    complex_load_half(Chi[1][1],PSIIMMH[1][1][0],Chi_mu,memory,Convert1,Convert2,Mask);
	    complex_load_half(Chi[1][2],PSIIMMH[1][2][0],Chi_mu,memory,Convert1,Convert2,Mask);
	  }

	  jump(lab_su3_mu);
	} else {
	  jump(lab_skip_mu);
	}
	
     ////////////////////////////////////////////////////////
     // Interior points are not already projected.
     ////////////////////////////////////////////////////////
	make_inst(DIRECTIVE,Target,lab_proj_mu); 

	if ( do_interior ) { 

	  complex_load(Chimu[2][0],PSIIMM[2][0][0],Chi_mu,FourSpinType);
	  complex_load(Chimu[3][1],PSIIMM[3][1][0],Chi_mu,FourSpinType);
	  complex_load(Chimu[2][2],PSIIMM[2][2][0],Chi_mu,FourSpinType);

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
	    int pmy = pm;
	    if ( sgny ) pmy=1-pm;
	    /*Gy
	     *  0 0  0  -1  [0] -+ [3]
	     *  0 0  1  0   [1] +- [2]
	     *  0 1  0  0
	     * -1 0  0  0
	     */
	    if ( pmy ==0 ) {
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
	for(int j=0;j<3;j++){
	  
	  complex_load(U[j][0],GIMM[j][0][0],U_p_s,GaugeType);
	  complex_load(U[j][1],GIMM[j][1][0],U_p_s,GaugeType);
	  complex_load(U[j][2],GIMM[j][2][0],U_p_s,GaugeType);

	  complex_twospin_color_dot(U[j][0],U[j][1],U[j][2],
				    Chi[0][0],Chi[0][1],Chi[0][2],
				    Chi[1][0],Chi[1][1],Chi[1][2],
				    UChi[0][j],UChi[1][j]);
	}

	/*****************************************************************
	 * Reconstruct / accumulate into Psi
	 *****************************************************************
	 */
	
	{
	  int pmy=pm;
	  if (sgny) pmy=1-pm;
	  
	  /*Upper two components are simple, by virtue of our choice of projection*/      
	  for(co=0;co<3;co++) complex_add(Psi[0][co],Psi[0][co],UChi[0][co]);
	  for(co=0;co<3;co++) complex_add(Psi[1][co],Psi[1][co],UChi[1][co]);
	  
	  /*Lower two components are complex, by virtue of our choice of projection*/      
	  if ( (mu==0) && (pm==0) ) for(co=0;co<3;co++) complex_AmiB(Psi[2][co],Psi[2][co],UChi[1][co]);
	  if ( (mu==0) && (pm==1) ) for(co=0;co<3;co++) complex_ApiB(Psi[2][co],Psi[2][co],UChi[1][co]);
	  if ( (mu==1) && (pmy==0) ) for(co=0;co<3;co++) complex_add (Psi[2][co],Psi[2][co],UChi[1][co]);
	  if ( (mu==1) && (pmy==1) ) for(co=0;co<3;co++) complex_sub (Psi[2][co],Psi[2][co],UChi[1][co]);
	  if ( (mu==2) && (pm==0) ) for(co=0;co<3;co++) complex_AmiB(Psi[2][co],Psi[2][co],UChi[0][co]);
	  if ( (mu==2) && (pm==1) ) for(co=0;co<3;co++) complex_ApiB(Psi[2][co],Psi[2][co],UChi[0][co]);
	  if ( (mu==3) && (pm==0) ) for(co=0;co<3;co++) complex_add (Psi[2][co],Psi[2][co],UChi[0][co]);
	  if ( (mu==3) && (pm==1) ) for(co=0;co<3;co++) complex_sub (Psi[2][co],Psi[2][co],UChi[0][co]);
	  
	  if ( (mu==0) && (pm==0) ) for(co=0;co<3;co++) complex_AmiB(Psi[3][co],Psi[3][co],UChi[0][co]);
	  if ( (mu==0) && (pm==1) ) for(co=0;co<3;co++) complex_ApiB(Psi[3][co],Psi[3][co],UChi[0][co]);
	  if ( (mu==1) && (pmy==0) ) for(co=0;co<3;co++) complex_sub (Psi[3][co],Psi[3][co],UChi[0][co]);
	  if ( (mu==1) && (pmy==1) ) for(co=0;co<3;co++) complex_add (Psi[3][co],Psi[3][co],UChi[0][co]);
	  if ( (mu==2) && (pm==0) ) for(co=0;co<3;co++) complex_ApiB(Psi[3][co],Psi[3][co],UChi[1][co]);
	  if ( (mu==2) && (pm==1) ) for(co=0;co<3;co++) complex_AmiB(Psi[3][co],Psi[3][co],UChi[1][co]);
	  if ( (mu==3) && (pm==0) ) for(co=0;co<3;co++) complex_add (Psi[3][co],Psi[3][co],UChi[1][co]);
	  if ( (mu==3) && (pm==1) ) for(co=0;co<3;co++) complex_sub (Psi[3][co],Psi[3][co],UChi[1][co]);
	}
	if ( (do_exterior==0) || (do_interior==0) ){
	  make_inst(DIRECTIVE,Target,lab_skip_mu);
	}
	queue_load_addr(U_p_s,MAT_IMM,U_p_s); 
	
      }      /*END MU/PM loops*/
    }
  }

  { 
    // Temporary allocation for epilogue fo loop
    int ptr    =temp3; // live only after direction loop

    if(scaleoutput==1) {
      int B = U[2][1];
      queue_iadd_imm(ptr,Complex_i,B_offset);
      complex_load(B,ZERO_IMM,ptr,ConstantType);
      for(sp=0;sp<4;sp++){
	for(co=0;co<3;co++){
	  simd_mul(Psi[sp][co],B,Psi[sp][co]);
	}
      }
    }
    /* Store result Psi & Pointer update*/    
    queue_iload(outptr,SPILL_IMM[SpillOutptr],PROC->StackPointer);
    queue_iload(Psi_p,SPILL_IMM[SpillPsiP],PROC->StackPointer);

    for(sp=0;sp<4;sp++){
      for(co=0;co<3;co++){
	//	simd_madd(nreg,Psi[sp][co],Psi[sp][co],nreg);
	complex_store(Psi[sp][co],PSIIMM[sp][co][0],outptr,FourSpinType);
      }
    }

    queue_iadd_imm(Psi_p,Psi_p,PSI_IMM);
    queue_istore(Psi_p,SPILL_IMM[SpillPsiP],PROC->StackPointer);
    queue_iadd_imm(outptr,outptr,PSI_IMM);
    queue_istore(outptr,SPILL_IMM[SpillOutptr],PROC->StackPointer);
    if(addto) touch(Psi_p,0);


    queue_iadd_imm(s_offset,s_offset,CHI_IMM);
    queue_iadd_imm(ptr,Complex_i,norm_offset);
    complex_store(nreg,ZERO_IMM,ptr,ConstantType);
  }

  stop_loop(branchls,s); 

  if ( ! do_interior ) { 
    //        int lab_hop = get_target_label(); 
    //        jump(lab_hop);
    //        make_inst(DIRECTIVE,Target,lab_nowork);
    //        queue_imul_imm(s_offset,Ls,PSI_IMM);
    //        queue_iadd_imm(Psi_p,Psi_p,s_offset);
    //        make_inst(DIRECTIVE,Target,lab_hop);
  }

  if ( LOCK_GAUGE ){
  if ( do_interior ) {
  queue_load_addr(s_offset,ZERO_IMM,U_p); // unlock from s_offset
  if ( GaugeType == Double ) {
  for(int i=0;i<6;i++) { // 9*8 32 bytes => 36 clines => 6*6
    unlock(s_offset,0);  // Unlock our old U's.
    unlock(s_offset,1);  
    unlock(s_offset,2);  
    unlock(s_offset,3);  
    unlock(s_offset,4);  
    unlock(s_offset,5); 
    queue_load_addr(s_offset,PSI_IMM,s_offset); 
  }
  } else { 
  for(int i=0;i<6;i++) { // 9*8 16 bytes => 18 clines => 6*3
    unlock(s_offset,0);  // Unlock our old U's.
    unlock(s_offset,1);  
    unlock(s_offset,2);  
    queue_load_addr(s_offset,PSI_IMM,s_offset); 
  }
  }
  }
  unlock(tab,0); 
  }
  queue_iadd_imm(U_p,U_p,GAUGE_SITE_IMM);
  queue_load_addr(tab,TAB_IMM[16],tab);
  touch(U_p,0);

  queue_iload(length,SPILL_IMM[SpillLength],PROC->StackPointer);
  stop_loop(branchsite,length);
  make_inst(DIRECTIVE,Target,retno);

                  /*
		   *
		   * EPILOGUE
		   *
		   */

  queue_isub_imm(PROC->StackPointer,PROC->StackPointer,hbias);
  restore_regs();
  free_stack();
  make_inst(DIRECTIVE,Exit_Routine,name);

  return;
}



