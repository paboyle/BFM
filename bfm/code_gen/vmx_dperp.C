/*
 *
 *  Copyright UKQCD Collaboration, March 2007.
 *  Written by Peter Boyle.
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

void dwf_dperp( char *);

/*Options flags*/
int human = 0;
char name[80];
int dagger = 0;

void do_writehint(int,int);
#define do_WH(A,B) ;

int addto = 0;
Datum FourSpinType=Double;
Datum GaugeType=Double;
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
    case 'a': addto = 1; exit(0); break;
    case 's': FourSpinType= GaugeType = Single; break;
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
  dwf_dperp(name);

  /*Filter through an virtual out of order processor */
  schedule_for_proc();

  /*Dump the resulting code*/
  dump_instruction_queue();

  return(0);
}

void dwf_dperp( char *name)
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
  reg_array_2d(Chi ,Cregs,4,3);// 2 spionr      - 6 regs 
  reg_array_2d(Psi ,Cregs,4,3); /*Output 4-spinor - 12 regs*/
  alreg(creg,Cregs);
  alreg(coeffm,Cregs);
  alreg(coeffp,Cregs);

  offset_3d(CHIIMM,FourSpinType,2,3,2*nsimd());
  offset_1d(COEFF,Double,4*nsimd());

  int t;

  /*
   * Integer registers
   */
  alreg(Psi_p,Iregs);   /*Pointer to current cpt in/out PSI field      */
  alreg(Chi_p,Iregs);   /*Pointer to current cpt in/out PSI field      */
  alreg(length,Iregs);  /*number of sites*/
  alreg(Complex_i,Iregs);  /*Point to (0,1)x Nsimd*/
  alreg(Ls,Iregs);        
  alreg(tab,Iregs);        

  alreg(site,Iregs);  /*number of sites*/
  alreg(sloop,Iregs);        
  alreg(args,Iregs);        
  alreg(tabent,Iregs);        
  alreg(s,Iregs);        
  alreg(s_m,Iregs);        
  alreg(s_p,Iregs);        
  alreg(s_pw,Iregs);        // wrapped versions
  alreg(s_mw,Iregs);        

  alreg(ptr,Iregs);        
  alreg(baserel,Iregs);        
  alreg(tmp,Iregs);        

  alreg(coeff_p,Iregs);        
  alreg(pref,Iregs);

  /*Useful integer immediate constants, in units of Fsize*/
  def_off( ZERO_IMM,Byte,0);
  def_off( CHI_IMM, FourSpinType,12*nsimd());
  def_off( PSI_IMM, FourSpinType,24*nsimd());

  int Isize      = def_offset(PROC->I_size,Byte,"Isize");
  int negone     = def_offset(-1,Byte,"negone");
  int plusone    = def_offset(1,Byte,"plusone");

  /*
   * Local variables for C++ code
   */
  int i,sp,co,j,k,nxt,ri;

  /*********************************************************************
   * Start of the "pseudo assembler proper.
   *********************************************************************/
  make_inst(DIRECTIVE,Enter_Routine,name);
  grab_stack(0);
  save_regs();

  /*
   * Define our arguments 
   */
  getarg(args); /*Pointer to arg list*/

  pragma(DCBT_SPACE,0);

  queue_iload(Chi_p, ZERO_IMM,args);  queue_load_addr(args,Isize,args);  
  queue_iload(Psi_p, ZERO_IMM,args);  queue_load_addr(args,Isize,args);  
  queue_iload(length,ZERO_IMM,args);  queue_load_addr(args,Isize,args);  
  queue_iload(Ls,    ZERO_IMM,args);  queue_load_addr(args,Isize,args);  
  queue_iload(tab,   ZERO_IMM,args);  queue_load_addr(args,Isize,args);  
  queue_iload(Complex_i,ZERO_IMM,args);    queue_load_addr(args,Isize,args);  
  queue_iload(coeff_p,ZERO_IMM,args);   

  /*Bugger 12 int regs for offsets as indexed loads only*/
  for ( i =0; i<6; i++ ) { 
    need_constant(i*2*SizeofDatum(FourSpinType)*nsimd());
  }
  int   two = args;
  int   negtwomass = coeff_p;

  queue_load_addr(two,COEFF[0],coeff_p);  
  queue_load_addr(negtwomass,COEFF[nsimd()*2],coeff_p);  

  /*
   * Site loop
   */

  queue_iload_imm(site,ZERO_IMM);
  retno = get_target_label(); /*Branch to exit if length <1*/
  check_iterations(length,retno); 
  branchsite = start_loop(length);

    /*
     * Here parity is NOT easy to compute as a function of site use a lookup table
     */
    queue_iload(tabent,ZERO_IMM,tab);
    queue_load_addr(tab,Isize,tab);

    /*Fifth dimension loop*/
    queue_iadd_imm(sloop,Ls,ZERO_IMM);
    queue_iload_imm(s,ZERO_IMM);


      /*checkerboarded site hopping arithmetic
       *result_cb is parity par of result
       *If result_cb == cb4d then neighbours on source_cb are s-1,s
       *If result_cb != cb4d then neighbours on source_cb are s,s+1
       */
#define POINTER_ARITHMETIC() \
      if ( par == 0 ) {\
	make_inst(IALUPIPE,IOR,ptr,two,two);   \
	queue_iadd_imm(s_p,s,ZERO_IMM); \
	queue_iadd_imm(s_m,s,negone); \
	queue_iadd    (s_mw,s_m,Ls); \
 \
	queue_cmovlt(s_m,negtwomass,ptr);\
	queue_cmovlt(s_m,s_mw,s_m);\
 \
	complex_load(coeffp,ZERO_IMM,two,Double);\
	complex_load(coeffm,ZERO_IMM,ptr,Double);\
      } else { \
	make_inst(IALUPIPE,IOR,ptr,negtwomass,negtwomass);   \
	queue_iadd_imm(s_m,s,ZERO_IMM); \
	queue_iadd_imm(s_p,s,plusone); \
	queue_isub(s_pw,s_p,Ls);\
	queue_cmovlt(s_pw,two,ptr);\
	queue_cmovlt(s_pw,s_p,s_pw);\
	make_inst(IALUPIPE,IOR,s_p,s_pw,s_pw);\
	complex_load(coeffm,ZERO_IMM,two,Double);\
	complex_load(coeffp,ZERO_IMM,ptr,Double);\
      } \
      queue_imul(baserel,site,Ls);    \
      queue_imul_imm(baserel,baserel,PSI_IMM);    \
      queue_iadd(ptr,baserel,Chi_p);    \
      queue_imul_imm(s_p,s_p,PSI_IMM);  \
      queue_iadd(s_p,s_p,ptr);    \
      queue_imul_imm(s_m,s_m,PSI_IMM);    \
      queue_iadd(s_m,s_m,ptr);    \
      if ( !dagger ) queue_iadd_imm(s_p,s_p,CHI_IMM); \
      if ( dagger )  queue_iadd_imm(s_m,s_m,CHI_IMM); 

    /*
      if ( !dagger ) { \
	do_prefetch(s_p,3);\
	do_prefetch(s_p,4);\
	do_prefetch(s_p,5);\
	do_prefetch(s_m,0);\
	do_prefetch(s_m,1);\
	do_prefetch(s_m,2);\
      } else {\
	do_prefetch(s_p,0);\
	do_prefetch(s_p,1);\
	do_prefetch(s_p,2);\
	do_prefetch(s_m,3);\
	do_prefetch(s_m,4);\
	do_prefetch(s_m,5);\
      }
    */
    int par_body[2];
    par_body[0] = get_target_label();
    par_body[1] = get_target_label();
    int finish  = get_target_label();

    conditional_branch_cmpzero ( BRANCH_GT, tabent, par_body[1] );

      // Code two cases based on tabent
    for(int par=0;par<2;par++){

      make_inst(DIRECTIVE,CmovTarget,par_body[par]);


      POINTER_ARITHMETIC();    

      int branchls = start_loop(sloop);

      for(int sp=0;sp<2;sp++){
	for(int co=0;co<3;co++){
	  if ( !dagger )
	    complex_load(Chi[sp+2][co],CHIIMM[sp][co][0],s_p,FourSpinType);
	  else
 	    complex_load(Chi[sp][co],CHIIMM[sp][co][0],s_p,FourSpinType);
	}
      }

      for(int sp=0;sp<2;sp++){
	for(int co=0;co<3;co++){
	  if ( dagger )
	    complex_load(Chi[sp+2][co],CHIIMM[sp][co][0],s_m,FourSpinType);
	  else
	    complex_load(Chi[sp][co],CHIIMM[sp][co][0],s_m,FourSpinType);
	}
      }

      queue_iadd(ptr,baserel,Psi_p);    // Offset from the base pointers
      queue_imul_imm(tmp,s,PSI_IMM);
      queue_iadd(ptr,tmp,ptr);  
      queue_iadd_imm(tmp,ptr,ZERO_IMM);   
      for(int osp=0;osp<2;osp++){
      for(int isp=0;isp<2;isp++){
	int sp=osp*2+isp;
	for(int co=0;co<3;co++){
	    if ( dagger ) { 
	      if ( sp>=2 ) { 
		complex_real_mul  (Psi[sp][co],coeffm,Chi[sp][co]);
	      } else { 
		complex_real_mul  (Psi[sp][co],coeffp,Chi[sp][co]);
	      }
	    } else { 
	      if ( sp>=2 ) { 
		complex_real_mul  (Psi[sp][co],coeffp,Chi[sp][co]);
	      } else { 
		complex_real_mul  (Psi[sp][co],coeffm,Chi[sp][co]);
	      }
	    }
	}
      }
      queue_iadd_imm(ptr,ptr,CHI_IMM);
      }
	pragma(STORE_LIM,1);
      for(int osp=0;osp<2;osp++){
      for(int isp=0;isp<2;isp++){
	int sp=osp*2+isp;
	for(int co=0;co<3;co++){
	  complex_store(Psi[sp][co],CHIIMM[isp][co][0],tmp,FourSpinType);
	}
      }
      queue_iadd_imm(tmp,tmp,CHI_IMM);
      }

      queue_iadd_imm(s,s,plusone);
      POINTER_ARITHMETIC();

    stop_loop(branchls,sloop);
    make_inst(BRCHPIPE,BRANCH,finish);
    }

    make_inst(DIRECTIVE,Target,finish);
    queue_iadd_imm(site,site,plusone);
  /*
   * TERMINATION point of the loop
   */
  stop_loop(branchsite,length);
  make_inst(DIRECTIVE,Target,retno);

  /*
   * EPILOGUE
   */
  restore_regs();
  free_stack();
  make_inst(DIRECTIVE,Exit_Routine,name);

  return;
}
  
