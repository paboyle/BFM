
/*
 *
 *  Copyright UKQCD Collaboration, Feb 2002.
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

#include "registers.h"

extern struct processor *PROC;


void vmx_gather_proj( char *);

void need_cache_line(int line);
void do_prefetch ( int ptr,int line );

int permute = 0;
int extract = 0;
int extr_hi  = 0;
int half_precision=0;
char name[80]="";
char procname[80] = "UNSPECIFIED"; /*Default processor type*/
Datum FourSpinType = Double;
Datum TwoSpinType = Double;

int main ( int argc, char *argv[])
{
  int arg;
  char *c;

  /*Process command line*/
  while ( ( arg = getopt(argc,argv,"ehcn:P:psS")) != EOF){
    switch(arg){
    case 'n': if (strlen(optarg)<30) strcpy(name,optarg); break;
    case 'P': if (strlen(optarg)<30) strcpy(procname,optarg); break;
    case 'p': permute = 1; break;
    case 'e': extract = 1; break;
    case 'c': half_precision=1; break;
    case 'h': extr_hi = 1; break;
    case 'S': TwoSpinType = Single; break;
    case 's': FourSpinType = TwoSpinType= Single; break;
    default: fprintf(stderr,"Usage: %s -[hs] -n routine_name -P proc\t",argv[0]); 
             exit (1); break;
    }
  }
  if ( half_precision ) { 
    TwoSpinType = Half;
  }
  /*Control output according to user set up args*/
  set_processor_optarg(procname);

  /*queue the naive asm code*/
  vmx_gather_proj(name);
  /*Filter through an virtual out of order processor */
  //  schedule_for_proc();
  /*Dump the resulting code*/
  dump_instruction_queue();
  return(0);
}

/*
void vmx_gather_proj(integer *sites,
                     Float *psi
                     integer comm_off, 
                     Float *comm_buf,
                     integer nface,
                    );
*/

void vmx_gather_proj( char *name)
{
  int dum = defargcount(1);

  /*Integer register usage*/
  alreg(args,Iregs);
  alreg(sites,Iregs);
  alreg(comm_offset,Iregs);
  alreg(nface,Iregs);
  alreg(psi_p,Iregs);
  alreg(comm_p,Iregs);
  alreg(Ls,Iregs);
  alreg(in_p,Iregs);
  alreg(Complex_i,Iregs);
  alreg(dir,Iregs);
  alreg(commr_p,Iregs);

  /*Floating register usage*/
  reg_array_2d(Psi,Cregs,4,3);
  reg_array_2d(Chi,Cregs,2,3);
  alreg(creg,Cregs);

  def_off(ZERO_IMM  ,FourSpinType,0);
  def_off(PSI_IMM  ,FourSpinType,24*nsimd());
  def_off(CHI_IMM  ,TwoSpinType,12*nsimd());
  offset_3d(PSIIMM,FourSpinType,4,3,2*nsimd());
  offset_3d(CHIIMM,TwoSpinType,2,3,2*nsimd());
  def_off(m1  ,Byte,-1);
  def_off(bits16,Byte,0xFFFF);
  def_off(thirtytwo,Byte,32);
  def_off(sixteen,Byte,16);

  int Isize      = def_offset(PROC->I_size,Byte,"Isize");

  int brchf,brchs,retno;       /*Branch target handles*/
  int i;

  make_inst(DIRECTIVE,Enter_Routine,name);
  int bias = grab_stack(64);
  save_regs();
  queue_iadd_imm(PROC->StackPointer,PROC->StackPointer,bias);

  getarg(args); /*Pointer to arg list*/

  queue_iload(sites , ZERO_IMM,args);  queue_load_addr(args,Isize,args);  
  queue_iload(psi_p , ZERO_IMM,args);  queue_load_addr(args,Isize,args);  
  queue_iload(comm_offset, ZERO_IMM,args);  queue_load_addr(args,Isize,args);  
  queue_iload(comm_p, ZERO_IMM,args);  queue_load_addr(args,Isize,args);  
  queue_iload(nface , ZERO_IMM,args);  queue_load_addr(args,Isize,args);  
  queue_iload(Ls    , ZERO_IMM,args);  queue_load_addr(args,Isize,args);  
  queue_iload(dir   , ZERO_IMM,args);  queue_load_addr(args,Isize,args);  
  queue_iload(Complex_i   , ZERO_IMM,args);  queue_load_addr(args,Isize,args);  
  if(extract){
    queue_iload(commr_p , ZERO_IMM,args);  queue_load_addr(args,Isize,args);  
  }

  alreg(s,Iregs);
  alreg(permreg,Cregs);
  alreg(Mask,Iregs);
  int  Convert1 = args; // Reuse the "args" register since it is no longer needed
  //  alreg(Convert1,Iregs);
  int Convert2;
  if ( half_precision ) {
    Convert2 = allocate_reg(Iregs,"Convert2"); //  Fixme: this will break on KNC due to register file size
  }
  int memory = PROC->StackPointer;

  if ( extract && permute ) {
    printf("cannot extract and permute\n");
    exit(-1);
  }

  for( i = 0;i<12;i++) { 
    need_constant(i*2*SizeofDatum(FourSpinType)*nsimd());
  }
  if ( ! half_precision ) {
    for( i = 0;i<6;i++) { 
      need_constant(i*2*SizeofDatum(TwoSpinType)*nsimd());
    }
  }
  complex_constants_prepare(creg,Complex_i);

  /*Branch to stack restore if length <1*/
  retno = get_target_label();
  check_iterations(nface,retno); 

  if ( half_precision ) {
    queue_iload_imm(Mask,ZERO_IMM);
    queue_ori(Mask,Mask,bits16);
    queue_lshift(Mask,Mask,thirtytwo);
    queue_ori(Mask,Mask,bits16);
    queue_lshift(Mask,Mask,sixteen);
  }

  int permute_mu=3;
  if ( permute ) {
    complex_simd_init(permreg);
    complex_simd_permute(permute_mu,Chi[0][0],Chi[0][0]);// Forces the permreg to get loaded
  }
  if ( extract ) {
    complex_simd_init(permreg);
    // Dropping the simd_permute call forces the load of permute bits inside innermost loop
  }
  // 
  if ( extract ) {
    queue_imul_imm(comm_offset,comm_offset,CHIIMM[1][0][0]);
  } else {
    queue_imul_imm(comm_offset,comm_offset,CHI_IMM);
  }
  queue_imul(comm_offset,comm_offset,Ls);
  queue_iadd(comm_p,comm_offset,comm_p);
  queue_iadd(commr_p,comm_offset,commr_p);

  PROC->latency[FLODPIPE] = 4;
  pragma(DCBT_SPACE,0);

  /**********************************************************************
   * Branch to the loop for THIS DIR based on value in dir register
   **********************************************************************
   */
  int branches[8];
  int branch_continue = get_target_label();
  queue_iadd_imm(s,dir,ZERO_IMM);
  for(int i=0;i<8;i++ ){ 
    branches[i] = get_target_label();
    conditional_branch_cmpzero ( BRANCH_EQ, s, branches[i]) ;
    queue_iadd_imm(s,s,m1);
  }
  make_inst(BRCHPIPE,BRANCH,branch_continue); 
  
  /**********************************************************************
   * Generate the 8 possible innermost loops
   **********************************************************************
   */
  for(int mu=0;mu<4;mu++) {
    for(int pm=0;pm<2;pm++){

      permute_mu=mu;

      int thisdir = 2*mu+pm;
      make_inst(DIRECTIVE,Target,branches[thisdir]);

      /* for i=0, i < nface */
      brchf = start_loop(nface);

      // Gather address calculation
      queue_iload(in_p,ZERO_IMM,sites);
      queue_imul(in_p,in_p,Ls);
      queue_imul_imm(in_p,in_p,PSI_IMM);
      queue_iadd(in_p,in_p,psi_p);


      queue_load_addr(sites,Isize,sites);

      queue_iadd_imm(s,Ls,ZERO_IMM);
      brchs = start_loop(s);


      complex_load(Psi[0][0],PSIIMM[0][0][0],in_p,FourSpinType); 
      complex_load(Psi[0][2],PSIIMM[0][2][0],in_p,FourSpinType); 
      complex_load(Psi[1][1],PSIIMM[1][1][0],in_p,FourSpinType); 
      complex_load(Psi[2][0],PSIIMM[2][0][0],in_p,FourSpinType); 
      complex_load(Psi[2][2],PSIIMM[2][2][0],in_p,FourSpinType); 
      complex_load(Psi[3][1],PSIIMM[3][1][0],in_p,FourSpinType); 

      complex_load(Psi[0][1],PSIIMM[0][1][0],in_p,FourSpinType); 
      complex_load(Psi[1][0],PSIIMM[1][0][0],in_p,FourSpinType); 
      complex_load(Psi[1][2],PSIIMM[1][2][0],in_p,FourSpinType); 
      complex_load(Psi[2][1],PSIIMM[2][1][0],in_p,FourSpinType); 
      complex_load(Psi[3][0],PSIIMM[3][0][0],in_p,FourSpinType); 
      complex_load(Psi[3][2],PSIIMM[3][2][0],in_p,FourSpinType); 

      queue_iadd_imm(in_p,in_p,PSI_IMM);

      if ( mu == 0 ) {
	
	/* Gx
	 *  0 0  0  i    [0]+-i[3]
	 *  0 0  i  0    [1]+-i[2]
	 *  0 -i 0  0
	 * -i 0  0  0
	 *
	 */
	
        if ( pm ==0 ) {
          for(int co=0;co<3;co++) complex_ApiB(Chi[1][co],Psi[1][co],Psi[2][co]);
          for(int co=0;co<3;co++) complex_ApiB(Chi[0][co],Psi[0][co],Psi[3][co]);
        } else {
          for(int co=0;co<3;co++) complex_AmiB(Chi[1][co],Psi[1][co],Psi[2][co]);
          for(int co=0;co<3;co++) complex_AmiB(Chi[0][co],Psi[0][co],Psi[3][co]);
        }
	
      } else if ( mu == 1 ) {
	
	/*Gy
	 *  0 0  0  -1  [0] -+ [3]
	 *  0 0  1  0   [1] +- [2]
	 *  0 1  0  0
	 * -1 0  0  0
	 */
	
        if ( pm ==0 ) {
          for(int co=0;co<3;co++) complex_add(Chi[1][co],Psi[1][co],Psi[2][co]);
          for(int co=0;co<3;co++) complex_sub(Chi[0][co],Psi[0][co],Psi[3][co]);
        } else {
          for(int co=0;co<3;co++) complex_sub(Chi[1][co],Psi[1][co],Psi[2][co]);
          for(int co=0;co<3;co++) complex_add(Chi[0][co],Psi[0][co],Psi[3][co]);
        }
	
      } else if ( mu == 2 ) {
	
	/*Gz
	 *  0 0  i  0   [0]+-i[2]
	 *  0 0  0 -i   [1]-+i[3]
	 * -i 0  0  0
	 *  0 i  0  0
	 */
        if ( pm ==0 ) {
          for(int co=0;co<3;co++) complex_ApiB(Chi[0][co],Psi[0][co],Psi[2][co]);
	  for(int co=0;co<3;co++) complex_AmiB(Chi[1][co],Psi[1][co],Psi[3][co]);
        } else {
          for(int co=0;co<3;co++) complex_AmiB(Chi[0][co],Psi[0][co],Psi[2][co]);
          for(int co=0;co<3;co++) complex_ApiB(Chi[1][co],Psi[1][co],Psi[3][co]);
        }
	
      } else if ( mu == 3 ) {
	
	/*Gt
	 *  0 0  1  0 [0]+-[2]
	 *  0 0  0  1 [1]+-[3]
	 *  1 0  0  0
	 *  0 1  0  0
	 */
        if ( pm ==0 ) {
          for(int co=0;co<3;co++) complex_add(Chi[0][co],Psi[0][co],Psi[2][co]);
	  for(int co=0;co<3;co++) complex_add(Chi[1][co],Psi[1][co],Psi[3][co]);
        } else {
          for(int co=0;co<3;co++) complex_sub(Chi[0][co],Psi[0][co],Psi[2][co]); 
          for(int co=0;co<3;co++) complex_sub(Chi[1][co],Psi[1][co],Psi[3][co]);
        }

      } 

      if ( extract ) {
	int hilo = 0;
	if ( extr_hi ) hilo = 1;

	complex_simd_extract (hilo ,permute_mu,Psi[0][0],Chi[0][0],Chi[0][1]);
	complex_simd_extract (hilo ,permute_mu,Psi[0][1],Chi[0][2],Chi[1][0]);
	complex_simd_extract (hilo ,permute_mu,Psi[0][2],Chi[1][1],Chi[1][2]);
	if ( half_precision ) { 
	  complex_store_half(Psi[0][0],CHIIMM[0][0][0],comm_p,memory,Convert1,Convert2,Mask);
	  complex_store_half(Psi[0][1],CHIIMM[0][1][0],comm_p,memory,Convert1,Convert2,Mask);
	  complex_store_half(Psi[0][2],CHIIMM[0][2][0],comm_p,memory,Convert1,Convert2,Mask);
	} else {
	  complex_store(Psi[0][0],CHIIMM[0][0][0],comm_p,TwoSpinType);
	  complex_store(Psi[0][1],CHIIMM[0][1][0],comm_p,TwoSpinType);
	  complex_store(Psi[0][2],CHIIMM[0][2][0],comm_p,TwoSpinType);
	}

	hilo=1-hilo;
	complex_simd_extract (hilo ,permute_mu,Psi[0][0],Chi[0][0],Chi[0][1]);
	complex_simd_extract (hilo ,permute_mu,Psi[0][1],Chi[0][2],Chi[1][0]);
	complex_simd_extract (hilo ,permute_mu,Psi[0][2],Chi[1][1],Chi[1][2]);
	if ( half_precision ) { 
	  complex_store_half(Psi[0][0],CHIIMM[0][0][0],commr_p,memory,Convert1,Convert2,Mask);
	  complex_store_half(Psi[0][1],CHIIMM[0][1][0],commr_p,memory,Convert1,Convert2,Mask);
	  complex_store_half(Psi[0][2],CHIIMM[0][2][0],commr_p,memory,Convert1,Convert2,Mask);
	} else {
	  complex_store(Psi[0][0],CHIIMM[0][0][0],commr_p,TwoSpinType);
	  complex_store(Psi[0][1],CHIIMM[0][1][0],commr_p,TwoSpinType);
	  complex_store(Psi[0][2],CHIIMM[0][2][0],commr_p,TwoSpinType);
	}
	
	queue_iadd_imm(comm_p,comm_p,CHIIMM[1][0][0]);
	queue_iadd_imm(commr_p,commr_p,CHIIMM[1][0][0]);

      } else {

	for(int spin=0;spin<2;spin++) {
	for(int co=0;co<3;co++) {
	    if ( permute ) { 
	     complex_simd_permute(permute_mu,Chi[spin][co],Chi[spin][co]);
            }
	    if ( half_precision ) { 
	      complex_store_half(Chi[spin][co],CHIIMM[spin][co][0],comm_p,memory,Convert1,Convert2,Mask);
	    } else { 
	      complex_store(Chi[spin][co],CHIIMM[spin][co][0],comm_p,TwoSpinType);
	    }
	}
	}
	queue_iadd_imm(comm_p,comm_p,CHI_IMM);
      }

      // end loop
      stop_loop(brchs,s);
      stop_loop(brchf,nface);
      make_inst(BRCHPIPE,BRANCH,branch_continue); 

    }
  }
  make_inst(DIRECTIVE,Target,branch_continue);

  // where to jump once loop is over
  make_inst(DIRECTIVE,Target,retno);

  // Cleanup
  queue_isub_imm(PROC->StackPointer,PROC->StackPointer,bias);
  restore_regs();
  free_stack();

  // finish
  make_inst(DIRECTIVE,Exit_Routine,name);

  return;
}








