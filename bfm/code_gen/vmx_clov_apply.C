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

void clov_apply( char *);

/*Options flags*/
int human = 0;
char name[80];
int dagger = 0;

Datum FourSpinType=Double;
Datum CloverType=Double;
int addto = 0;
char procname[80]="UNSPECIFIED";

int main ( int argc, char *argv[])
{
    struct processor *myproc;
    int arg;
    char *c;

    name[0] = '\0';
    /*Process command line*/
    while ( ( arg = getopt(argc,argv,"aRn:P:ds")) != EOF)
    {
        switch(arg)
        {
        case 'R':
            human = 1;
            break;
        case 'n':
            if (strlen(optarg)<30) strcpy(name,optarg);
            break;
        case 'P':
            if (strlen(optarg)<30) strcpy(procname,optarg);
            break;
        case 'd':
            dagger = 1;
            break;
        case 'a':
            addto = 1;
            break;
        case 's':
            FourSpinType=Single;
            CloverType=Single;
            break;
        default:
            fprintf(stderr,"Usage: %s -[Rd] -n routine_name\n",argv[0]);
            fprintf(stderr,"\tR -> human readable .S format \n");
            fprintf(stderr,"\td -> dagger decompose\n");
            exit (1);
            break;
        }
    }

    /*Control output according to user set up args*/
    set_processor_optarg(procname);
    set_human_readable(human);
    setup_cmadds();

    /*Write the naive asm code*/
    clov_apply(name);

    /*Filter through an virtual out of order processor */
    schedule_for_proc();

    /*Dump the resulting code*/
    dump_instruction_queue();

    return(0);
}


void clov_apply( char *name)
{
    /*
     * This marks the argument registers as defined by ABI as off limits
     * to us until they are freed by "getarg()";
     */
    int dum = defargcount(1);

    int retno;
    int branchsite;

    /*-------------------------------------------------------------------------------
     * registers used
     *-------------------------------------------------------------------------------
     */

    reg_array_1d(Chi ,Cregs,6); // 2 spionr      - 6 regs
    reg_array_1d(Psi ,Cregs,6); // 2 spinors     - 6 regs
    reg_array_1d(Clo ,Cregs,15); 
        
    offset_2d(PSIIMM,FourSpinType,6,2*nsimd());
    offset_2d(CLOIMM,CloverType,21,2*nsimd());

    /*
     * Integer registers
     */
    alreg(Clo_p,Iregs);    /*Pointer to the current cpt of gauge field    */
    alreg(Chi_p,Iregs);    /*Pointer to the input four spinor             */
    alreg(Psi_p,Iregs);    /*Pointer to current cpt output PSI field      */
    alreg(length,Iregs);   /*number of sites*/
    //alreg(tab,Iregs);      /*Pointer to current entry in offset table*/

    alreg(args,Iregs);

    /*Useful integer immediate constants, in units of Fsize*/
    def_off(ZERO_IMM,Byte,0);
    def_off(SPINOR, FourSpinType, 12*nsimd());
    def_off(CLOVER, CloverType, 42*nsimd());
    
    int Isize      = def_offset(PROC->I_size,Byte,"Isize");
    
    /*--------------------------------------------------------------------
     * Start of the "pseudo assembler proper.
     *--------------------------------------------------------------------
     */
    make_inst(DIRECTIVE,Enter_Routine,name);

    grab_stack(0);
    save_regs();

    /*     * Define our arguments     */
    getarg(args);  /*Pointer to arg list*/

    queue_iload(Chi_p, ZERO_IMM,args);
    queue_load_addr(args,Isize,args);
    queue_iload(Psi_p, ZERO_IMM,args);
    queue_load_addr(args,Isize,args);
    queue_iload(Clo_p,   ZERO_IMM,args);
    queue_load_addr(args,Isize,args);
    queue_iload(length,ZERO_IMM,args);

    for (int i =0; i<6; i++ )
    {
        need_constant(i*2*SizeofDatum(FourSpinType)*nsimd());
    }
    for (int i =0; i<21; i++ )
    {
        need_constant(i*2*SizeofDatum(CloverType)*nsimd());
    }
  
    retno = get_target_label(); /*Branch to exit if length <1*/
    check_iterations(length,retno);

    /*     * Site loop              */
    branchsite = start_loop(length);


    for(int pp=0; pp<6; pp++)
    {
        complex_load(Psi[pp],PSIIMM[pp][0],Psi_p,FourSpinType);
    } 
    for(int qq=0; qq<6; qq++)
    {
        complex_load(Clo[qq],CLOIMM[qq][0],Clo_p,CloverType);
    }

    for(int rr=0; rr<6; rr++)
    {
        complex_mul     (Chi[rr], Clo[rr],  Psi[rr]);
    }

    for(int ss=0; ss<15; ss++)
    {
        complex_load(Clo[ss],CLOIMM[ss+6][0],Clo_p,CloverType);
    }    

    complex_conjmadd(Chi[0], Clo[0],  Psi[1]);
    complex_conjmadd(Chi[0], Clo[1],  Psi[2]);
    complex_conjmadd(Chi[0], Clo[3],  Psi[3]);
    complex_conjmadd(Chi[0], Clo[6],  Psi[4]);
    complex_conjmadd(Chi[0], Clo[10], Psi[5]);

    complex_madd    (Chi[1], Clo[0],  Psi[0]);
    complex_conjmadd(Chi[1], Clo[2],  Psi[2]);
    complex_conjmadd(Chi[1], Clo[4],  Psi[3]);
    complex_conjmadd(Chi[1], Clo[7],  Psi[4]);
    complex_conjmadd(Chi[1], Clo[11], Psi[5]);

    complex_madd    (Chi[2], Clo[1],  Psi[0]);
    complex_madd    (Chi[2], Clo[2],  Psi[1]);
    complex_conjmadd(Chi[2], Clo[5], Psi[3]);
    complex_conjmadd(Chi[2], Clo[8], Psi[4]);
    complex_conjmadd(Chi[2], Clo[12], Psi[5]);

    complex_madd    (Chi[3], Clo[3],  Psi[0]);
    complex_madd    (Chi[3], Clo[4], Psi[1]);
    complex_madd    (Chi[3], Clo[5], Psi[2]);
    complex_conjmadd(Chi[3], Clo[9], Psi[4]);
    complex_conjmadd(Chi[3], Clo[13], Psi[5]);

    complex_madd    (Chi[4], Clo[6], Psi[0]);
    complex_madd    (Chi[4], Clo[7], Psi[1]);
    complex_madd    (Chi[4], Clo[8], Psi[2]);
    complex_madd    (Chi[4], Clo[9], Psi[3]);
    complex_conjmadd(Chi[4], Clo[14], Psi[5]);

    complex_madd    (Chi[5], Clo[10], Psi[0]);
    complex_madd    (Chi[5], Clo[11], Psi[1]);
    complex_madd    (Chi[5], Clo[12], Psi[2]);
    complex_madd    (Chi[5], Clo[13], Psi[3]);
    complex_madd    (Chi[5], Clo[14], Psi[4]);
  

    for(int sp=0; sp<6; sp++)
    {
         complex_store(Chi[sp],PSIIMM[sp][0],Chi_p,FourSpinType);
    }
    
    queue_iadd_imm(Psi_p,Psi_p,SPINOR);
    queue_iadd_imm(Chi_p,Chi_p,SPINOR);
    queue_iadd_imm(Clo_p,Clo_p,CLOVER);

    /* TERMINATION point of the loop*/
    stop_loop(branchsite,length);


    make_inst(DIRECTIVE,Target,retno);

    /*	     EPILOGUE               */

    restore_regs();
    free_stack();
    make_inst(DIRECTIVE,Exit_Routine,name);

    return;
}


