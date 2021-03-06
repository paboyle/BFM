/**
 * \file
 * \brief Provides support for Intel MIC Knights Corner.
 * \author Bernhard Mendl
 *
 *  Copyright Peter Boyle and Glasgow University, 2000.
 *
 *  It is provided under the GNU pubic License V2
 *  It is provided as is and is not guaranteed fit for any purpose.
 */


#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdint.h>
#include "processor.h"
#include "mic.h"


/**
 * If 1, all registers are saved on stack when entering a routine. If 0, only callee-save registers (64-bit Linux) are saved.
 */
//#define SAVE_ALL_REGS 1


extern struct processor *PROC;

static void reset_regs_to_save()
{
    int i;
    for (i = 0; i < PROC->iregs->regfile; i++) PROC->iregs->save_on_stack[i] = 0;
    for (i = 0; i < PROC->fregs->regfile; i++) PROC->fregs->save_on_stack[i] = 0;
}

struct processor *create_knc()
{
    struct regs *iregs = new regs;
    struct regs *fregs = new regs;
    struct processor *proc = new processor;
    struct pipeline *pipes = new pipeline[2];  // 2 hardware pipelines
    struct instructionset *iset = new instructionset;
    int i;

    iregs->regfile = 16;
    fregs->regfile = 32;  // NB: We use 512-bit vector registers here!
    iregs->fmt = "%%r%d";    // Ugly x86 reg names: %rax, %rbx etc. are handled in dump.C::make_reg_name()
    fregs->fmt = "%%zmm%d";  // NB: 512-bit vector regs!

    for (i = 0; i < iregs->regfile; i++)
    {
        iregs->allocated[i] = 0;
        iregs->mnemonic[i] = NULL;
        iregs->save_on_stack[i] = 0;
    }
    for (i = 0; i < fregs->regfile; i++)
    {
        fregs->allocated[i] = 0;
        fregs->mnemonic[i] = NULL;
        fregs->save_on_stack[i] = 0;
    }

    proc->iregs = iregs;
    proc->fregs = fregs;
    proc->pipes = pipes;
    proc->nsimd = 4;  /* Number of COMPLEX elements in a SIMD operation */
    proc->iset = iset;
    proc->FP_size = 8;
    proc->FSP_size = 4;
    proc->I_size = 8;
    proc->IS_size = 4;
    proc->usenops = 1;  // TODO: are there nops?! ICC uses things like "orl %eax, %eax". We set usenops to 1; let's see what happens
    proc->delayslot = 0;  // QUESTION: what does the stuff below mean?!?!
    proc->RegWindows = 0;
    proc->Bias = 0;  // first cache line on a stack frame is used to store the previous, potentially unaligned stack pointer. Bias defines the offset from that on stack is usable for registers...
    proc->lmq_space = 1;
    proc->l1_locking = 0;
    proc->load_reorder = 1;
    proc->cross_negative = 0;
    proc->nmicmask = 8;
    proc->IREGS_SIZE = 8;
    proc->FREGS_SIZE = 64;
    proc->MAX_ARGS_BY_REGS = 6;  // the 7th and any further function arg is passed on the stack (according to x86-64 Linux ABI)


    /*{ ABI definition ... */
    proc->ABI_argbase = 0;  // How many registers are skipped until they are used for passing arguments.
    /* Allocate dedicated registers */
    proc->StackPointer = 7;  // index of stack pointer in regfile
    proc->ReturnPointer = -1;

    /* ABI defines a number of callee save registers */
    proc->Isave = 0;
#if (STACK_POLICY & STACK_POLICY_ALL_IREGS)
    for (i = 0; i < proc->iregs->regfile; i++)
    {
        proc->Isave |= (1<<i);
    }
#elif (STACK_POLICY & STACK_POLICY_LINUX_ABI)
    // save callee-save registers (64-bit Linux ABI)
    proc->Isave |= (1<<8) | (1<<9);  // registers %rbx and %rbp
    proc->Isave |= (1<<12) | (1<<13) | (1<<14) | (1<<15);  // regs %r12..%r15
#endif
    // Some Iregs don't need to be saved:
    proc->Isave &= ~(1<<10);  // %r10 is used as scratch reg in Enter_Routine and is already pushed on stack at beginning of Enter_Routine. No need to save again. See directive_print_knc() and make_reg_name().
    proc->Isave &= ~(1<<(proc->StackPointer));  // SP always used as SP; no need to save on stack

    //fprintf(stderr, "Debug: Isave = %x\n", proc->Isave);


    proc->Fsave = 0;
#if (STACK_POLICY & STACK_POLICY_ALL_FREGS)
    for (i = 0; i < proc->fregs->regfile; i++)
    {
        proc->Fsave = proc->Fsave | (1<<i);
    }
#endif
    /* ...END ABI }*/


    // We mark %rsp as allocated. Then "the programme" cannot use and clobber it.
    proc->iregs->allocated[proc->StackPointer] = 1;
    proc->iregs->mnemonic[proc->StackPointer] = "SP";


    proc->CacheLine = 64;  /* 64 byte cache line */

    proc->npipes = 2;
    for (i = 0; i < N_PIPE_TYPES; i++)
    {
        proc->switchargs[i] = 1;
    }


    pipes[0].issuegroups = SIMDMICMASK | SIMDMICSWIZMASK | SIMDPERMUTEMASK  /* instructions for MIC VPU */
        | FMACMASK | FADDMASK | FMULMASK   /* because we use the VPU for scalar floating-point! */
        | FMOVMASK ;
    pipes[1].issuegroups = IALUMASK | BRCHMASK
        | LOADMASK | STORMASK | CACHMASK | FLODMASK | FSTRMASK;

    proc->latency[FMACPIPE] = 4;
    proc->latency[FADDPIPE] = 4;
    proc->latency[FMULPIPE] = 4;
    proc->latency[FMOVPIPE] = 3;
    proc->latency[FLODPIPE] = 3;
    proc->latency[FSTRPIPE] = 2;
    proc->latency[IALUPIPE] = 2;
    proc->latency[LOADPIPE] = 2;
    proc->latency[STORPIPE] = 2;
    proc->latency[SIMDPERMUTEPIPE] = 6;
    proc->latency[SIMDMICPIPE]     = 4;
    proc->latency[SIMDMICSWIZPIPE] = 6;

    create_knc_instructions(proc);
    return proc;
}


void create_knc_instructions(struct processor *proc)
{
    int i, j;

    struct instruction *fmacs = new instruction [nfmacs];
    struct instruction *fadds = new instruction [nfadds];
    struct instruction *fmuls = new instruction [nfmuls];
    struct instruction *fmovs = new instruction [nfmovs];
    struct instruction *floads= new instruction [nfloads];
    struct instruction *fsaves= new instruction [nfsaves];

    struct instruction *loads = new instruction [nloads];
    struct instruction *cache = new instruction [ncache];
    struct instruction *adds  = new instruction [nadds];
    struct instruction *branch   = new instruction [nbrch];
    struct instruction *directive= new instruction [ndirective];
    struct instruction *stores   = new instruction [nstore];
    struct instruction *mics     = new instruction [nsimdmic];
    struct instruction *micswizs = new instruction [nsimdmicswiz];
    struct instruction *permutes = new instruction [nsimdpermute];
    struct instruction *weird    = new instruction [nweirdops];

    // NB: Instructions' arguments are ordered as in Intel MIC ISA reference manual, i.e. leftmost arg is destination and is instruction->args[0].

    /* Define the FMACS. Just use vector unit. Make scalar store pull out the first element. */
    fmacs[FMADD].mnemonic  = (char *)NULL;  // MIC has 3 variants of FMADD etc. We use 2 of them below.
    fmacs[FMSUB].mnemonic  = (char *)NULL;
    fmacs[FNMADD].mnemonic = (char *)NULL;
    fmacs[FNMSUB].mnemonic = (char *)NULL;
    fmacs[FMADD213].mnemonic  = (char *)"vfmadd213pd  ";  // 1 = +2*1 + 3
    fmacs[FMADD231].mnemonic  = (char *)"vfmadd231pd  ";  // 1 = +2*3 + 1
    fmacs[FNMADD213].mnemonic = (char *)"vfnmadd213pd ";  // 1 = -2*1 + 3
    fmacs[FNMADD231].mnemonic = (char *)"vfnmadd231pd ";  // 1 = -2*3 + 1
    fmacs[FMSUB213].mnemonic  = (char *)"vfmsub213pd  ";  // 1 = +2*1 - 3
    fmacs[FMSUB231].mnemonic  = (char *)"vfmsub231pd  ";  // 1 = +2*3 - 1
    for (i = 0; i < nfmacs; i++) {
        fmacs[i].pipe = FMACPIPE;
        fmacs[i].code = i;
        fmacs[i].narg = 4;  // 3 vectors + 1 mask
        for (j = 0; j < fmacs[i].narg; j++) {
            fmacs[i].argtypes[j] = Fregs;
            fmacs[i].inout[j]    = inarg;
        }
        fmacs[i].argtypes[1] = MICmask;
        fmacs[i].inout[0] = outarg;
        fmacs[i].args[1] = MIC_MASK_SCALFLOP;  // We set vector mask here for convenience.  Then we don't need to pass it to make_inst() every time.  Need to hack make_inst() such that args[1] is NOT overwritten by make_inst()!
    }


    /* Define the FADDS -- Use the vector unit for scalar ops */
    fadds[FADD].mnemonic = (char *)"vaddpd       ";
    fadds[FSUB].mnemonic = (char *)"vsubpd       ";
    for (i = 0; i < nfadds; i++)
    {
        fadds[i].pipe = FADDPIPE;
        fadds[i].code = i;
        fadds[i].narg = 4;  // 3 vectors + 1 mask
        for (j = 0; j < fadds[i].narg; j++) {
            fadds[i].argtypes[j] = Fregs;
            fadds[i].inout[j]    = inarg;
        }
        fadds[i].argtypes[1] = MICmask;
        fadds[i].inout[0] = outarg;
        fadds[i].args[1] = MIC_MASK_SCALFLOP;
    }


    /* Define FMUL -- Use the vector unit for scalar ops */
    fmuls[FMUL].mnemonic = (char *)"vmulpd       ";
    fmuls[FMUL].pipe = FMULPIPE;
    fmuls[FMUL].code = FMUL;
    fmuls[FMUL].narg = 4;  // 3 vec + 1 mask
    for (j = 0; j < fmuls[FMUL].narg; j++) {
        fmuls[FMUL].argtypes[j] = Fregs;
        fmuls[FMUL].inout[j]    = inarg;
    }
    fmuls[FMUL].argtypes[1] = MICmask;
    fmuls[FMUL].inout[0] = outarg;
    fmuls[FMUL].args[1] = MIC_MASK_SCALFLOP;


    /* Define the FMOV */
    fmovs[FMOV].mnemonic = (char *)"vmovapd      ";
    fmovs[FSPLAT].mnemonic = NULL;
    fmovs[FNEG].mnemonic = NULL;
    for (i = 0; i < nfmovs; i++) {
        fmovs[i].pipe = FMOVPIPE;
        fmovs[i].code = i;
        fmovs[i].narg = 3;
        for (j = 0; j < fmovs[i].narg; j++) {
            fmovs[i].argtypes[j] = Fregs;
            fmovs[i].inout[j]    = inarg;
        }
    }
    fmovs[FMOV].inout[0] = outarg;
    fmovs[FMOV].argtypes[1] = MICmask;
    fmovs[FMOV].args[1] = MIC_MASK_SCALFLOP;


    /* Floating point loads */
    floads[FLOAD_EA].mnemonic        = (char *)"vbroadcastsd ";  // load single float64 and broadcast 1to8
    floads[FLOAD_REG_EA].mnemonic    = (char *)"vmovapd      ";
    floads[FSLOAD_EA].mnemonic      = (char *)"vbroadcastss ";  // load single float32 and broadcast 1to16
    floads[CLOAD_INDEXED].mnemonic   = (char *)"vmovapd      ";  // load full float64 vector
    floads[CSLOAD_INDEXED].mnemonic = (char *)"vmovaps      ";  // load full float32 vector

    for (i = 0; i < nfloads; i++) {
        floads[i].pipe = FLODPIPE;
        floads[i].code = i;
        floads[i].narg = 3;
        floads[i].argtypes[0] = Fregs;
        floads[i].argtypes[1] = SIMM;
        floads[i].argtypes[2] = Iregs;
        floads[i].inout[0]  = outarg;
        floads[i].inout[1]  = inarg;
        floads[i].inout[2]  = inarg;
    }
    

    fsaves[FSAVE_EA].mnemonic        = (char *)"vpackstorelpd";  // save single float64
    fsaves[FSAVE_REG_EA].mnemonic    = (char *)"vmovapd      ";  // save full register, no mask!
    fsaves[FSSAVE_EA].mnemonic      = (char *)"vpackstorelps";  // save single float32
    fsaves[CSAVE_INDEXED].mnemonic   = (char *)"vmovnrngoapd ";  // store full vector  NOTE: using no-read hint & non-globally ordered storess
    fsaves[CSSAVE_INDEXED].mnemonic = (char *)"vmovnrngoaps ";
    for (i = 0; i < nfsaves; i++) {
        fsaves[i].pipe = FSTRPIPE;
        fsaves[i].code = i;
        fsaves[i].narg = 4;
        fsaves[i].argtypes[0] = Fregs;    // the vector that contains our scalar float/double
        fsaves[i].argtypes[1] = SIMM;     // offset
        fsaves[i].argtypes[2] = Iregs;    // address register
        fsaves[i].argtypes[3] = MICmask;  // mask
        for (j = 0; j < fsaves[i].narg; j++) {
            fsaves[i].inout[j] = inarg;
        }
        fsaves[i].args[3] = MIC_MASK_SCALFLOP;  // NB: last arg (mask) needn't be set in make_inst().
    }
    fsaves[CSAVE_INDEXED].narg   = 3;  // no mask applied
    fsaves[CSSAVE_INDEXED].narg = 3;
    fsaves[FSAVE_REG_EA].narg = 3;  // no mask, save full reg




    /* Define the permute pipe */
    permutes[VPERM].mnemonic    = (char *)"vpermd       ";  // shuffles arbitrarily according to index register
    permutes[VPERMIMM].mnemonic = (char *)"vpshufd      ";  // shuffle within lanes according to immediate pattern
    permutes[VPERMCOARSE].mnemonic = (char *)"vpermf32x4   ";  // permute across lanes
    for (i = 0; i < nsimdpermute; i++)
    {
        permutes[i].pipe = SIMDPERMUTEPIPE;
        permutes[i].code = i;
        permutes[i].narg = 3;
        for (j = 0; j < permutes[i].narg; j++)
        {
            permutes[i].argtypes[j] = Fregs;
            permutes[i].inout[j] = inarg;
        }
        permutes[i].inout[0] = outarg;
    }
    permutes[VPERMIMM].argtypes[2] = SIMM;
    permutes[VPERMCOARSE].narg = 3;
    permutes[VPERMCOARSE].inout[2] = inarg;
    permutes[VPERMCOARSE].argtypes[1] = SIMM;
    permutes[VPERMCOARSE].argtypes[2] = Fregs;


    /* Define the MIC support */
    /* unswizzled */
    mics[MIC_MOV].mnemonic       = (char *)"vmovapd      ";  // XXX: VPORQ or VMOVAPD for a simple reg-reg copy? What about latencies? Don't forget to adjust narg etc.
    mics[MIC_ADD].mnemonic       = (char *)"vaddpd       ";
    mics[MIC_SUB].mnemonic       = (char *)"vsubpd       ";
    mics[MIC_MUL].mnemonic       = (char *)"vmulpd       ";
    mics[MIC_FMADD213].mnemonic  = (char *)"vfmadd213pd  ";
    mics[MIC_FMADD231].mnemonic  = (char *)"vfmadd231pd  ";
    mics[MIC_FNMADD213].mnemonic = (char *)"vfnmadd213pd ";
    mics[MIC_FNMADD231].mnemonic = (char *)"vfnmadd231pd ";
    mics[MIC_FMSUB213].mnemonic  = (char *)"vfmsub213pd  ";
    mics[MIC_FMSUB231].mnemonic  = (char *)"vfmsub231pd  ";
    mics[MIC_AND_MEMOP].mnemonic = (char *)"vpandd       ";
    mics[MIC_ADD_MEMOP].mnemonic = (char *)"vaddpd       ";
    mics[MIC_SUB_MEMOP].mnemonic = (char *)"vsubpd       ";
    mics[MIC_MUL_MEMOP].mnemonic = (char *)"vmulpd       ";
    mics[MIC_FMADD231_MEMOP].mnemonic  = (char *)"vfmadd231pd  ";
    mics[MIC_XOR].mnemonic       = (char *)"vpxorq       ";  // int64 XOR
    for (i = 0; i < nsimdmic; i++)
    {
        mics[i].pipe = SIMDMICPIPE;
        mics[i].code = i;
        mics[i].narg = 3;
        for (j = 0; j < mics[i].narg; j++)
        {
            mics[i].argtypes[j] = Fregs;
            mics[i].inout[j]    = inarg;
        }
        mics[i].inout[0] = outarg;
    }
    mics[MIC_MOV].narg = 2;

    for (i=MIC_AND_MEMOP;i<=MIC_FMADD231_MEMOP;i++){
        mics[i].narg = 4;
	mics[i].argtypes[0]=Fregs; 
	mics[i].argtypes[1]=Fregs; 
	mics[i].argtypes[2]=SIMM;     // Mem operand, base+immediate offset.
	mics[i].argtypes[3]=Iregs;    // Mem operand, base+immediate offset.
    }


    /* swizzled, unmasked */
    micswizs[MIC_ADD_SWIZ].mnemonic = (char *)"vaddpd       ";
    micswizs[MIC_SUB_SWIZ].mnemonic = (char *)"vsubpd       ";
    micswizs[MIC_MUL_SWIZ].mnemonic = (char *)"vmulpd       ";
    micswizs[MIC_FMADD231_SWIZ].mnemonic = (char *)"vfmadd231pd       ";
    for (i = 0; i < MIC_ADD_MASK_SWIZ; i++)  // walk through all unmasked ops
    {
        micswizs[i].pipe = SIMDMICSWIZPIPE;
        micswizs[i].code = i;
        micswizs[i].narg = 4;
        for (j = 0; j < micswizs[i].narg; j++)
        {
            micswizs[i].argtypes[j] = Fregs;
            micswizs[i].inout[j]    = inarg;
        }
        micswizs[i].inout[0] = outarg;
        micswizs[i].argtypes[3] = MICswiz;
    }

    /* swizzled, masked */
    micswizs[MIC_ADD_MASK_SWIZ].mnemonic  = (char *)"vaddpd       ";
    micswizs[MIC_SUB_MASK_SWIZ].mnemonic  = (char *)"vsubpd       ";
    micswizs[MIC_SUBR_MASK_SWIZ].mnemonic = (char *)"vsubrpd      ";
    micswizs[MIC_BLEND_MASK_SWIZ].mnemonic = (char *)"vblendmpd    ";
    micswizs[MIC_VPERMF32X4].mnemonic = (char *)"vpermf32x4    ";
    micswizs[CLOAD_SPLAT].mnemonic = (char *)"vbroadcastsd    ";
    micswizs[CSLOAD_SPLAT].mnemonic = (char *)"vbroadcastss   ";

    for (i = MIC_ADD_MASK_SWIZ; i < nsimdmicswiz; i++)  // walk through all masked ops
    {
        micswizs[i].pipe = SIMDMICSWIZPIPE;
        micswizs[i].code = i;
        micswizs[i].narg = 5;
        for (j = 0; j < micswizs[i].narg; j++)
        {
            micswizs[i].argtypes[j] = Fregs;
            micswizs[i].inout[j]    = inarg;
        }
        micswizs[i].inout[0]    = outarg;
        micswizs[i].argtypes[1] = MICmask;
        micswizs[i].argtypes[4] = MICswiz;
    }
    micswizs[MIC_VPERMF32X4].narg=4;
    micswizs[MIC_VPERMF32X4].argtypes[3]=SIMM;

    micswizs[CLOAD_SPLAT].narg=4;
    micswizs[CLOAD_SPLAT].argtypes[0]=Fregs;
    micswizs[CLOAD_SPLAT].argtypes[1]=MICmask;
    micswizs[CLOAD_SPLAT].argtypes[2]=SIMM;
    micswizs[CLOAD_SPLAT].argtypes[3]=Iregs;

    micswizs[CSLOAD_SPLAT].narg=4;
    micswizs[CSLOAD_SPLAT].argtypes[0]=Fregs;
    micswizs[CSLOAD_SPLAT].argtypes[1]=MICmask;
    micswizs[CSLOAD_SPLAT].argtypes[2]=SIMM;
    micswizs[CSLOAD_SPLAT].argtypes[3]=Iregs;


    /* Define loads */
    loads[ILD_EA].mnemonic       = (char *)"movq         ";
    loads[LOAD_REG_EA].mnemonic  = (char *)"movq         ";
    loads[LOAD_ADDR].mnemonic    = (char *)"addq         ";
    loads[LOAD_IMM].mnemonic     = (char *)"movq         ";
    loads[ILD_EA_SHORT].mnemonic = NULL;  // See fix_short_integer()
    for (i = 0; i < nloads; i++)
    {
        loads[i].narg = 2;
        loads[i].argtypes[0] = Iregs;
        loads[i].argtypes[1] = SIMM;
        loads[i].inout[0] = outarg;
        loads[i].inout[1] = inarg;
        loads[i].pipe  = LOADPIPE;
        loads[i].code  = i;
    }
    loads[LOAD_REG_EA].narg = 3;
    loads[LOAD_REG_EA].argtypes[2] = Iregs;
    loads[LOAD_REG_EA].inout[2] = inarg;
    loads[ILD_EA].narg = 3;
    loads[ILD_EA].argtypes[2] = Iregs;
    loads[ILD_EA].inout[2]    = inarg;




    /*Define STORE - really the register save and restore*/
    /*Since these must be full register width we abstract them*/
    /* differently from the floating and integer which can be a lower precision*/
    stores[SAVE_REG_EA].mnemonic = (char *)"movq         ";
    stores[IST_EA].mnemonic      = (char *)"movq         ";
    for (i = 0; i < nstore; i++)
    {
        stores[i].narg = 3;
        stores[i].argtypes[0] = Iregs;
        stores[i].argtypes[1] = SIMM;
        stores[i].argtypes[2] = Iregs;
        stores[i].inout[0] = inarg;
        stores[i].inout[1] = inarg;
        stores[i].inout[2] = inarg;
        stores[i].pipe  = STORPIPE;
        stores[i].code  = i;
    }


    /* Define cache control codes */
    // TODO: use exclusive or NT prefetches?
    // TODO: use clflush instead of clevictn?
    cache[PREF_IMM].mnemonic   = (char *)"vprefetch0   ";  // prefetch into MIC L1d (MIC ISA ref man p.641)
    cache[WRITE_HINT].mnemonic = (char *)"clevict0     ";  // Write hint, evict from L1d? (MIC ISA ref man p.624)  // XXX: or use exclusive prefetch here??!
    cache[FLUSH].mnemonic = (char *)"clevict0 ";  // evict from L2 (MIC ISA ref man p.626)
    cache[FLUSHL2].mnemonic = (char *)"clevict1 ";  // evict from L2 (MIC ISA ref man p.626)
    //cache[TOUCHLOCKSET].mnemonic = "dcbtls ";
    //cache[TOUCHLOCKCLEAR].mnemonic = "dcblc ";
    cache[TOUCHL2].mnemonic = (char *)"vprefetch1 ";  // prefetch into MIC L2 (MIC ISA ref man p.643)
    //cache[COUNTER_HINT].mnemonic = "mtctr ";

    cache[PREF0].mnemonic = (char *)"vprefetch0   ";
    cache[PREF1].mnemonic = (char *)"vprefetch1   ";
    cache[PREF0NT].mnemonic = (char *)"vprefetchnta ";
    cache[PREF1NT].mnemonic = (char *)"vprefetch2   ";
    for (i = 0; i < ncache; i++)
    {
        cache[i].narg        = 2;
        cache[i].argtypes[0] = SIMM;
        cache[i].argtypes[1] = Iregs;
        cache[i].inout[0]    = inarg;
        cache[i].inout[1]    = inarg;
        cache[i].pipe  = CACHPIPE;
        cache[i].code  = i;
    }


    /* Integer arithmetic pipe */
    // NB: checked against Intel docs [Intel 64 arch software devel manual, 2013-03, 253665-046US] and [Intel Xeon Phi ISA ref manual, 2012-09, 327364-001en]
    adds[IADD].mnemonic   = (char *)"addq         ";
    adds[ISUB].mnemonic   = (char *)"subq         ";
    adds[IOR].mnemonic    = (char *)"orq          ";
    adds[CMOVGE].mnemonic = (char *)"cmovae";
    adds[CMOVGT].mnemonic = (char *)"cmova";
    adds[CMOVLE].mnemonic = (char *)"cmovbe";
    adds[CMOVLT].mnemonic = (char *)"cmovb";
    adds[ISUB_PREDICATE].mnemonic = NULL;  // QUESTION: never used?
    adds[IOR_PREDICATE].mnemonic  = (char *)"orq          ";
    adds[IMUL].mnemonic     = (char *)"imulq        ";
    adds[IMUL_IMM].mnemonic = (char *)"imulq        ";
    adds[IADD_IMM].mnemonic = (char *)"addq         ";
    adds[IAND_IMM].mnemonic = (char *)"andq         ";
    adds[MOV].mnemonic      = (char *)"movq         ";
    adds[INEG].mnemonic     = (char *)"negq         ";
    adds[RSHIFT].mnemonic   = (char *)"shrq         ";  // logical / unsigned shift
    for (i = 0; i < nadds; i++) {
        adds[i].narg = 2;
        for (j = 0; j < adds[i].narg; j++) {
            adds[i].argtypes[j] = Iregs;
            adds[i].inout[j] = inarg;
        }
        adds[i].inout[0] = outarg;
        adds[i].pipe = IALUPIPE;
        adds[i].code = i;
    }
    adds[INEG].narg = 1;
    adds[IMUL_IMM].narg = 3;
    adds[IMUL_IMM].argtypes[1] = SIMM;
    adds[IMUL_IMM].argtypes[2] = Iregs;
    adds[IMUL_IMM].inout[2] = inarg;

    adds[IADD_IMM].argtypes[1] = SIMM;
    adds[IAND_IMM].argtypes[1] = SIMM;
    adds[RSHIFT].argtypes[1] = SIMM;




    /*Branch instructions*/
    branch[BRANCH_GT].mnemonic = (char *)"jg    ";
    branch[BRANCH_EQ].mnemonic = (char *)"je    ";
    branch[BRANCH_LT].mnemonic = (char *)"jl    ";
    branch[BRANCH].mnemonic    = (char *)"jmp   ";
    branch[BRANCH_GE].mnemonic = (char *)"jge   ";
    branch[BRANCH_LE].mnemonic = (char *)"jle   ";
    branch[BRANCH_CTR].mnemonic= (char *)NULL;
    branch[BRANCH_RET].mnemonic= (char *)"ret   ";
    for (i = 0; i < nbrch; i++)
    {
        branch[i].narg = 1;
        branch[i].argtypes[0] = Label;
        branch[i].inout[0] = inarg;
        branch[i].pipe = BRCHPIPE;
        branch[i].code = i;
    }

    /*Assembler directives*/
    directive[Target].mnemonic = (char *)"";
    directive[CmovTarget].mnemonic = (char *)"";
    directive[Cache_Align].mnemonic = (char *)"";
    directive[Enter_Routine].mnemonic = (char *)"";
    directive[Exit_Routine].mnemonic = (char *)"";
    directive[Pipe_Flush].mnemonic = (char *)"";
    directive[NOP].mnemonic  = (char *)"nop   ";
    directive[FNOP].mnemonic = (char *)"fnop  ";
    directive[LS_BARRIER].mnemonic = (char *)"";
    directive[COMMENT].mnemonic    = (char *)"\t# ";

    for (i = 0; i < ndirective; i++)
    {
        directive[i].narg = 1;
        directive[i].argtypes[0] = Label;
        directive[i].pipe  = DIRECTIVE;
        directive[i].code  = i;
    }
    weird[WARSTALL].mnemonic     = (char *)"Write after read";
    weird[LSPIPESPARSE].mnemonic = (char *)"L/S pipe busy";
    weird[CACHESTALL].mnemonic   = (char *)"Stalled by Cache unit";
    for (i = 0; i < nweirdops; i++)
    {
        weird[i].narg = 0;
        weird[i].argtypes[0] = Iregs;
        weird[i].argtypes[1] = Iregs;
        weird[i].pipe  = WEIRDRULES;
        weird[i].code  = i;
    }

    /*Now set these up as the instruction set*/
    proc->iset->instructions[SIMDHUMMERPIPE] = NULL;
    proc->iset->instructions[SIMDVMXPIPE] = NULL;
    proc->iset->instructions[SIMDMICPIPE] = mics;
    proc->iset->instructions[SIMDMICSWIZPIPE] = micswizs;
    proc->iset->instructions[SIMDPERMUTEPIPE] = permutes;
    proc->iset->instructions[FMACPIPE] = fmacs;
    proc->iset->instructions[FADDPIPE] = fadds;
    proc->iset->instructions[FMULPIPE] = fmuls;
    proc->iset->instructions[CACHPIPE] = cache;
    proc->iset->instructions[IALUPIPE] = adds;
    proc->iset->instructions[LOADPIPE] = loads;
    proc->iset->instructions[FLODPIPE] = floads;
    proc->iset->instructions[FSTRPIPE] = fsaves;
    proc->iset->instructions[BRCHPIPE] = branch;
    proc->iset->instructions[DIRECTIVE] = directive;
    proc->iset->instructions[STORPIPE] = stores;
    proc->iset->instructions[FMOVPIPE] = fmovs;
    proc->iset->instructions[WEIRDRULES] = weird;
}


static int Hdecrement;


int start_loop_knc(int counter)
{
    (void)counter;  // suppresses "unused parameter" compiler warning
    int lab = get_target_label();               /* Target label of head of loop */
    Hdecrement = def_offset(-1,Byte,"minus1");  /* Defines counter decrement value */
    make_inst(DIRECTIVE, Target, lab);
    return lab;
}


void stop_loop_knc(int branchno, int counter)
{
  queue_load_addr(counter,Hdecrement,counter);        /* Decrements counter by Hdecrement */
  make_inst(IALUPIPE,IOR_PREDICATE,counter,counter);  /* Sets the predicate reg*/
  make_inst(BRCHPIPE,BRANCH_GT,branchno);
}


void check_iterations_knc(int cntreg, int retno)
{
    make_inst(IALUPIPE, IOR_PREDICATE,cntreg,cntreg); /*Sets the predicate reg*/
    make_inst(BRCHPIPE, BRANCH_LE,retno);   /*Branch if cntreg <= 0*/
}


void directive_print_knc(struct instruction *i)
{
    uint64_t mask_patterns;
    const char *const SCRATCH = "r10";  // the register used for doing stack alignment and MIC mask preparation. MUST NOT clobber any input register, should not use be %rax or %rdx, which might be hold return values.

    switch (i->code)
    {
    case Target:
        printf("%s:\n", make_label(i->args[0]));
        break;
    case CmovTarget:
        printf("%s:\n", make_label(i->args[0]));
        break;
    case Cache_Align:
        printf(".align 4\n");  // Note: 4 ok? Doesn't seem to be important anyway because this directive is apparently never issued on MIC (and BG/L, BG/Q).
        break;


        // We must align the stack to 64 Bytes! Otherwise MIC vector regs (fregs) will cause bus errors.
        // We must store the former value of stack pointer. In order not to loose a valuable register, we store it on stack. Take care for the offset!
    case Enter_Routine:
        printf(".text\n");
        printf(".align 16,0x90\n");
        printf(".globl %s\n", i->mnemonic);
        printf("%s:\n", i->mnemonic);

        printf("\tpushq    %%%s\n", SCRATCH);           // save scratch reg on stack
        printf("\tmovq     %%rsp, %%%s\n", SCRATCH);    // save SP to scratch reg
        printf("\tandq     $-64, %%rsp\n");             // align stack pointer
        printf("\tmovq     %%%s, (%%rsp)\n", SCRATCH);  // push previous SP onto stack

        // define / set any mask registers we ever use in bagel...
        mask_patterns = (0xF0ULL<<48)|(0xCCULL<<32) | (0x1<<16UL) | (0xAAUL);
        printf("\tmovq     $0x%lx, %%%s\n", mask_patterns, SCRATCH);  // NB: this register holds up to 4 mask patterns . Use %rax to not clobber any input registers!
        printf("\tkextract $0, %%%s, %%k%d   # even elements\n", SCRATCH, MIC_MASK_EVEN);      // mask is 170 = 0xAA
        printf("\tkextract $2, %%%s, %%k%d   # even pairs\n", SCRATCH, MIC_MASK_EVEN_PAIRS);  // mask is 204 = 0xCC
        printf("\tkextract $1, %%%s, %%k%d   # scalar flop\n", SCRATCH, MIC_MASK_SCALFLOP);  // mask is 1
        printf("\tkextract $3, %%%s, %%k%d   # scalar flop\n", SCRATCH, MIC_MASK_MU3);  // mask is 0xF0
        printf("\tknot     %%k%d,  %%k%d       # odd elements\n", MIC_MASK_EVEN, MIC_MASK_ODD);  // mask is 85  = 0x55; NB: we negate MIC_MASK_EVEN to save space in %rax
        printf("\tknot     %%k%d,  %%k%d       # odd pairs\n", MIC_MASK_EVEN_PAIRS, MIC_MASK_ODD_PAIRS);

	//        mask_patterns = (0xcc<<32) | (0x1<<16) | (0xaa);
	//        printf("\tmovq     $%lu, %%%s\n", mask_patterns, SCRATCH);  // NB: this register holds up to 4 mask patterns . Use %rax to not clobber any input registers!

        printf("\t# ----- end of Enter_Routine -----\n\n");
        break;
    case Exit_Routine:
        printf("\n\t# ----- Exit_Routine -----\n");
        printf("\tmovq     (%%rsp), %%rsp\n");  // restore previous stack pointer
        printf("\tpopq     %%%s\n", SCRATCH);           // restore scratch reg
        printf("\tret\n");
        printf(".align 16,0x90\n");
        printf(".type %s,@function\n", i->mnemonic);
        printf(".size %s,.-%s\n", i->mnemonic, i->mnemonic);
        printf(".data\n");
        reset_regs_to_save();
        break;
    case Pipe_Flush:
        printf("# Instruction order barrier inserted\n");
        break;
    case NOP:
        printf("  nop\n");
        break;
    case FNOP:
        printf("  fnop\n");
        break;
    case COMMENT:
        printf("%s%s\n", i->mnemonic, (char *)(i->args[0]));
        break;
    default:
        printf("Error: Unknown assembler directive requested %d\n", i->code);
        exit(0);
        break;
    }
}


#if 1
void queue_cmovge_knc(int cond,int src,int dst)
{
  int skipto = get_target_label();

  make_inst(IALUPIPE,IOR_PREDICATE,cond,cond,cond); /*Sets the predicate reg*/
  make_inst(BRCHPIPE,BRANCH_LT,skipto);             /*Branch if cntreg < 0*/
  make_inst(IALUPIPE,MOV,dst,src);              /*Copy if cond >= 0 */
  make_inst(DIRECTIVE,CmovTarget,skipto);
}

void queue_cmovgt_knc(int cond,int src,int dst)
{
  int skipto = get_target_label();

  make_inst(IALUPIPE,IOR_PREDICATE,cond,cond,cond); /*Sets the predicate reg*/
  make_inst(BRCHPIPE,BRANCH_LE,skipto);            /*Branch if cntreg <= 0*/
  make_inst(IALUPIPE,MOV,dst,src);             /*Copy if cond > 0 */
  make_inst(DIRECTIVE,CmovTarget,skipto);
}

void queue_cmovlt_knc(int cond,int src,int dst)
{
  int skipto = get_target_label();

  make_inst(IALUPIPE,IOR_PREDICATE,cond,cond,cond); /*Sets the predicate reg*/
  make_inst(BRCHPIPE,BRANCH_GE,skipto);             /*Branch if cntreg >= 0*/
  make_inst(IALUPIPE,MOV,dst,src);              /*Copy if cond < 0 */
  make_inst(DIRECTIVE,CmovTarget,skipto);
}

void queue_cmovle_knc(int cond,int src,int dst)
{
  int skipto = get_target_label();

  make_inst(IALUPIPE,IOR_PREDICATE,cond,cond,cond); /*Sets the predicate reg*/
  make_inst(BRCHPIPE,BRANCH_GT,skipto);             /*Branch if cond > 0*/
  make_inst(IALUPIPE,MOV,dst,src);              /*Copy if cond <= 0 */
  make_inst(DIRECTIVE,CmovTarget,skipto);
}
#else

void queue_cmovgt_knc(int cond, int src, int dst)
{
    make_inst(IALUPIPE,IOR_PREDICATE,cond,cond);
    make_inst(IALUPIPE,CMOVGT,dst,src);
}


void queue_cmovlt_knc(int cond, int src, int dst)
{
    make_inst(IALUPIPE,IOR_PREDICATE,cond,cond);
    make_inst(IALUPIPE,CMOVLT,dst,src);
}


void queue_cmovle_knc(int cond, int src, int dst)
{
    make_inst(IALUPIPE,IOR_PREDICATE,cond,cond);
    make_inst(IALUPIPE,CMOVLE,dst,src);
}


void queue_cmovge_knc(int cond, int src, int dst)
{
    make_inst(IALUPIPE,IOR_PREDICATE,cond,cond);
    make_inst(IALUPIPE,CMOVLT,dst,src);
}
#endif

void cache_print_knc(struct instruction *i)
{
    switch (i->code)
    {
    case TOUCHL2:  // FIXME
        printf("\t%s 2, \t", i->mnemonic);
        arg_print(i, 0);
        printf(",");
        arg_print(i, 1);
        printf("\n");
        break;
    case FLUSH:  // FIXME
    case FLUSHL2:  // FIXME
        printf("\t%s \t", i->mnemonic);
        arg_print(i, 0);
        printf("(");
        arg_print(i, 1);
        printf(")\n");
        break;
    case TOUCHLOCKSET:
    case TOUCHLOCKCLEAR:
    case WRITE_HINT:
    case PREF_IMM:
        printf("\t%s ", i->mnemonic);
        arg_print(i, 0);
        printf("(");
        arg_print(i, 1);
        printf(")\n");
        break;
    case COUNTER_HINT:
        printf("\t%s\t", i->mnemonic);
        arg_print(i, 0);
        printf("\n");
        break;

    case PREF0:
    case PREF1:
    case PREF0NT:
    case PREF1NT:
        printf("\t%s ", i->mnemonic);
        arg_print(i, 0);
        printf("(");
        arg_print(i, 1);
        printf(")\n");
        break;

    default:
        printf("Error: Unknown cache hint requested\n");
        exit(0);
        break;
    }
}



