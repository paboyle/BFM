#ifndef _PROCESSOR_H_
#define _PROCESSOR_H_
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "prefetch.h"

/*
*
*  Copyright Peter Boyle, University of Glasgow, 2000.
*  This software is provided for NON-COMMERCIAL use only,
*  It is provided under the GNU pubic License V2
*  It is provided as is and is not guaranteed fit for any purpose.
*
*/

/*Some limits for convenience for multi-dimensional arrays*/
#define MAX_WIDTH 8   /* QUESTION: maximum number of physical pipelines */
#define MAX_CYCLES (2048*64)  // QUESTION: arbitrary number? What is the typ. number of cycles?
#define MAX_REGS 128

#define MAX_INST_ARG 7  /* maximum number of arguments of instructions */
#define MAX_FUNC_ARG 10  /* maximum number of args of functions (not instructions) */

enum Datum { Byte, Double, Single, Integer, ShortInteger, Half };
int SizeofDatum(enum Datum t);
void test_proc(void);
void initialise_streams(void);
void streams_reg_alloc(int val);

/*Define the types of instruction - later called issue groups*/
/*Philosophy is to define all types that are likely*/
/*to map to execution units */
enum EXECUTIONPIPES
{
    FMACPIPE = 0,
    FADDPIPE,//1
    FMULPIPE,//2
    FMOVPIPE,//3
    FLODPIPE,//4
    FSTRPIPE,//5
    IALUPIPE,//6
    LOADPIPE,//7
    STORPIPE,//8
    BRCHPIPE,//9
    CACHPIPE,//10
    DIRECTIVE,//11
    WEIRDRULES,//12
    PRAGMA, /*Allows optimisation hints to be passed to the scheduler*/ //13

/*BlueGene/L,P support*/
    SIMDHUMMERPIPE,//14
    SIMDPERMUTEPIPE,//15

/*Altivec support*/
    SIMDVMXPIPE,//16

    /* Intel MIC support: swizzled & unswizzled issuegroups allow easier handling of instructions / operands. Also easier to handle different latencies! */
    SIMDMICPIPE,//17    no swizzles
    SIMDMICSWIZPIPE,  // with swizzles

    N_PIPE_TYPES
};

/*Bit mask to check if we can issue to a pipe*/
#define FMACMASK (1<<FMACPIPE)
#define FADDMASK (1<<FADDPIPE)
#define FMULMASK (1<<FMULPIPE)
#define FMOVMASK (1<<FMOVPIPE)
#define IALUMASK (1<<IALUPIPE)
#define LOADMASK (1<<LOADPIPE)
#define BRCHMASK (1<<BRCHPIPE)
#define STORMASK (1<<STORPIPE)
#define FLODMASK (1<<FLODPIPE)
#define FSTRMASK (1<<FSTRPIPE)
#define CACHMASK (1<<CACHPIPE)
#define SIMDHUMMERMASK (1<<SIMDHUMMERPIPE)
#define SIMDPERMUTEMASK (1<<SIMDPERMUTEPIPE)
#define SIMDVMXMASK (1<<SIMDVMXPIPE)
#define SIMDMICMASK (1<<SIMDMICPIPE)
#define SIMDMICSWIZMASK (1<<SIMDMICSWIZPIPE)

/*Define the instructions we may want to use*/

enum FMAC_OPS
{
    FMADD,
    FMSUB,
    FNMADD,
    FNMSUB,
    FMADD213,  // using MIC's VPU for scalar floating-point ops
    FMADD231,
    FNMADD213,
    FNMADD231,
    FMSUB213,
    FMSUB231,
    nfmacs
} ;
enum FMOV_OPS
{
    FMOV,
    FNEG,
    FSPLAT,
    nfmovs
};
enum FADD_OPS
{
    FADD,
    FSUB,
    nfadds
};

enum FMUL_OPS
{
    FMUL,
    nfmuls
};

enum LOAD_OPS{
    ILD_EA,
    LOAD_REG_EA,  // load reg from stack
    LOAD_ADDR,
    LOAD_IMM,
    ILD_EA_SHORT,
    nloads
};

enum BRANCH_OPS{
    BRANCH_GT,
    BRANCH_EQ,
    BRANCH_LT,
    BRANCH,
    BRANCH_CTR, /*PowerPC specific - try to use as platform independent*/
    /* a paradigm as possible unless there really *is* a performance hit*/
    BRANCH_GE,
    BRANCH_LE,
    BRANCH_RET,
    nbrch
};

enum STORE_OPS{
    SAVE_REG_EA,  // store reg on stack
    IST_EA,
    nstore
};
enum FLOAD_OPS{
    FLOAD_EA,     /*Scalar load - splat to all elements in SIMD case*/
    FLOAD_REG_EA,  // load reg from stack
    CLOAD_INDEXED,/*Complex load (or SIMD load)*/
    FSLOAD_EA, /*Single precision load stores for bandwidth*/
    CSLOAD_INDEXED,
    VLOAD_INDEXED,/*VMX*/
    CLOAD_SPLAT_INDEXED,/*Complex load (or SIMD load)*/
    CSLOAD_SPLAT_INDEXED,/*Complex load (or SIMD load)*/
    nfloads
};
enum FSAVE_OPS{
    FSAVE_EA,
    FSAVE_REG_EA,  // store reg on stack
    CSAVE_INDEXED,
    FSSAVE_EA,
    CSSAVE_INDEXED,
    VSAVE_INDEXED,/*VMX*/
    nfsaves
};
enum WEIRD_OPS{ /*Only used if we want to print dependencies */
    /*These are fictitious instructions which denote stall type*/
    /*when the stall is not caused by a single other instruction*/
    /*Only supported on PPC 440 currently                       */
    WARSTALL, LSPIPESPARSE, CACHESTALL, nweirdops
} ;
enum PRAGMA_OPS{
    /*Scheduler hints*/
    DCBT_SPACE,  /*Space between DCBT's or equivalent  */
    DCBT_PRE,    /*No L/S use for N-cycles before DCBT */
    DCBT_POST,   /*No L/S use for N-cycles after  DCBT */
    LOAD_LIM,    /*Insert a LS pipe bubble if N loads in a row */
    STORE_LIM,   /*Insert a LS pipe bubble if N stores in a row */
    LS_LIM,      /*Insert a LS pipe bubble if N loads or stores in a row */
    LOAD_SPACE,
    STORE_SPACE,
    STORE_INORDER,   /*Preserve order of stores */
    npragmas
};

/*CACHE CONTROL INSTRUCTIONS - NEED TO THINK ABOUT PRECISION*/
/*AND CACHELINE INFORMATION     */
/*I added the mtctr ppc op here */
/*This is really becoming a 1 argument category with no effect on "architectural" regs*/
/*Should really change the PIPE category name to reflect this*/
enum CACHE_OPS{
    PREF_IMM,
    WRITE_HINT,
    COUNTER_HINT,
    FLUSH,
    FLUSHL2,
    TOUCHLOCKSET,
    TOUCHLOCKCLEAR,
    TOUCHL2,
    PREF0,  // prefetch into L1
    PREF1,  // prefetch into L2
    PREF0NT,  // pref. into L1 with non-temporal hint
    PREF1NT,  //    -"-     L2  -"-
    ncache
};

enum IALU_OPS{
    IADD,
    ISUB,
    IMUL,
    IOR, /*NB OR RD,RS,RS == MOV RD,RS*/
    CMOVGE,/*NOT PRESENT ON THE POWER PC - use branch with hint*/
    CMOVLE,/*NOT PRESENT ON THE POWER PC - use branch with hint*/
    CMOVGT,/*NOT PRESENT ON THE POWER PC - use branch with hint*/
    CMOVLT,/*NOT PRESENT ON THE POWER PC - use branch with hint*/
    ISUB_PREDICATE, /*Predicate recording subtraction PPC only*/
    /*Since we count down in loops and subtract in*/
    /*bounds checking, I think this is the only predicate*/
    /*Recording instruction needed....*/
    IOR_PREDICATE,  /*Predicate recording or PPC only*/
    IMUL_IMM,       /*Multiply by immediate value - needed to workout offsets*/
    /*from site no and atom immediate*/
    IADD_IMM,
    IAND_IMM,
    MOV,  // "OR rd,rs,rs" is not avail. on x86-based MIC
    INEG,
    RSHIFT,
    LSHIFT,
    AND,
    ANDC,
    IOR_IMM,
    nadds
};

/*Assembler Directives (like branch target, align, etc)*/
enum ASM_OPS{
    Target,
    CmovTarget,
    Cache_Align,
    Enter_Routine,
    Exit_Routine,
    Pipe_Flush,
    FNOP,
    NOP,
    PUSH_WINDOW,
    POP_WINDOW,
    LS_BARRIER,
    DEBUG_MSG_INT,
    DEBUG_MSG_FLOAT,
    DEBUG_MSG_COMPLEX,
    COMMENT,
    ndirective
};
enum SIMDHUMMER_OPS {
/*Needed for faster linalg (Dot, VAXPY, VAXMY etc...*/
    VMADD,
    VNMADD,
    VMSUB,
    VNMSUB,
    VMUL,
    VNMUL,
    VADD,
    VSUB,
    /*CMUL A+iB C+iD
    * MUL     RR ,RI,
    * MADD   -II , IR
    *
    *CMADD A+iB C+iD
    * MADD    RR ,RI
    * MADD   -II , IR
    *
    *CONJMUL A-iB C+iD
    *  MUL  RR,RI
    *  MADD II,-IR
    *
    *CONJMADD A-iB C+iD
     *  MADD  RR,RI
     *  MADD  II,-IR
     */
     VMUL_RR_RI,
     VMADD_RR_RI,
     VMADD_MII_IR,
     VMADD_II_MIR,
     VMUL_MII_IR,
     VNEG,
    VZERO,
     nsimdhummer
};
enum SIMDPERMUTE_OPS {
    VPERM,
    VPERMIMM,
    VPERMCOARSE,
    VPERMFINE,
    nsimdpermute
};
/*VMX/Altivec instructions*/
enum SIMDVMX_OPS {
    VMX_MADD,
    VMX_NMSUB,
    VMX_VADD,
    VMX_VSUB,
    nsimdvmx
};


/* Intel MIC instructions */
/*
 NOTE: This issuegroup MUST NOT contain MIC instructions which use mask registers in an "uncurly" way, e.g. KEXTRACT. If an instruction uses a mask register the mask reg must be enclosed by curly braces in the assembly code. If "uncurly" MIC instructions are required a new issuegroup must be introduced or dump.C::arg_print() etc. must be modified appropriately!

how to implement swizzles and masks?
define following abstract MIC instructions:
MIC_ADD,      narg=3, argtypes={vreg, vreg, vreg}  (NB: dest. is left most arg.)
MIC_ADD_MASK, narg=4, argtypes={vreg, mask, vreg, vreg}
MIC_ADD_SWIZ, narg=4, argtypes={vreg, vreg, vreg, swiz}
MIC_ADD_MASK_SWIZ, narg=5, argtypes={vreg, mask, vreg, vreg, swiz}
*/
enum SIMDMIC_OPS
{
    MIC_MOV,
    MIC_ADD,
    MIC_SUB,
    MIC_MUL,
    MIC_FMADD213,
    MIC_FMADD231,
    MIC_FNMADD213,
    MIC_FNMADD231,
    MIC_FMSUB213,
    MIC_FMSUB231,
    MIC_AND_MEMOP,
    MIC_ADD_MEMOP,
    MIC_SUB_MEMOP,
    MIC_MUL_MEMOP,
    MIC_FMADD231_MEMOP,
    MIC_XOR,
    nsimdmic
};
enum SIMDMICSWIZ_OPS
{
    // NOTE: first list unmasked ops, then list masked ops
    MIC_ADD_SWIZ = 0,
    MIC_MUL_SWIZ,
    MIC_SUB_SWIZ,
    MIC_FMADD231_SWIZ,
    MIC_ADD_MASK_SWIZ,  // NOTE: this must be the first masked op. Important as it is used as limit in create_knc_instructions()!
    MIC_SUB_MASK_SWIZ,
    MIC_SUBR_MASK_SWIZ,
    MIC_BLEND_MASK_SWIZ,
    MIC_VPERMF32X4,
    CLOAD_SPLAT,        /*Complex load (or SIMD load)*/    
    CSLOAD_SPLAT,        /*Complex load (or SIMD load)*/
    nsimdmicswiz
};


/**
 * Define all types of argument to mnemonics
 *Presently just Fregs, Iregs, Effective Addresses = SIMM(Ireg), and immediate val SIMM,Label
 *Condition register - an ireg on alpha, and a predicate on PPC
 */
enum ARGUMENTS
{
    Fregs,
    Iregs,
    SIMM,  // signed immediate
    Label,
    Pregs,  /*Predicate register - unused at present*/
    Cycles, /*Cycle count for scheduler pragmas*/
    Cregs,
    MICmask,  /* Intel MIC 16-bit mask registers */
    MICswiz,  /* Intel MIC swizzle patterns */
    UIMM,  // unsigned immediate
};


// pass strings to make_inst()? Then replace int by corresponding string in make_inst()
enum MIC_SWIZ  // NB: order MUST correspond to mic_swizzles!
{
    MIC_SWIZ_DCBA,
    MIC_SWIZ_CDAB,
    MIC_SWIZ_BADC,
    MIC_SWIZ_DACB,
    MIC_SWIZ_AAAA,
    MIC_SWIZ_BBBB,
    MIC_SWIZ_CCCC,
    MIC_SWIZ_DDDD,
    nmicswiz
};


extern const char *mic_swizzles[];  // see processor.C




enum ARGDIRS  // QUESTION: what about args which are both in&out?
{
    inarg    = 1,
    outarg   = 2,
    inoutarg = 3   // NB: (inarg|outarg)
};

/*Define a register file - should really have size as a parameter but */
/*Easy to sort later as required I think*/
struct regs
{
    int regfile;  // number of registers
    const char *mnemonic[MAX_REGS];
    int allocated[MAX_REGS];
    const char *fmt;  // NB: "format"
    int save_on_stack[MAX_REGS];  // regs alloc'd during a function execution get a flag set to 1
};

/*
* Each execution unit gets a pipeline - slightly
* clunky when issue width is narrower than the number of units.
* but then would need to encode issue rules.
* Should be alright just doing it this way since we want to keep the FPU
* busy at all times, so just pretend the combination of other all units
* constitutes (a) separate pipe(s)
*/
struct pipeline
{
    int issuegroups;  /* Bit mask flagging which kinds of instruction are allowed to issue */
};

struct instruction
{
    char *mnemonic;
    int narg;  ///< number of arguments
    int pipe;  ///< issue group (enum EXECUTIONPIPES), e.g., IALUPIPE or FMACPIPE
    int code;  ///< the instr's index within its issue group, e.g., IADD or FMUL
    int argtypes[MAX_INST_ARG];  ///< argument types (enum ARGUMENTS), e.g., Iregs or SIMM
    int inout[MAX_INST_ARG];     /* read or write - use enum names inarg and outarg*/
    int args[MAX_INST_ARG];      ///< The values of the arguments, i.e., index of iregs/fregs or offset handle
    struct instruction *staller;
    //int latency;  // QUESTION: is this member ever used?
};

struct instructionset
{
    struct instruction *instructions[N_PIPE_TYPES];
};


/*
* This is where we say how the different instructions can issue
* on a given processor.
*/
struct processor
{
    int npipes; ///< Number of issues allowed simultaneously
    int id;  ///< processor id (enum procids)
    int nsimd;///< SIMD vector width ... this represents the number of COMPLEX elements .. so for BGQ it is 2, for BGL/P it is 1
    struct pipeline *pipes;  ///< hardware pipelines
    struct regs *iregs;  ///< int register file
    struct regs *fregs;  ///< float/vector register file
    struct instructionset *iset;
    int FP_size;    ///< size of double precision floating point
    int FSP_size;   ///< size of single precision floating point
    int I_size;     ///< size of int
    int IS_size;    ///< size of short int
    int CacheLine;  ///< size of cache line
    int ABI_argbase;
    int StackPointer;  ///< iregs index of stack pointer register
    int ReturnPointer;
    long Isave;  ///< bit mask: which iregs are to be saved on stack
    long Fsave;  ///< bit mask: fregs to be saved on stack
    int switchargs[N_PIPE_TYPES];
    int latency[N_PIPE_TYPES];
    int usenops;
    int delayslot;
    int RegWindows;  // see stack.C: grab_stack()
    int Bias;  ///< offset between SP and actual stack address (SPARC v9)      // see stack.C: grab_stack()
    int lmq_space;
    int l1_locking;
    int load_reorder;
    int cross_negative;
    int nmicmask;  ///< number of Intel MIC mask registers
    int ImmediateOffsets;  ///< have immediate offsets? [1|0]
    int IREGS_SIZE;
    int FREGS_SIZE;  ///< Size of one Freg in bytes
    int MAX_ARGS_BY_REGS;  ///< Maximum number of function arguments passed via registers. The (MAX_ARGS_BY_REGS + 1)th (and any further) argument is passed on the stack.
    processor()
    {
        nsimd = 1; lmq_space = 0; l1_locking = 0;
        ImmediateOffsets=1;
        cross_negative = 0; load_reorder = 1; nmicmask = 0;
        IREGS_SIZE = 8; FREGS_SIZE = 8;
        MAX_ARGS_BY_REGS = 1000;  // NB: virtually all arguments are passed via registers by default
    };
};

/*
* Dependency types
*/
enum DEPENDENCIES
{
    NO_DEPENDENCY=0,
    READ_A_WRITE , /*True dependency - must wait instruction latency*/
    WRITE_A_READ , /*Can be fixed with a temporary reg*/
    BARRIER        /*raising messes up semantics */
              /*i .e. cant move instructions across a branch etc...*/
};


/*Linked list of instructions that gets reordered by the scheduler*/
struct inst_queue_entry{
    struct instruction *inst;
    struct inst_queue_entry *prev;
    struct inst_queue_entry *next;
};


/*Linked list to contain all registered immediate arguments*/
struct offset {
    int offst;
    int units;
    const char *name;
    struct offset *next;
};


/*
* Defines prefetch and write hint streams used by auto managing
* stream functions
*/

enum PATTERNS {
    LINEAR,
    STRIDED,
    LOOKUP,
    STREAM_IN,
    STREAM_OUT,
    LOOKUP_MULTI,
    POINTER
};

struct stream {
    int inout;   /*Direction*/
    int pointer; /*Register for cur pointer to */
    int basereg; /*Register for base pointer to matrix*/
    int preg;    /*Register for pointer to prefetch from*/
    int atom; /* handle for the immediate value that defines the bytes per site*/
    enum PATTERNS pattern; /*Name of the access pattern (previous enum)*/
    int stride;  /*handle for immediate value of stride (equals atom if LINEAR)*/
    int lines;
    int minline;
    int counter;
    int lookupreg; /*A pointer to a look up shift table*/
    int way;       /*Adding multi-way lookup tables*/
    int nway;
    int l1dist;  ///< Handle for L1 prefetch distance (in Bytes)
    int l2dist;  ///< Handle for L2 prefetch distance (in Bytes)
};


/*
* A software rotating register ... hope to also be able to
* map to IA64 rotating regs and Hitachi VPP reg windows
*/
struct rotating_reg {
    int kind;
    int *regs;
    int last;
    int num;
};

/*
*
* Define functions which are common to all processors
* this is really the user interface or proto-language definition
*
* Though it is stretching things a bit to call this a language
*
*/

/*Set the processor of choice*/
void set_processor(struct processor *);
void set_abi_caller_save_all(void);
void set_label_prefix(const char *pfx);

/* Register management routines*/
int  defargcount(int narg);
void getarg(int ireg);
int allocate_reg( int kind,const char *name); /*kind = Iregs or Fregs*/
int allocate_tmp_reg(enum ARGUMENTS kind);
void free_reg( int kind,int num);

/*
* Processor query routines
*/
int have_hummer(void);
int have_permute(void);
int have_mic();
int nsimd(void);


/**
 * Checks whether processor is based on x86. This auxiliary function is called in 3-operand scalar operations in order to decide whether we have to perform an additional register copy. This copy might be necessary on x86 as it does not support non-destructive 3-operand scalar operations.
 */
int is_x86();


// Some shortcuts for convenience...
void enter_routine(const char *const name, int nbyte);
void exit_routine();
void target(int label);


/*
* rotating register support
*/
struct rotating_reg *create_rotating_reg(int kind,int len,const char *name);
int get_rotating_register(struct rotating_reg *rp);
void reset_rotating_register(struct rotating_reg *rp);

/*Immediate argument management*/

int def_offset(int units, enum Datum type,const char * name);
int get_offset(int handle);
int get_offset_handle(int off, enum Datum size);
const char *get_offset_name( int handle );
struct offset *get_offset_struct(int handle);

/*Loop management routines*/
int get_target_label(void);

/*Generic instruction insertion*/
void make_inst(int PIPE, int OP, ...);

void arg_print(struct instruction *i, int j);

/*Stream optimisation set up routine*/
extern int N_Ahead; /*Streaming distance ahead - 3 times around outer loop*/
struct stream * create_stream(int atom,int base_reg,int countreg,
    enum PATTERNS dir, enum PATTERNS kind,...);
void iterate_stream(struct stream *);
void iterate_table(struct stream *); /*Introduced for share multi stream tables*/
void set_lookup_nxt_pointer (struct stream *mystream,int nxt);

/*Simplified and multi instruction insertion routines*/
void queue_prefetch  (struct stream *mystream);
void queue_write_hint(struct stream *mystream);
void queue_fload     (int dest,int Himm,int basereg,enum Datum type = Double);
void queue_fstore    (int src ,int Himm,int basereg,enum Datum type = Double);
void queue_iload      (int dest,int Himm,int basereg);
void queue_iload_short(int dest,int Himm,int basereg);
void queue_load_addr (int dest,int Himm,int basereg );
void queue_instruction(struct instruction *inst);
void queue_iadd( int dest, int a, int b );
void queue_isub( int dest, int a, int b );
void queue_iadd_imm( int dest, int a, int himm );
void queue_isub_imm( int dest, int a, int himm );
void queue_iand_imm( int dest, int a, int himm );
void queue_rshift(int dest, int src, int shifthandle);
void queue_lshift(int dst,int src,int himm);
void queue_ori   (int dst,int src,int huimm);
void queue_and   (int dst,int src1,int src2);
void queue_andc   (int dst,int src1,int src2);
void queue_imul( int dest, int a, int b );
void queue_imul_imm( int dest, int a, int himm );
void queue_fadd( int dest, int a, int b );
void queue_ficonv( int dest, int src );
void queue_fmul( int dest, int a, int b );
void queue_fsub( int dest, int a, int b );
void queue_fneg(int dest,int src);
void queue_fmov(int dest,int src);
void queue_fmac(int acct,int a,int b); /*d += a*b*/
void queue_fnmac(int acc,int a,int b); /*d -= a*b*/
void queue_fmadd(int dest,int a,int b,int c); /*d = c+ a*b */
void queue_fmsub(int d, int a, int b, int c);  /* d = a*b - c */
void queue_fnmsub(int dest,int a,int b,int c);/*d = c- a*b */
void queue_istore ( int reg, int offhandle, int basereg );

void pragma(int,int);
void debugC(int); /*Inserts a debug print of cmplx vec for noarch*/
void debugF(int); /*Inserts a debug print of float for noarch*/
void debugI(int); /*Inserts a debug print of int for noarch*/

/*My crude optimiser*/
void schedule_for_proc(void);
void optimise_copies();

/*Output routines*/
void set_human_readable(int h);
void dump_instruction_queue(void);
char *make_label(int i);
char *make_reg_name(int kind, int reg);
void dump_regs(int);
void inst_print(struct instruction *i);

/*These routines depend on whether or not you have an FMADD instruction*/
void setup_cmadds(void);
/*
 * If we have no FMADD pipe we need to allocate temporary
 *  registers to chain multiplies and adds
 */


void queue_bounds_check(int work, int dst);
void queue_barrier(void);
void queue_mbar(void);
void queue_pref_set_predicate(int counter);
void queue_iload_imm(int dest,int imm);


/*Complex support*/
void complex_load(int dest,int Himm,int basereg,enum Datum type = Double);
void complex_rsplat(int dest,int Himm,int basereg,enum Datum type = Double);
void complex_csplat(int dest,int Himm,int basereg,enum Datum type = Double);
void complex_store(int src ,int Himm,int basereg,enum Datum type = Double);
void complex_load_half(int dest,int Himm,int basereg,int membuf,int t1,int t2,int mask);
void complex_store_half(int src ,int Himm,int basereg,int membuf,int t1,int t2,int mask);
void complex_simd_init(int permreg);

/*Mult by I support: ea_reg should be pointer to (0.0,1.0)*/
void complex_constants_prepare(int creg, int i_reg);
enum { cmplx_i, cmplx_neg_i, cmplx_zero };
void complex_constants_assert( int constant);

/*Linear combination support*/
void complex_add(int e, int a, int c);
void complex_sub(int e, int a, int c);
void complex_ApiB(int e, int a, int c);/*e = a + i B*/
void complex_AmiB(int e, int a, int c);/*e = a - i B*/

void complex_real_madd(int e, int real_reg, int comp_reg,int addto);
void complex_real_madd(int e, int real_reg, int comp_reg);
void complex_real_mul(int e, int real_reg, int comp_reg);
// for consistency with complex_madd[to]()...
#define complex_real_maddto(e, r, c, g) complex_real_madd((e), (r), (c), (g))

void simd_mul(int e, int a, int b);
void simd_madd(int e, int a, int b,int g);
void simd_nmadd(int e, int a, int b,int g);
void simd_add(int e, int a, int b);
void simd_msub(int e, int a, int c, int g);
void simd_zero(int a);
void simd_reduce_add(int basereg, int offhandle, int v, enum Datum type);


/*Complex register variants*/


void complex_mul(int e, int a, int c);
void complex_madd(int e, int a, int c);
void complex_madd(int e, int a, int c,int d);
void complex_conjmul(int e, int a, int c);
void complex_conjmadd(int e, int a, int c);
void complex_maddto(int h,int e,int a,int c);
void complex_conjmaddto(int h, int e,int a,int c);

/*
* Provided as general primitives for
*
*/
void complex_simd_permute(int mu, int dest, int src);
void complex_simd_extract(int hi, int mu, int dest , int src1 , int src2);
void complex_simd_merge(int hi, int mu, int dest1, int dest2, int src );
void complex_simd_permute_internal(int dest, int src1, int src2,int perm);


/*Fload register variants --- backwards compat*/
void queue_cmul    (int e, int f, int a, int b, int c, int d);
void queue_cmadd   (int e, int f, int a, int b, int c, int d);
void queue_conjmul (int e, int f, int a, int b, int c, int d);
void queue_conjmadd(int e, int f, int a, int b, int c, int d);
void queue_cmaddto (int h,int g,int e,int f,int a,int b,int c,int d);
void queue_conjmaddto(int h,int g,int e,int f,int a,int b,int c,int d);


void complex_twospin_color_dot(int u1,int u2,int u3,
			       int in_h1,int in_h2,int in_h3,
			       int in_l1,int in_l2,int in_l3,
			       int out_h,int out_l);


void complex_six_cmuls(int e1,int a1,int c1,
    int e2,int a2,int c2,
    int e3,int a3,int c3,
    int e4,int a4,int c4,
    int e5,int a5,int c5,
    int e6,int a6,int c6);
void complex_six_cmadds(int e1,int a1,int c1,
    int e2,int a2,int c2,
    int e3,int a3,int c3,
    int e4,int a4,int c4,
    int e5,int a5,int c5,
    int e6,int a6,int c6);
void complex_six_cmaddtos(
		       int h1, int e1,int a1,int c1,
                       int h2, int e2,int a2,int c2,
                       int h3, int e3,int a3,int c3,
		       int h4, int e4,int a4,int c4,
                       int h5, int e5,int a5,int c5,
                       int h6, int e6,int a6,int c6);


void complex_three_cmuls(int e1,int a1,int c1,
    int e2,int a2,int c2,
    int e3,int a3,int c3);
void complex_three_cmadds(int e1,int a1,int c1,
    int e2,int a2,int c2,
    int e3,int a3,int c3);
void complex_three_cmaddtos(int h1,int e1,int a1,int c1,
    int h2,int e2,int a2,int c2,
    int h3,int e3,int a3,int c3);

void complex_three_conjmuls(int e1,int a1,int c1,
    int e2,int a2,int c2,
    int e3,int a3,int c3);
void complex_three_conjmadds(int e1,int a1,int c1,
    int e2,int a2,int c2,
    int e3,int a3,int c3);
void complex_three_conjmaddtos(int h1,int e1,int a1,int c1,
    int h2,int e2,int a2,int c2,
    int h3,int e3,int a3,int c3);
void complex_zero(int reg);


void queue_three_cmuls(int e1,int f1,int a1,int b1,int c1,int d1,
    int e2,int f2,int a2,int b2,int c2,int d2,
    int e3,int f3,int a3,int b3,int c3,int d3);


void queue_three_cmadds(int e1,int f1,int a1,int b1,int c1,int d1,
    int e2,int f2,int a2,int b2,int c2,int d2,
    int e3,int f3,int a3,int b3,int c3,int d3);

void queue_three_cmaddtos(int h1,int g1,int e1,int f1,int a1,int b1, int c1, int d1,
    int h2,int g2,int e2,int f2,int a2,int b2, int c2, int d2,
    int h3,int g3,int e3,int f3,int a3,int b3, int c3, int d3);


void queue_three_conjmuls(int e1,int f1,int a1,int b1,int c1,int d1,
    int e2,int f2,int a2,int b2,int c2,int d2,
    int e3,int f3,int a3,int b3,int c3,int d3);

void queue_three_conjmadds(int e1,int f1,int a1,int b1,int c1,int d1,
    int e2,int f2,int a2,int b2,int c2,int d2,
    int e3,int f3,int a3,int b3,int c3,int d3);

void queue_three_conjmaddtos(int h1,int g1,int e1,int f1,int a1,int b1, int c1, int d1,
    int h2,int g2,int e2,int f2,int a2,int b2, int c2, int d2,
    int h3,int g3,int e3,int f3,int a3,int b3, int c3, int d3);



/*
* Save and restore registers as defined by the ABI - so probably define in
* processor specific c code
*/

int grab_stack ( int );
void free_stack ( void );
void save_regs(void);
void restore_regs(void);

#ifdef alreg
#error
#endif
#ifdef def_off
#error
#endif

/*Macros for defining registers and immediate arguments*/
void *arralloc(size_t size, int ndim, ...);

#define def_off( A , T, B ) int A = def_offset ( B, T,(char *)#A )
#define def_dp_off( A , B ) int A = def_offset ( B, Double,(char *)#A )
#define def_sp_off( A , B ) int A = def_offset ( B, Single,(char *)#A )
#define def_byte_off( A , B ) int A = def_offset ( B, Byte,(char *)#A )
#define alreg(A,T)  int A = allocate_reg(T,(char *)#A)

#define reg_array_1d( C, T , d )        register_array_1d C((char *)#C,T,d)
#define reg_array_2d( C, T , d1,d2 )    register_array_2d C((char *)#C,T,d1,d2)
#define reg_array_3d( C, T , d1,d2,d3 ) register_array_3d C( (char *) #C,T,d1,d2,d3)

#define offset_1d( C, T, d )        offset_array_1d C((char *)#C,d,0,T)
#define offset_2d( C, T, d1,d2 )    offset_array_2d C((char *)#C,d1,d2,0,T)
#define offset_3d( C, T, d1,d2,d3 ) offset_array_3d C((char *)#C,d1,d2,d3,0,T)

#define offset_2d_dp( C,d1,d2) offset_2d(C,Double,d1,d2)

#define offset_1d_shft( C, T, d,s )        offset_array_1d C((char *)#C,d,s,T)
#define offset_2d_shft( C, T, d1,d2,s )    offset_array_2d C((char *)#C,d1,d2,s,T)
#define offset_3d_shft( C, T, d1,d2,d3,s ) offset_array_3d C((char *)#C,d1,d2,d3,s,T)

/*Define the available processors*/
enum procids {
    ALPHA264,
    ALPHA264_SINGLE,
    ALPHA164,
    ALPHA164_SINGLE,
    ALPHA064,
    ALPHA064_SINGLE,
    KNC,         /* Intel MIC Knights Corner */
    KNC_SINGLE,  /* Intel MIC Knights Corner, single precision */
    PPC440,
    PPC440_SINGLE,
    PPC440BGL,
    PPCBGQ,
    POWERIII,
    POWERIII_SINGLE,
    USPARCII,
    USPARCII_SINGLE,
    UNKNOWN,  // NOTE: UNKNOWN ids must be at the end. See stack.C:save_regs().
    UNKNOWN_SINGLE,
    UNKNOWN_BGQ,
    UNKNOWN_VMX,
    nprocids
};
extern const char *procargs[nprocids];


/*Processor dependent functions -> use set of function pointers*/
void set_scheduler_weird(void);
void set_scheduler_bgq(void );
struct processor * create_processor(int);
int start_loop(int counter);
void stop_loop(int branchno,int counter);
void directive_print(struct instruction *i);
void directive_print_unknown(struct instruction *i);
void queue_cmovge(int cond,int src,int dst);
void queue_cmovgt(int cond,int src,int dst);
void queue_cmovle(int cond,int src,int dst);
void queue_cmovlt(int cond,int src,int dst);
void queue_mov(int src,int dst);

void cache_print(struct instruction *i);
void check_iterations(int cntreg,int retno);
void conditional_branch_cmpzero( int branch_type, int condreg, int target  ) ;
void conditional_branch_cmpimm( int branch_type, int condreg, int himm, int target ) ;


void need_cache_line(int line);
void need_constant (int);
int  get_constant(int val);
int  get_constant_reg(int myreg);
void do_prefetch ( int ptr,int line );
void do_flush    ( int ptr,int line );
void prefetch ( int ptr,int line );
void flush    ( int ptr,int line );
void l1_lock ( int ptr,int line );
void l1_unlock( int ptr,int line );
void l2_touch ( int ptr,int line );


void set_processor_optarg(char *opt);
void print_eff_addr(struct instruction *i);
int get_I_size(void);
#include "alpha.h"
#include "alpha264.h"
#include "mic.h"
#include "powerpc.h"
#include "unknown.h"
#include "usparc.h"
#endif
