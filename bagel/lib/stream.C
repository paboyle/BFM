/*
 *
 *  Copyright Peter Boyle and Glasgow University 2000.
 *
 *  It is provided under the GNU pubic License V2
 *  It is provided as is and is not guaranteed fit for any purpose.
 *
 */

#include "processor.h"

extern struct processor *PROC;
void do_writehint ( int ptr,int line );

/* {...Cache Streaming related functions*/
#define MAX_PLINES 24
#define MAX_LOOKUP 4
#define MAX_MULTIWAY 5

int lookupno =0;
const char *gnames[MAX_LOOKUP]={"LOOKUP0","LOOKUP1","LOOKUP2","LOOKUP3"};

/*
 * Support for integer constant registers needed for indexed loads
 * and PPC cache touch in particular which doesn't support immediate opcode encoded constants.
 * Bastards. Also needed for BG/L HUMMER loads, guess they were short of opcode space.
 * Perhaps should start using a load&update form to reduce the number of constants required
 */
const char *pnames[MAX_PLINES]={"CONST0","CONST1","CONST2","CONST3","CONST4",
                          "CONST5","CONST6","CONST7","CONST8","CONST9",
                          "CONST10","CONST11","CONST12","CONST13","CONST14",
                          "CONST15","CONST16","CONST17","CONST18","CONST19",
                          "CONST20","CONST21","CONST22","CONST23"
                         };
const char *prnames[MAX_PLINES]={"CONSTR0","CONSTR1","CONSTR2","CONSTR3","CONSTR4",
                           "CONSTR5","CONSTR6","CONSTR7","CONSTR8","CONSTR9",
                           "CONSTR10","CONSTR11","CONSTR12","CONSTR13","CONSTR14",
                           "CONSTR15","CONSTR16","CONSTR17","CONSTR18","CONSTR19",
                           "CONSTR20","CONSTR21","CONSTR22","CONSTR23"
                         };
int pregs[MAX_PLINES]={-1,-1,-1,-1,-1,
                       -1,-1,-1,-1,-1,
                       -1,-1,-1,-1,-1,
                       -1,-1,-1,-1,-1,
                       -1,-1,-1,-1
                      };
int pimms[MAX_PLINES]={-1,-1,-1,-1,-1,
                       -1,-1,-1,-1,-1,
                       -1,-1,-1,-1,-1,
                       -1,-1,-1,-1,-1,
                       -1,-1,-1,-1
                      };
int Prefetch;
int cmovreg=-1;  ///< Iregs
int scratch=-1;  ///< Iregs
int ahead_imm=-1;
int shift_ahead_imm=-1;
int lookupmulti_imms[MAX_MULTIWAY+1]= {
  -1,-1,-1,-1,-1,-1
};

const char *lookupmulti_immstr[MAX_MULTIWAY+1]= {
  "zeromulti","onemulti","twomulti","threemulti","fourmulti","fivemulti"
};

int Bdy = -1; /*Immediate number for iterations check        */
              /*Switch to prefetching from stack when within */
              /*Bdy sites of inner loop                      */

int Isize = -1;
int minusone = -1;

int  need_pref_regs ( void );

void set_lookup_pointer (struct stream *mystream);
int stream_reg_alloc=1;

void streams_reg_alloc(int val)
{
  stream_reg_alloc=val;
}

void initialise_streams(void)
{
  int i;

  test_proc();

  Prefetch   = def_offset( PROC->CacheLine,Byte,"CACHELINE");
  Isize      = def_offset(PROC->I_size,Byte,"Isize");
  Bdy        = def_offset(-(N_Ahead+2),Byte,"Bdy");
  ahead_imm  = def_offset(N_Ahead,Byte,"Nahead");
  shift_ahead_imm  = def_offset(N_Ahead*PROC->I_size,Byte,"SHFT_PREF_IMM");

  for(i=0;i<=MAX_MULTIWAY;i++){
   lookupmulti_imms[i]=
        def_offset(i*PROC->I_size,Byte,lookupmulti_immstr[i]);
  }
  if (stream_reg_alloc){
    cmovreg = allocate_reg(Iregs,"predicate");
    scratch = allocate_reg(Iregs,"prefscratch");
  }
  minusone = def_offset(-1,Byte,"minusone");
}


/**
 * For \p kind == \c LINEAR expected variable arguments are HANDLES for L1 prefetch distance and L2 prefetch distance in units of Bytes.
 */
struct stream * create_stream(int atom,int base_reg,int countreg,
                               enum PATTERNS dir, enum PATTERNS kind,...)
{
  struct stream *ret = new stream;
  int size = get_offset(atom);
  int i,line;
  va_list ap;

  va_start ( ap,kind);

  ret->inout = dir;
  ret->pointer = base_reg;
  ret->basereg = base_reg;
  ret->atom = atom;
  ret->counter = countreg;
  ret->pattern = kind;
        /*no of cache Lines covering one atom*/
  ret->lines = 1  + ((size-1 )/(PROC->CacheLine));  // XXX: This assumes data to be aligned at CL!!

  if ( dir != STREAM_IN && dir != STREAM_OUT ) {
    puts("Create stream: direction not specified");
    exit(1);
  }

  switch ( kind ){
  case LINEAR:
    ret->stride = ret->atom;
    ret->lookupreg = -1;
    ret->minline = (N_Ahead*size)/PROC->CacheLine;
    /*Check we have grabbed the offsets required   */
    /*On powerpc these have to be stored in Iregs  */
    /*since Iops are at a premium and the dcbt inst*/
    /*doesn't have and immediate offset            */
    for ( i = 0 ; i < ret->lines ; i++ ){
      line =  i + (N_Ahead*size)/PROC->CacheLine ;
      need_cache_line(line);
    }

    // for Bernhard's prefetching...
    ret->l1dist = va_arg(ap, int);
    ret->l2dist = va_arg(ap, int);
    break;
  case STRIDED:
    ret->stride = va_arg(ap,int); /*Stride should be an offset handle, */
                                  /*not the number itself*/
    ret->lookupreg = -1;
    ret->minline = 0;
    for ( i = 0 ; i < ret->lines ; i++ ){
      need_cache_line(i);
    }

    // for Bernhard's prefetching...
    ret->l1dist = va_arg(ap, int);
    ret->l2dist = va_arg(ap, int);
    break;
  case POINTER:
    ret->basereg = -1;
    ret->stride = -1;
    ret->lookupreg = va_arg(ap,int); /*This is the table pointer*/
    for ( i = 0 ; i < ret->lines ; i++ ){
      need_cache_line(i);
    }
    set_lookup_pointer(ret);
    break;
  case LOOKUP:
    /*We need to keep a copy of the start of the array*/
    if ( lookupno >= MAX_LOOKUP ) {
       puts("Error out of indirection regs");
       exit(0);
    }
    ret->basereg = allocate_reg(Iregs,gnames[lookupno++]);
    queue_mov(base_reg, ret->basereg);
    ret->stride = -1;
    ret->lookupreg = va_arg(ap,int);
    ret->minline = 0;
    for ( i = 0 ; i < ret->lines ; i++ ){
      need_cache_line(i);
    }
    need_cache_line(1);
    set_lookup_pointer(ret);
    break;

  case LOOKUP_MULTI: /*A multi way look up table*/

    /*We need to keep a copy of the start of the array*/
    if ( lookupno >= MAX_LOOKUP ) {
       puts("Error out of indirection regs");
       exit(0);
    }
    ret->basereg = allocate_reg(Iregs,gnames[lookupno++]);
    queue_mov(base_reg, ret->basereg);
    ret->stride = -1;
    ret->way  = va_arg(ap,int); /*Lookup field way in */
    ret->lookupreg = va_arg(ap,int);
    ret->nway = va_arg(ap,int); /*a table with nway ints per site*/
    ret->minline = 0;
    for ( i = 0 ; i < ret->lines ; i++ ){
      need_cache_line(i);
    }
    need_cache_line(1);
    set_lookup_pointer(ret);
    break;

  default:
    puts("create_stream: Unkown access pattern");
    exit(1);
  }
  va_end(ap);
  return(ret);
}

void need_cache_line(int line)
{
  need_constant(PROC->CacheLine*line);
}
int  get_constant_reg(int myreg)
{
  for (int reg = 0; reg < MAX_PLINES; reg++) {
    if ( pregs[reg] == myreg ) {
      return get_offset(pimms[reg]);
    }
  }
  printf("Oops couldn't locate constant register %d\n",myreg);
  exit(-1);
}

int  get_constant(int val)
{
  for (int reg = 0; reg < MAX_PLINES; reg++) {
    if ( pimms[reg] >= 0 ) {
      if ( get_offset(pimms[reg]) == val ) {
        return reg;
      }
    }
  }
  return -1;
}
void need_constant(int num)
{

  if (need_pref_regs () ) {

    /*Find free line*/
    if ( get_constant(num) != -1 ) {
      return;
    }
    for (int reg = 0; reg < MAX_PLINES; reg++) {
      if ( pimms[reg] == -1 ) {
    pimms[reg] = def_offset(num,Byte,pnames[reg]);
    /*Power doesn't provide immediate offsets for prefetches & Hummer loads  */
    /*So stick the value in a register to save recalculation of*/
    /*offset value                                             */
    if ( pregs[reg]== -1 ){
      pregs[reg] = allocate_reg(Iregs,prnames[reg]);
      make_inst(LOADPIPE,LOAD_IMM,pregs[reg],pimms[reg]);
    }
    return;
      }
    }
    printf("Couldn't allocate constant %d\n",num);
    exit(-1);
  }
}

int need_pref_regs ( void )
{
  if ( PROC->ImmediateOffsets == 1 ) return 0;
  else return 1;
}

void internal_prefetch ( int op, int ptr,int cline );

void l1_lock ( int ptr,int line )
{
  if (PROC->l1_locking) internal_prefetch(TOUCHLOCKSET,ptr,line);
  else internal_prefetch(PREF_IMM,ptr,line);
}
void l2_touch ( int ptr,int line )
{
  if (PROC->l1_locking) internal_prefetch(TOUCHL2,ptr,line);
  else internal_prefetch(PREF_IMM,ptr,line);
}
void l1_unlock( int ptr,int line )
{
  if (PROC->l1_locking) internal_prefetch(TOUCHLOCKCLEAR,ptr,line);
}
void prefetch ( int ptr,int cline )
{
  do_prefetch(ptr,cline);
}
void do_prefetch ( int ptr,int cline )
{
  internal_prefetch(PREF_IMM,ptr,cline);
}

void internal_prefetch ( int op, int ptr,int cline )
{
//  printf("op=%d,ptr=%d,cline=%d\n", op, ptr, cline);  // DEBUGGING
  if ( need_pref_regs() ){
    int constant = PROC->CacheLine * cline;
    int reg = get_constant(constant);
    if (pregs[reg]<0) {printf("Error: unallocated reg constant %d\n",constant); exit(-1);}
    make_inst(CACHPIPE,op,pregs[reg],ptr);
  }
  else if (PROC->id == KNC || PROC->id == KNC_SINGLE) {
    fprintf(stderr, "Error: Prefetching must be ported to KNC!\n");
    //    abort();
  }
  else {
    make_inst(CACHPIPE,op,cline,ptr);
  }
}
void flush( int ptr, int cline ) {
  do_flush(ptr,cline);
}
void do_flush( int ptr, int cline )
{
  int constant = PROC->CacheLine * cline;
  int reg = get_constant(constant);
  if ( need_pref_regs() ){
    if (pregs[reg]<0) {puts("Error: unallocated pref line"); exit(-1);}
    make_inst(CACHPIPE,FLUSH,pregs[reg],ptr);
  }
  else {
    /* Not supported on anything other than PPC yet */
    /*       make_inst(CACHPIPE,PREF_IMM,pimms[reg],ptr); */
    return;
  }
}
void do_writehint ( int ptr,int cline )
{
  /*
   *
   */
  int constant = PROC->CacheLine * cline;
  int reg = get_constant(constant);
  if ( need_pref_regs() ){
    if (pregs[reg]<0) {puts("Error: unallocated pref line"); exit(-1);}
    make_inst(CACHPIPE,WRITE_HINT,pregs[reg],ptr);
    if ( pregs[reg] == 0 ) {
      puts("do_writehint: Error reg zero is weird");
      exit(0);
    }
  } else {
    queue_load_addr(scratch,pimms[reg],ptr);
    make_inst(CACHPIPE,WRITE_HINT,scratch);
  }
}

void queue_prefetch  (struct stream *mystream)
{
  int i;
  int line;

  switch ( mystream->pattern ){
  case LINEAR:
    if (PROC->id == KNC || PROC->id == KNC_SINGLE) {
        // use non-temporal prefetches for streams
        prefetch_data_nt(L2L1, get_offset(mystream->atom), mystream->basereg, get_offset(mystream->l2dist), get_offset(mystream->l1dist));
    }
    else {
      for(i=0;i<mystream->lines;i++){
        line = mystream->minline + i;
        do_prefetch(mystream->basereg,line);
      }
    }
    break;
  case STRIDED:
    if (is_x86()) {
      // Why not multiply immediates at "compile-time"? We can do this on x86.
      int off = get_offset(mystream->stride) * get_offset(ahead_imm);
      int off_handle = def_offset(off, Byte, "strided_offset");
      queue_iadd_imm(scratch, mystream->basereg, off_handle);
      prefetch_data_nt(L2, get_offset(mystream->atom), scratch, get_offset_handle(0, Byte));  // XXX:
    } else {
      queue_iload_imm(scratch,mystream->stride);
      queue_imul_imm(scratch,scratch,ahead_imm);
      queue_iadd(scratch,scratch,mystream->basereg);
      for(i=0;i<mystream->lines;i++){
        do_prefetch(scratch,i);
      }
    }
    break;
  case LOOKUP:
    /*Load 3 ahead in shift table => scratch */
    queue_iload(scratch,shift_ahead_imm,mystream->lookupreg);
    /*  times atom to get byte offset into array*/
    queue_imul_imm(scratch,scratch,mystream->atom);
    /* add to base register */
    queue_iadd(scratch,mystream->basereg,scratch);
    /* Cmov to stack pointer if off limits in array*/
    queue_bounds_check ( mystream->counter,scratch);

    for(i=0;i<mystream->lines;i++){
      do_prefetch(scratch,i);
    }
    do_prefetch(mystream->lookupreg,1);
    break;
  case LOOKUP_MULTI:
    set_lookup_nxt_pointer(mystream,scratch);
    for(i=0;i<mystream->lines;i++){
      do_prefetch(scratch,i);
    }
  }
}

void set_lookup_pointer (struct stream *mystream)
{
  int ZeroImm = get_offset_handle(0,Byte);
  if ( mystream -> pattern == LOOKUP ){

    queue_iload(scratch,ZeroImm,mystream->lookupreg);
    queue_imul_imm(scratch,scratch,mystream->atom);
    /* add to base register and store in pointer*/
    queue_iadd(mystream->pointer, mystream->basereg, scratch);
  }  else if ( mystream-> pattern == LOOKUP_MULTI ){

   /*Load our way in shift table => scratch */
   queue_iload(scratch,lookupmulti_imms[mystream->way],mystream->lookupreg);

  }  else if ( mystream-> pattern == POINTER ){

    queue_iload(mystream->pointer,ZeroImm,mystream->lookupreg);

  } else {

    puts("Error: set_lookup_basereg on non-lookup stream");
    exit(0);

  }
}



void set_lookup_nxt_pointer (struct stream *mystream,int nxt)
{
  if ( mystream -> pattern == LOOKUP ){

    queue_iload(scratch,Isize,mystream->lookupreg);
    queue_imul_imm(scratch,scratch,mystream->atom);
    queue_iadd(nxt, mystream->basereg, scratch);
    queue_load_addr(scratch,minusone,mystream->counter);
    queue_cmovle(scratch,PROC->StackPointer,nxt);

  } else if ( mystream -> pattern == LOOKUP_MULTI ){

    queue_load_addr(scratch,lookupmulti_imms[mystream->nway],mystream->lookupreg);
    queue_iload(nxt,lookupmulti_imms[mystream->way],scratch);
    queue_load_addr(scratch,minusone,mystream->counter);
    queue_cmovle(scratch,PROC->StackPointer,nxt);

  } else {
    puts("Error: set_lookup_basereg on non-lookup stream");
    exit(0);
  }

}

void iterate_table(struct stream * mystream)
{
  if ( mystream -> pattern == LOOKUP_MULTI ){
    queue_load_addr(mystream->lookupreg,
        lookupmulti_imms[mystream->nway],mystream->lookupreg);
    do_prefetch(mystream->lookupreg,1);
  } else {
    puts("iterate_table: pooh"); exit(0);
  }
}

void iterate_stream(struct stream * mystream)
{

  if ( mystream -> pattern == LOOKUP ){
    queue_load_addr(mystream->lookupreg,Isize,mystream->lookupreg);
    set_lookup_pointer(mystream);
  } else if ( mystream -> pattern == POINTER ) {
    queue_load_addr(mystream->lookupreg,Isize,mystream->lookupreg);
    set_lookup_pointer(mystream);
  } else if ( mystream -> pattern == LOOKUP_MULTI ) {
    set_lookup_pointer(mystream);
  } else if ( mystream -> pattern == STRIDED ){
    queue_load_addr(mystream->basereg,mystream->stride,mystream->basereg);
  } else if (mystream->pattern == LINEAR) {
    queue_load_addr(mystream->basereg,mystream->atom,mystream->basereg);
  }
  else {
    fprintf(stderr, "Error: Invalid stream pattern.\n");
    abort();
  }
}

void queue_pref_set_predicate( int counter )
{
  queue_load_addr(cmovreg,Bdy,counter);
  if ( PROC->iset->instructions[IALUPIPE][IOR_PREDICATE].mnemonic !=NULL ){
    if (is_x86()) {
        make_inst(IALUPIPE, IOR_PREDICATE, cmovreg, cmovreg);
    }
    else {
        make_inst(IALUPIPE,IOR_PREDICATE,cmovreg,cmovreg,cmovreg);
    }
  }
}

void queue_bounds_check ( int counter, int reg )
{
  (void)counter;  // suppresses "unused parameter" compiler warning
  queue_cmovlt(cmovreg,PROC->StackPointer,reg);
}

void queue_write_hint(struct stream *mystream)
{
  int i;

  if ( mystream ->pattern == LINEAR ){
    int reg = mystream->minline*PROC->CacheLine;
    queue_load_addr(scratch,pimms[reg],mystream->basereg);
    queue_bounds_check(mystream->counter,scratch);

    for(i=0;i<mystream->lines;i++){
      make_inst(CACHPIPE,WRITE_HINT,scratch);
      if ( i+1 <mystream->lines )
        queue_load_addr(scratch,Prefetch,scratch);
    }
  }
  else {
    queue_prefetch(mystream);
  }
}
