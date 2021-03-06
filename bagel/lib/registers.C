/*
 *
 *  Copyright Peter Boyle and Glasgow University, 2000.
 *
 *  It is provided under the GNU pubic License V2
 *  It is provided as is and is not guaranteed fit for any purpose.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "processor.h"

extern struct processor *PROC;
extern struct inst_queue_entry *HEAD;

#define REG_FREED 3

/* { ... REGISTER management routines */
int allocate_reg( int kind,const char *name)
{
  struct regs *set;
  int i;
  int nregs;
  test_proc();
  switch ( kind ){
    case Fregs:
      set = PROC->fregs;
      nregs = 1;
      break;
    case Iregs:
      set = PROC->iregs;
      nregs = 1;
      break;
    case Cregs:
      set = PROC->fregs;
      if ( have_hummer() || have_mic() ) {
    nregs = 1;
      }
      else nregs = 2*nsimd();
      break;
    default: puts("Bad Argument: Allocate reg"); exit(0);
  }

  for ( i =0 ; i < set->regfile; i++){
    if ( set->allocated[i]==0 ){
      if ( nregs==1 ) {
        set->mnemonic[i] = name;
        set->allocated[i] = 1;
        set->save_on_stack[i] = 1;
	//	if ( kind == Cregs )
	  //	  fprintf(stderr,"Using SIMD reg for vector : (p%d,s%d)\n",i,i);
      } else if ( nregs>1 ) {

	for (int v=0;v<nregs;v++){

	  char *reg = (char *)malloc(strlen(name)+6);
	  sprintf(reg,"%s_%d",name,v);
	  if ( set->allocated[i+v]==0 ){
	    set->mnemonic[i+v] = reg;
	    set->allocated[i+v] = 1;
	    set->save_on_stack[i+v] = 1;
	  } else {
	    puts("Could not alloc contiguous registers"); exit(0);
	  }
        }
      }
      //      fprintf(stderr, "Allocated reg %d of type %d.\n", i, kind);
      return(i);
    }
  }

  for ( i =0 ; i < set->regfile; i++){
    if ( set->allocated[i] == REG_FREED ){
      //      fprintf(stderr,"Warning reusing register %s\n",make_reg_name(kind,i));
      set->mnemonic[i] = name;
      set->allocated[i] = 1;
      set->save_on_stack[i] = 1;
      return(i);
      if ( nregs > 1 ) {
    puts("Error: Could not reuse register");
    exit(0);
      }
    }
  }
  puts("Error: Could not allocate register");
  puts(name);
  exit(0);
}


/**
 * \brief Frees a register.
 * \param[in] kind The register type: Iregs, Fregs, Cregs.
 * \param[in] num The register's index in its register file.
 */
void free_reg( int kind,int num)
{
  struct regs *set;
  int nregs;
  test_proc();
  switch ( kind ){
    case Fregs:
      set = PROC->fregs;
      nregs = 1;
      break;
    case Iregs:
      set = PROC->iregs;
      nregs = 1;
      break;
    case Cregs:
      set = PROC->fregs;
      if ( have_hummer() || have_mic() ) nregs = 1;
      else nregs = 2;
      break;
    default: puts("Bad Argument: Free reg: Unknown reg type."); exit(0);
  }

  if ( num > set -> regfile || num < 0 ) {
     puts("Bad Argument: Free reg"); exit(0);
  }

  if ( !set->allocated[num] ){
    puts("Error: freeing unallocated reg");
    printf("Kind = %d, num = %d",kind,num);
  }
  set->allocated[num] = REG_FREED; /*Tag for reuse but keep mnemonic around unless actually reused*/
  if ( nregs == 2 ) {
    set->allocated[num+1] = REG_FREED; /*Tag for reuse but keep mnemonic around unless actually reused*/
  }
  if ( nregs > 2 ) {
    exit(0);
  }
  return;
}

/* ...}*/


struct rotating_reg *create_rotating_reg(int kind,int len,const char *name)
{
  struct rotating_reg *rp = new rotating_reg;
  int i;
  char *rnam;
  int namlen;
  rp -> kind = kind;
  rp -> last = 0;
  rp -> num = len;
  rp -> regs = (int *)malloc(len * sizeof(int));
  namlen = strlen(name)+4; /*4 since may port to IA64 with 3 digit regs*/
  for ( i = 0; i<len; i++){
    rnam = (char *)malloc(namlen);
    sprintf(rnam,"%s%d",name,i);
    rp->regs[i] = allocate_reg(kind,rnam);
  }
  return(rp);
}

int get_rotating_register(struct rotating_reg *rp)
{
  int r = rp->regs[rp->last];

  rp->last ++;
  if ( rp->last >= rp->num ){
    rp->last = 0;
  }
  return(r);
}
void reset_rotating_register(struct rotating_reg *rp)
{
  rp->last = 0;
}


#include <assert.h>
#define STR_BUF_SIZE 64

/**
 * \brief Allocates a temporary register of \p kind (e.g., Fregs, Iregs).
 * \details
 * This function does the same as allocate_reg(), but automatically chooses a unique name for the temporary register. The idea is to make human-readable assembly unambiguous.
 * Registers are freed with free_reg().
 *
 * \see allocate_reg()
 * \see free_reg()
 */
int allocate_tmp_reg(enum ARGUMENTS kind)
{
    static int counter = 0;
    char *str;
    int reg, n;

    str = (char *)malloc(STR_BUF_SIZE);
    assert(str != NULL);
    n = snprintf(str, STR_BUF_SIZE, "tmp%d_%d", (int)kind, counter);
    assert(n > 0 && n < STR_BUF_SIZE);
    ++counter;
    reg = allocate_reg(kind, str);
    free(str);
    return reg;
}






