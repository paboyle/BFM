

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
#include "processor.h"

#ifndef MAX
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))
#endif

extern struct processor * PROC;

int dependency ( struct instruction * i1, struct instruction *i2 );
int issue_early ( int min_dep,int max_dep,
                  struct inst_queue_entry *thisent,
                  struct instruction *execution_map[][MAX_WIDTH],
                  struct processor *proc);
void  execution_report(struct instruction * execution_map[][MAX_WIDTH]);
int deps_could_issue ( int min_dep, int max_dep,struct instruction *inst,
                  struct instruction *execution_map[][MAX_WIDTH],
                  struct processor *proc );
struct inst_queue_entry *find_next_dependent(struct inst_queue_entry *thisent);
struct inst_queue_entry *find_prev_dependent(struct inst_queue_entry *thisent);
struct inst_queue_entry *HEAD = NULL;
struct inst_queue_entry *HEAD1 = NULL;

void report_section( FILE *fptr, int start, int stop,
 int Fadds,int Fmuls, int Fmadds,
 int loads, int stores, int Iops,
 int nops, int fnops, int Islots, int Fslots
		     );
int WeirdShit = 0;
int BGQloads  = 0;
int dcbt_space = 3;   /*Space them 4 apart*/
int dcbt_pre = 0;     /*Assume dcbt pipelines perfectly by default*/
int dcbt_post = 0;
int store_lim = 1000; /*In practice infinite defaults*/
int load_lim = 1000;
int ls_lim = 1000;
int store_space = 0;
int load_space = 0;
int store_noreorder = 0;

void pragma_action(struct instruction *);
void pragma_action(struct instruction *pragma)
{
  switch(pragma->code ) {
  case DCBT_SPACE :
    dcbt_space = pragma->args[0];
    break;
  case DCBT_PRE :
    dcbt_pre = pragma->args[0];
    break;
  case DCBT_POST :
    dcbt_post = pragma->args[0];
    break;
  case STORE_LIM :
    store_lim = pragma->args[0];
    break;
  case STORE_SPACE :
    store_space = pragma->args[0];
    break;
  case STORE_INORDER :
    store_noreorder = pragma->args[0];
    break;
  case LOAD_LIM :
    load_lim = pragma->args[0];
    break;
  case LS_LIM :
    ls_lim = pragma->args[0];
    break;
  case LOAD_SPACE :
    load_space = pragma->args[0];
    fprintf (stderr,"Load space %d\n",load_space);fflush(stderr);
    break;
  default:
    printf ("bad pragma %d\n",pragma->code);
    exit(0);
    break;
  }
}

void set_scheduler_weird(void )
{
  WeirdShit = 1;
}

void set_scheduler_bgq(void )
{
  BGQloads = 1;
}

/* {...OPTIMISATION ROUTINES - NOT REALLY SMART, BUT EFFECTIVE*/
/*
 *
 * Take an instruction stream and create and execution record for our modelled
 * processor - this may not be completely accurate, but then we just need good enough
 *
 */


/* XXX: experimental {..... */

/*
 * After queuing all instructions and before re-scheduling, a new optimisation step could be inserted.
 * There are quite a lot of MIC vector MOVs which are not really necessary but costly. We could check within each function if a MIC_MOV is really necessary. If not the MIC_MOV is removed. Registers need to be renamed........
 *
 * NOTE: This function is HIGHLY EXPERIMENTAL and not complete yet!!!
 * TODO: WHAT ABOUT LOOPS?!??
 */
void optimise_copies()
{
	int i, src, dest, can_rewrite, line, innerline;
	struct inst_queue_entry *walk, *innerwalk;
	struct instruction *inst;

	if (!have_mic()) {  
		fprintf(stderr, "Note: Copy-Optimiser is only available for MIC at the moment.\n");
		return;
	}

	fprintf(stderr, "!!!\n!!! WARNING: Using a highly experimental optimiser!\n!!!\n");

	walk = HEAD;
	HEAD = NULL;

	line = 0;
	while (walk != NULL) {  // walk through instruction queue
		++line;  // for book-keeping

		if (walk->inst->pipe == SIMDMICPIPE && walk->inst->code == MIC_MOV) {  // search for MIC_MOV
			src  = walk->inst->args[1];
			dest = walk->inst->args[0];
			fprintf(stderr, "Info: Optimiser: MIC vector reg-reg copy found: %2d->%2d.\n", src, dest);

			// Now look for re-use of source reg

			can_rewrite = 1;
			innerwalk = walk->next;
			while (can_rewrite && innerwalk != NULL &&
				!(innerwalk->inst->pipe == DIRECTIVE && innerwalk->inst->code == Exit_Routine))  // don't cross function boundaries
			{
				inst = innerwalk->inst;

				if (inst->pipe == DIRECTIVE && inst->code == Enter_Routine) {
					fprintf(stderr, "Error: Optimiser: Enter_Routine found; something's badly screwed up!\n");
					exit(-1);
				}

				// "continue" for pipelines that do NOT operate on MIC vectors... We can skip our "expensive" analysis for these pipes.
				if (inst->pipe == IALUPIPE || inst->pipe == BRCHPIPE
					|| inst->pipe == DIRECTIVE || inst->pipe == WEIRDRULES
					|| inst->pipe == PRAGMA)
				{
					innerwalk = innerwalk->next;
					continue;
				}

				// analyse for potential re-use of "src"
				for (i = 0; i < inst->narg; i++) {
					if (inst->args[i] == src && inst->argtypes[i] == Fregs) {
						if (inst->inout[i] == inarg) {
							fprintf(stderr, "Info: Optimiser: Source reg %2d is re-used. pipe=%d,op=%d.\n", src, inst->pipe, inst->code);
							can_rewrite = 0;
							break;
						}
						else {
							if (!(inst->pipe == FLODPIPE && inst->code == FLOAD_REG_EA)) {  // can_rewrite = 1, if restoring from stack. But in that case we are not allowed to replace reg in restore!
								fprintf(stderr, "Info: Optimiser: reg %2d overwritten.\n", src);
								can_rewrite = 0;
							}
							// TODO........
							// FIXME: what about FMADD where out is also in?!
							// At the moment we optimise only if "src" is not re-used at all (or at stack-restore at best).
							// What about loops?
						}
					}
				}

				innerwalk = innerwalk->next;
			}

			if (can_rewrite) {
				// We found an unnecessary MIC_MOV. Now fix it!

				// (1) rename regs in subsequent instrs
				innerwalk = walk->next;
				innerline = 0;
				while (innerwalk != NULL &&
					!(innerwalk->inst->pipe == DIRECTIVE && innerwalk->inst->code == Exit_Routine))  // don't cross function boundaries
				{
					++innerline;
					inst = innerwalk->inst;

					if (inst->pipe == FLODPIPE && inst->code == FLOAD_REG_EA) {
						// don't rename if in stack restoration ...
						innerwalk = innerwalk->next;
						continue;
					}

					for (i = 0; i < inst->narg; i++) {
						if (inst->args[i] == dest && inst->argtypes[i] == Fregs) {
							inst->args[i] = src;  // XXX: is it ok to replace EVERYWHERE??
							fprintf(stderr, "Info: Optimiser: Replaced %2d by %2d in instr #%d. pipe=%d, op=%d.\n",
								dest, src, line+innerline, inst->pipe, inst->code);
						}
					}
					innerwalk = innerwalk->next;
				}

				// (2) remove the MIC_MOV instruction from queue
				struct inst_queue_entry *del = walk;
				walk = walk->next;  // skip MIC_MOV instr.
				delete del;  // NB: to avoid mem leak
			}
		}

		queue_instruction(walk->inst);
		walk = walk->next;
	}
}

/* XXX: experimental ....}*/



void schedule_for_proc(void)
{
/*Maximum 8 way super scalar*/
  struct inst_queue_entry *entry;
  struct instruction * execution_map[MAX_CYCLES][MAX_WIDTH];
  struct instruction *inst ;

  int max_dep; /*This is the last cycle previously issued to*/
  int min_dep; /*This is the minimum cycle from to which we should greedy issue*/
               /*Also the earliest cycle from which we consider dependencies*/
               /*... So max_dep - min_dep is "sortakinda" the peephole size in compiler speak*/
  int cycle;
  int pipe;

  if ( PROC->id >= UNKNOWN ){
    return;
  }

  /*Fill record with nops*/
  for(cycle = 0; cycle<MAX_CYCLES;cycle ++ ) {
    for(pipe = 0; pipe<MAX_WIDTH;pipe++){
      execution_map[cycle][pipe] = (struct instruction *)NULL;
    }
  }


  /*Loop through instructions trying to issue as early as allowed*/

  entry = HEAD;
  cycle = 0;
  max_dep=0;
  min_dep=0;

  while (entry != NULL ) {

    if ( entry->inst->pipe == PRAGMA ) { /*Set up the scheduler and discard the insn*/
      fprintf(stderr,"Pragma %d %s\n",entry->inst->code,entry->inst->mnemonic);
      pragma_action(entry->inst);

    } else {

      cycle = issue_early ( min_dep,max_dep,entry,execution_map,PROC );
      if (cycle > max_dep ) max_dep = cycle ;

    }

    entry = entry -> next;
  }

  /*Could write a file with execution map here*/
  execution_report(execution_map);


  HEAD1 = HEAD;
  HEAD = NULL;
  /*Now generate a instruction stream with the reordered execution map here*/
  for(cycle = 0; cycle<MAX_CYCLES;cycle ++ ) {
    for(pipe = 0; pipe<PROC->npipes;pipe++){

      inst = execution_map[cycle][pipe];

      if ( inst  != (struct instruction *)NULL ){
        queue_instruction(execution_map[cycle][pipe]);
        if ( inst->pipe==DIRECTIVE && inst->code == Target
                                   && inst->argtypes[0] == Label ){
          for(pipe++;pipe<PROC->npipes;pipe++)
            if ( execution_map[cycle][pipe] != NULL ) {
              printf("Error issuing align"); exit(0);
	    }
        }
      }
    }
  }
}

/*
 * Delay issuing instructions whose next dependent inst will stall
 * Crude attempt at forcing critical paths to be interleaved
 */

struct inst_queue_entry *find_next_dependent(struct inst_queue_entry *thisent)
{
  struct inst_queue_entry *dep;

  dep = thisent;
  do{
    dep = dep ->next;
  }  while ( dep != NULL &&
             (dependency(dep->inst,thisent->inst) == NO_DEPENDENCY) );
  return(dep);
}

struct inst_queue_entry *find_prev_dependent(struct inst_queue_entry *thisent)
{
  struct inst_queue_entry *dep;

  dep = thisent;
  do {
    dep = dep ->prev;
  }  while ( dep != NULL &&
	     (dependency(thisent->inst,dep->inst) == NO_DEPENDENCY ) );
  return(dep);
}

int pref_space=3;
int pref_base=3;
int pref_eo=0;
int fstore_space = 1 ;
int space_stores = 0;
void set_space_stores(void )
{
  space_stores = 1;
  fstore_space = 2;
}

int deps_could_issue ( int min_dep, int max_dep,struct instruction *inst,
                  struct instruction *execution_map[][MAX_WIDTH],
                  struct processor *proc)
{
  int cycle;
  int issue_cycle;
  int pipe;
  struct instruction *depop;

  inst->staller = NULL;
  issue_cycle = min_dep;
  for (cycle = max_dep; cycle >=min_dep; cycle -- ){
    for (pipe = proc->npipes -1; pipe >=0 ; pipe -- ){

      if (execution_map[cycle][pipe]!= NULL ) {
	depop = execution_map[cycle][pipe];
        switch(dependency(inst,execution_map[cycle][pipe])){

        case READ_A_WRITE:
          if ( cycle + proc->latency[depop->pipe] > issue_cycle ){
            issue_cycle = cycle + proc->latency[depop->pipe] ;
	    inst->staller = execution_map[cycle][pipe];
	  }
          break;

        case BARRIER:
          if ( cycle + 1 > issue_cycle){
            issue_cycle = cycle + 1 ;
	    inst->staller = execution_map[cycle][pipe];
          }
          break;
        case WRITE_A_READ:
          if ( WeirdShit ) {
           if ( cycle + 2 > issue_cycle){
            issue_cycle = cycle + 2 ;
	    inst->staller = execution_map[cycle][pipe];
           }
	  } else {
           if ( cycle + 1 > issue_cycle){
	    inst->staller = execution_map[cycle][pipe];
            issue_cycle = cycle + 1 ;
           }
	  }
          break;
        case NO_DEPENDENCY:

          break;
	}
      }
    }
  }

  return(issue_cycle);
}


void  execution_report(struct instruction * execution_map[][MAX_WIDTH])
{
  int cycle;
  int start=0;
  int stop=0;
  int pipe;

  int Fadds=0;
  int Fmuls=0;
  int Fmadds=0;
  int loads=0;
  int stores=0;
  int Fsubs=0;
  int Fslots=0;
  int Islots=0;
  int Iops=0;
  int nops=0;
  int fnops=0;
  int lab=-1;
  int done = 0;
  struct instruction *inst;

  FILE *fptr = fopen("execution.map","w");

  for (pipe = 0;pipe<PROC->npipes; pipe ++){
    if (BRCHMASK & PROC->pipes[pipe].issuegroups)      fprintf(fptr,"Ibox \t");
    else if (FMACMASK & PROC->pipes[pipe].issuegroups) fprintf(fptr,"FMAC \t");
    else if (FMULMASK & PROC->pipes[pipe].issuegroups) fprintf(fptr,"FMUL \t");
    else if (FADDMASK & PROC->pipes[pipe].issuegroups) fprintf(fptr,"FADD \t");
  }

  fprintf(fptr,"\n");
  fprintf(fptr,"------------------------------------------------\n");
  for ( cycle= 0; done==0; cycle ++){
    for ( pipe = 0; pipe<PROC->npipes && done ==0; pipe ++){

      if ( (inst = execution_map[cycle][pipe]) ){

       if ( inst->pipe==DIRECTIVE && inst->code == Target
                                  && inst->argtypes[0]== Label
          ){

          lab = inst->args[0];
          stop = cycle;
          report_section( fptr,start,stop,
                          Fadds,Fmuls, Fmadds,
			    loads, stores, Iops,
                          nops, fnops, Islots, Fslots );
          start = cycle;
          Fadds=0;          Fmuls=0;          Fmadds=0;         loads=0;
          stores=0;         Fsubs=0;          Fslots=0;         Islots=0;
          Iops=0;           nops=0;           fnops=0;
fprintf(fptr,"%s ------------------------------------------------\n",
	make_label(lab));
       }

       if (inst->pipe==BRCHPIPE && (inst->code==BRANCH_GT || inst->code==BRANCH_CTR))
	 {
          lab = inst->args[1];
          stop = cycle;
          report_section( fptr,start,stop,
                          Fadds,Fmuls, Fmadds,
			    loads, stores, Iops,
                          nops, fnops, Islots, Fslots );
          start = cycle;
          Fadds=0;          Fmuls=0;          Fmadds=0;         loads=0;
          stores=0;         Fsubs=0;          Fslots=0;         Islots=0;
          Iops=0;           nops=0;           fnops=0;
fprintf(fptr,"=>%s ---------------------------------------------\n",
	make_label(lab));

	 }

        if ( inst->pipe==DIRECTIVE && inst->code == Exit_Routine ){
          done = 1;
        }

        if ( inst->pipe != DIRECTIVE){
          fprintf(fptr,"%s\t",inst->mnemonic);
	}
        else fprintf(fptr,"  -  \t");

        if ( inst->pipe == FADDPIPE ) Fadds ++;
        else if ( inst->pipe == FMULPIPE ) Fmuls ++;
        else if ( inst->pipe == FMACPIPE ) Fmadds ++;
        else if ( inst->pipe == FLODPIPE ) loads ++;
        else if ( inst->pipe == FSTRPIPE ) stores ++;

        // MIC pipelines {...
        else if (inst->pipe == SIMDMICPIPE && (inst->code == MIC_ADD || inst->code == MIC_SUB)) Fadds += (PROC->nsimd*2);
        else if (inst->pipe == SIMDMICPIPE && (inst->code == MIC_MUL)) Fmuls += (PROC->nsimd*2);
        else if (inst->pipe == SIMDMICPIPE && (inst->code == MIC_FMADD213 || inst->code == MIC_FMADD231 || inst->code == MIC_FNMADD213 || inst->code == MIC_FNMADD231)) {
			Fadds += (PROC->nsimd*2);
			Fmuls += (PROC->nsimd*2);
		}
		else if (inst->pipe == SIMDMICSWIZPIPE && (inst->code == MIC_ADD_MASK_SWIZ || inst->code == MIC_SUB_MASK_SWIZ || inst->code == MIC_SUBR_MASK_SWIZ)) Fadds += (PROC->nsimd*2);  // FIXME: mask does reduce actual number of ops!! We neglect that for now...
        else if (inst->pipe == SIMDMICSWIZPIPE && (inst->code == MIC_MUL_SWIZ)) Fmuls += (PROC->nsimd*2);
        // FIXME: SIMD screwed up utilisation report of floating pipe!
        // ...} MIC

        else Iops ++;


        if (BRCHMASK & PROC->pipes[pipe].issuegroups) Islots ++;
        else if (FMULMASK & PROC->pipes[pipe].issuegroups) Fslots ++;
        else if (FADDMASK & PROC->pipes[pipe].issuegroups) Fslots ++;
        else if (FMACMASK & PROC->pipes[pipe].issuegroups) Fslots ++;
        else if (SIMDMICMASK & PROC->pipes[pipe].issuegroups) Fslots += (2*PROC->nsimd);  // XXX: buggy: MIC_MOV does not contribute
        else if (SIMDMICSWIZMASK & PROC->pipes[pipe].issuegroups) Fslots += (2*PROC->nsimd);
      }
      else {
        fprintf(fptr,"  -  \t");
          if (BRCHMASK & PROC->pipes[pipe].issuegroups) { nops ++; Islots++;}
          else if (FMULMASK & PROC->pipes[pipe].issuegroups) {fnops ++; Fslots++;}
          else if (FADDMASK & PROC->pipes[pipe].issuegroups) {fnops ++; Fslots++;}
          else if (FMACMASK & PROC->pipes[pipe].issuegroups) {fnops ++; Fslots++;}
          else fprintf(fptr,"WARNING: unidentified nop\n");
      }
    }
    fprintf(fptr,"\n");
  }

  fclose(fptr);

}

void report_section( FILE *fptr,int start, int stop,
  int Fadds,int Fmuls, int Fmadds,
  int loads, int stores, int Iops,
 int nops, int fnops, int Islots, int Fslots
                   )
{
  fprintf(fptr,"\n################### Section Analysis ################\n");
  fprintf(fptr,"%d Floating operations\n",Fadds+Fmuls+2*Fmadds);
  fprintf(fptr,"%d Multiplies, %d Adds, %d Madds\n",Fmuls,Fadds,Fmadds);
  fprintf(fptr,"%d Floads, %d Fstores => %d bytes in, %d bytes out\n",
               loads,stores,loads*PROC->FP_size,stores*PROC->FP_size);
  fprintf(fptr,"%d Ibox ops, %d Fnops, %d Unops\n",Iops,fnops,nops);

  if ( stop != start )
  fprintf(fptr,"Expected IPC = %.1f\n",
                       ((double)Fslots+Islots - nops -fnops) / (stop-start));
  if ( Fslots != 0 )
  fprintf(fptr,"%d Fslots => %.1f %% utilisation of Floating pipe\n",Fslots,
                           ((double) Fadds+Fmuls+Fmadds)*100.0/Fslots);
  fprintf(fptr,"################### Section Analysis ################\n\n");

  return;

}



int issue_early ( int min_dep, int max_dep,struct inst_queue_entry *thisent,
                  struct instruction *execution_map[][MAX_WIDTH],
                  struct processor *proc )
{
  int cycle;
  int issue_cycle;
  int pipe;
  int issued;
  int dont_issue;
  int conflict;
  int prev_cycle;
  int sparse_stall;
  struct instruction *inst;

  /*
   * Find the earliest cycle to which we can issue
   * based upon instruction dependencies
   *
   */
  inst = thisent->inst;
  issue_cycle = deps_could_issue(min_dep,max_dep,inst,execution_map,proc);

  /*
   * have to slot the instruction wherever possible
   * on or after issue cycle
   */
  cycle = issue_cycle ;
  issued = 0;

  while ( !issued ) {

    for(pipe = 0; pipe < proc->npipes && (!issued); pipe++ ){

      if ( execution_map[cycle][pipe] == NULL ){

        if (  (proc->pipes[pipe].issuegroups & (1 << inst->pipe ))
           || (inst->pipe == DIRECTIVE ) ){
	  dont_issue = 0;

	  if ( cycle > 2 ){

	    for (conflict = proc->npipes -1; conflict >=0 ; conflict -- ){


	      /*
	       *Check for store-store spacing > store_space
	       */
	      for ( prev_cycle = cycle-store_space;
		    prev_cycle < cycle ; prev_cycle++) {

 	        if (execution_map[prev_cycle][conflict] != NULL ) {

		  if (  inst->pipe == FSTRPIPE 	&&
			execution_map[prev_cycle][conflict]->pipe==FSTRPIPE ){

  		    dont_issue = 1;
		    inst->staller = &
		      PROC->iset->instructions[WEIRDRULES][CACHESTALL];
		  }
		}
	      }

	      /*
	       *Check for load-load spacing > load_space
	       */
	      for ( prev_cycle = cycle-load_space;
		    prev_cycle < cycle ; prev_cycle++) {

 	        if (execution_map[prev_cycle][conflict] != NULL ) {

		  if (  ((inst->pipe == FLODPIPE) || (inst->pipe == CACHPIPE)) &&
			((execution_map[prev_cycle][conflict]->pipe==FLODPIPE)
		       ||(execution_map[prev_cycle][conflict]->pipe==CACHPIPE))
                         ){

		    fprintf(stderr,"Stalling due to load-space = %d\n",load_space );
  		    dont_issue = 1;
		    inst->staller = &
		      PROC->iset->instructions[WEIRDRULES][CACHESTALL];
		  }
		}
	      }


	      /*
	       *Check for dcbt-dcbt spacing > dcbt_space
	       */
	      for ( prev_cycle = cycle-dcbt_space;
		    prev_cycle < cycle ; prev_cycle++) {

 	        if (execution_map[prev_cycle][conflict] != NULL ) {

		    if (  (inst->pipe == CACHPIPE )
			  &&  (inst->code == execution_map[prev_cycle][conflict]->code)
			  &&  (execution_map[prev_cycle][conflict]->pipe==CACHPIPE )
			  //			  &&  (execution_map[prev_cycle][conflict]->code == PREF_IMM)
			){

  		    dont_issue = 1;
		    inst->staller = &
		      PROC->iset->instructions[WEIRDRULES][CACHESTALL];
		  }
		}
	      }

	      /*
	       *Check for l/s-dcbt spacing > dcbt_pre when issuing dcbt
	       */
	      for ( prev_cycle = cycle-dcbt_pre; prev_cycle < cycle ; prev_cycle++) {

 	        if (execution_map[prev_cycle][conflict] != NULL ) {

		  if (  inst->pipe == CACHPIPE
 		     && ( execution_map[prev_cycle][conflict]->pipe==FSTRPIPE
                        ||execution_map[prev_cycle][conflict]->pipe==FLODPIPE
                        ||execution_map[prev_cycle][conflict]->pipe==LOADPIPE
                        ||execution_map[prev_cycle][conflict]->pipe==STORPIPE
			) ) {

  		    dont_issue = 1;
		    inst->staller = &
		      PROC->iset->instructions[WEIRDRULES][CACHESTALL];
		  }
		}

	      }

	      /*
	       *Check for l/s-dcbt spacing > dcbt_pre when issuing l/s
	       */
	      for ( prev_cycle = cycle+1; prev_cycle <= cycle+dcbt_pre ; prev_cycle++) {

		if (execution_map[prev_cycle][conflict] != NULL ) {

		  if ( ( inst->pipe == FSTRPIPE
                      || inst->pipe == FLODPIPE
                      || inst->pipe == LOADPIPE
                      || inst->pipe == STORPIPE
			 )
		       && execution_map[prev_cycle][conflict]->pipe==CACHPIPE ){


  		    dont_issue = 1;
		    inst->staller = &
		      PROC->iset->instructions[WEIRDRULES][CACHESTALL];
		  }
		}

	      }




	      /*
	       *Check for dcbt-ls spacing > dcbt_post when issuing dcbt
	       */
	      for ( prev_cycle = cycle+1;
		    prev_cycle <= cycle+dcbt_post ; prev_cycle++) {

 	        if (execution_map[prev_cycle][conflict] != NULL ) {

		  if (  ( inst->pipe == CACHPIPE && inst->code == PREF_IMM )
 		     && ( execution_map[prev_cycle][conflict]->pipe==FSTRPIPE
                        ||execution_map[prev_cycle][conflict]->pipe==FLODPIPE
                        ||execution_map[prev_cycle][conflict]->pipe==LOADPIPE
                        ||execution_map[prev_cycle][conflict]->pipe==STORPIPE
			  ) ) {


  		    dont_issue = 1;
		    inst->staller = &
		      PROC->iset->instructions[WEIRDRULES][CACHESTALL];
		  }
		}
	      }


	      /*
	       *Check for dcbt-l/s spacing > dcbt_pre when issuing l/s
	       */
	      for ( prev_cycle = cycle-dcbt_post;
		    prev_cycle < cycle ; prev_cycle++) {

		if (execution_map[prev_cycle][conflict] != NULL ) {

		  if ( ( inst->pipe == FSTRPIPE
		      || inst->pipe == FLODPIPE
		      || inst->pipe == LOADPIPE
		      || inst->pipe == STORPIPE
			 )
		       && (
			   execution_map[prev_cycle][conflict]->pipe==CACHPIPE
   		       && execution_map[prev_cycle][conflict]->code ==PREF_IMM)
			   ){


  		    dont_issue = 1;
		    inst->staller = &
		      PROC->iset->instructions[WEIRDRULES][CACHESTALL];
		  }
		}

	      }


	      /*
	       *Check for LMQ conflicts
	      for ( prev_cycle = cycle-PROC->lmq_space;
		    prev_cycle < cycle ; prev_cycle++) {

		if (execution_map[prev_cycle][conflict] != NULL ) {
		  // Apply checks between all indexed loads
		  if ( (inst->pipe == FLODPIPE)
		       && (execution_map[prev_cycle][conflict]->pipe==FLODPIPE)
		       && (execution_map[prev_cycle][conflict]->argtypes[1]==Iregs)
		       && (inst->pipe==FLODPIPE)
		       && (inst->argtypes[1] == Iregs)

			   ){


		    int index_A = get_constant_reg(execution_map[prev_cycle][conflict]->args[1]);
		    int index_B = get_constant_reg(inst->args[1]);
		    if ( (index_A & 0xFFC0) == (index_B & 0xFFC0) ) {
		      dont_issue = 1;
		      inst->staller = &PROC->iset->instructions[WEIRDRULES][CACHESTALL];
		    }

		  }
		}

	      }
	       */




	    }




	    /*
	     *Check for store-store runs > store_lim
	     */
            if ( inst ->pipe == FSTRPIPE ) {

 	     for (conflict = proc->npipes -1; conflict >=0 ; conflict -- ){
	      sparse_stall = 1 ;
	      for ( prev_cycle = MAX(cycle-store_lim,0) ;
                    prev_cycle < cycle ; prev_cycle ++ ) {

		if (execution_map[prev_cycle][pipe]!= NULL ){

		  if ( execution_map[prev_cycle][pipe]->pipe!=FSTRPIPE ) {
  	           sparse_stall = 0 ;
		  }
		} else {
		  sparse_stall = 0 ;
		}

	      }
	     }
	    /*
	     *Check for l/s-store runs > ls_lim
	     */
	     if ( sparse_stall ) dont_issue = 1;

 	     for (conflict = proc->npipes -1; conflict >=0 ; conflict -- ){
	      sparse_stall = 1 ;
	      for ( prev_cycle = MAX(cycle-ls_lim,0) ;
		    prev_cycle < cycle ; prev_cycle ++ ) {

		if (execution_map[prev_cycle][pipe]!= NULL ){

		  if (  execution_map[prev_cycle][pipe]->pipe!=FSTRPIPE
                     && execution_map[prev_cycle][pipe]->pipe!=FLODPIPE ) {
		    sparse_stall = 0 ;
		  }
		} else {
		  sparse_stall = 0 ;
		}

	      }
	     }
	     if ( sparse_stall ) dont_issue = 1;

	    }

	    /*
	     *Check for load-load runs > ls_lim
	     */
            if ( inst ->pipe == FLODPIPE ) {

 	     for (conflict = proc->npipes -1; conflict >=0 ; conflict -- ){
	      sparse_stall = 1 ;
	      for ( prev_cycle = MAX(cycle-load_lim,0) ;
		    prev_cycle < cycle ; prev_cycle ++ ) {

		if (execution_map[prev_cycle][pipe]!= NULL ){

		  if ( execution_map[prev_cycle][pipe]->pipe!=FLODPIPE ) {
  	           sparse_stall = 0 ;
		  }
		} else {
		  sparse_stall = 0 ;
		}

	      }
	     }
	     if ( sparse_stall ) dont_issue = 1;

	    /*
	     *Check for ls-load runs > ls_lim
	     */
 	     for (conflict = proc->npipes -1; conflict >=0 ; conflict -- ){
	      sparse_stall = 1 ;
	      for ( prev_cycle = MAX(cycle-ls_lim,0) ;
		    prev_cycle < cycle ; prev_cycle ++ ) {

		if (execution_map[prev_cycle][pipe]!= NULL ){

		  if (  execution_map[prev_cycle][pipe]->pipe!=FSTRPIPE
                     && execution_map[prev_cycle][pipe]->pipe!=FLODPIPE ) {
  	           sparse_stall = 0 ;
		  }
		} else {
		  sparse_stall = 0 ;
		}

	      }
	     }
	     if ( sparse_stall ) dont_issue = 1;

	    }

	  }
	  if ( ! dont_issue ){
	    execution_map[cycle][pipe] = inst;
	    issued = 1;
	  }

	}
      }
    }
    if ( !issued )  cycle ++;
    if ( cycle >= MAX_CYCLES ) {
       puts(" Can't issue : hit max cycles ");
       inst_print(inst);
       printf("%d, %d\n",inst->pipe,CACHPIPE);
       exit(0);
    }
  }
  return (cycle);
}

extern int pregs[];
extern int pimms[];
int find_offset(int off_reg)
{
  int idx = -1;
  for(int i=0;i<24;i++) if ( pregs[i] == off_reg ) idx = i;
  if ( idx == -1 ) exit(0);
  return get_offset(pimms[idx]);
}

int dependency ( struct instruction * i1, struct instruction *i2 )
{
  int dep;
  int i,j;

  dep = NO_DEPENDENCY;

  /*Never reorder two condition code updates*/
  if ( ( i1->pipe == IALUPIPE && i1->code == IOR_PREDICATE )
   &&  ( i2->pipe == IALUPIPE && i2->code == IOR_PREDICATE ) )
    return (BARRIER);
  if ( ( i1->pipe == IALUPIPE && i1->code == IOR_PREDICATE )
   &&  ( i2->pipe == IALUPIPE && i2->code == IAND_IMM ) )
    return (BARRIER);
  if ( ( i2->pipe == IALUPIPE && i2->code == IOR_PREDICATE )
   &&  ( i1->pipe == IALUPIPE && i1->code == IAND_IMM ) )
    return (BARRIER);
  if ( ( i2->pipe == IALUPIPE && i2->code == IAND_IMM )
   &&  ( i1->pipe == IALUPIPE && i1->code == IAND_IMM ) )
    return (BARRIER);

  if (PROC->id == KNC || PROC->id == KNC_SINGLE) {
	/* Don't reorder condition code update with IADD_IMM or IADD, which modifies %rflags */
	if ( ( i1->pipe == IALUPIPE && i1->code == IOR_PREDICATE )
	 &&  ( i2->pipe == IALUPIPE && i2->code == IADD_IMM ) )
		return (BARRIER);
    if ( ( i1->pipe == IALUPIPE && i1->code == IADD_IMM )
     &&  ( i2->pipe == IALUPIPE && i2->code == IOR_PREDICATE ) )
		return (BARRIER);

	if ( ( i1->pipe == IALUPIPE && i1->code == IOR_PREDICATE )
	 &&  ( i2->pipe == LOADPIPE && i2->code == LOAD_ADDR ) )  // NB: implemented as ADDQ on MIC
		return (BARRIER);
    if ( ( i1->pipe == LOADPIPE && i1->code == LOAD_ADDR )
     &&  ( i2->pipe == IALUPIPE && i2->code == IOR_PREDICATE ) )
		return (BARRIER);
  }

  // BGQ can suppress reorder of loads as the LMQ is a bitch
  if ( PROC->load_reorder == 0 ) {
    if ( i1->pipe == FLODPIPE && i2->pipe == FLODPIPE ) {
      return(BARRIER);
    }
    if ( i1->pipe == FLODPIPE && i2->pipe == CACHPIPE ) {
      return(BARRIER);
    }
    if ( i1->pipe == CACHPIPE && i2->pipe == FLODPIPE ) {
      return(BARRIER);
    }
    if ( i1->pipe == CACHPIPE && i2->pipe == CACHPIPE ) {
      return(BARRIER);
    }
  }


  /*Never reorder two condition code updates*/
  int excl1 =0;
  int excl2 =0;

  if ( ( i1->pipe == IALUPIPE && i1->code == IOR_PREDICATE )
     ||( i1->pipe == IALUPIPE && i1->code == CMOVGE )
     ||( i1->pipe == IALUPIPE && i1->code == CMOVGT )
     ||( i1->pipe == IALUPIPE && i1->code == CMOVLE )
     ||( i1->pipe == IALUPIPE && i1->code == CMOVLT ) ) excl1 = 1;

  if ( ( i2->pipe == IALUPIPE && i2->code == IOR_PREDICATE )
     ||( i2->pipe == IALUPIPE && i2->code == CMOVGE )
     ||( i2->pipe == IALUPIPE && i2->code == CMOVGT )
     ||( i2->pipe == IALUPIPE && i2->code == CMOVLE )
     ||( i2->pipe == IALUPIPE && i2->code == CMOVLT ) ) excl2 = 1;

  if ( excl1 && excl2 ) return BARRIER;

  if ( i1->pipe == DIRECTIVE && i1->code == LS_BARRIER ){
    if ( i2->pipe == FLODPIPE || i2->pipe == FSTRPIPE ) {
      return (BARRIER);
    }
    if ( i2->pipe == LOADPIPE || i2->pipe == STORPIPE ) {
      return (BARRIER);
    }
    return(NO_DEPENDENCY);
  }
  if ( i2->pipe == DIRECTIVE && i2->code == LS_BARRIER ){
    if ( i1->pipe == FLODPIPE || i1->pipe == FSTRPIPE ) {
      return (BARRIER);
    }
    if ( i1->pipe == LOADPIPE || i1->pipe == STORPIPE ) {
      return (BARRIER);
    }
    return(NO_DEPENDENCY);
  }
  if ( i1->pipe == FSTRPIPE && i2->pipe == FSTRPIPE ) {
    if ( store_noreorder ) return (BARRIER);
  }


  /*Check for ASM directives - never raise*/
  if ( i1->pipe == DIRECTIVE || i2->pipe == DIRECTIVE) {
    dep = BARRIER;
    return(dep);
  }

  /*Check for CACHEOPS directives - reorder*/
  if ( i1->pipe == CACHPIPE && i2->pipe == CACHPIPE) {
    dep = BARRIER;
    return(dep);
  }

  /*Check for Branch instructions - never raise*/
  if ( i1->pipe == BRCHPIPE || i2->pipe == BRCHPIPE) {
    dep = BARRIER;
    return(dep);
  }

  if ( i1->narg <=0 || i2->narg <= 0 ) return ( NO_DEPENDENCY );

  /*
   * Preserve load ordering
   */
  if ( BGQloads ) {
    if ( ( i1->pipe == FLODPIPE) && ( i2->pipe == FLODPIPE) ) {
      dep = READ_A_WRITE;
      return dep;
    }
  }

  /*BGQ load spacing to avoid LMQ collision
  if ( BGQloads ) {
    if ( ( i1->pipe == FLODPIPE)
      && ( (i1->code == CLOAD_INDEXED)||(i1->code == CSLOAD_INDEXED))
      && ( i2->pipe == FLODPIPE)
      && ( (i2->code == CLOAD_INDEXED)||(i2->code == CSLOAD_INDEXED))
      && ( i1->args[2] == i2->args[2] )
      && ( (find_offset(i1->args[1])&0xFFC0) ==
	   (find_offset(i2->args[1])&0xFFC0) )
	 ) {
      fprintf(stderr,"Stalling QPX load %d <- %d(%d)\n",i1->args[0],find_offset(i1->args[1]),i1->args[2]);
      fprintf(stderr,"Staller           %d <- %d(%d)\n",i2->args[0],find_offset(i2->args[1]),i2->args[2]);
      fprintf(stderr,"Staller           %lx %lx\n",i1,i2);
      dep = READ_A_WRITE;
      return(dep);
      }
  }
  */


  /*Check for READ after WRITE */
  for(i=0;i<i2->narg;i++){
    if ( i2->inout[i] == outarg &&
         (i2->argtypes[i] == Fregs || i2->argtypes[i] == Iregs ) ){

      for(j=0;j<i1->narg;j++){
        if ( i1->argtypes[j] == i2->argtypes[i] && i1->args[j]==i2->args[i]
             && i1->inout[j] == inarg){
          dep = READ_A_WRITE;
          return(dep);
	}
      }
    }
  }


  /*Check for WRITE after READ or WRITE */
  for(i=0;i<i2->narg;i++){
    if ( i2->argtypes[i] == Fregs || i2->argtypes[i] == Iregs  ){

      for(j=0;j<i1->narg;j++){
        if ( i1->argtypes[j] == i2->argtypes[i] && i1->args[j]==i2->args[i]
             && i1->inout[j] == outarg){
          dep = WRITE_A_READ;
          return(dep);
	}
      }
    }
  }



  return(dep);
}

void queue_instruction(struct instruction *inst)
{
  struct inst_queue_entry *ent ;
  struct inst_queue_entry *tail;


  ent = new inst_queue_entry;
  if ( !ent ){ perror("Malloc failed"); exit(0);}

  ent -> inst = inst;
  ent -> prev = NULL;
  ent -> next = NULL;

  /*If first call start the list*/
  if ( HEAD == NULL ) {
    HEAD = ent;
    return;
  }
  /*Otherwise */
  tail = HEAD;
  while ( tail->next != NULL )
    {
      tail = tail -> next;
    };
  tail->next = ent;
  ent->prev = tail;
}
/* ...END OPTIMISATION ROUTINES}*/





