#include "bfm.h"
#include <stdio.h>
#include <stdlib.h>

int main (int argc,char **argv )
{
  if ( argc != 5 ) { 
   printf("Usage: %s lx ly lz lt \n All must be even\n",argv[0]);
   printf("argc is %d\n",argc);
    for ( int i=0;i<argc;i++)
      printf("%d %s\n",i,argv[i]);
    exit(-1);
  }

  int lx = atoi(argv[1]);
  int ly = atoi(argv[2]);
  int lz = atoi(argv[3]);
  int lt = atoi(argv[4]);

  int nrow[4];
  nrow[0] = lx;
  nrow[1] = ly;
  nrow[2] = lz;
  nrow[3] = lt;

  bfmarg dwfa;
  dwfa.solver = WilsonFermion;

  bfm    dwf;

  dwfa.node_latt[0]  = lx;
  dwfa.node_latt[1]  = ly;
  dwfa.node_latt[2]  = lz;
  dwfa.node_latt[3]  = lt;

  dwfa.local_comm[0]  = 1;
  dwfa.local_comm[1]  = 1;
  dwfa.local_comm[2]  = 1;
  dwfa.local_comm[3]  = 1;

  dwfa.Ls = 1;
  dwfa.mass = 0.0;
  dwfa.Csw  = 0.0;

  printf("Initialising bfm operator\n");
  dwf.init(dwfa);

  Fermion_t psi_h = dwf.allocFermion();
  Fermion_t chi_h = dwf.allocFermion();

  // Naive Dslash

  // cb is cb of result, 1-cb is cb of input field
  for(int cb=0;cb<2;cb++){

    /*Import this checkerboard of QDP fields to bagel*/
    printf("Importing psi field cb %d\n",cb);

    // Fill the other checkerboard.
    for(int dag=0;dag<2;dag++){

      printf("Checking cb=%d dag=%d\n",cb,dag); 

      dwf.dslash(psi_h,
                 chi_h, 
		 cb,dag);

      printf("Done\n"); 

    }
  }

}







