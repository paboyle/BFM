#ifndef BFM_GCR_H
#define BFM_GCR_H

#include "bfm.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

template <class Float>
int bfmbase<Float>::GCR_prec(Fermion_t psi, Fermion_t src)
{

}
template <class Float>
int bfmbase<Float>::GCR_prec_restart(Fermion_t psi, Fermion_t src,int k)
{
  Fermion_t p[k];
  Fermion_t mp;
  Fermion_t mr;
  Fermion_t r;
  complex<double> a;
  double d;
  double b[k];

  r  =threadedAllocFermion(mem_fast); 
  mp =threadedAllocFermion(mem_fast); 
  mr =threadedAllocFermion(mem_fast); 
  for(int i=0;i<k;i++) {
    p[i] =threadedAllocFermion(mem_fast); 
  }


  // Restarts
  int restarts=30;
  
  for(int o=0;o<restarts;o++) {
    

    //Initial guess residual r0 = src - A psi
    //                       p0 = r0
    Mprec(psi,mp,tmp,DaggerNo);
    axpy (r, mp, src,-1.0);
    axpy (p[0], mp, src,-1.0);

    for(int i=0;i<k;i++) { 
    
      d=Mprec(p[i],mp,tmp,DaggerNo);

      a = inner(r,mp)/ d;

      caxpy(psi,psi,p[i],a);
      caxpy(r,r,mp,-a);

      d=Mprec(r,mr,tmp,DaggerNo);


      axpy(p[i+1],r,r,0.0);
      for(int j=0;j<=i;j++) {
	b[j] = -inner(mr,mp)/d;
	axpy(p[i+1],p[i+1],p[j],b[j]);
      }

      if ( ) { 

	threadedFreeFermion(r);
	threadedFreeFermion(mp);
	threadedFreeFermion(x);
	for(int j=0;j<k;j++) {
	  threadedFreeFermion(p[i]);
	}
      }

    }
  }

}


#endif
