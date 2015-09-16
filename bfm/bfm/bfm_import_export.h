#ifndef BFM_IMP_EX_H
#define BFM_IMP_EX_H
/*Peter Boyle April 2007*/
#include "bfm.h"
#include <omp.h>

#include <stdio.h>
#include <stdlib.h>

extern "C" {
void   qpx_single_to_double (int N, float *  in, double *out);
void   qpx_double_to_single (int N, double *  in, float *out);
}

/////////////////////////////
// Ls=16 -> Ls=8; 
// depth = 3, block = 5
// 3+3+5+5 -> 3+4+1+1 multigrid
/////////////////////////////
/*
template <class Float>
void bfmbase<Float>::ProjReduceLs (int Block,std::vector<double> &Weights, Fermion_t in, Fermion_t out, Fermion_t tmp)
{
  integer sLs = Ls-2*Block+2;
  integer bLs = Ls;
  int depth = (bLs - 2*Block)/2; 
  double nn;
  int me = this->thread_barrier();

  // Verbatim copy near surface
  fill(tmp,0.0);
  for(int s=0;s<depth;s++){
    axpby_ssp(tmp,0.0,in,1.0,in,s,s);
    axpby_ssp(tmp,0.0,in,1.0,in,sLs-1-s,bLs-1-s);
  }

  //
  // Take example: Ls=16, block=5, sLs=8
  //
  // depth = 3;
  // Middle vectors are (0,0,0, 1,1,1,1,1, 0,0,0,0,0,  0,0,0 )/sqrt(5)   [normalised]
  //                    (0,0,0, 0,0,0,0,0, 1,1,1,1,1,  0,0,0 )/sqrt(5)
  //
  axpby_ssp(tmp,0.0,in,0.0,in,depth,depth);
  axpby_ssp(tmp,0.0,in,0.0,in,depth+1,depth+1);
  for(int s=0;s<Block;s++){
    double scale;

    //0,1,2  3         4              5,6,7
    //0,1,2  3,4,5,6,7 8,9,10,11,12  13,14,15

    scale = Weights[s];
    axpby_ssp(tmp,1.0,tmp,scale,in,depth,depth+s);

    scale= Weights[Block-1-s];
    axpby_ssp(tmp,1.0,tmp,scale,in,depth+1,depth+s+Block);   
  }                                                     

  // Must re-layout for shrunken Ls

  int thrlen, throff;
  int work =this->node_cbvol;
  this->thread_work(work,me,thrlen,throff);
  int Nspinco=12;
  Float *o_p = (Float *)out;
  Float *t_p = (Float *)tmp;
  for (int site=throff;site<throff+thrlen;site++ ) { 
    for(int s=0;s<sLs;s++){
      for ( int co=0;co<Nspinco;co++ ) { 
      for ( int reim=0;reim<2;reim++ ) { 
	int s_idx = this->bagel_idx5d_tmp (sLs,site,s,reim,co,Nspinco) ;
	int b_idx = this->bagel_idx5d_tmp (bLs,site,s,reim,co,Nspinco) ;
	o_p[s_idx]=t_p[b_idx];
      }}
    }
  }
  this->thread_barrier();
}
template <class Float>
void bfmbase<Float>::PromIncreaseLs (int Block,std::vector<double> &Weights, Fermion_t in, Fermion_t out, Fermion_t tmp)
{
  int sLs = Ls-2*Block+2;
  int bLs = Ls;
  int depth = (bLs - 2*Block)/2; 
  double nn;

  // Must re-layout for enlarged Ls
  int thrlen, throff;
  int work =this->node_cbvol;
  Float *i_p = (Float *)in;
  Float *t_p = (Float *)tmp;
  int me;
  this->fill(tmp,0.0);
  this->fill(out,0.0);
  this->thread_work(work,me,thrlen,throff);
  int Nspinco=12;
  for (int site=throff;site<throff+thrlen;site++ ) { 
    for(int s=0;s<sLs;s++){
      for ( int co=0;co<Nspinco;co++ ) { 
      for ( int reim=0;reim<2;reim++ ) { 
	int s_idx = bagel_idx5d_tmp (sLs,site,s,reim,co,Nspinco) ;
	int b_idx = bagel_idx5d_tmp (bLs,site,s,reim,co,Nspinco) ;
	t_p[b_idx]=i_p[s_idx];
      }}
    }
  }
  this->thread_barrier();

  // Verbatim copy near surface
  for(int s=0;s<depth;s++){
    axpby_ssp(out,0.0,tmp,1.0,tmp,s,s);
    axpby_ssp(out,0.0,tmp,1.0,tmp,bLs-1-s,sLs-1-s);
  }
  // Middle vectors are (0,0,0, 1,1,1,1,1, 0,0,0,0,0,  0,0,0 )/sqrt(5)   [normalised]
  //                    (0,0,0, 0,0,0,0,0, 1,1,1,1,1,  0,0,0 )/sqrt(5)
  for(int s=0;s<Block;s++){
    double scale;
    scale= Weights[s];
    axpby_ssp(out,0.0,tmp,scale,tmp,depth+s,depth);           //0,1,2  3         4              5,6,7

    scale= Weights[Block-1-s];
    axpby_ssp(out,0.0,tmp,scale,tmp,depth+s+Block,depth+1);   //0,1,2  3,4,5,6,7 8,9,10,11,12  13,14,15
  }                                                     
}
*/

template <class Float>
void      bfmbase<Float>::precisionChange (Fermion_t in, Fermion_t out, int control, int cb)
{
  int thrlen, throff;
  int work =this->node_cbvol*this->cbLs;
  int me;
  
  this->thread_work(work,me,thrlen,throff);

#ifdef BGQ
#error
  // argument is number of complex simd vectors
  int idx=24*throff;
  int swords = (12*thrlen)/this->nsimd;//simdwords
  double *dp;
  float  *fp;
  if ( control == SingleToDouble ) {
    dp = (double *)out;
    fp = (float  *)in;
    qpx_single_to_double (swords,&fp[idx],&dp[idx]);
  } else { 
    dp = (double *)in;
    fp = (float  *)out;
    qpx_double_to_single (swords,&dp[idx],&fp[idx]);
    }
#else
  for (int site=throff;site<thrlen+throff;site++ ) { 
    int idx=24*site;
    double *dp;
    float  *fp;
    if ( control == SingleToDouble ) {
      dp = (double *)out;
      fp = (float  *)in;
      for(int i=0;i<24;i++)      dp[idx+i] = fp[idx+i];
    } else { 
      dp = (double *)in;
      fp = (float  *)out;
      for(int i=0;i<24;i++)      fp[idx+i] = dp[idx+i];
    }
  }
#endif 
  thread_barrier();
}

template <class Float>
void      bfmbase<Float>::importFermion  (QDPdouble *psi, Fermion_t handle,int cb)
{
  impexFermion(psi,handle,1,cb,0);
}
template <class Float>
void      bfmbase<Float>::exportFermion  (QDPdouble *psi, Fermion_t handle,int cb)
{
  impexFermion(psi,handle,0,cb,0);
}

template <class Float>
integer bfmbase<Float>::chroma_idx(int x[4],int reim,int i, int i_size)
{
  int ccb   = ((x[0]+x[1]+x[2]+x[3])&0x1); /*Work out local checkerboard of site*/ 
                         /*FIXME Here need to worry about base_parity
                          *as chroma indexes with global and not local parity
                          */
  int chroma_cbsite[4];
  chroma_cbsite[0] = x[0];
  chroma_cbsite[1] = x[1];
  chroma_cbsite[2] = x[2];
  chroma_cbsite[3] = x[3];

  int csite= x[0] + node_latt[0]*(x[1]
                  +       node_latt[1]*(x[2]
                  +             node_latt[2]*x[3]));
  csite = csite/2;
  int cbvol = (node_latt[0]*node_latt[1]*node_latt[2]*node_latt[3])/2;
  int idx = (ccb*cbvol+csite)*i_size*2 + i*2 + reim;

  return idx;  
}

// Double store the gauge field.
// Plan : do the cshift & adj in QDP, not in my code!
template <class Float>
void      bfmbase<Float>::importGauge(QDPdouble *gauge, int dir)
{
  int Ndircoco=72; /*8 directions stored*/
  int Ncoco=9;
  omp_set_num_threads(this->nthread);
#pragma omp parallel 
  {

#pragma omp for 
  for (int site=0;site<node_latt[0]*node_latt[1]*node_latt[2]*node_latt[3];site++ ) { 
    
    int x[4] ;
    int s=site;
    x[0]=s%node_latt[0];    s=s/node_latt[0];
    x[1]=s%node_latt[1];    s=s/node_latt[1];
    x[2]=s%node_latt[2];    s=s/node_latt[2];
    x[3]=s%node_latt[3];
    
      for ( int coco=0;coco<9;coco++ ) { 
      for ( int reim=0;reim<2;reim++ ) { 

        Float * bagel = this->u;
	int bbase = dir*9;
        int cidx = chroma_idx(x,reim,coco,Ncoco);
        int bidx = bagel_idx(x,reim,coco+bbase,Ndircoco,0);
        bagel[bidx] = gauge[cidx];


      }}
  }
  }
}



template <class Float>
void      bfmbase<Float>::exportGauge(QDPdouble *gauge, int dir)
{
  int Ndircoco=72; /*8 directions stored*/
  int Ncoco=9;
  omp_set_num_threads(this->nthread);
#pragma omp parallel 
  {

#pragma omp for 
  for (int site=0;site<node_latt[0]*node_latt[1]*node_latt[2]*node_latt[3];site++ ) { 
    
    int x[4] ;
    int s=site;
    x[0]=s%node_latt[0];    s=s/node_latt[0];
    x[1]=s%node_latt[1];    s=s/node_latt[1];
    x[2]=s%node_latt[2];    s=s/node_latt[2];
    x[3]=s%node_latt[3];
    
      for ( int coco=0;coco<9;coco++ ) { 
      for ( int reim=0;reim<2;reim++ ) { 

        Float * bagel = this->u;
	int bbase = dir*9;
        int cidx = chroma_idx(x,reim,coco,Ncoco);
        int bidx = bagel_idx(x,reim,coco+bbase,Ndircoco,0);
        gauge[cidx] = bagel[bidx];

      }}
  }
  }
}


template <class Float>
void      bfmbase<Float>::exportForce(Matrix_t handle,QDPdouble *force, int dir,int cb) 
{
  int Ndircoco=36; /*4 directions stored*/
  int Ncoco=9;
  omp_set_num_threads(this->nthread);
#pragma omp parallel 
  {

#pragma omp for 
  for (int site=0;site<node_latt[0]*node_latt[1]*node_latt[2]*node_latt[3];site++ ) { 
    
    int x[4] ;
    int s=site;
    x[0]=s%node_latt[0];    s=s/node_latt[0];
    x[1]=s%node_latt[1];    s=s/node_latt[1];
    x[2]=s%node_latt[2];    s=s/node_latt[2];
    x[3]=s%node_latt[3];
    
    if ( ((x[0]+x[1]+x[2]+x[3])&0x1) == cb ) {      

      for ( int coco=0;coco<9;coco++ ) { 
      for ( int reim=0;reim<2;reim++ ) { 

        Float * bagel = (Float *)handle;
	int bbase = dir*9;
        int cidx = chroma_idx(x,reim,coco,Ncoco);
        int bidx = bagel_idx (x,reim,coco+bbase,Ndircoco,1);
        force[cidx] = bagel[bidx];

      }}
    }
  }
  }
}


template <class Float>
void      bfmbase<Float>::impexFermion  (QDPdouble *psi, Fermion_t handle,int doimport,int cb,int s)
{
  int Nspinco=12;
  omp_set_num_threads(this->nthread);

#pragma omp parallel
  {
#pragma omp for 
  for (int site=0;site<node_latt[0]*node_latt[1]*node_latt[2]*node_latt[3];site++ ) { 
    
    int x[4] ;
    int si=site;
    x[0]=si%node_latt[0];    si=si/node_latt[0];
    x[1]=si%node_latt[1];    si=si/node_latt[1];
    x[2]=si%node_latt[2];    si=si/node_latt[2];
    x[3]=si%node_latt[3];
    
    int sp;
    if ( precon_5d ) sp = s;
    else sp = 0;

    if ( ((x[0]+x[1]+x[2]+x[3]+sp)&0x1) == cb ) {      
      for ( int co=0;co<Nspinco;co++ ) { 
      for ( int reim=0;reim<2;reim++ ) { 
        Float * bagel = (Float *)handle;
        int bidx = bagel_idx5d(x,s,reim,co,Nspinco,1);
        int cidx = chroma_idx(x,reim,co,Nspinco);
	if ( doimport ) bagel[bidx] = psi[cidx]; 
	else psi[cidx] = bagel[bidx] ;
      }}
    }      

  }
  }
}

template <class Float>
Fermion_t bfmbase<Float>::threadedAllocFermion   (int mem_type)
{
  Fermion_t ret;
  int me = thread_barrier();
  if ( me == 0 ) {
    ret = this->allocFermion(mem_type);
  }
  ret = this->thread_bcast(me,ret);
  thread_barrier();
  return ret;
}
template <class Float>
void      bfmbase<Float>::threadedFreeFermion    (Fermion_t handle)
{
  int words = 24 * nsimd * (simd_cbvol);
  int bytes = words*cbLs*sizeof(Float);
  int me = thread_barrier();
  if ( me == 0 ) { 
    Float * ferm = ( Float * ) handle;
    bfm_free(ferm,bytes);
  }
  thread_barrier();
}

template <class Float>
Fermion_t bfmbase<Float>::allocFermion   (int mem_type)
{
  int words = 24 * nsimd * simd_cbvol *cbLs;
  int bytes = words*sizeof(Float);
  Float *ferm = (Float *)bfm_alloc(bytes,mem_type);

  if(ferm == 0){
    this->Error("bfmbase::allocFermion\n");
    fflush(stdout);
    exit(-1);
  }
  return (Fermion_t)ferm;
}
template <class Float>
void      bfmbase<Float>::freeFermion    (Fermion_t handle)
{
  int words = 24 * nsimd * (simd_cbvol);
  int bytes = words*cbLs*sizeof(Float);
  Float * ferm = ( Float * ) handle;
  bfm_free(ferm,bytes);
}

	template <class Float>
Fermion_t bfmbase<Float>::threadedAllocCompactFermion   (int mem_type)
{
  Fermion_t ret;
  int me = this->thread_barrier();
  if ( me == 0 ) {
    ret = this->allocCompactFermion(mem_type);
  }
  ret = this->thread_bcast(me,ret);
  this->thread_barrier();
  return ret;
}

template <class Float>
Fermion_t bfmbase<Float>::allocCompactFermion   (int mem_type)
{
  return allocFermion(mem_type);
}


template <class Float>
void bfmbase<Float>::zeroMatrix(Matrix_t mat)
{
  int ncocodir = 9*4;
  Float *x = (Float *) mat;
  Float *xx;
  Float nrm = 0.0;

  int me,thrlen,throff;
  thread_work(simd_cbvol*nsimd,me,thrlen,throff);

  for(int s=0;s<thrlen;s++){

    xx = &x[(s+throff)*ncocodir*2 ];

    for(int co=0;co<ncocodir;co++){
    for(int reim=0;reim<2;reim++){
      int idx = reim + co*2;
      xx[idx] = 0.0 ;
    }}

  }  
  thread_barrier();
}

template <class Float>
Matrix_t bfmbase<Float>::allocMatrix   (void)
{
  int bytes = 18*4 * nsimd * (simd_cbvol)*sizeof(Float);
  Float *mat = (Float *)bfm_alloc(bytes,mem_slow);

  if(mat == 0){
    this->Error("bfmbase::allocMatrix\n");
    fflush(stdout);
    exit(-1);
  }
  return (Matrix_t)mat;
}
template <class Float>
void      bfmbase<Float>::freeMatrix    (Matrix_t handle)
{

  int bytes = 18*4 * nsimd * (simd_cbvol)*sizeof(Float);
  Float * ferm = ( Float * ) handle;
  bfm_free(ferm,bytes);
}


template <class Float> template<typename FloatEXT>
void bfmbase<Float>::cps_importGauge(FloatEXT *importme)
{
  multi1d<LatticeColorMatrix> U(Nd);

  omp_set_num_threads(this->nthread);

  int Ndircoco=72;
  int Ncoco = 9;
  QDPdouble *U_p;

  int vol4d =
    this->node_latt[0] *
    this->node_latt[1] *
    this->node_latt[2] *
    this->node_latt[3];

  for (int mu=0;mu<Nd;mu++) {
    U_p = (QDPdouble *)&(U[mu].elem(0).elem());

#pragma omp parallel for 
    for (int site=0;site<vol4d;site++ ) {
      int x[4];
      int s=site;
      x[0]=s%this->node_latt[0];    s=s/this->node_latt[0];
      x[1]=s%this->node_latt[1];    s=s/this->node_latt[1];
      x[2]=s%this->node_latt[2];    s=s/this->node_latt[2];
      x[3]=s%this->node_latt[3];
      
      int qidx_base = this->chroma_idx(x, 0, 0, Ncoco);

      for(int coco = 0; coco < Ncoco; ++coco) {
        for ( int reim = 0; reim < 2; ++reim) {
          // int qidx = this->chroma_idx(x,reim,coco,Ncoco);
          int qidx = qidx_base + reim + coco * 2;

          int siteoff = mu + Nd * site;
          int cidx = reim + 2 * (coco + Ncoco * siteoff);
          U_p[qidx] = importme[cidx];
        }} // reim,coco
    } // x
  }//mu

  // to bfm
  this->importGauge(U);
}



#ifdef BFM_QDP_BINDING
//begin karthee-clover
template <class Float>
void      bfmbase<Float>::importClover(multi1d< RScalar <REAL> > &clover_diag, multi1d< RComplex <REAL> > &clover_offdiag, Float *bagel_clover)
{
    int Nelem=42; /*6x6 at each site ; number of complex*/
    int Nspinco=12;
    int Nc=3;
    int vol = this->node_latt[0]*this->node_latt[1]*this->node_latt[2]*this->node_latt[3];

#pragma omp parallel
    {
#pragma omp for
        for (int site=0; site<node_latt[0]*node_latt[1]*node_latt[2]*node_latt[3]; site++ )
        {

            int x[4] ;
            int s=site;
            x[0]=s%node_latt[0];
            s=s/node_latt[0];
            x[1]=s%node_latt[1];
            s=s/node_latt[1];
            x[2]=s%node_latt[2];
            s=s/node_latt[2];
            x[3]=s%node_latt[3];

	    int ccb   = ((x[0]+x[1]+x[2]+x[3])&0x1); /*Work out local checkerboard of site*/
            int cbvol = vol/2;
            int csite = cbvol*ccb+site/2;

	    for(int jj = 0; jj < 2; jj++)
	    {
               for(int ii = 0; ii < 2*Nc; ii++)
               {
		  for ( int reim=0; reim<2; reim++ )
                  {
                      Float * bagel;
                      bagel = bagel_clover;       
                      int shift = jj*21+ii;                                              
                      int bidx = bagel_idx(x,reim,shift,Nelem,0);
                      if(reim == 0)
			bagel[bidx] = clover_diag[csite*12 + jj*6 + ii].elem();
                      else
                         bagel[bidx] = 0.0;                      
                    
                  }
               }
            }
            
            for(int jj = 0; jj < 2; jj++)
            {
                 for(int ii = 0; ii < 2*Nc*Nc-Nc; ii++)
                 {
  	             for ( int reim=0; reim<2; reim++ ) 
                     {
                        Float * bagel;
                        bagel = bagel_clover;    
                        int shift = jj*21+2*Nc+ii;                                         
                        int bidx = bagel_idx(x,reim,shift,Nelem,0);
                       if(reim == 0)
                 	 bagel[bidx] = clover_offdiag[csite*30 + jj*15 + ii].real(); 
                       else
                         bagel[bidx] = clover_offdiag[csite*30 + jj*15 + ii].imag();                                     	        
                     }
                  }
               
            }
            
        }
    }
}
//end karthee-clover
#endif

#endif

