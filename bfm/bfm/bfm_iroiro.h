#ifndef BFM_IROIRO_H
#define BFM_IROIRO_H


template <class Float> 
integer iroiro_gauge_idx(int row,int column,int mu,int site[4]) 
{
  int vol = node_latt[0]*node_latt[1]*node_latt[2]*node_latt[3];
  int lex = site[0]+node_latt[0]*(site[1]+node_latt[1]*(site[2]+node_latt[2]*site[3]));

  integer idx = row*6+column*2+18*lex + vol*18*mu;
  return idx;
}

template <class Float> 
integer iroiro_ferm_idx(int col,int spin,int site[4])
{
  int lex = site[0]+node_latt[0]*(site[1]+node_latt[1]*(site[2]+node_latt[2]*site[3]));
  integer idx = 2*col+6*spin+lex*12;
}


template <class Float> 
void      bfmbase<Float>::importGaugeIroIro(QDPdouble *gauge, int dir)
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

#endif
