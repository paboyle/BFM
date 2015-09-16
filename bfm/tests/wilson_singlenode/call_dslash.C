#include <bfm.h>
#include <stdio.h>
#include <stdlib.h>


#undef GENERATE

void * thr_main(void *p);

Fermion_t psi_h;
Fermion_t chi_h;
Fermion_t check;
Fermion_t diff ;

int cb, dag;
double delta;
double n1;
double n2;

//extern double src[];
extern double cb0dag0[];
extern double cb0dag1[];
extern double cb1dag0[];
extern double cb1dag1[];

char * files[4] = { "cb0dag0.C",
                    "cb0dag1.C",
                    "cb1dag0.C",
                    "cb1dag1.C" };

char * array_names[4] = {"cb0dag0",
                         "cb0dag1",
                         "cb1dag0",
                         "cb1dag1" };

double * arrays[4] = { 
  cb0dag0,
  cb0dag1,
  cb1dag0,
  cb1dag1
};

typedef double Float;
class singletest : public bfm {
public:

  void unitGauge(void)
  {
    int vol = node_latt[0]
            * node_latt[1]
            * node_latt[2]
            * node_latt[3];
    int words = 18 * vol;

    int bytes = sizeof(QDPdouble)*words;

    QDPdouble *gauge = (QDPdouble *)malloc(bytes);
    bzero(gauge,bytes);

    // No easy way to support bagel-2      
    // without CHROMA/QDP as I used QDP to do the cshift & adjoint
    // for double stored links. just use unit gauge.
    for(int v=0;v<vol;v++ ) {
      gauge[v*18]    = 1.0;
      gauge[v*18+8]  = 1.0;
      gauge[v*18+16] = 1.0;
    }
    int ndir=4;
    if ( doubleStored() ) { 
      ndir = 8;
    }
    for(int mu=0;mu<ndir;mu++) {
      importGauge(gauge, mu); 
    }

#ifdef GENERATE
    dump(this->u,words*ndir,"gauge.C","gauge");
#endif

 }
  void randFermion(Fermion_t x_t){

    int vol = node_latt[0]
            * node_latt[1]
            * node_latt[2]
            * node_latt[3];
    int words = 24 * vol;
    int bytes = sizeof(QDPdouble)*words;

    QDPdouble *psi = (QDPdouble *)malloc(bytes);

    for(int w=0;w<words;w++ ) {
      psi[w] = drand48();
    }

    importFermion(psi, x_t, 0);

  };

void dump(Fermion_t psi, char *file, char *array_name)
{
  int vol = node_latt[0]
          * node_latt[1]
          * node_latt[2]
          * node_latt[3];

  int words = 24 * vol;

  int bytes = sizeof(QDPdouble)*words;

  double * psi_p = (QDPdouble *)malloc(bytes);
  bzero(psi_p, bytes);

  //  exportFermion(psi_p,psi,0);
  //  dump(psi_p, words, file, array_name);
  dump((double *)psi, words, file, array_name);

}
void dump(double *psi_p, int words, char *file, char *array_name)
{

  FILE * fp = fopen(file,"w");
  fprintf(fp,"double %s[] __attribute__((aligned(32))) = {\n",array_name);
  for(int w=0;w<words;w++){
    fprintf(fp,"\t%le,\n",(double)*psi_p++);
  }
  fprintf(fp,"0.0\n");
  fprintf(fp,"};\n");
  fclose(fp);
};

};

singletest    dwf;


int main (int argc,char **argv )
{

  int generate = 0;
  if ( argc > 1 ) generate = 1;

  int lx = 4;
  int ly = 4;
  int lz = 4;
  int lt = 4;

  int nrow[4];
  nrow[0] = lx;
  nrow[1] = ly;
  nrow[2] = lz;
  nrow[3] = lt;

  bfmarg dwfa;
  dwfa.solver = WilsonFermion;
  dwfa.threads = 1;

  dwfa.node_latt[0]  = lx;
  dwfa.node_latt[1]  = ly;
  dwfa.node_latt[2]  = lz;
  dwfa.node_latt[3]  = lt;

  dwfa.local_comm[0]  = 0;
  dwfa.local_comm[1]  = 0;
  dwfa.local_comm[2]  = 0;
  dwfa.local_comm[3]  = 0;

  dwfa.Ls = 1;
  dwfa.mass = 0.0;
  dwfa.Csw  = 0.0;

  printf("Initialising bfm operator\n");
  printf("drand48 seed = 0\n");
  srand48(0);

  dwf.init(dwfa);

  psi_h = dwf.allocFermion();
  chi_h = dwf.allocFermion();
  check = dwf.allocFermion();
  diff  = dwf.allocFermion();
  
  dwf.randFermion(psi_h);
#ifdef GENERATE
  dwf.dump(psi_h,"src.C","src");
#endif

  //  dwf.importFermion(src,psi_h,0);
  
  dwf.unitGauge();

  printf("cb0dag0 is %lx\n",(unsigned long)cb0dag0);
  printf("cb0dag1 is %lx\n",(unsigned long)cb0dag1);
  printf("cb1dag0 is %lx\n",(unsigned long)cb1dag0);
  printf("cb1dag1 is %lx\n",(unsigned long)cb1dag1);
  // Naive Dslash

  // cb is cb of result, 1-cb is cb of input field

  int idx=0;
  for(cb=0;cb<2;cb++){
    
    /*Import this checkerboard of QDP fields to bagel*/

    // Fill the other checkerboard.
    for(dag=0;dag<2;dag++){

      
      printf("Checking cb=%d dag=%d %lx \n",cb,dag,
	     (unsigned long)arrays[idx]);

      dwf.importFermion(arrays[idx],check,0);

      thr_main(NULL);

#ifdef GENERATE
      dwf.dump(chi_h,files[idx],array_names[idx]);
#else
      printf("Norm of difference is %le\n",delta);
      //printf("Norm result %le\n",n2);
      //printf("Norm check  %le\n",n1);
#endif
      idx++;
    }
  }
  printf("Done\n"); 
}

void * thr_main(void *p)
{

  dwf.dslash(psi_h,
	     chi_h, 
	     cb,dag);
#ifndef GENERATE
  dwf.axpy(diff,check,chi_h,-1.0);
  delta = dwf.norm(diff);
  n1    = dwf.norm(check);
  n2    = dwf.norm(chi_h);
#endif

}





