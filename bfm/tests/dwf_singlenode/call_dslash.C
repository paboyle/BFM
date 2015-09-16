#include <bfm.h>
#include <stdio.h>
#include <stdlib.h>


#define GENERATE

void * thr_main(void *p);

Fermion_t psi_h;
Fermion_t chi_h;
Fermion_t check;
Fermion_t diff ;

int cb, dag;
double delta;
double n1;
double n2;

#ifdef GENERATE
double cb0dag0[1];
double cb0dag1[1];
double cb1dag0[1];
double cb1dag1[1];
#else
//extern double src[];
extern double cb0dag0[];
extern double cb0dag1[];
extern double cb1dag0[];
extern double cb1dag1[];
#endif

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
          * node_latt[3]
          * Ls;

  int words = 24 * vol;

  int bytes = sizeof(QDPdouble)*words;

  double * psi_p = (QDPdouble *)malloc(bytes);
  bzero(psi_p, bytes);

  //exportFermion(psi_p,psi,0);
  //dump(psi_p, words, file, array_name);
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

void dump_integers(FILE * fp,char * name,integer * array, int words)
{
  fprintf(fp,"integer %s[] = { \n",name);
  for(int w=0;w<words;w++){
    fprintf(fp,"\t%lu,\n",(unsigned long)array[w]);
  }
  fprintf (fp,"0\n};\n");
}
void pointers_dump(void)
{
  FILE *fp;
  int mu;
  int cb, pm;
  int table_size;

  fp = fopen("pointers.C","w");
  fprintf(fp,"#include <bagel_int.h>\n");

  table_size = NMinusPlus*(this->simd_cbvol+1)*Nmu;
  dump_integers(fp,"shift_table_0",this->shift_table[0],table_size);
  dump_integers(fp,"shift_table_1",this->shift_table[1],table_size);

  table_size = (this->simd_cbvol+8);
  dump_integers(fp,"cb_table_0",this->cb_table[0],table_size);
  dump_integers(fp,"cb_table_1",this->cb_table[1],table_size);

  table_size =this->simd_allbound;
  dump_integers(fp,"face_table_0_0_0",this->face_table[0][0][0],table_size);
  dump_integers(fp,"face_table_0_0_1",this->face_table[0][0][1],table_size);
  dump_integers(fp,"face_table_0_0_2",this->face_table[0][0][2],table_size);
  dump_integers(fp,"face_table_0_0_3",this->face_table[0][0][3],table_size);
  dump_integers(fp,"face_table_0_1_0",this->face_table[0][1][0],table_size);
  dump_integers(fp,"face_table_0_1_1",this->face_table[0][1][1],table_size);
  dump_integers(fp,"face_table_0_1_2",this->face_table[0][1][2],table_size);
  dump_integers(fp,"face_table_0_1_3",this->face_table[0][1][3],table_size);
  dump_integers(fp,"face_table_1_0_0",this->face_table[1][0][0],table_size);
  dump_integers(fp,"face_table_1_0_1",this->face_table[1][0][1],table_size);
  dump_integers(fp,"face_table_1_0_2",this->face_table[1][0][2],table_size);
  dump_integers(fp,"face_table_1_0_3",this->face_table[1][0][3],table_size);
  dump_integers(fp,"face_table_1_1_0",this->face_table[1][1][0],table_size);
  dump_integers(fp,"face_table_1_1_1",this->face_table[1][1][1],table_size);
  dump_integers(fp,"face_table_1_1_2",this->face_table[1][1][2],table_size);
  dump_integers(fp,"face_table_1_1_3",this->face_table[1][1][3],table_size);

}

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
  int ls = 4;

  bfmarg dwfa;
  dwfa.solver = DWF;
  dwfa.threads = 1;

  dwfa.node_latt[0]  = lx;
  dwfa.node_latt[1]  = ly;
  dwfa.node_latt[2]  = lz;
  dwfa.node_latt[3]  = lt;

  dwfa.local_comm[0]  = 1;
  dwfa.local_comm[1]  = 1;
  dwfa.local_comm[2]  = 1;
  dwfa.local_comm[3]  = 1;

  dwfa.Ls = ls;
  dwfa.mass = 0.0;
  dwfa.Csw  = 0.0;
  dwfa.precon_5d = 1;

  printf("Initialising bfm operator\n");
  printf("drand48 seed = 0\n");
  srand48(0);

  dwf.init(dwfa);
  dwf.pointers_dump();

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

#ifndef GENERATE
      dwf.importFermion(arrays[idx],check,0);
#endif

      thr_main(NULL);

#ifdef GENERATE
      dwf.dump(chi_h,files[idx],array_names[idx]);
#else
      printf("Norm of difference is %le\n",delta);
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
	     cb,dag,0); // No dperp for now
#ifndef GENERATE
  dwf.axpy(diff,check,chi_h,-1.0);
  delta = dwf.norm(diff);
  n1    = dwf.norm(check);
  n2    = dwf.norm(chi_h);
#endif

}





