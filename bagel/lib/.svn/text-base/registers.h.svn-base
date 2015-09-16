
/*
 *
 *  Copyright Peter Boyle, University of Glasgow, 2000.
 *  This software is provided for NON-COMMERCIAL use only,
 *  It is provided under the GNU pubic License V2
 *  It is provided as is and is not guaranteed fit for any purpose.
 *
 */


#ifndef _REGISTERS_H_
#define _REGISTERS_H_
class register_array_1d {
  int *regs;
 public:
  register_array_1d(char *Var,int type, int d1)
  {
    regs = (int *)arralloc(sizeof(int),1,d1);
    int i;
    char *rstr;
    for(i=0;i<d1;i++) {
      rstr = (char *) malloc(80);
      sprintf(rstr,"%s%d",Var,i);
      regs[i] = allocate_reg(type,rstr);
    }
  };
  int operator[](int index)
  {  
    return(regs[index]);
  };
};

class register_array_2d {
  int **regs;

 public:
  register_array_2d(char *Var, int type,int d1, int d2)
  {
    regs = (int **)arralloc(sizeof(int),2,d1,d2);
    int i,j;
    char *rstr;
    for(i=0;i<d1;i++) {
    for(j=0;j<d2;j++) {
      rstr = (char *)malloc(80);
      sprintf(rstr,"%s%d%d",Var,i,j);
      regs[i][j] = allocate_reg(type,rstr);
    }
    }
  };
  int *operator[](int index)
  {  
    return(regs[index]);
  };
};

class register_array_3d {
  int ***regs;
 public:
  register_array_3d(char *Var, int type,int d1,int d2,int d3)
  {
    regs = (int ***) arralloc(sizeof(int),3,d1,d2,d3);
    int i,j,k;
    char *rstr;
    for(i=0;i<d1;i++) {
    for(j=0;j<d2;j++) {
    for(k=0;k<d3;k++) {
      rstr = (char *)malloc(80);
      sprintf(rstr,"%s%d%d%d",Var,i,j,k);
      regs[i][j][k] = allocate_reg (type,rstr);
    }
    }
    }
  };
  int **operator[](int index)
  {  
    return(regs[index]);
  };
};


class offset_array_1d {
  int *offs;
 public:
  offset_array_1d(char *Var, int d1,int shft=0,enum Datum t=Double)
  {
    offs = (int *)arralloc(sizeof(int),1,d1);
    int i;
    char *rstr;
    for(i=0;i<d1;i++) {
      rstr = (char *)malloc(80);
      sprintf(rstr,"%s%d",Var,i);
      offs[i] = def_offset(i+shft,t,rstr);
    }
  };
  int operator[](int index)
  {  
    return(offs[index]);
  };
};


class offset_array_2d {
  int **offs;
 public:
  offset_array_2d(char *Var, int d1, int d2,int shft=0,Datum type=Double)
  {
    offs = (int **)arralloc(sizeof(int),2,d1,d2);
    int i,j;
    char *rstr;
    for(i=0;i<d1;i++) {
    for(j=0;j<d2;j++) {
      rstr = (char *)malloc(80);
      sprintf(rstr,"%s%d%d",Var,i,j);
      offs[i][j] = def_offset(i*d2+j+shft,type,rstr);
    }
    }
  };
  int *operator[](int index)
  {  
    return(offs[index]);
  };
};

class offset_array_3d {
  int ***offs;
 public:
  offset_array_3d(char *Var, int d1,int d2,int d3,int shft=0,Datum type=Double)
  {
    offs = (int ***) arralloc(sizeof(int),3,d1,d2,d3);
    int i,j,k;
    char *rstr;
    for(i=0;i<d1;i++) {
    for(j=0;j<d2;j++) {
    for(k=0;k<d3;k++) {
      rstr = (char *)malloc(80);
      sprintf(rstr,"%s%d%d%d",Var,i,j,k);
      offs[i][j][k] = def_offset(i*d2*d3+j*d3+k+shft,type,rstr);
    }
    }
    }
  };
  int **operator[](int index)
  {  
    return(offs[index]);
  };
};

#endif
