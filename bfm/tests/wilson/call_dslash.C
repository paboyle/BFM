#include <chroma.h>
 
typedef double Float;

#define Printf if ( QMP_is_primary_node() ) printf
//#define Printf printf

#include <bfm.h>
#include <bfm_qdp_chroma_linop.h>

typedef LatticeFermion T;
typedef LatticeFermion Phi;
typedef multi1d<LatticeColorMatrix> U;
typedef multi1d<LatticeColorMatrix> P;
typedef multi1d<LatticeColorMatrix> Q;
typedef bfm bfm_t;

using namespace Chroma;
void compare_result (LatticeFermion A, LatticeFermion B, bfm_t & dwf);

void      testFermion  (LatticeFermion &psi,  bfm_t   & dwf)
{
  multi1d<int> x(4);

  int Nspinco=12;
  Float *Psi_p = (Float *) &psi.elem(0).elem(0).elem(0).real();

 QDPIO::cout << "testFermion : Checking I get the CHROMA array ordering correct by comparing to peekSite\n";

  for ( x[3]=0; x[3]<dwf.node_latt[3];x[3]++ ) { 
  for ( x[2]=0; x[2]<dwf.node_latt[2];x[2]++ ) { 
  for ( x[1]=0; x[1]<dwf.node_latt[1];x[1]++ ) { 
  for ( x[0]=0; x[0]<dwf.node_latt[0];x[0]++ ) { 
    int xx[4];
    xx[0]=x[0];
    xx[1]=x[1];
    xx[2]=x[2];
    xx[3]=x[3];

    int print = 0;

    for ( int co=0;co<3;co++ ) { 
    for ( int sp=0;sp<4;sp++ ) { 

        int spco = co+3*sp;

	int reim=0;
        int cidx = dwf.chroma_idx(xx,reim,spco,Nspinco);
        Fermion ferm = peekSite(psi,x); 
        ColorVector cv = peekSpin(ferm,sp);
        Complex cp = peekColor(cv,co);
	if ( Psi_p[cidx] != toDouble(real(cp) ) ) {
	  print = 1;			  
        }

	reim=1;
        cidx = dwf.chroma_idx(xx,reim,spco,Nspinco);
        ferm = peekSite(psi,x); 
        cv = peekSpin(ferm,sp);
        cp = peekColor(cv,co);
	if ( Psi_p[cidx] != toDouble(imag(cp) ) ) {

          print = 1;					  

	}

      }}

      if ( print ) 
       QDPIO::cout <<x[0]<<" "<<x[1]<<" "<<x[2]<<" "<<x[3]<< endl;

     for ( int co=0;co<3;co++ ) { 
     for ( int sp=0;sp<4;sp++ ) { 

        int spco = co+3*sp;

	int reim=0;
        int cidx = dwf.chroma_idx(xx,reim,spco,Nspinco);
        Fermion ferm = peekSite(psi,x); 
        ColorVector cv = peekSpin(ferm,sp);
        Complex cp = peekColor(cv,co);
	if ( print ) {
	 QDPIO::cout << "("<<Psi_p[cidx] << "," <<Psi_p[cidx+1] << ") " << cp<<endl;
        }

     }}      

  }}}}

}

int main (int argc,char **argv )
{
  Chroma::initialize(&argc,&argv);

  if ( argc != 5 ) { 
   Printf("Usage: %s lx ly lz lt \n All must be even\n",argv[0]);
   Printf("argc is %d\n",argc);
    for ( int i=0;i<argc;i++)
      Printf("%d %s\n",i,argv[i]);
    exit(-1);
  }


  /********************************************************
   * Setup QDP
   ********************************************************
   */
  multi1d<int> nrow(Nd);
  nrow[0] = atoi(argv[1]);
  nrow[1] = atoi(argv[2]);
  nrow[2] = atoi(argv[3]);
  nrow[3] = atoi(argv[4]);

  Layout::setLattSize(nrow);
  Layout::create();

  int lx = QDP::Layout::subgridLattSize()[0];
  int ly = QDP::Layout::subgridLattSize()[1];
  int lz = QDP::Layout::subgridLattSize()[2];
  int lt = QDP::Layout::subgridLattSize()[3];

  /********************************************************
   * Setup DWF operator
   ********************************************************
   */
  bfmarg dwfa;
  bfm_t  dwf;
  dwfa.solver = WilsonFermion;
  dwfa.node_latt[0]  = lx;
  dwfa.node_latt[1]  = ly;
  dwfa.node_latt[2]  = lz;
  dwfa.node_latt[3]  = lt;

  multi1d<int> procs = QDP::Layout::logicalSize();
  Printf("%d dim machine\n\t", procs.size());
  for(int mu=0;mu<4;mu++){
    if ( procs[mu]>1 ) {
      dwfa.local_comm[mu] = 0;
    } else { 
      dwfa.local_comm[mu] = 1;
    }
  }

  // Make up a random gauge field.
  Complex      cc=Real(1.0);
  ColorMatrix cm = zero;
  pokeColor(cm,cc,0,0);
  pokeColor(cm,cc*10,1,1);
  pokeColor(cm,cc*100,2,2);

  multi1d<LatticeColorMatrix> u(Nd);
  HotSt(u);
  /*
    u[0] = cm;
    u[1] = cm*2.0;
    u[2] = cm*3.0;
    u[3] = cm*4.0;
  */
  
  //  u[0]=zero;
  //  u[1]=zero;
  //  u[2] = zero;
  //  u[3] = zero;

  //  for(int mu=0;mu<4;mu++){
  //u[mu]=where(Layout::latticeCoordinate(0)==7,LatticeColorMatrix(zero),u[mu]);
  //u[mu]=where(Layout::latticeCoordinate(0)==3,LatticeColorMatrix(zero),u[mu]);
  //  }

  // Make up a gaussian source and a zero result vector
  LatticeFermion psi; 
#if 1
  gaussian(psi);

#else
  psi=zero;
  ColorVector cv=zero; pokeColor(cv,cc,0);
  Fermion     fv=zero; pokeSpin(fv,cv,0);

  multi1d<int> mysite(4);
  mysite[0]=3;           
  mysite[1]=3;
  mysite[2]=3;
  mysite[3]=0;

  // Works: 1,0,0,0   3,0,0,0, 3,3,3,0
  // Broke: 0,0,0,1
  pokeSite(psi,fv,mysite);
#endif  

  LatticeFermion chi;
  chi = zero;

  dwfa.Ls = 1;
  dwfa.mass = 0.1;
  dwfa.Csw = 0.0;

  Printf("Initialising bfm operator\n");
  dwf.init(dwfa);
  
  Printf("bfm operator using vector length %d\n",dwf.simd());
  Fermion_t psi_h = dwf.allocFermion();
  Fermion_t chi_h = dwf.allocFermion();

  dwf.importGauge(u);

  //testFermion  (psi,dwf);

  // Naive Dslash
#ifdef CHROMA_3
  multi1d<int> bcs(Nd); bcs[0] = bcs[1] = bcs[2] = bcs[3] = 1;
  Handle<FermBC<T,U,U> > fbc(new SimpleFermBC< T, U, U >(bcs));
  Handle<CreateFermState<T,U,U> > cfs( new CreateSimpleFermState<T,U,U>(fbc));
  Handle<FermState<T,U,U> > fs ((*cfs)(u));
  QDPWilsonDslash D(fs);
#else
  QDPWilsonDslash D(u);
#endif
  LatticeFermion chi_qdp;
  LatticeFermion noise;
  gaussian(noise);

  // cb is cb of result, 1-cb is cb of input field
  for(int cb=0;cb<2;cb++){

    /*Import this checkerboard of QDP fields to bagel*/
    Printf("Importing psi field cb %d\n",cb);

    // Fill the other checkerboard.
    chi = zero;
    chi_qdp = zero;

    Float *psi_p = (Float *) &psi.elem(0).elem(0).elem(0).real();
    Float *chi_p = (Float *) &chi.elem(0).elem(0).elem(0).real();
    dwf.importFermion(psi_p,psi_h,1-cb);

    double tn = dwf.norm(psi_h);
    printf("Norm of bagel input vector %le\n",tn);

    for(int dag=0;dag<2;dag++){

      Printf("Checking cb=%d dag=%d\n",cb,dag); 

      dwf.dslash(psi_h,
                 chi_h, 
		 cb,dag);

      dwf.exportFermion(chi_p,chi_h,cb);

      // Check the result
      PlusMinus pm;
      if ( dag ) pm = MINUS;
      else pm = PLUS;

#define DEBUG
#ifdef DEBUG1
        LatticeFermion tmp;

        {
         QDPIO::cout << "U[mu]"<<endl;
	  int c1,c2,mu;
          for(mu=0;mu<4;mu++) {
	  for(c1=0;c1<3;c1++) { for(c2=0;c2<3;c2++) {
	   QDPIO::cout << peekColor(u[mu].elem(0),c1,c2) ;
            }
           QDPIO::cout << endl;
          }
           QDPIO::cout << endl;
	  }
         QDPIO::cout << "Udag[mu]"<<endl;
          for(mu=0;mu<4;mu++) {
          LatticeColorMatrix udag = adj(shift(u[mu],BACKWARD,mu));
	  for(c1=0;c1<3;c1++) { for(c2=0;c2<3;c2++) {
	   QDPIO::cout << peekColor(udag.elem(0),c1,c2) ;
            }
           QDPIO::cout << endl;
          }
           QDPIO::cout << endl;
	  }
	}

	LatticeFermion fm;
	fm = zero;
        tmp = spinReconstructDir0Minus(u[0] * shift(spinProjectDir0Minus(psi), FORWARD, 0));
       QDPIO::cout << "Uchip0"<<endl;
       QDPIO::cout << tmp.elem(0) << endl;
	fm.elem(0) = tmp.elem(0);
        tmp = spinReconstructDir0Plus(shift(adj(u[0]) * spinProjectDir0Plus(psi), BACKWARD, 0));
       QDPIO::cout << "Uchim0"<<endl;
	cout << tmp.elem(0) << endl;
	fm.elem(0) = fm.elem(0) +  tmp.elem(0);
	cout << "Recon running sum "<<fm.elem(0) <<endl;

        tmp = spinReconstructDir1Minus(u[1] * shift(spinProjectDir1Minus(psi), FORWARD, 1));
       QDPIO::cout << "Uchip1"<<endl;
       QDPIO::cout << tmp.elem(0) << endl;
	fm.elem(0) = fm.elem(0) +  tmp.elem(0);
	
        tmp = spinReconstructDir1Plus(shift(adj(u[1]) * spinProjectDir1Plus(psi), BACKWARD, 1));
       QDPIO::cout << "Uchim1"<<endl;
	cout << tmp.elem(0) << endl;
	fm.elem(0) = fm.elem(0) +  tmp.elem(0);
	cout << "Recon running sum "<<fm.elem(0) <<endl;
    
        tmp = spinReconstructDir2Minus(u[2] * shift(spinProjectDir2Minus(psi), FORWARD, 2));
       QDPIO::cout << "Uchip2"<<endl;
       QDPIO::cout << tmp.elem(0) << endl;
	fm.elem(0) = fm.elem(0) +  tmp.elem(0);
	
        tmp = spinReconstructDir2Plus(shift(adj(u[2]) * spinProjectDir2Plus(psi), BACKWARD, 2));
       QDPIO::cout << "Uchim2"<<endl;
	cout << tmp.elem(0) << endl;
	fm.elem(0) = fm.elem(0) +  tmp.elem(0);
	cout << "Recon running sum "<<fm.elem(0) <<endl;
    

        tmp = spinReconstructDir3Minus(u[3] * shift(spinProjectDir3Minus(psi), FORWARD, 3));
       QDPIO::cout << "Uchip3"<<endl;
       QDPIO::cout << tmp.elem(0) << endl;
	fm.elem(0) = fm.elem(0) +  tmp.elem(0);
	
        tmp = spinReconstructDir3Plus(shift(adj(u[3]) * spinProjectDir3Plus(psi), BACKWARD, 3));
       QDPIO::cout << "Uchim3"<<endl;
	cout << tmp.elem(0) << endl;
	fm.elem(0) = fm.elem(0) +  tmp.elem(0);
	cout << "Recon running sum "<<fm.elem(0) <<endl;
    
 
        for(int mu=0;mu<4;mu++) {
          tmp = shift(psi, FORWARD, mu);
  	 QDPIO::cout << "psi x+"<<mu<< endl;
  	 QDPIO::cout << tmp.elem(0) << endl;
        }
        for(int mu=0;mu<4;mu++) {
          tmp = shift(psi, BACKWARD, mu);
  	 QDPIO::cout << "psi x-"<<mu<< endl;
  	 QDPIO::cout << tmp.elem(0) << endl;
        }

#endif

      D.apply(chi_qdp, psi, pm, cb);

#ifdef DEBUG
     QDPIO::cout << "QDP gets"<< endl;
      Float *chi_q = (Float *) &chi_qdp.elem(0).elem(0).elem(0).real();
      FILE * fp = fopen("foo","w");
      for(int s=0 ; s < lx*ly*lz*lt/2 ; s++ ) { 
      for(int i=0 ; i < 24 ; i++ ) { 
	fprintf(fp,"%le %le\n",chi_p[i+24*s],chi_q[i+s*24]);
      }
      fprintf(fp,"\n");
      }
     QDPIO::cout << chi_qdp.elem(0) << endl;
#endif

    double bn = dwf.norm(psi_h);
    double qn = toDouble(norm2(psi));

    cout << "bn,qn "<<bn<< " "<<qn<<endl;

      Double n2 = norm2(chi-chi_qdp);
      
     QDPIO::cout << "|| Bagel - QDP || = "<< n2 << endl;

     n2 = norm2(chi);
     QDPIO::cout << "|| Bagel || = "<< n2 << endl;
     n2 = norm2(chi_qdp);
     QDPIO::cout << "|| QDP || = "<< n2 << endl;
#ifdef DEBUG
  compare_result(chi,chi_qdp,dwf);
  exit(0);
#endif
    }
  }


  Printf("Done\n"); 

  multi1d<int> site(4);
  LatticeFermion Gauss; 
  LatticeFermion X; 
  LatticeFermion Y;
  Fermion frm;
  
  gaussian(X);
  X=0.1*X;
  gaussian(Y);

  multi1d<LatticeColorMatrix> F(4);
  multi1d<LatticeColorMatrix> Fb(4);

  {
    void * force= dwf.allocMatrix();
    
    for( int dag=0;dag<2;dag++){
      
      if (dag == 0 ) 
	D.deriv(F,X,Y,PLUS);
      else
	D.deriv(F,X,Y,MINUS);
      
      for(int cb=0;cb<2;cb++){
	dwf.importFermion(X,psi_h,cb);
	dwf.importFermion(Y,chi_h,1-cb);
	dwf.zeroMatrix(force);
	dwf.MeoDeriv(psi_h,
		     chi_h, 
		     force,
		     cb,dag);
	dwf.exportForce(force,Fb,cb);
	
      }

      for(int mu=0;mu<4;mu++){
	Real h(-0.5);
	QDPIO::cout << "dag "<< dag<<"mu "<<mu<<" n2diff " << norm2(Fb[mu]-h*F[mu])<<endl;
	QDPIO::cout << "dag "<< dag<<"mu "<<mu<<" n2diff " << norm2(Fb[mu])<<endl;
	QDPIO::cout << "dag "<< dag<<"mu "<<mu<<" n2diff " << norm2(h*F[mu])<<endl;
      }
    }
    dwf.freeMatrix(force);
  }
  EvenOddPrecWilsonFermAct FA(cfs,dwfa.mass);
  Handle< DiffLinearOperator<Phi,P,Q> > M(FA.linOp(fs));


  {
    void * force[2];
    force[0] = dwf.allocMatrix();
    force[1] = dwf.allocMatrix();

    for( int dag=0;dag<2;dag++){
      
      Printf("Calling mprec deriv\n");
      int cb=1; // Odd Odd component here
      dwf.zeroMatrix(force[0]);
      dwf.zeroMatrix(force[1]);

      dwf.importFermion(X,psi_h,cb);
      dwf.importFermion(Y,chi_h,cb);
      dwf.MprecDeriv(psi_h,chi_h,force,dag);
      Fb=zero;
      for(int cb=0;cb<2;cb++)
	dwf.exportForce(force[cb],Fb,cb);
    
      Printf("Calling S_f deriv\n");
      if( dag==1 ) 
	M->deriv(F, X, Y, MINUS);
      else
	M->deriv(F, X, Y, PLUS);
      
      for(int mu=0;mu<4;mu++){
	QDPIO::cout << "dag "<< dag<<"mu "<<mu<<" n2diff " << norm2(Fb[mu]-F[mu])<<endl;
	QDPIO::cout << "dag "<< dag<<"mu "<<mu<<" n2diff " << norm2(Fb[mu])<<endl;
	QDPIO::cout << "dag "<< dag<<"mu "<<mu<<" n2diff " << norm2(F[mu])<<endl;
      }
    }
    dwf.freeMatrix(force[0]);
    dwf.freeMatrix(force[1]);


  }

  exit(0);
}


void compare_result (LatticeFermion A, LatticeFermion B, bfm_t & dwf)
{

  multi1d<int> x(4);

  int Nspinco=12;

  Printf("Result vectors \n");
  for ( x[3]=0; x[3]<dwf.node_latt[3];x[3]++ ) { 
  for ( x[2]=0; x[2]<dwf.node_latt[2];x[2]++ ) { 
  for ( x[1]=0; x[1]<dwf.node_latt[1];x[1]++ ) { 
  for ( x[0]=0; x[0]<dwf.node_latt[0];x[0]++ ) { 

    int xx[4];
    xx[0]=x[0];
    xx[1]=x[1];
    xx[2]=x[2];
    xx[3]=x[3];
    Printf("site %d %d %d %d\n",x[0],x[1],x[2],x[3]);
    double nna=0;
    double nnb=0;
    for ( int sp=0;sp<4;sp++ ) { 
      for ( int co=0;co<3;co++ ) { 

        int spco = co+3*sp;

	int reim=0;
        Fermion ferm = peekSite(A,x); 
        ColorVector cv = peekSpin(ferm,sp);
        Complex ca = peekColor(cv,co);

        ferm = peekSite(B,x); 
        cv = peekSpin(ferm,sp);
        Complex cb = peekColor(cv,co);
	nna=nna+toDouble(real(ca*conj(ca)));
	nnb=nnb+toDouble(real(cb*conj(cb)));
	if(abs(toDouble(real(ca))-(toDouble(real(cb)))) > 1.0e-9 ){
	  Printf("s %d c %d r %le %le\n",sp,co,toDouble(real(ca)),toDouble(real(cb)));
	}
	if(abs(toDouble(imag(ca))-(toDouble(imag(cb)))) > 1.0e-9 ){
	  Printf("s %d c %d i %le %le\n",sp,co,toDouble(imag(ca)),toDouble(imag(cb)));
	}

    }}
    if(nna+nnb>0){
      Printf("\t\t nna %le / nnb %le\n",nna,nnb);
    }

  }}}}


}






