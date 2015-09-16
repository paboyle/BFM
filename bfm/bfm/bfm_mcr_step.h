void MCR_step(Fermion_t psi,Fermion_t r,Fermion_t p,Fermion_t Ap,Fermion_t tmp, Fermion_t mtmp,double &c,double &d)
{
  double a;
  double b;
  double cp;

  int k = this->iter;

  //c = real( dot(r,Ar) );		//c = rAr
  a = c/d;
       axpy(psi,p,psi,a);		        //x = x + ap
  cp = axpy_norm(r,Ap,r,-a);		//r = r - aAp

    if ( ((k%10 == 0) && (verbose!=0)) || (verbose > 10)){
      if ( isBoss() && !me ) {
	printf("bfmbase::MCR_prec: k= %d r^2= %le %le %lx\n",k,cp,sqrt(cp/ssq),this);
      }
    }

  d = Mprec(r,mp,tmp,DaggerNo); //(r,Ar)
      Mprec(mp,Ar,tmp,DaggerYes);
  //d = real( dot(r,Ar) );
  b = d/c;
  c = d;
  axpy(p,p,r,b);
  d = axpy_norm(Ap,Ap,Ar,b);


}
