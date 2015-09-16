#ifndef __bfm_hdcg_wrapper_h__
#define __bfm_hdcg_wrapper_h__
#include <BfmHDCG.h>
using namespace std;

// Heirarchically deflated conjugate gradient
class HDCG_wrapper {
  
 public:

    typedef LatticeFermion T4;
    typedef LatticeColorMatrix U;
    typedef multi1d<LatticeColorMatrix> Q;
    typedef float rFloat;
    BfmHDCGParams Parms;
    BfmHDCG<double> *ldop_d;

    bfm_qdp<rFloat> dop_r;
    bfm_qdp<double> dop_d;
    bfm_qdp<float>  dop_f;

  void HDCG_subspace_rotate(multi1d<Real> &phases_by_pi,multi1d<int> &phases_dir, int sgn)
  {

    /* TwistedBC does the following
    for(int i=0; i < Nd-1; i++) {
       int direction = phases_dir[i]; 
       Real arg = phases_by_pi[i]*onepi / Real( nrow[ direction ] );
       Real re = cos(arg);
       Real im = sin(arg);
       Complex phase = cmplx( re, im );
       u[ direction ] *= phase;
    }
     * Suppose 
     *    D psi = 0
     * Want to modify subspace vector psi so that
     *    D^tw psi^tw = 0
     * If 
     *    psi^tw(x) = e^{-i x arg} psi(x)
     * Then 
     *    U_\mu^tw psi^tw(x+mu) = e^{i x arg} U_\mu psi(x+mu) 
    */
    dop_d.BossLog("HDCG Rotating subspace for twisted BCs\n");

    LatticeReal theta = zero;
    for(int i=0; i < Nd-1; i++) {
      Real onepi = M_PI;
      int mu = phases_dir[i];
      theta = theta + Layout::latticeCoordinate(mu) * phases_by_pi[i]* onepi / Layout::lattSize()[mu] ;
    }

    LatticeComplex eitheta = cmplx(cos(theta),-sgn*sin(theta));

    int Ls=dop_d.Ls;
    multi1d<LatticeFermion> copy(Ls);
    for(int v=0;v<ldop_d->Nvec;v++){
      dop_d.exportFermion(copy,ldop_d->subspace_d[v],1);
      for(int s=0;s<Ls;s++){
	copy[s]=eitheta*copy[s];
      }
      dop_d.importFermion(copy,ldop_d->subspace_d[v],1);
    }
  }

  void HDCG_set_mass(double mass)
  {
    dop_d.BossLog("HDCG setting mass %le\n",mass);
    dop_d.GeneralisedFiveDimEnd();
    dop_d.mass = mass;
    dop_d.GeneralisedFiveDimInit();

    dop_f.GeneralisedFiveDimEnd();
    dop_f.mass = mass;
    dop_f.GeneralisedFiveDimInit();

    dop_r.GeneralisedFiveDimEnd();
    dop_r.mass = mass;
    dop_r.GeneralisedFiveDimInit();

  }

  void HDCG_init(BfmHDCGParams & Parms_,
		 bfmActionParams &BAP) {

    Parms = Parms_;

    // Complicated init sequence
    bfmarg bfma;
    bfmActionParams *bfmap = (bfmActionParams *) &bfma;

    // Physica params
    *bfmap = BAP;
    
    // Algorithm & code control
    bfma.Threads(omp_get_max_threads());
    bfma.Verbose(0);
    bfma.time_report_iter=-100;
    bfma.max_iter     = 10000;
    bfma.residual     = 1.0e-8;
  
    //Geometry
    bfma.node_latt[0] =  QDP::Layout::subgridLattSize()[0];
    bfma.node_latt[1] =  QDP::Layout::subgridLattSize()[1];
    bfma.node_latt[2] =  QDP::Layout::subgridLattSize()[2];
    bfma.node_latt[3] =  QDP::Layout::subgridLattSize()[3];
    
    multi1d<int> procs = QDP::Layout::logicalSize();
    for(int mu=0;mu<4;mu++){
      if (procs[mu]>1) bfma.local_comm[mu] = 0;
      else             bfma.local_comm[mu] = 1;
    }


    dop_f.init(bfma);
    dop_f.BossLog("HDCG inititialised single prec linop\n");
    dop_d.init(bfma);
    dop_d.BossLog("HDCG inititialised double prec linop\n");
    bfma.Ls   = Parms.SubspaceRationalLs;
    bfma.mass = Parms.SubspaceRationalMass;
    dop_r.BossLog("HDCG inititialising subspace generation Ls=%d mass %le\n",bfma.Ls,bfma.mass);
    dop_r.init(bfma);
    dop_r.BossLog("HDCG inititialised subspace generation\n");

    // Merge the Parms into constructor.
    ldop_d = new BfmHDCG<double>(Parms.Ls,Parms.NumberSubspace,Parms.Block,Parms.Block,&dop_d,&dop_f);
    ldop_d->SetParams(Parms);

  }

  void HDCG_end(void){
    dop_d.end();
    dop_f.end();
    dop_r.end();
    ldop_d->end();
    HDCG_subspace_free();
    delete ldop_d;
  }


  template<typename FloatEXT>
  void HDCG_gauge_import_cps(FloatEXT *gauge){
    dop_d.cps_importGauge<FloatEXT>(gauge);
    dop_f.cps_importGauge<FloatEXT>(gauge);
    dop_r.cps_importGauge<FloatEXT>(gauge);
  }
  
  void HDCG_gauge_import(multi1d<LatticeColorMatrix> &gauge){
    dop_d.importGauge(gauge);
    dop_f.importGauge(gauge);
    dop_r.importGauge(gauge);
  }
  void HDCG_subspace_init(){
    //    int sloppy_relax=Parms.SubspaceRationalRefine;
    int sloppy_relax=1;

    uint32_t seed = ldop_d->ncoor[0]*17
                  + ldop_d->ncoor[1]*251;
                  + ldop_d->ncoor[2]*1023;
                  + ldop_d->ncoor[3]*8191;
    srand48(seed);

    if ( sloppy_relax ) { 
      ldop_d->RelaxSubspace<rFloat>(&dop_r);
    } else {
      ldop_d->RelaxSubspace<double>(&dop_d);
    }
    ldop_d->SinglePrecSubspace();
  }
  void HDCG_subspace_refine(){

    int Nvec = Parms.NumberSubspace;

    Fermion_t src = dop_d.allocFermion();
    Fermion_t sol = dop_d.allocFermion();
    Fermion_t tmp = dop_d.allocFermion();

    int threads = dop_d.threads;
    BfmHDCG<float>  *ldop_f;
    ldop_f = new BfmHDCG<float>(Parms.Ls,Parms.NumberSubspace,Parms.Block,Parms.Block,&dop_d,&dop_f);
    ldop_f->CloneSubspace<double>(*ldop_d);
    ldop_f->SetParams(Parms);
    //    ldop_f->PcgSingleShift = 0;
    ldop_f->PcgSingleShift = Parms.SubspaceRationalRefineLo;
    dop_d.residual = Parms.SubspaceRationalRefineResidual;

    ldop_f->PcgType        = PcgAdef2fSingleShift;

#pragma omp parallel 
    {
#pragma omp for 
      for(int i=0;i<threads;i++) {
	for(int v=0;v<Nvec;v++){

	  double nn=dop_d.norm(ldop_d->subspace_r[v]);
	  dop_d.ThreadBossLog("HDCG Refining vector[%d] of norm %le\n",v,nn);

	  dop_d.copy(src,ldop_d->subspace_r[v]);
	  dop_d.set_zero(sol);
	  ldop_f->Pcg(ldop_d,sol,src,tmp);
	  dop_d.copy(ldop_d->subspace_r[v],sol);

	  nn=dop_d.norm(ldop_d->subspace_r[v]);
	  //	  dop_d.ThreadBossLog(" Refined vector[%d] new norm %le\n",v,nn);
	}
	for(int v=0;v<Nvec;v++){
	  dop_d.copy(ldop_d->subspace_d[v],ldop_d->subspace_r[v]);
	}
      }
    }
    ldop_f->end();
    delete ldop_f;

    dop_d.freeFermion(tmp);
    dop_d.freeFermion(src);
    dop_d.freeFermion(sol);
    
    ldop_d->FreeSingleSubspace();
    ldop_d->OrthogonaliseSubspace();
    ldop_d->SinglePrecSubspace();

  }
  void HDCG_subspace_free(void){
    ldop_d->FreeDoubleSubspace();
    ldop_d->FreeSingleSubspace();
  }
  void HDCG_subspace_compute(int sloppy){

    int LdopDeflVecs = Parms.LdopDeflVecs;
    int threads = dop_d.threads;
    int Nvec = Parms.NumberSubspace;

    if ( sloppy ) { 
      //ldop_d->ComputeLittleMatrixColored<float>(&dop_f,ldop_d->subspace_f);
      //      ldop_d->ComputeLittleMatrixColoredOneHopHermitian<float>(&dop_f,ldop_d->subspace_f);
      ldop_d->ComputeLittleMatrixColoredOneHopHermitian<float>(&dop_r,ldop_d->subspace_ro);
      ldop_d->LdopDeflationBasisSize=0;
      if ( LdopDeflVecs )
	ldop_d->LdopDeflationBasisTrivial();

    } else { 
      //ldop_d->ComputeLittleMatrixColored<float>(&dop_f,ldop_d->subspace_f);
      //      ldop_d->ComputeLittleMatrixColored<double>(&dop_d,ldop_d->subspace_d);
      ldop_d->ComputeLittleMatrixColoredOneHopHermitian<float>(&dop_f,ldop_d->subspace_f);
      //ldop_d->ComputeLittleMatrixColoredOneHopHermitian<float>(&dop_f,ldop_d->subspace_f);

      // NB due to BG/Q QPX have hard coded DeflationBasisSize must be multiple of 4.
      // This may become even more on future machines.
#if 1
      if( Parms.SubspaceSurfaceDepth < Parms.Ls/2) {
	dop_d.BossMessage("Zeroing s=[%d,%d]\n",Parms.SubspaceSurfaceDepth,Parms.Ls-Parms.SubspaceSurfaceDepth-1);
#pragma omp parallel 
	{
#pragma omp for 
	  for(int i=0;i<threads;i++) {
	    for(int v=0;v<Nvec;v++){
	      double nno = dop_d.norm(ldop_d->subspace_d[v]);
	      for(int s=Parms.SubspaceSurfaceDepth; s<Parms.Ls-Parms.SubspaceSurfaceDepth;s++){
		dop_d.axpby_ssp(ldop_d->subspace_d[v],0.0,ldop_d->subspace_d[v],0.0,ldop_d->subspace_d[v],s,s);
	      }
	      double nnr = dop_d.norm(ldop_d->subspace_d[v]);
	      dop_d.ThreadBossLog("Restricted vec[%d] %le -> %le \n",v,nno,nnr);
	    }
	  }
	}
      }
#endif

      uint32_t seed = ldop_d->ncoor[0]*11
                  + ldop_d->ncoor[1]*161;
                  + ldop_d->ncoor[2]*357;
                  + ldop_d->ncoor[3]*1203;
      srand48(seed);

      int LdopDeflVecs = Parms.LdopDeflVecs;

      ldop_d->LdopDeflationBasisSize=0;
      if ( LdopDeflVecs > 0 ) {
	//      ldop_d->LdopDeflationBasisTrivial();
	int min =16;
	if ( ldop_d->LdopDeflationBasisSize> min )  min = ldop_d->LdopDeflationBasisSize;
	for (int v = min ; v<=LdopDeflVecs; v+=4){
	  ldop_d->LdopDeflationBasisInit(v);
	}
	ldop_d->LdopDeflationBasisDiagonalise(LdopDeflVecs);
      }
    }
  }

  void HDCG_evl_invert(Fermion_t solution, Fermion_t source, double residual, int maxit){

    dop_d.residual = residual;
    dop_d.max_iter = maxit;

    Fermion_t tmp = dop_d.allocFermion();

    BfmHDCG<float>  *ldop_f;
    ldop_f = new BfmHDCG<float>(Parms.Ls,Parms.NumberSubspace,Parms.Block,Parms.Block,&dop_d,&dop_f);
    ldop_f->CloneSubspace<double>(*ldop_d);
    ldop_f->SetParams(Parms);
    ldop_f->PcgType=PcgAdef2f;
    dop_f.BossLog("LittleDopSolvers %d %d \n",ldop_f->LittleDopSolver,ldop_d->LittleDopSolver);
    if ( ldop_d->LdopDeflationBasisSize==0 ) {
      ldop_f->LittleDopSolver = LittleDopSolverCG;
      ldop_d->LittleDopSolver = LittleDopSolverCG;
    } else { 
      ldop_f->LittleDopSolver = LittleDopSolverDeflCG;
      ldop_d->LittleDopSolver = LittleDopSolverDeflCG;
    }

    int threads = dop_d.threads;
#pragma omp parallel
    {
#pragma omp for
      for(int i=0;i<threads;i++) {

	int me = dop_d.thread_barrier();

	// checkboardsource
	ldop_f->AnalyseSpectrum=0;
	if ( Parms.Flexible == 1 ) { 
	  dop_f.ThreadBossLog("Solving with single precision Ldop, fPcg\n");
	  ldop_f->fPcg(ldop_d,solution,source,tmp);
	} else if ( Parms.Flexible == 2 ) {  
	  dop_d.ThreadBossLog("Solving with double precision Ldop, fPcg\n");
	  ldop_d->fPcg(ldop_d,solution,source,tmp);
	} else { 
	  dop_f.ThreadBossLog("Solving with single precision Ldop, Pcg\n");
	  ldop_f->Pcg(ldop_d,solution,source,tmp);
	}

      }
    }

    dop_d.freeFermion(tmp);

    ldop_f->end();
    delete ldop_f;

  }
  void HDCG_invert(Fermion_t solution[2], Fermion_t source[2], double residual, int maxit){

    dop_d.residual = residual;
    dop_d.max_iter = maxit;

    Fermion_t src = dop_d.allocFermion();
    Fermion_t tmp = dop_d.allocFermion();
    Fermion_t Mtmp= dop_d.allocFermion();
    Fermion_t resid= dop_d.allocFermion();

    BfmHDCG<float>  *ldop_f;
    ldop_f = new BfmHDCG<float>(Parms.Ls,Parms.NumberSubspace,Parms.Block,Parms.Block,&dop_d,&dop_f);
    ldop_f->CloneSubspace<double>(*ldop_d);
    ldop_f->SetParams(Parms);
    ldop_f->PcgType=PcgAdef2f;

    int threads = dop_d.threads;
#pragma omp parallel
    {
#pragma omp for
      for(int i=0;i<threads;i++) {

	int me = dop_d.thread_barrier();

	// checkboardsource
	dop_d.MooeeInv(source[Even],tmp,DaggerNo,Even);
	dop_d.Meo     (tmp,src,Odd,DaggerNo);
	dop_d.axpy    (tmp,src,source[Odd],-1.0);
	dop_d.Mprec(tmp,src,Mtmp,DaggerYes);

	ldop_f->AnalyseSpectrum=0;
	if ( Parms.Flexible == 1 ) { 
	  dop_f.ThreadBossLog("Solving with single precision Ldop, fPcg\n");
	  ldop_f->fPcg(ldop_d,solution[Odd],src,tmp);
	} else if ( Parms.Flexible == 2 ) {  
	  dop_d.ThreadBossLog("Solving with double precision Ldop, fPcg\n");
	  ldop_d->fPcg(ldop_d,solution[Odd],src,tmp);
	} else { 
	  dop_f.ThreadBossLog("Solving with single precision Ldop, Pcg\n");
	  ldop_f->Pcg(ldop_d,solution[Odd],src,tmp);
	}

	// sol_e = M_ee^-1 * ( src_e - Meo sol_o )...
	dop_d.Meo(solution[Odd],tmp,Even,DaggerNo);
	dop_d.axpy(src,tmp,source[Even],-1.0);
	dop_d.MooeeInv(src,solution[Even],DaggerNo,Even);

	// Even defect
	dop_d.Meo(solution[Odd],tmp,Even,DaggerNo);
	dop_d.Mooee(solution[Even],Mtmp,DaggerNo,Even);
	dop_d.axpy(Mtmp,Mtmp,tmp,1.0);
	dop_d.axpy(Mtmp,Mtmp,source[Even],-1.0);

	double rf=dop_d.norm(Mtmp);
	// Odd defect
	dop_d.Meo(solution[Even],tmp,Odd,DaggerNo);
	dop_d.Mooee(solution[Odd],Mtmp,DaggerNo,Odd);
	dop_d.axpy(Mtmp,Mtmp,tmp,1.0);
	dop_d.axpy(Mtmp,Mtmp,source[Odd],-1.0);
	rf+=dop_d.norm(Mtmp);

	double f = dop_d.norm(source[Odd]);
	f       += dop_d.norm(source[Even]);

	if (dop_d.isBoss() && !me ) printf("bfm_hdcg_wrapper: unprec sol true residual is %le\n",sqrt(rf/f));

      }
    }

    dop_d.freeFermion(src);
    dop_d.freeFermion(tmp);
    dop_d.freeFermion(Mtmp);
    dop_d.freeFermion(resid);

    ldop_f->end();
    delete ldop_f;

  }    


  void HDCG_invert(multi1d<T4> &qdp_sol,multi1d<T4> &qdp_src,double residual, int maxit){

    dop_d.residual = residual;
    dop_d.max_iter = maxit;

    Fermion_t solution[2];
    Fermion_t source[2];
    for(int cb=0;cb<2;cb++) {
      solution[cb] = dop_d.allocFermion();
      source[cb] = dop_d.allocFermion();
      dop_d.importFermion(qdp_sol,solution[cb],cb);
      dop_d.importFermion(qdp_src,source[cb],cb);
    }

    Fermion_t src = dop_d.allocFermion(); 
    Fermion_t tmp = dop_d.allocFermion(); 
    Fermion_t Mtmp= dop_d.allocFermion(); 
    Fermion_t resid= dop_d.allocFermion(); 

    BfmHDCG<float>  *ldop_f;
    ldop_f = new BfmHDCG<float>(Parms.Ls,Parms.NumberSubspace,Parms.Block,Parms.Block,&dop_d,&dop_f); 
    ldop_f->CloneSubspace<double>(*ldop_d);
    ldop_f->SetParams(Parms);
    dop_d.BossLog("Cloned subspace for single prec ldop : LdopM1iter =%d \n",ldop_f->LdopM1iter);
    dop_d.BossLog("                    double prec ldop : LdopM1iter =%d \n",ldop_d->LdopM1iter);

    int threads = dop_d.threads;

#pragma omp parallel 
  {
#pragma omp for 
    for(int i=0;i<threads;i++) {

      int me = dop_d.thread_barrier();

      // checkboardsource
      dop_d.MooeeInv(source[Even],tmp,DaggerNo,Even);
      dop_d.Meo     (tmp,src,Odd,DaggerNo);
      dop_d.axpy    (tmp,src,source[Odd],-1.0);
      dop_d.Mprec(tmp,src,Mtmp,DaggerYes);  
      dop_d.ThreadBossLog("Solving for flexible = %d \n",Parms.Flexible);
      ldop_f->AnalyseSpectrum=0;
      if ( Parms.Flexible == 1 ) { 
	dop_f.ThreadBossLog("Solving with single precision Ldop, fPcg\n");
	ldop_f->fPcg(ldop_d,solution[Odd],src,tmp);
      } else if ( Parms.Flexible == 2 ) {  
	dop_d.ThreadBossLog("Solving with double precision Ldop, fPcg\n");
	ldop_d->fPcg(ldop_d,solution[Odd],src,tmp);
      } else { 
	dop_f.ThreadBossLog("Solving with single precision Ldop, Pcg\n");
	ldop_f->Pcg(ldop_d,solution[Odd],src,tmp);
      }

      // sol_e = M_ee^-1 * ( src_e - Meo sol_o )...
      dop_d.Meo(solution[Odd],tmp,Even,DaggerNo);
      dop_d.axpy(src,tmp,source[Even],-1.0);
      dop_d.MooeeInv(src,solution[Even],DaggerNo,Even);

      // Even defect
      dop_d.Meo(solution[Odd],tmp,Even,DaggerNo);
      dop_d.Mooee(solution[Even],Mtmp,DaggerNo,Even);
      dop_d.axpy(Mtmp,Mtmp,tmp,1.0);
      dop_d.axpy(Mtmp,Mtmp,source[Even],-1.0);

      double rf=dop_d.norm(Mtmp);
      // Odd defect
      dop_d.Meo(solution[Even],tmp,Odd,DaggerNo);
      dop_d.Mooee(solution[Odd],Mtmp,DaggerNo,Odd);
      dop_d.axpy(Mtmp,Mtmp,tmp,1.0);
      dop_d.axpy(Mtmp,Mtmp,source[Odd],-1.0);
      rf+=dop_d.norm(Mtmp);

      double f = dop_d.norm(source[Odd]);
      f       += dop_d.norm(source[Even]);
      
      dop_d.ThreadBossLog("HDCG unprec sol true residual is %le\n",sqrt(rf/f));

    }
  }


  for(int cb=0;cb<2;cb++) {
    dop_d.exportFermion(qdp_sol,solution[cb],cb);
  }

  dop_d.freeFermion(solution[0]);
  dop_d.freeFermion(solution[1]);
  dop_d.freeFermion(source[0]);
  dop_d.freeFermion(source[1]);
  dop_d.freeFermion(src);
  dop_d.freeFermion(tmp);
  dop_d.freeFermion(Mtmp);
  dop_d.freeFermion(resid);

  ldop_f->end();
  delete ldop_f;

}    

};

#endif
