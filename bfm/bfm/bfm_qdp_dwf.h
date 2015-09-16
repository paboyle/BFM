#ifndef _BFM_QDP_DWF_H_
#define _BFM_QDP_DWF_H_

typedef T4 T;
#include <bfm_cg_mixed_prec.h>
#include <omp.h>

/*****************************************************************************
 * Domain Wall Fermion interface
 *****************************************************************************
 */
/////////////////////////////////////////////////////////////
// Inverts M, used for propagators
// May need an interface for CPS to return 
// 5d observables such as J5q, rather than "output_stem"
/////////////////////////////////////////////////////////////
template <class Float_f,class Float_d>
int dwf_CG(LatticePropagator &sol, 
	   LatticePropagator &src,
	   multi1d<LatticeColorMatrix> &U,
	   std::string output_stem,
	   int Ls,
	   int DodwfConservedVector,
	   Real mass,
	   Real M5,int Ncg, Real resid[],int max_iter[]);

template <class Float_f,class Float_d>
int dwf_CG(LatticePropagator &sol, 
	   LatticePropagator &src,
	   multi1d<multi1d<LatticeFermion> > &evecs,
	   multi1d<Complex> &evals,
	   multi1d<LatticeColorMatrix> &U,
	   std::string output_stem,
	   int Ls,
	   int DodwfConservedVector,
	   Real mass,
	   Real M5,int Ncg, Real resid[],int max_iter[]);

//////////////////////////////////////////////////////
// Does same as the above for a single component source
// Maps (perhaps) to FMatInv in CPS
//////////////////////////////////////////////////////
template <class Float_f,class Float_d>
int dwf_CG(LatticeFermion &sol, 
	   LatticeFermion &src,
	   multi1d<LatticeColorMatrix> &U,
	   std::string output_stem,
	   int Ls,
	   Real mass,
	   Real M5,int Ncg, Real resid[],int max_iter[]);

template <class Float_f,class Float_d>
int dwf_CG(LatticeFermion &sol, 
	   LatticeFermion &src,
	   multi1d<multi1d<LatticeFermion> > &evecs,
	   multi1d<Complex> &evals,
	   multi1d<LatticeColorMatrix> &U,
	   std::string output_stem,
	   int Ls,
	   Real mass,
	   Real M5,int Ncg, Real resid[],int max_iter[]);

//////////////////////////////////////////////////////
// Inverts the odd parity of MatPCdagMatPC.
// Maps to FMatEvlInv in CPS
//////////////////////////////////////////////////////
template <class Float>
int dwf_CG_precMdagM_oo(multi1d<LatticeFermion> &sol, 
			multi1d<LatticeFermion> &src,
			multi1d<LatticeColorMatrix> &U,
			std::string output_stem,
			int Ls,
			Real mass,
			Real M5, Real residual, int max_iter);

//////////////////////////////////////////////////////////
// Restarted mixed prec inverter of any of the above.
//////////////////////////////////////////////////////////
template <class Float_f,class Float_d>
int dwf_restarted_invert(inv_type_t inv_type,
			 multi1d<LatticeFermion> &sol, 
			 multi1d<LatticeFermion> &src,
	   		 multi1d<multi1d<LatticeFermion> > &evecs,
	   		 multi1d<Complex> &evals,
			 multi1d<LatticeColorMatrix> &U,
			 std::string output_stem,
			 int Ls,
			 Real mass,
			 Real M5,
			 int Ncg,
			 Real residual[],
			 int max_iter[]);

template <class Float_f,class Float_d>
int dwf_restarted_invert(inv_type_t inv_type,
			 multi1d<LatticeFermion> &sol, 
			 multi1d<LatticeFermion> &src,
			 multi1d<LatticeColorMatrix> &U,
			 std::string output_stem,
			 int Ls,
			 Real mass,
			 Real M5,
			 int Ncg,
			 Real residual[],
			 int max_iter[]);


// Internal
template <class Float>
int dwf_invert(inv_type_t inv_type,
	       multi1d<LatticeFermion> &sol, 
	       multi1d<LatticeFermion> &src,
	       multi1d<multi1d<LatticeFermion> > &evecs,
	       multi1d<Complex> &evals,
	       multi1d<LatticeColorMatrix> &U,
	       std::string output_stem,
	       int Ls,
	       Real mass,
	       Real M5, Real residual, int max_iter);

/**********************************************************************
 * DOMAIN WALL SOLVERS
 **********************************************************************
 */

template <class Float_f,class Float_d>
int dwf_CG(LatticeFermion &psi, 
	   LatticeFermion &chi,
	   multi1d<LatticeColorMatrix> &U,
	   std::string output_stem,
	   int Ls,
	   Real mass,
	   Real M5,int Ncg,
	   Real residual[],
	   int max_iter[])
{
 multi1d<multi1d<LatticeFermion> > evecs(0);
 multi1d<Complex> evals(0);
 int iter = dwf_CG<Float_f,Float_d>(psi,chi,evecs,evals,U,output_stem,
	   Ls,mass,M5,Ncg,residual,max_iter);
  return iter;
}

template <class Float_f,class Float_d>
int dwf_CG(LatticeFermion &psi, 
	   LatticeFermion &chi,
	   multi1d<multi1d<LatticeFermion> > &evecs,
	   multi1d<Complex> &evals,
	   multi1d<LatticeColorMatrix> &U,
	   std::string output_stem,
	   int Ls,
	   Real mass,
	   Real M5,int Ncg,
	   Real residual[],
	   int max_iter[])
{
  multi1d<LatticeFermion> chi_5(Ls), psi_5(Ls);

  chi_5[0]   = chiralProjectPlus(chi);
  for(int s=1;s<Ls-1;s++) chi_5[s] = zero;
  chi_5[Ls-1]= chiralProjectMinus(chi);
	
  for(int s=0;s<Ls;s++) psi_5[s]=zero;

  int iter =  dwf_restarted_invert<Float_f,Float_d>(CG_PREC_M,
						    psi_5,chi_5,evecs,evals,U,
						    output_stem,
						    Ls,mass,M5,
						    Ncg,residual,max_iter);

  QDPIO::cout << "psi_5[0] "<< norm2(psi_5[0])<<endl;
  QDPIO::cout << "psi_5[Ls-1] "<< norm2(psi_5[Ls-1])<<endl;

  psi = chiralProjectMinus(psi_5[0]);
  psi += chiralProjectPlus(psi_5[Ls-1]);

  if ( bfmarg::onepluskappanorm ) {
    Real rescale(5.0-M5);
    psi = psi * rescale;
  }

  return iter;
}

// Internal interfaces
template <class Float>
int dwf_CG_precMdagM_oo(multi1d<LatticeFermion> &sol, 
			multi1d<LatticeFermion> &src,
			multi1d<LatticeColorMatrix> &U,
			std::string output_stem,
			int Ls,
			Real mass,
			Real M5, Real residual, int max_iter)
{
  int iter =  dwf_invert<Float>(CG_PREC_MDAGM, sol,src,
				  output_stem,
				  Ls,mass,M5,
				  residual,max_iter);

  if ( bfmarg::onepluskappanorm ) {
    Real rescale(5.0-M5);
    for(int s=0;s<Ls;s++)  sol[s] = sol[s] * rescale * rescale;
  }
  return iter;
}

template <class Float_f,class Float_d>
int dwf_CG(LatticePropagator &sol, 
	   LatticePropagator &src,
	   multi1d<LatticeColorMatrix> &U,
	   std::string output_stem,
	   int Ls,
	   int do_PImunu,
	   Real mass,
	   Real M5,int Ncg,
	   Real residual[],
	   int max_iter[])
{
  multi1d<multi1d<LatticeFermion> > evecs(0);
  multi1d<Complex> evals(0);
  int iter = dwf_CG<Float_f,Float_d>(
			sol,src,
			evecs,evals,U,output_stem,
	   		Ls,do_PImunu,mass,M5,
			Ncg,residual,max_iter);
  return iter;
}
 
template <class Float_f,class Float_d>
int dwf_CG(LatticePropagator &sol, 
	   LatticePropagator &src,
	   multi1d<multi1d<LatticeFermion> > &evecs,
	   multi1d<Complex> &evals,
	   multi1d<LatticeColorMatrix> &U,
	   std::string output_stem,
	   int Ls,
	   int do_PImunu,
	   Real mass,
	   Real M5,int Ncg,
	   Real residual[],
	   int max_iter[]
	   //	   LatticePropagator &midpoint
	   )
{
  int doXML=0;
  std::string output    = output_stem + ".xml";
  std::string outputbin = output_stem + ".VV.bin";
 
  LatticePropagator midpoint;
  LatticeFermion chi,psi;

  /***************************************/
  /* Collect some useful 5d correlators  */
  /***************************************/
  LatticeComplex PA = zero;
  LatticeComplex PP = zero;
  LatticeComplex PAconsv = zero;
  LatticeComplex PJ5q    = zero;

  multi1d<int> iterations(Nc*Ns);
  multi1d<Real> residuals(Nc*Ns);

  // Optionally do the vacuum polarisation insertions
  multi1d<LatticePropagator> sol_5;
  if ( do_PImunu ) {
    sol_5.resize(Ls);
    for(int s=0;s<Ls;s++) sol_5[s]=zero;
  }
  
  for(int color_source = 0; color_source < Nc; ++color_source) {
    for(int spin_source = 0; spin_source < Ns; ++spin_source) {
	
      int idx = spin_source*3+color_source;

      // Extract a fermion source
      BfmPropToFerm(src, chi, color_source, spin_source);

      /**************************************************************
       * Zero source implies zero solution 
       ***************************************************************
       */
      Real nrm = sqrt(norm2(chi));
      if (toFloat(nrm) == 0.0) {

	psi = zero;
	iterations[idx] = 0;
        BfmFermToProp(psi, sol, color_source, spin_source);

      } else {

	multi1d<LatticeFermion> chi_5(Ls), psi_5(Ls);


	/*******************************************************************
	 * Do the CG :
	 *
	 * Normalize the source in case it is really huge or small -
	 * a trick to avoid overflows or underflows.
	 *
	 * Splat the chiralities of the source to the walls
	 *******************************************************************/

	Real fact = 1.0/nrm;
	chi *= fact;
	chi_5[0]   = chiralProjectPlus(chi);
	for(int s=1;s<Ls-1;s++) chi_5[s] = zero;
	chi_5[Ls-1]= chiralProjectMinus(chi);
	
	for(int s=0;s<Ls;s++) psi_5[s]=zero;

	// Call the the CG
	int iter = 
	  dwf_restarted_invert<Float_f,Float_d>(CG_PREC_M,
						psi_5,chi_5,evecs,evals,U,
						output_stem,
						Ls,mass,M5,
						Ncg,residual,max_iter);
	iterations[idx] = iter;

        fact = Real(1) / fact;
	if ( bfmarg::onepluskappanorm ) {
	  Real rescale(5.0-M5);
	  fact = fact *rescale;
	}

        for(int s=0;s<Ls;s++) psi_5[s] *= fact;

	/*******************************************************************
	 * Accumulate important 5d matrix elements
	 *******************************************************************
	 */
	      /****************************/  
	      /* Get the midpoint beastie */
	      /****************************/  

	psi = chiralProjectPlus(psi_5[Ls/2-1]);
	psi += chiralProjectMinus(psi_5[Ls/2]);
	BfmFermToProp(psi,midpoint,color_source,spin_source);
	PJ5q+= localInnerProduct(psi,psi);
	
	      /****************************/  
	      /* Get the conserved axial  */
	      /****************************/  
	LatticeComplex C;
	LatticeFermion    p5d;
	LatticeFermion us_p5d;
	int mu=Nd-1;
	int g5=Ns*Ns-1;
	int gt=8;

	for(int s=0;s<Ls;s++){

	  p5d = psi_5[Ls-1-s] ;
	  us_p5d = U[mu]*shift(psi_5[s],FORWARD,mu);

	  C = 0.5*localInnerProduct(Gamma(g5)*p5d,Gamma(gt)*us_p5d);
	  C-= 0.5*localInnerProduct(Gamma(g5)*p5d,us_p5d);

	  p5d = psi_5[s];
	  us_p5d = U[mu]*shift(psi_5[Ls-1-s],FORWARD,mu);

	  C+= 0.5*localInnerProduct(Gamma(g5)*us_p5d,Gamma(gt)*p5d);
	  C+= 0.5*localInnerProduct(Gamma(g5)*us_p5d,p5d);

	  if (s < Ls/2) PAconsv -= C;
	  else          PAconsv += C;

	}


	/****************************/  
	/* Get the full VV AA       */
	/****************************/  

	if ( do_PImunu ) {
	  for(int s=0;s<Ls;s++){
	    BfmFermToProp(psi_5[s], sol_5[s], color_source, spin_source);
	  }
	}

	/********************************************************
	 * Merge the chiralities from the walls to reconstruct 4d prop
	 ********************************************************/
	psi = chiralProjectMinus(psi_5[0]);
	psi += chiralProjectPlus(psi_5[Ls-1]);

        BfmFermToProp(psi, sol, color_source, spin_source);

	/*******************************************************************
	 * Accumulate important 4d matrix elements
	 *******************************************************************
	 */
	PA += localInnerProduct(psi,Gamma(1<<Nd-1)*psi);
	PP += localInnerProduct(psi,psi);

      }
    } /* end loop over spin_source */
  } /* end loop over color_source */

  
#ifdef BFM_CHROMA_BINDING
  /************************************************************************/
  /* Contract PImunu for V and A                                          */
  /************************************************************************/
  if ( do_PImunu ) {
   BinaryFileWriter WR(outputbin);

   // 3d FT
   SftMom Phases(2,false,Nd-1);
   int length = Phases.numSubsets();
   int nmom   = Phases.numMom();
   

    for(int mu=0;mu<Nd;mu++){ // Sink direction (of conserved current)
    for(int nu=0;nu<Nd;nu++){ // source direction (of local current)

      // Collect the P channel
      LatticeComplex PAc = zero;
      LatticeComplex CP;
      // AA
      LatticeComplex AAc = zero;
      LatticeComplex CA;
      // VV
      LatticeComplex VVc = zero;
      LatticeComplex CV;

      LatticePropagator p5d;
      LatticePropagator us_p5d;

      SpinMatrix G5 =1.0; G5  = G5*Gamma(15);
      SpinMatrix Gmu=1.0; Gmu = Gmu*Gamma(1<<mu);
      SpinMatrix Gnu=1.0; Gnu = Gnu*Gamma(1<<nu);

      for(int s=0;s<Ls;s++){
 
        p5d = sol_5[Ls-1-s] ;
        us_p5d = U[mu]*shift(sol_5[s],FORWARD,mu);

        CV = 0.5*trace(Gnu*G5*adj(p5d)*G5*Gmu*us_p5d);
        CV-= 0.5*trace(Gnu*G5*adj(p5d)*G5*us_p5d);

        CA = 0.5*trace(Gnu*adj(p5d)*G5*Gmu*us_p5d);
        CA-= 0.5*trace(Gnu*adj(p5d)*G5*us_p5d);

        CP = 0.5*trace(adj(p5d)*G5*Gmu*us_p5d);
        CP-= 0.5*trace(adj(p5d)*G5*us_p5d);

        p5d = sol_5[s];
        us_p5d = U[mu]*shift(sol_5[Ls-1-s],FORWARD,mu);

        CV+= 0.5*trace(Gnu*G5*adj(us_p5d)*G5*Gmu*p5d);
        CV+= 0.5*trace(Gnu*G5*adj(us_p5d)*G5*p5d);
        CA+= 0.5*trace(Gnu*adj(us_p5d)*G5*Gmu*p5d);
        CA+= 0.5*trace(Gnu*adj(us_p5d)*G5*p5d);
        CP+= 0.5*trace(adj(us_p5d)*G5*Gmu*p5d);
        CP+= 0.5*trace(adj(us_p5d)*G5*p5d);

        VVc += CV;
        if (s < Ls/2) {
	  PAc -= CP;
	  AAc -= CA;
	} else {
	  PAc += CP;
	  AAc += CA;
	}
	
      }// s Loop

      ////////////////////////////////////
      // Do the fourier transform
      ////////////////////////////////////
      multi2d<DComplex> corr(nmom,length);

      corr = Phases.sft(PAc);
      write(WR,corr);
      
      corr = Phases.sft(AAc);
      write(WR,corr);

      corr = Phases.sft(VVc);
      write(WR,corr);

    }
    }//munu

    multi1d<multi1d<int> > momlist(nmom);
    for(int m=0;m<nmom;m++) momlist[m] = Phases.numToMom(m);
    write(WR,momlist);

  }
#endif


  /************************************************************************/
  /* Spatially sum the correlation functions                              */
  /************************************************************************/
  Set &tslice = GetTimeslice();
  int length = tslice.numSubsets();
  multi1d<DComplex> corr(length);
  multi1d<Real>    rcorr(length);

  /***********************************************************************/
  /* Write the XM-HELL                                                   */
  /***********************************************************************/
  if ( doXML ) {
    XMLFileWriter xml(output);
    push(xml,"DWF_prop");
    write(xml,"iterations",iterations);
    write(xml,"residuals",residuals);
    push(xml,"DWF_observables");

      corr = sumMulti(PP,tslice);  
      for(int t=0;t<length;t++)rcorr[t]=real(corr[t]);
      write(xml,"PP",rcorr);

      corr = sumMulti(PA,tslice);  
      for(int t=0;t<length;t++)rcorr[t]=real(corr[t]);
      write(xml,"PA",rcorr);

      corr = sumMulti(PAconsv,tslice);  
      for(int t=0;t<length;t++)rcorr[t]=real(corr[t]);
      write(xml,"PAconsv",rcorr);

      corr = sumMulti(PJ5q,tslice);  
      for(int t=0;t<length;t++)rcorr[t]=real(corr[t]);
      write(xml,"PJ5q",rcorr);

    pop(xml);     
    pop(xml);      
  } else {
      corr = sumMulti(PP,tslice);  
      for(int t=0;t<length;t++) {
	rcorr[t]=real(corr[t]);
	QDPIO::cout << "PP["<< t<<"] = " << rcorr[t]<<endl;
      }

      corr = sumMulti(PA,tslice);  
      for(int t=0;t<length;t++) { 
	rcorr[t]=real(corr[t]);
	QDPIO::cout << "PA["<< t<<"] = " << rcorr[t]<<endl;
      }

      corr = sumMulti(PAconsv,tslice);  
      for(int t=0;t<length;t++) { 
	rcorr[t]=real(corr[t]);
	QDPIO::cout << "PAconsv["<< t<<"] = " << rcorr[t]<<endl;
      }

      corr = sumMulti(PJ5q,tslice);  
      for(int t=0;t<length;t++) {
	rcorr[t]=real(corr[t]);
	QDPIO::cout << "PJ5q["<< t<<"] = " << rcorr[t]<<endl;
      }

  }
  return 0; 
}

template <class Float_f,class Float_d>
int dwf_restarted_invert(inv_type_t inv_type,
			 multi1d<LatticeFermion> &sol, 
			 multi1d<LatticeFermion> &src,
			 multi1d<LatticeColorMatrix> &U,
			 std::string output_stem,
			 int Ls,
			 Real mass,
			 Real M5,
			 int Ncg,
			 Real residual[],
			 int max_iter[])
{
  multi1d<multi1d<LatticeFermion> > evecs(0);
  multi1d<Complex> evals(0);

  return dwf_restarted_invert<Float_f,Float_d>(inv_type,sol,src,
					       evecs,evals,
					       U,output_stem,Ls,mass,M5,Ncg,residual,max_iter);
}

// Internal interfaces
template <class Float_f,class Float_d>
int dwf_restarted_invert(inv_type_t inv_type,
			 multi1d<LatticeFermion> &sol, 
			 multi1d<LatticeFermion> &src,
	       		 multi1d<multi1d<LatticeFermion> > &evecs,
	       		 multi1d<Complex> &evals,
			 multi1d<LatticeColorMatrix> &U,
			 std::string output_stem,
			 int Ls,
			 Real mass,
			 Real M5,
			 int Ncg,
			 Real residual[],
			 int max_iter[])
{
  multi1d<int> iter(Ncg);

  multi1d<LatticeFermion> sol_f(Ls);
  std::string stem_f = output_stem + ".prop.f";
  std::string stem_d = output_stem + ".prop.d";

  for(int s=0;s<Ls;s++) sol_f[s] = sol[s];
  for(int cg=0;cg<Ncg-1;cg++){
    iter[cg] = 
      dwf_invert<Float_f>(inv_type,
			  sol_f,src,
			  evecs, evals,
			  U,stem_f,
			  Ls,mass,M5,
			  residual[cg],max_iter[cg]);
  }
  int iter_all=0;
  for(int s=0;s<Ls;s++) sol[s] = sol_f[s];

  {
    int cg = Ncg-1;
    iter[cg] = dwf_invert<Float_d>(inv_type,sol,src,evecs,evals,U,stem_d,
				   Ls,mass,M5,
				   residual[cg],max_iter[cg]);
    iter_all+=iter[cg];
  }

  for(int cg=0;cg<Ncg;cg++){
    QDPIO::cout <<"bfm_qdp::"<<iter[cg]<<" pass " << cg<< " iterations"<<endl;
  }
  return iter_all;
}

// Internal interfaces
template <class Float>
int dwf_invert(inv_type_t inv_type,multi1d<LatticeFermion> &sol, 
	       multi1d<LatticeFermion> &src,
	       multi1d<multi1d<LatticeFermion> > &evecs,
	       multi1d<Complex> &evals,
	       multi1d<LatticeColorMatrix> &U,
	       std::string output_stem,
	       int Ls,
	       Real mass,
	       Real M5, Real residual, int max_iter)
{
  // Set up BAGEL object
  int lx = QDP::Layout::subgridLattSize()[0];
  int ly = QDP::Layout::subgridLattSize()[1];
  int lz = QDP::Layout::subgridLattSize()[2];
  int lt = QDP::Layout::subgridLattSize()[3];

  multi1d<int> procs = QDP::Layout::logicalSize();

  bfmarg dwfa;

  //Physics parameters
  dwfa.solver       = HtCayleyTanh;
  dwfa.Ls           = Ls;
  dwfa.M5           = toDouble(M5);
  dwfa.mass         = toDouble(mass);
  dwfa.precon_5d    = 0;
  dwfa.max_iter     = max_iter;
  dwfa.residual     = toDouble(residual);

  //Geometry
  dwfa.node_latt[0] = lx;
  dwfa.node_latt[1] = ly;
  dwfa.node_latt[2] = lz;
  dwfa.node_latt[3] = lt;
  for(int mu=0;mu<4;mu++){
    if (procs[mu]>1) dwfa.local_comm[mu] = 0;
    else             dwfa.local_comm[mu] = 1;
  }

  QDPIO::cout << "bfm_qdp:: Initialising BAGEL-2 solver "<<endl;

  bfm_qdp<Float> dwf; dwf.init(dwfa);

  dwf.importGauge(U);

  for(int cb=0;cb<2;cb++){
    dwf.qdp_psi_h[cb] = dwf.allocFermion();
    dwf.qdp_chi_h[cb] = dwf.allocFermion();
    dwf.importFermion(src,dwf.qdp_psi_h[cb],cb);
    dwf.importFermion(sol,dwf.qdp_chi_h[cb],cb);
  }
  multi1d<bfm_fermion> bfm_evecs( evecs.size() );
  multi1d<double> bfm_evals( evals.size() );
  for(int n = 0; n< evecs.size(); n++){
    dwf.qdp_evals[n] = toDouble( real(evals[n]) );
    for(int cb=0;cb<2;cb++){
      dwf.qdp_evecs[n][cb] = dwf.allocFermion();
      dwf.importFermion(evecs[n],bfm_evecs[n][cb],cb);
    }
  }

  QDPIO::cout << "bfm_qdp:: Spawning "<<dwf.threads<<"threads for CG "<<endl;
  int iter;
  dwf.inv_type = inv_type;
  bfm_spawn_cg(dwf);
  QDPIO::cout << "bfm_qdp:: collecting results from threads "<<endl;

  for(int n = 0; n< evecs.size(); n++){
    for(int cb=0;cb<2;cb++){
      dwf.freeFermion(dwf.qdp_evecs[n][cb]);
    }
  }

  for(int cb=0;cb<2;cb++){
    dwf.exportFermion(sol,dwf.qdp_chi_h[cb],cb);
    dwf.freeFermion(dwf.qdp_psi_h[cb]);
    dwf.freeFermion(dwf.qdp_chi_h[cb]);
  }
  dwf.end();

#ifdef BFM_CHROMA_BINDING
  /*******************************************************/
  /* verify the solution                                 */
  /*******************************************************/
  if ( dwf.inv_type ==  CG_PREC_M ) {

    QDPIO::cout << "bfm_qdp:: Verifying solution" <<endl;

    /*
     * Compute the residual according to CHROMA
     */
    multi1d<LatticeFermion> regress(Ls);
    multi1d<int> bcs(Nd);
    bcs[0] = bcs[1] = bcs[2] = bcs[3] = 1;

    Handle< FermBC<T,G,G> > fbc(new SimpleFermBC< T, G, G >(bcs));
    Handle<CreateFermState<T,G,G> > cfs( new CreateSimpleFermState<T,G,G>(fbc));
    
    UnprecDWFermActArray  S_f(cfs, M5, mass, Ls);
    Handle< FermState<T,G,G> > fs( S_f.createState(U) );

    Handle< UnprecLinearOperatorArray<T,G,G> > M(S_f.unprecLinOp(fs,mass));
  
    // Check the result
    PlusMinus pm = PLUS;

    (*M)(regress,sol,pm);

    regress = regress - src;
    QDPIO::cout << "bfm_qdp:: QDP regression check :  |M sol - src| = " << norm2(regress) << endl;
    if ( toDouble( norm2(regress) / norm2(src) ) > 1.0e-5 ) { 
      QDPIO::cout << "bfm_qdp:: QDP regression check : This is worryingly large - PANIC"<< endl;
      exit(-1);
    }
  }
#endif
  return iter;

}

int dwf_mixed_precision_CG(LatticeFermion &sol, 
			   LatticeFermion &src,
			   multi1d<LatticeColorMatrix> &U,		     
			   int Ls, Real mass, Real M5,double residual,int max_iter);

int dwf_mixed_precision_CG(multi1d<LatticeFermion> &sol, 
			   multi1d<LatticeFermion> &src,
			   multi1d<LatticeColorMatrix> &U,		     
			   int Ls, Real mass, Real M5,double residual,int max_iter);

#endif
