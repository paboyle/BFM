#ifndef _BFM_QDP_DWF_GMINUS2_H_
#define _BFM_QDP_DWF_GMINUS2_H_
template <class Float_f,class Float_d>
int dwf_gminus2(LatticePropagator &src, // expect point source or Z2Wall for this
		multi1d<LatticeColorMatrix> &Uperi,
		multi1d<LatticeColorMatrix> &Utwist,
		std::string output_stem,
		int Ls,
		Real mass,
		Real M5,int Ncg, Real resid[],int max_iter[]);

int dwf_gminus2(LatticePropagator &src, // expect point source or Z2Wall for this
		multi1d<LatticeColorMatrix> &Uperi,
		multi1d<LatticeColorMatrix> &Utwist,
		std::string output_stem,
		int Ls,
		Real mass,
		  Real M5,int Ncg, Real resid[],int max_iter[])
{
  multi1d<multi1d<LatticeFermion> > evecs(0);
  multi1d<Complex> evals(0);
  
  int doXML=0;
  std::string output    = output_stem + ".xml";
  std::string outputbin = output_stem + ".VV.bin";
  
  /***************************************/
  /* Collect some useful 5d correlators  */
  /***************************************/
  LatticePropagator midpoint_twist, midpoint_peri;
  LatticeComplex PA = zero;
  LatticeComplex PP = zero;
  LatticeComplex PAconsv = zero;
  LatticeComplex PJ5q    = zero;

  LatticeFermion src_4,ferm_twist,ferm_peri;
  LatticePropagator prop_peri, prop_twist;

  multi1d<LatticeFermion> src_5(Ls), ferm_5_peri(Ls), ferm_5_twist(Ls);
  multi1d<LatticePropagator> prop_5_peri(Ls), prop_5_twist(Ls);

  prop_5_peri = zero;
  prop_5_twist= zero;

  for(int color_source = 0; color_source < Nc; ++color_source) {
    for(int spin_source = 0; spin_source < Ns; ++spin_source) {
      
      int idx = spin_source*3+color_source;
      
      // Extract a fermion source
      BfmPropToFerm(src, src_4, color_source, spin_source);
      Real nrm = sqrt(norm2(src_4));

      if (toFloat(nrm) != 0.0) {

	/////////////////////////////////////////////////
	// Build normalised source in src_5
	/////////////////////////////////////////////////
	Real fact = 1.0/nrm;
	src_4 *= fact;
	src_5[0]   = chiralProjectPlus(src_4);
	for(int s=1;s<Ls-1;s++) src_5[s] = zero;
	src_5[Ls-1]= chiralProjectMinus(src_4);
	
	/////////////////////////////////////////////////
	// Solve for both twisted and untisted BC's
	/////////////////////////////////////////////////
	for(int s=0;s<Ls;s++) ferm_5_peri [s]=zero;
	for(int s=0;s<Ls;s++) ferm_5_twist[s]=zero;

	dwf_restarted_invert<Float_f,Float_d>(CG_PREC_M,
					      ferm_5_peri,src_5,evecs,evals,Uperi,
					      output_stem,
					      Ls,mass,M5,
					      Ncg,residual,max_iter);

	dwf_restarted_invert<Float_f,Float_d>(CG_PREC_M,
					      ferm_5_twist,src_5,evecs,evals,Utwist,
					      output_stem,
					      Ls,mass,M5,
					      Ncg,residual,max_iter);

        fact = Real(1) / fact;
	if ( bfmarg::onepluskappanorm ) {
	  Real rescale(5.0-M5);
	  fact = fact *rescale;
	}
        for(int s=0;s<Ls;s++) ferm_5_twist[s] *= fact;
        for(int s=0;s<Ls;s++) ferm_5_peri[s] *= fact;

	/*******************************************************************
	 * Accumulate important 5d matrix elements
	 *******************************************************************
	 */
	      /****************************/  
	      /* Get the midpoint beastie */
	      /****************************/  

	LatticeFermion psi_twist,psi_peri;
	psi_twist = chiralProjectPlus (ferm_5_twist[Ls/2-1]);
	psi_twist+= chiralProjectMinus(ferm_5_twist[Ls/2]);
	psi_peri  = chiralProjectPlus (ferm_5_peri [Ls/2-1]);
	psi_peri += chiralProjectMinus(ferm_5_peri [Ls/2]);

	BfmFermToProp(psi_twist,midpoint_twist,color_source,spin_source);
	BfmFermToProp(psi_peri, midpoint_peri ,color_source,spin_source);

	PJ5q+= localInnerProduct(psi_twist,psi_peri);
	
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

	  p5d = ferm_5_twist[Ls-1-s] ;
	  us_p5d = Uperi[mu]*shift(ferm_5_peri[s],FORWARD,mu);

	  C = 0.5*localInnerProduct(Gamma(g5)*p5d,Gamma(gt)*us_p5d);
	  C-= 0.5*localInnerProduct(Gamma(g5)*p5d,us_p5d);

	  // Is this correct?
	  p5d = ferm_5_peri[s];
	  us_p5d = Utwist[mu]*shift(ferm_5_twist[Ls-1-s],FORWARD,mu);

	  C+= 0.5*localInnerProduct(Gamma(g5)*us_p5d,Gamma(gt)*p5d);
	  C+= 0.5*localInnerProduct(Gamma(g5)*us_p5d,p5d);

	  if (s < Ls/2) PAconsv -= C;
	  else          PAconsv += C;

	}


	/****************************/  
	/* Get the full VV AA       */
	/****************************/  

	for(int s=0;s<Ls;s++){
	  BfmFermToProp(ferm_5_twist[s], prop_5_twist[s], color_source, spin_source);
	  BfmFermToProp(ferm_5_peri [s], prop_5_peri [s], color_source, spin_source);
	}

	/********************************************************
	 * Merge the chiralities from the walls to reconstruct 4d prop
	 ********************************************************/
	ferm_twist = chiralProjectMinus(ferm_5_twist[0]);
	ferm_twist+= chiralProjectPlus (ferm_5_twist[Ls-1]);

	ferm_peri = chiralProjectMinus(ferm_5_peri[0]);
	ferm_peri+= chiralProjectPlus (ferm_5_peri[Ls-1]);

        BfmFermToProp(ferm_peri , prop_peri, color_source, spin_source);
        BfmFermToProp(ferm_twist, prop_twist, color_source, spin_source);

	/*******************************************************************
	 * Accumulate important 4d matrix elements
	 *******************************************************************
	 */
	PA += localInnerProduct(ferm_twist,Gamma(1<<Nd-1)*ferm_peri);
	PP += localInnerProduct(ferm_twist,ferm_peri);

      }
    } /* end loop over spin_source */
  } /* end loop over color_source */


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
 
        p5d = prop_5_twist[Ls-1-s] ;
        us_p5d = Uperi[mu]*shift(prop_5_peri[s],FORWARD,mu);

        CV = 0.5*trace(Gnu*G5*adj(p5d)*G5*Gmu*us_p5d);
        CV-= 0.5*trace(Gnu*G5*adj(p5d)*G5*us_p5d);

        CA = 0.5*trace(Gnu*adj(p5d)*G5*Gmu*us_p5d);
        CA-= 0.5*trace(Gnu*adj(p5d)*G5*us_p5d);

        CP = 0.5*trace(adj(p5d)*G5*Gmu*us_p5d);
        CP-= 0.5*trace(adj(p5d)*G5*us_p5d);

        p5d = prop_5_peri[s];
        us_p5d = Utwist[mu]*shift(prop_5_twist[Ls-1-s],FORWARD,mu);

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

      LatticeComplex CAL = zero;
      LatticeComplex CVL = zero;
 
      CAL = -trace (Gnu*adj(prop_twist)*Gmu*prop_peri) ;
      CVL =  trace (Gnu*G5*adj(prop_twist)*G5*Gmu*prop_peri ) ;

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

      corr = Phases.sft(CAL);
      write(WR,corr);

      corr = Phases.sft(CVL);
      write(WR,corr);
    }
    }//munu

    


    multi1d<multi1d<int> > momlist(nmom);
    for(int m=0;m<nmom;m++) momlist[m] = Phases.numToMom(m);
    write(WR,momlist);

  }

#endif
