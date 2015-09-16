#ifndef _BFM_QDP_WILSON_H_
#define _BFM_QDP_WILSON_H_

typedef T4 T;

/*****************************************************************************
 * Wilson Fermion interface
 *****************************************************************************
 */
/////////////////////////////////////////////////////////////
// Inverts M, used for propagators
// May need an interface for CPS to return 
/////////////////////////////////////////////////////////////
template <class Float_f,class Float_d>
int wilson_CG(LatticePropagator &sol, 
  	      LatticePropagator &src,
	      multi1d<LatticeColorMatrix> &U,
	      Real mass,
	      int Ncg, Real resid[],int max_iter[]);

//////////////////////////////////////////////////////
// Does same as the above for a single component source
// Maps (perhaps) to FMatInv in CPS
//////////////////////////////////////////////////////
template <class Float_f,class Float_d>
int wilson_CG(LatticeFermion &sol, 
	      LatticeFermion &src,
	      multi1d<LatticeColorMatrix> &U,
	      Real mass,
	      int Ncg, Real resid[],int max_iter[]);

//////////////////////////////////////////////////////
// Inverts the odd parity of MatPCdagMatPC.
// Maps to FMatEvlInv in CPS
//////////////////////////////////////////////////////
template <class Float>
int wilson_CG_precMdagM_oo(LatticeFermion &sol, 
			   LatticeFermion &src,
			multi1d<LatticeColorMatrix> &U,
			Real mass,
			Real residual, int max_iter);


//////////////////////////////////////////////////////////
// Restarted mixed prec inverter of any of the above.
//////////////////////////////////////////////////////////
template <class Float_f,class Float_d>
int wilson_restarted_invert(inv_type_t inv_type,
			    LatticeFermion &sol, 
			    LatticeFermion &src,
			    multi1d<LatticeColorMatrix> &U,
			    Real mass,
			    int Ncg,
			    Real residual[],
			    int max_iter[]);

// Internal
template <class Float>
int wilson_invert(inv_type_t inv_type,
		  LatticeFermion &sol, 
		  LatticeFermion &src,
		  multi1d<LatticeColorMatrix> &U,
	          Real mass,
		  Real residual, int max_iter);



/**********************************************************************
 * WILSON SOLVERS
 **********************************************************************
 */
template <class Float_f,class Float_d>
int wilson_CG(LatticeFermion &psi, 
	      LatticeFermion &chi,
	      multi1d<LatticeColorMatrix> &U,
	      Real mass,
	      int Ncg,
	      Real residual[],
	      int max_iter[])
{
  int iter =  wilson_restarted_invert<Float_f,Float_d>(CG_PREC_M,
						    psi,chi,U,
						    mass,
						    Ncg,residual,max_iter);
  if ( bfmarg::onepluskappanorm ) {
    Real rescale(4.0+mass);
    psi = psi * rescale;
  }
  return iter;
}

// Internal interfaces
template <class Float>
int wilson_CG_precMdagM_oo(LatticeFermion &sol, 
			   LatticeFermion &src,
			   multi1d<LatticeColorMatrix> &U,
			   Real mass,
			   Real residual, int max_iter)
{
  int iter =  wilson_invert<Float>(CG_PREC_MDAGM, sol,src,
				   mass,
				   residual,max_iter);

  if ( bfmarg::onepluskappanorm ) {
    Real rescale(4.0+mass);
    sol = sol * rescale*rescale;
  }
  return iter;
}



template <class Float_f,class Float_d>
int wilson_CG(LatticePropagator &sol, 
	   LatticePropagator &src,
	   multi1d<LatticeColorMatrix> &U,
	   Real mass,
	   int Ncg,
	   Real residual[],
	   int max_iter[])
{

  multi1d<int> iterations(Nc*Ns);
  multi1d<Real> residuals(Nc*Ns);
  LatticeFermion chi,psi;
  
  for(int color_source = 0; color_source < Nc; ++color_source) {
    for(int spin_source = 0; spin_source < Ns; ++spin_source) {
	
      int idx = spin_source*3+color_source;

      // Extract a fermion source
      BfmPropToFerm(src, chi, color_source, spin_source);

      /**************************************************************
       * Zero source implies zero solution 
       ***************************************************************
       */
      psi=zero;
      Real nrm = sqrt(norm2(chi));
      if (toFloat(nrm) == 0.0) {

	iterations[idx] = 0;
        BfmFermToProp(psi, sol, color_source, spin_source);

      } else {

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

	// Call the the CG
	int iter = 
	  wilson_restarted_invert<Float_f,Float_d>(CG_PREC_M,
						psi,chi,U,
						mass,
						Ncg,residual,max_iter);
	iterations[idx] = iter;

        fact = Real(1) / fact;
	if ( bfmarg::onepluskappanorm ) {
	  Real rescale(4.0+mass);
	  fact = fact * rescale;
	}
        psi *= fact;

        BfmFermToProp(psi, sol, color_source, spin_source);

      }
    } /* end loop over spin_source */
  } /* end loop over color_source */

  return 0; 
}


// Internal interfaces
template <class Float_f,class Float_d>
int wilson_restarted_invert(inv_type_t inv_type,
			    LatticeFermion &sol, 
			    LatticeFermion &src,
			    multi1d<LatticeColorMatrix> &U,
			    Real mass,
			    int Ncg,
			    Real residual[],
			    int max_iter[])
{
  multi1d<int> iter(Ncg);


  for(int cg=0;cg<Ncg-1;cg++){
    iter[cg] = 
      wilson_invert<Float_f>(inv_type,
			     sol,src,U,
			     mass,
			     residual[cg],max_iter[cg]);
  }
  int iter_all=0;
  {
    int cg = Ncg-1;
    iter[cg] = wilson_invert<Float_d>(inv_type,sol,src,U,
				   mass,residual[cg],max_iter[cg]);
    iter_all+=iter[cg];
  }

  for(int cg=0;cg<Ncg;cg++){
    QDPIO::cout <<"bfm_qdp::"<<iter[cg]<<" pass " << cg<< " iterations"<<endl;
  }
  return iter_all;
}

// Internal interfaces
template <class Float>
int wilson_invert(inv_type_t inv_type,LatticeFermion &sol, 
 		  LatticeFermion &src,
		  multi1d<LatticeColorMatrix> &U,
		  Real mass, Real residual, int max_iter)
{
  // Set up BAGEL object
  int lx = QDP::Layout::subgridLattSize()[0];
  int ly = QDP::Layout::subgridLattSize()[1];
  int lz = QDP::Layout::subgridLattSize()[2];
  int lt = QDP::Layout::subgridLattSize()[3];

  multi1d<int> procs = QDP::Layout::logicalSize();

  bfmarg wilsa;

  //Physics parameters
  wilsa.solver       = WilsonFermion;
  wilsa.Ls           = 1;
  wilsa.M5           = 0.0;
  wilsa.mass         = toDouble(mass);
  wilsa.precon_5d    = 0;
  wilsa.max_iter     = max_iter;
  wilsa.residual     = toDouble(residual);

  //Geometry
  wilsa.node_latt[0] = lx;
  wilsa.node_latt[1] = ly;
  wilsa.node_latt[2] = lz;
  wilsa.node_latt[3] = lt;
  for(int mu=0;mu<4;mu++){
    if (procs[mu]>1) wilsa.local_comm[mu] = 0;
    else             wilsa.local_comm[mu] = 1;
  }

  QDPIO::cout << "bfm_qdp:: Initialising BAGEL-2 solver "<<endl;

  bfm_qdp<Float> wils; wils.init(wilsa);

  wils.importGauge(U);

  for(int cb=0;cb<2;cb++){
    wils.qdp_psi_h[cb] = wils.allocFermion();
    wils.qdp_chi_h[cb] = wils.allocFermion();
    wils.importFermion(src,wils.qdp_psi_h[cb],cb);
    wils.importFermion(sol,wils.qdp_chi_h[cb],cb);
  }

  QDPIO::cout << "bfm_qdp:: Spawning "<<wils.threads<<"threads for CG "<<endl;
  int iter;
  wils.inv_type = inv_type;
  bfm_spawn_cg(wils);
  QDPIO::cout << "bfm_qdp:: collecting results from threads "<<endl;

  for(int cb=0;cb<2;cb++){
    wils.exportFermion(sol,wils.qdp_chi_h[cb],cb);
    wils.freeFermion(wils.qdp_psi_h[cb]);
    wils.freeFermion(wils.qdp_chi_h[cb]);
  }
  wils.end();

#ifdef BFM_CHROMA_BINDING
  /*******************************************************/
  /* verify the solution                                 */
  /*******************************************************/
  if ( wils.inv_type ==  CG_PREC_M ) {

    QDPIO::cout << "bfm_qdp:: Verifying solution" <<endl;

    /*
     * Compute the residual according to CHROMA
     */
    LatticeFermion regress;
    multi1d<int> bcs(Nd);
    bcs[0] = bcs[1] = bcs[2] = bcs[3] = 1;

    Handle< FermBC<T,G,G> > fbc(new SimpleFermBC< T, G, G >(bcs));
    Handle<CreateFermState<T,G,G> > cfs( new CreateSimpleFermState<T,G,G>(fbc));
    UnprecWilsonFermAct  S_f(cfs, mass);
    Handle< FermState<T,G,G> > fs( S_f.createState(U) );
    Handle< UnprecLinearOperator<T,G,G> > M(S_f.linOp(fs));
  
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


#endif
