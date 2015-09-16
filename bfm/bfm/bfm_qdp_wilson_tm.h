#ifndef _BFM_QDP_WILSON_TM_H_
#define _BFM_QDP_WILSON_TM_H_

/*****************************************************************************
 * Twisted Mass Wilson Fermion interface. Only implement enough here that
 * I can support the RBC-UKQCD AuxDet HMC.
 *****************************************************************************
 */

//////////////////////////////////////////////////////
// Inverts the odd parity of MatPCdagMatPC.
// Maps to FMatEvlInv in CPS
//////////////////////////////////////////////////////
template <class Float>
int wilsontm_CG_precMdagM_oo(LatticeFermion &sol, 
			     LatticeFermion &src,
			     multi1d<LatticeColorMatrix> &U,
			     Real mass, Real mu,
			     Real residual, int max_iter);

// Internal
template <class Float>
int wilsontm_invert(inv_type_t inv_type,
		  LatticeFermion &sol, 
		  LatticeFermion &src,
		  multi1d<LatticeColorMatrix> &U,
	          Real mass,Real mu,
		  Real residual, int max_iter);



// Internal interfaces
template <class Float>
int wilsontm_CG_precMdagM_oo(LatticeFermion &sol, 
			     LatticeFermion &src,
			     multi1d<LatticeColorMatrix> &U,
			     Real mass,Real mu,
			     Real residual, int max_iter)
{
  int iter =  wilsontm_invert<Float>(CG_PREC_MDAGM, sol,src,U,
				     mass,mu,
				     residual,max_iter);

  if ( bfmarg::onepluskappanorm ) {
    Real rescale( (4.0+mass)*(4.0+mass) + mu*mu);
    sol = sol * rescale;
  }
  return iter;
}

// Internal interfaces
template <class Float>
int wilsontm_invert(inv_type_t inv_type,LatticeFermion &sol, 
		    LatticeFermion &src,
		    multi1d<LatticeColorMatrix> &U,
		    Real mass, Real mu, Real residual, int max_iter)
{
  // Set up BAGEL object
  int lx = QDP::Layout::subgridLattSize()[0];
  int ly = QDP::Layout::subgridLattSize()[1];
  int lz = QDP::Layout::subgridLattSize()[2];
  int lt = QDP::Layout::subgridLattSize()[3];

  multi1d<int> procs = QDP::Layout::logicalSize();

  bfmarg wilsa;

  //Physics parameters
  wilsa.solver       = WilsonTM;
  wilsa.Ls           = 1;
  wilsa.M5           = 0.0;
  wilsa.mass         = toDouble(mass);
  wilsa.twistedmass  = toDouble(mu);
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
  QDPIO::cout << "FIXME : Cannot verify solution as no reference code" << endl;
#endif
   return iter;

}


#endif
