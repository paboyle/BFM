#include <bfm.h>
#include <bfm_qdp.h>

Set &GetTimeslice(void)
{
  if ( !timeslice_is_made )
    timeslice.make(TSliceFunc(Nd-1));
  return timeslice;
}


// Explicit template instantiations.
// Internal interfaces
template int dwf_invert<float> (inv_type_t inv,
				multi1d<LatticeFermion> &sol, 
				multi1d<LatticeFermion> &src, 
	       			multi1d<multi1d<LatticeFermion> > &evecs,
	       			multi1d<Complex> &evals,
				multi1d<LatticeColorMatrix> &U,
				std::string output_stem, 
				int Ls, Real mass, Real M5, Real residual, int max_iter);
// Internal interfaces
template int dwf_invert<double> (inv_type_t inv,
				 multi1d<LatticeFermion> &sol, 
				 multi1d<LatticeFermion> &src, 
	       			 multi1d<multi1d<LatticeFermion> > &evecs,
	       			 multi1d<Complex> &evals,
				 multi1d<LatticeColorMatrix> &U,
				 std::string output_stem, 
				 int Ls, Real mass, Real M5, Real residual, int max_iter);

// Internal interfaces
template int dwf_restarted_invert<float,double> (inv_type_t inv,
						 multi1d<LatticeFermion> &sol, 
						 multi1d<LatticeFermion> &src, 
	       					 multi1d<multi1d<LatticeFermion> > &evecs,
	       					 multi1d<Complex> &evals,
						 multi1d<LatticeColorMatrix> &U,
						 std::string output_stem, 
						 int Ls, Real mass, Real M5, int Ncg, Real residual[], int max_iter[]);

// Internal interfaces
template int dwf_restarted_invert<float,float> (inv_type_t inv,
						multi1d<LatticeFermion> &sol, 
						multi1d<LatticeFermion> &src, 
	       					multi1d<multi1d<LatticeFermion> > &evecs,
	       					multi1d<Complex> &evals,
						multi1d<LatticeColorMatrix> &U,
						std::string output_stem, 
						int Ls, Real mass, Real M5, int Ncg, Real residual[], int max_iter[]);

// Internal interfaces
template int dwf_restarted_invert<double,double> (inv_type_t inv,
						  multi1d<LatticeFermion> &sol, 
						  multi1d<LatticeFermion> &src,
	       					  multi1d<multi1d<LatticeFermion> > &evecs,
	       					  multi1d<Complex> &evals, 
						  multi1d<LatticeColorMatrix> &U,
						  std::string output_stem, 
						  int Ls, Real mass, Real M5, int Ncg, Real residual[], int max_iter[]);

template 
int dwf_CG<float,double>
           (LatticePropagator &sol, 
	    LatticePropagator &src,
	    multi1d<LatticeColorMatrix> &U,
	    std::string output_stem,
	    int Ls,
	    int do5d,
	    Real mass,
	    Real M5,int Ncg,Real residual[],int max_iter[]);

template 
int dwf_CG<float,double>
           (LatticePropagator &sol, 
	    LatticePropagator &src,
	    multi1d<multi1d<LatticeFermion> > &evecs,
	    multi1d<Complex> &evals,
	    multi1d<LatticeColorMatrix> &U,
	    std::string output_stem,
	    int Ls,
	    int do5d,
	    Real mass,
	    Real M5,int Ncg,Real residual[],int max_iter[]);

template 
int dwf_CG<double,double>
           (LatticePropagator &sol, 
	    LatticePropagator &src,
	    multi1d<LatticeColorMatrix> &U,
	    std::string output_stem,
	    int Ls,
	    int do5d,
	    Real mass,
	    Real M5,int Ncg,Real residual[],int max_iter[]);

template 
int dwf_CG<double,double>
           (LatticePropagator &sol, 
	    LatticePropagator &src,
	    multi1d<multi1d<LatticeFermion> > &evecs,
	    multi1d<Complex> &evals,
	    multi1d<LatticeColorMatrix> &U,
	    std::string output_stem,
	    int Ls,
	    int do5d,
	    Real mass,
	    Real M5,int Ncg,Real residual[],int max_iter[]);

template 
int dwf_CG<float,float>
          (LatticePropagator &sol, 
	   LatticePropagator &src,
	   multi1d<LatticeColorMatrix> &U,
	   std::string output_stem,
	   int Ls,
	   int do5d,
	   Real mass,
	   Real M5,int Ncg,Real residual[],int max_iter[]);

template 
int dwf_CG<float,float>
          (LatticePropagator &sol, 
	   LatticePropagator &src,
	   multi1d<multi1d<LatticeFermion> > &evecs,
	   multi1d<Complex> &evals,
	   multi1d<LatticeColorMatrix> &U,
	   std::string output_stem,
	   int Ls,
	   int do5d,
	   Real mass,
	   Real M5,int Ncg,Real residual[],int max_iter[]);

/* FermToProp and PropToFerm lifted from CHROMA - sorry guys*/
/* You lifted my routines in what you call bagel-qdp anyway... touche*/


  //! Extract a LatticeFermion from a LatticePropagator
  /*!
   * \ingroup ferm
   *
   * \param a      Source propagator
   * \param b      Destination fermion
   * \param color_index  Color index
   * \param spin_index   Spin index
   */
  void BfmPropToFerm(const LatticePropagatorF& b, LatticeFermionF& a, 
		  int color_index, int spin_index)
  {
    for(int j = 0; j < Ns; ++j)
    {
      LatticeColorMatrixF bb = peekSpin(b, j, spin_index);
      LatticeColorVectorF aa = peekSpin(a, j);

      for(int i = 0; i < Nc; ++i)
	pokeColor(aa, peekColor(bb, i, color_index), i);

      pokeSpin(a, aa, j);
    }
  }

  //! Extract a LatticeFermion from a LatticePropagator
  /*!
   * \ingroup ferm
   *
   * \param a      Source propagator
   * \param b      Destination fermion
   * \param color_index  Color index
   * \param spin_index   Spin index
   */
  void BfmPropToFerm(const LatticePropagatorD& b, LatticeFermionD& a, 
		  int color_index, int spin_index)
  {
    for(int j = 0; j < Ns; ++j)
    {
      LatticeColorMatrixD bb = peekSpin(b, j, spin_index);
      LatticeColorVectorD aa = peekSpin(a, j);

      for(int i = 0; i < Nc; ++i)
	pokeColor(aa, peekColor(bb, i, color_index), i);

      pokeSpin(a, aa, j);
    }
  }


  void BfmFermToProp(const LatticeFermionF& a, LatticePropagatorF& b, 
		  int color_index, int spin_index)
  {
    for(int j = 0; j < Ns; ++j)
    {
      LatticeColorMatrixF bb = peekSpin(b, j, spin_index);
      LatticeColorVectorF aa = peekSpin(a, j);

      for(int i = 0; i < Nc; ++i)
	pokeColor(bb, peekColor(aa, i), i, color_index);

      pokeSpin(b, bb, j, spin_index);
    }
  }

  //! Insert a LatticeFermion into a LatticePropagator
  /*!
   * \ingroup ferm
   *
   * \param a      Source fermion
   * \param b      Destination propagator
   * \param color_index  Color index
   * \param spin_index   Spin index
   */
  void BfmFermToProp(const LatticeFermionD& a, LatticePropagatorD& b, 
		  int color_index, int spin_index)
  {
    for(int j = 0; j < Ns; ++j)
    {
      LatticeColorMatrixD bb = peekSpin(b, j, spin_index);
      LatticeColorVectorD aa = peekSpin(a, j);

      for(int i = 0; i < Nc; ++i)
	pokeColor(bb, peekColor(aa, i), i, color_index);

      pokeSpin(b, bb, j, spin_index);
    }
  }
