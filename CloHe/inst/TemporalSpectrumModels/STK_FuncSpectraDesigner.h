/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2016  Serge Iovleff, Université Lille 1, Inria

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2 of the
 License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this program; if not, write to the
 Free Software Foundation, Inc.,
 59 Temple Place,
 Suite 330,
 Boston, MA 02111-1307
 USA

 Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 */

/*
 * Project:  stkpp::CloHe
 * created on: 8 août 2016
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_FuncSpectraDesigner.h
 *  @brief In this file we define the STK_FuncSpectraDesigner.h class.
 **/

#ifndef STK_FUNCSPECTRADESIGNER_H
#define STK_FUNCSPECTRADESIGNER_H

#include <Regress.h>
#include "STK_BSplineWCoefficients.h"

namespace STK
{

template<int Size_>
class FuncSpectraDesigner
{
  public:
    // spectra container (CArrayPoint<double, size_>)
    typedef typename hidden::ClassifTraits<Size_ >::Spectra Spectra;

    // array of spectra (ArraySp)
    typedef typename hidden::ClassifTraits<Size_ >::ArraySp ArraySp;
    // array of values during time for a given pixel (CArrayVector<Spectra>)
    typedef typename hidden::ClassifTraits<Size_ >::SeriesSpectra SeriesSpectra;
    // 3D array of all values (Array1D<SeriesSpectra>)
    typedef typename hidden::ClassifTraits<Size_ >::ArraySeriesSpectra ArraySeriesSpectra;

    enum
    {
      size_ =  hidden::ClassifTraits<Size_ >::size_,
    };
    /** constructor.
     * @param times array with the sampling times of all the observations
     * @param labels array with the labels of all the observations
     * @param kernelName name of the kernel to use for the covariance matrix
     * @param h width to use in the kernel
     * @param dim number of basis to use
     * @param degree degree of the BSpline to use
     * @param posKnots way to put the knots
     **/
    FuncSpectraDesigner( Times const& times
                       , Labels const& labels
                       , String const& kernelName, Real const& h
                       , String const& basisName, CPointXi const& dim
                       , int degree, String const& posKnots);
    /** destructor. */
    ~FuncSpectraDesigner();

    /** @return array with the knots of each class */
    inline Array1D< VectorX > const&  knots() const { return knots_;}
    /** @return array with the sampling times for each classes */
    inline Array1D< VectorX > const&  samplingTimes() const { return samplingTimes_;}
    /** @return  vector of dimensions for each class */
    inline CPointXi const& dim() const { return dim_;}
    /** @return  array with the coefficients of each class */
    inline Array1D< ArrayXX > const& coefficients() const { return coefficients_;}
    /**  @return array with the first date of each class */
    inline Array1D< Real > const&  tmin() const { return tmin_;}
    /**  @return array with the last date of each class */
    inline Array1D< Real > const& tmax() const { return tmax_;}

    /** compute the design matrix for an individual (a pixel) using the
     *  sampling times of the spectra, the inverse square roots of the
     *  covariance matrix and weight the resulting matrices by the inverse
     *  square root of the covariance matrix.
     *  @param tt the vector with the date of the observation
     *  @param xt,yt the weighted design matrix xt = \Sigma_i^{-1/2} B^i and
     *  the weighted vector of the observations yt = \Sigma_i^{-1/2} y^i
     *  @param k number of the class
     *  @return determinant of the covariance matrix
     **/
    Real computeDesign(CPointX const& tt, ArrayXX& xt,  ArraySp& yt, int k) const;
    /** compute safely the design matrix for an individual (a pixel) using the
     *  sampling times of the spectra, the inverse square roots of the
     *  covariance matrix and weight the resulting matrices by the inverse
     *  square root of the covariance matrix.
     *  @param tt the vector with the date of the observation
     *  @param xt,yt the weighted design matrix xt = \Sigma_i^{-1/2} B^i and
     *  the weighted vector of the observations yt = \Sigma_i^{-1/2} y^i
     *  @param k class number
     *  @return determinant of the covariance matrix
     **/
    Real computeSafeDesign(CPointX const& tt, ArrayXX& xt,  ArraySp& yt, int k) const;
    /** compute the design matrix for an individual (a pixel) using the
     *  sampling times of the spectra
     *  @param tt the vector with the date of the individual
     *  @param k class number
     *  @return the design matrix for class k associated to dates tt
     **/
    ArrayXX getDesignMatrix(CPointX const& tt, int k) const;
    /** compute safely the design matrix for an individual (a pixel) using the
     *  sampling times of the spectra
     *  @param tt the vector with the date of the individual
     *  @param k class number
     *  @return the design matrix for class k associated to dates tt
     **/
    ArrayXX getSafeDesignMatrix(CPointX const& tt, int k) const;

  protected:
    /** constant reference on the sampled days */
    Times const& times_;
    /** constant reference on the labels */
    Labels const& labels_;
    // type of basis to use
    Basis::TypeBasisFunction typeBasis_;
    /** vector of dimensions for each class */
    CPointXi const& dim_;
    /** array with the coefficients for each classes */
    Array1D< VectorX > knots_;
    /** array with the sampling times for each classes */
    Array1D< VectorX > samplingTimes_;
    /** array with the coefficients for each classes */
    Array1D< ArrayXX > coefficients_;
    /** first time */
    Array1D< Real > tmin_;
    /** last time */
    Array1D< Real > tmax_;

  private:
    /** Vector with all kernel values (366 values) */
    CVectorX allKernel_;
    /** build the vector allKernel_ */
    void buildKernel(String const& kernelName, Real const& h);
    /** build the vector of sampling times */
    void buildSamplingTimes();
    /** build the vector of coefficients_ matrices */
    void buildCoefficients();
    /** build the vector of coefficients_ matrices using Bspline basis
     *  @param degree degree of the BSpline basis
     *  @param posKnots method for positioning the knots (density, periodic or uniform)
     **/
    void buildCoefficients( int degree, Regress::KnotsPosition posKnots);
};

/* Default constructor */
template<int Size_>
FuncSpectraDesigner<Size_>::FuncSpectraDesigner( Times const& times
                                               , Labels const& labels
                                               , String const& kernelName, Real const& h
                                               , String const& basisName, CPointXi const& dim
                                               , int degree, String const& posKnots
                                               )
                                               : times_(times)
                                               , labels_(labels)
                                               , typeBasis_(Basis::stringToTypeBasisFunction(basisName))
                                               , dim_(dim)
                                               , knots_(dim_.range()), coefficients_(dim_.range())
                                               , samplingTimes_(dim_.range())
                                               , tmin_(dim_.range()), tmax_(dim_.range())
                                               , allKernel_( Range(0,366) )
{
#ifdef STK_CLOHE_DEBUG
  stk_cout << _T("Entering FuncSpectraDesigner<Size_>::FuncSpectraDesigner\n");
#endif
  // create kernel
  buildKernel(kernelName, h);
#ifdef STK_CLOHE_DEBUG
  stk_cout << _T("build Kernel terminated\n");
  stk_cout << " allKernel_ =\n" << allKernel_.transpose();
#endif
  // compute sampling times
  buildSamplingTimes();
#ifdef STK_CLOHE_DEBUG
  stk_cout << _T("Sampling times computed\n");
  for (int k=samplingTimes_.begin(); k<samplingTimes_.end(); ++k)
  {
    stk_cout << _T("classe :") << k << _T("\n");
    stk_cout << samplingTimes_[k].transpose();
  }
  #endif
  // create design matrices
  if (typeBasis_ == Basis::bspline_)
  { buildCoefficients( degree, Regress::stringToKnotsPosition(posKnots));}
  else
  { buildCoefficients();}

#ifdef STK_CLOHE_DEBUG
  stk_cout << _T("FuncSpectraDesigner<Size_>::FuncSpectraDesigner constructor terminated\n");
#endif
}
/* destructor */
template<int Size_>
FuncSpectraDesigner<Size_>::~FuncSpectraDesigner(){}

/* compute the design matrix for an individual (a pixel) using the
 *  sampling times of the spectra, the inverse square roots of the
 *  covariance matrix and weight the resulting matrices by the inverse
 *  square root of the covariance matrix.
 *  @param tt the vector with the date of the observation
 *  @param xt,yt the weighted design matrix xt = \Sigma_i^{-1/2} B^i
 *  @return determinant of the covariance matrix
 **/
template<int Size_>
Real FuncSpectraDesigner<Size_>::computeDesign(CPointX const& tt, ArrayXX& xt, ArraySp& yt, int k) const
{
  // get design matrix xt
  xt.move(getDesignMatrix(tt, k));
  // build covariance matrix with kernel
  CSquareX c(tt.range());
  for(int t=tt.begin(); t< tt.end(); ++t)
  {
    c(t,t) = allKernel_[0];
    for(int u=t+1; u< tt.end(); ++u)
    { c(t,u) = allKernel_[(int)std::abs(tt[u] - tt[t])];}
  }
  /* compute the eigenvalues decomposition of c */
  lapack::SymEigen<CSquareX> decomp_(c);
  decomp_.run();

  xt = decomp_.ginvsqrt(c) * xt; // c will be overwritten by ginvsqrt()
  yt = c * yt;
  return decomp_.det();
}

/* compute the design matrix for an individual (a pixel) using the
 *  sampling times of the spectra, the inverse square roots of the
 *  covariance matrix and weight the resulting matrices by the inverse
 *  square root of the covariance matrix.
 *  @param tt the vector with the date of the observation
 *  @param xt,yt the weighted design matrix xt = \Sigma_i^{-1/2} B^i
 *  @return determinant of the covariance matrix
 **/
template<int Size_>
Real FuncSpectraDesigner<Size_>::computeSafeDesign(CPointX const& tt, ArrayXX& xt, ArraySp& yt, int k) const
{
  // get design matrix xt
  xt.move(getSafeDesignMatrix(tt, k));
  // build covariance matrix with kernel
  CSquareX c(tt.range());
  for(int t=tt.begin(); t< tt.end(); ++t)
  {
    c(t,t) = allKernel_[0];
    for(int u=t+1; u< tt.end(); ++u)
    { c(t,u) = allKernel_[(int)std::abs(tt[u] - tt[t])];}
  }
  /* compute the eigenvalues decomposition of c */
  lapack::SymEigen<CSquareX> decomp_(c);
  decomp_.run();

  xt = decomp_.ginvsqrt(c) * xt; // c will be overwritten by ginvsqrt()
  yt = c * yt;
  return decomp_.det();
}

/* compute the design matrix for an individual (a pixel) using the
 *  sampling times of the spectrum
 *  @param tt the vector with the date of the observation
 *  @param xt the design matrix
 *  @return the design matrix
 **/
template<int Size_>
ArrayXX FuncSpectraDesigner<Size_>::getDesignMatrix(CPointX const& tt, int k) const
{
  //stk_cout << "tt=" << tt;
  ArrayXX xt(tt.range(), coefficients_[k].cols());
  for (int t=tt.begin(); t<tt.end(); ++t)
  { xt.row(t) = coefficients_[k].row( (int)tt[t] );}
  return xt;
}

/* compute the design matrix for an individual (a pixel) using the
 *  sampling times of the spectrum
 *  @param tt the vector with the date of the observation
 *  @param xt the design matrix
 *  @return the design matrix
 **/
template<int Size_>
ArrayXX FuncSpectraDesigner<Size_>::getSafeDesignMatrix(CPointX const& tt, int k) const
{
  ArrayXX xt(tt.range(), coefficients_[k].cols());
  int imin = coefficients_[k].beginRows(), imax = coefficients_[k].endRows()-1;
  for (int t=tt.begin(); t<tt.end(); ++t)
  {
    int i = std::max( std::min(imax, (int)tt[t]), imin);
    xt.row(t) = coefficients_[k].row( i );
  }
  return xt;
}

/* build the vector allKernel_ */
template<int Size_>
void FuncSpectraDesigner<Size_>::buildKernel(String const& kernelName, Real const& h)
{
  // create kernel
  Kernel::kernelType kind = Kernel::stringToKernelType(kernelName);
  Kernel::IKernelBase<CPointX>* p_kernel;
  switch (kind)
  {
    case Kernel::gaussian_:
      p_kernel = new Kernel::Gaussian<CPointX>(h);
      break;
    case Kernel::laplace_:
      p_kernel = new Kernel::Laplace<CPointX>(h);
      break;
    case STK::Kernel::rationalQuadratic_:
      p_kernel = new Kernel::RationalQuadratic<CPointX>(h);
      break;
    default:
      STKRUNTIME_ERROR_1ARG(FuncSpectraDesigner,kernelName,invalid kernel name);
      break;
  }
  for (int d= allKernel_.begin(); d < allKernel_.end(); ++d)
  { allKernel_[d] = p_kernel->value((Real)d);}
  if (p_kernel) delete p_kernel;
}

/* build the vector of coefficients_ */
template<int Size_>
void FuncSpectraDesigner<Size_>::buildSamplingTimes()
{
  // count the sampled days
  CArrayXXi count( samplingTimes_.range(), Range(0,366), 0 );
  for (int i = times_.begin(); i < times_.end(); ++i)
  {
    for(int t = times_[i].begin(); t < times_[i].end(); ++t)
    { count(labels_[i], (int)times_[i][t] )++;}
  }
  // build sampling times
  for (int k=samplingTimes_.begin(); k<samplingTimes_.end(); ++k)
  {
    samplingTimes_[k].resize(Range(0,366));
    for (int i = count.lastIdxCols(); i>= count.beginCols(); --i)
    {
      if (count(k, i) == 0) { samplingTimes_[k].erase(i);}
      else                  { samplingTimes_[k][i] = i;}

    }
  }
}

/* build the vector of coefficients_ */
template<int Size_>
void FuncSpectraDesigner<Size_>::buildCoefficients()
{
  IBasis<VectorX>* p_basis;
  for (int k=dim_.begin(); k< dim_.end(); ++k)
  {
    switch (typeBasis_)
    {
      case Basis::chebyshev_:
        p_basis = new ChebyshevCoefficients<VectorX>(samplingTimes_[k], dim_[k]);
        break;
      case Basis::sines_:
        p_basis = new SinesCoefficients<VectorX>(samplingTimes_[k], dim_[k]);
        break;
      case Basis::cosines_:
        p_basis = new CosinesCoefficients<VectorX>(samplingTimes_[k], dim_[k]);
        break;
      case Basis::trigonometric_:
        p_basis = new TrigonometricCoefficients<VectorX>(samplingTimes_[k], dim_[k]);
        break;
      default:
        STKRUNTIME_ERROR_NO_ARG(FuncSpectraDesigne::buildCoefficients,invalid typeBasis_);
        break;
    }
    p_basis->run();
    coefficients_[k] = p_basis->coefficients();
    tmin_[k] = p_basis->minValue();
    tmax_[k] = p_basis->maxValue();
    delete p_basis;
  }
#ifdef STK_CLOHE_DEBUG
  stk_cout << _T("Coefficients computed\n");
  for (int k=coefficients_.begin(); k<coefficients_.end(); ++k)
  {
    stk_cout << _T("classe :") << k << _T("\n");
    stk_cout << coefficients_[k].row(0);
  }

#endif
}
/* build the vector of coefficients_ */
template<int Size_>
void FuncSpectraDesigner<Size_>::buildCoefficients( int degree, Regress::KnotsPosition posKnots
                                                  )
{
  for (int k=dim_.begin(); k< dim_.end(); ++k)
  {
    BSplineWCoefficients builder(times_, labels_, dim_[k], degree,  posKnots);
    builder.setClass(k);
    builder.run();
    knots_[k] = builder.knots();
    coefficients_[k] = builder.coefficients();
    tmin_[k] = builder.tmin();
    tmax_[k] = builder.tmax();
  }
}


} // namespace STK

#endif /* STK_FUNCSPECTRADESIGNER_H */
