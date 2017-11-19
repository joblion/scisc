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

/** @file STK_FuncSpectraDesign.h
 *  @brief In this file we define the STK_FuncSpectraDesign.h class.
 **/

#ifndef STK_FUNCSPECTRADESIGN_H
#define STK_FUNCSPECTRADESIGN_H

#include "STK_BSplineWCoefficients.h"

namespace STK
{

template<int Size_>
class FuncSpectraDesign
{
  public:
    // spectra container (CArrayPoint<double, size_>)
    typedef typename hidden::ClassifTraits<Size_ >::Spectra Spectra;

    // array of spectra (ArraySpectra)
    typedef typename hidden::ClassifTraits<Size_ >::ArraySpectra ArraySpectra;
    // array of values during time for a given pixel (CArrayVector<Spectra>)
    typedef typename hidden::ClassifTraits<Size_ >::SeriesSpectra SeriesSpectra;
    // 3D array of all values (Array1D<SeriesSpectra>)
    typedef typename hidden::ClassifTraits<Size_ >::ArraySeriesSpectra ArraySeriesSpectra;

    enum
    {
      size_ =  hidden::ClassifTraits<Size_ >::size_,
    };
    /** constructor.
     * @param t array with the sampling times of all the observations
     * @param kernelName name of the kernel to use for the covariance matrix
     * @param h width to use in the kernel
     * @param dim number of basis to use
     * @param degree degree of the BSPline to use
     * @param posKnots way to put the knots
     **/
    FuncSpectraDesign( Times const& t
                     , String const& kernelName, Real const& h
                     , int dim, int degree, String const& posKnots);
    /** constructor.
     * @param t array with the sampling times of all the observations
     * @param kernelName name of the kernel to use for the covariance matrix
     * @param h width to use in the kernel
     * @param dim number of basis to use
     * @param degree degree of the BSPline to use
     * @param posKnots way to put the knots
     **/
    FuncSpectraDesign( Times const& t
                     , String const& kernelName, Real const& h
                     , CPointXi const& dim, int degree, String const& posKnots);
    /** destructor. */
    ~FuncSpectraDesign();
    /** @return a reference on the array of knots */
    VectorX const& knots() const { return builder_.knots();}
    /** @return a reference on the array of coefficients */
    ArrayXX const& coefficients() const { return builder_.coefficients();}
    /** compute the design matrix for an individual (a pixel) using the
     *  sampling times of the spectra, the inverse square roots of the
     *  covariance matrix and weight the resulting matrices by the inverse
     *  square root of the covariance matrix.
     *  @param tt the vector with the date of the observation
     *  @param xt,yt the weighted design matrix xt = \Sigma_i^{-1/2} B^i and
     *  the weighted vector of the observations yt = \Sigma_i^{-1/2} y^i
     *  @return determinant of the covariance matrix
     **/
    Real computeDesign(CPointX const& tt, ArrayXX& xt,  ArraySpectra& yt) const;
    /** compute the design matrix for an individual (a pixel) using the
     *  sampling times of the spectra
     *  @param tt the vector with the date of the observation
     *  @param xt the design matrix
     *  @return the design matrix xt
     **/
    ArrayXX getDesignMatrix(CPointX const& tt) const;

  private:
    /** Vector with all kernel values (366 values) */
    CVectorX allKernel_;
    /** builder of the design matrix */
    BSplineWCoefficients builder_;
};

/* Default constructor */
template<int Size_>
FuncSpectraDesign<Size_>::FuncSpectraDesign( Times const& times
                                           , String const& kernelName, Real const& h
                                           , int dim, int degree, String const& posKnots
                                           )
                                           : allKernel_( Range(0,366) )
                                           , builder_(times, dim, degree,  Regress::stringToKnotsPosition(posKnots))
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
      STKRUNTIME_ERROR_1ARG(FuncSpectraDesign,kernelName,invalid kernel name);
      break;
  }
  for (int d= allKernel_.begin(); d < allKernel_.end(); ++d)
  { allKernel_[d] = p_kernel->value((Real)d);}
  if (p_kernel) delete p_kernel;

  // create BSpline full basis
  // sample all the year
  CPointX allT( 366 );
  for (int i = allT.begin(); i < allT.end(); ++i) { allT[i] = i;}

  // builder of the design matrix
  builder_.setData(allT);
  builder_.run();

#ifdef STK_CLOHE_DEBUG
  stk_cout << _T("FuncSpectraDesign<Size_>::FuncSpectraDesign constructor terminated\n");
  stk_cout << " allKernel_ =\n" << allKernel_.transpose();
  stk_cout << " builder_.coefficients() =\n" << builder_.coefficients();
#endif
}
/* Default constructor */
template<int Size_>
FuncSpectraDesign<Size_>::FuncSpectraDesign( Times const& times
                                           , String const& kernelName, Real const& h
                                           , CPointXi const& dim, int degree, String const& posKnots
                                           )
                                           : allKernel_( Range(0,366) )
                                           , builder_(times, dim.front(), degree,  Regress::stringToKnotsPosition(posKnots))
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
      STKRUNTIME_ERROR_1ARG(FuncSpectraDesign,kernelName,invalid kernel name);
      break;
  }
  for (int d= allKernel_.begin(); d < allKernel_.end(); ++d)
  { allKernel_[d] = p_kernel->value((Real)d);}
  if (p_kernel) delete p_kernel;

  // create BSpline full basis
  // sample all the year
  CPointX allT( 366 );
  for (int i = allT.begin(); i < allT.end(); ++i) { allT[i] = i;}

  // builder of the design matrix
  builder_.setData(allT);
  builder_.run();

#ifdef STK_CLOHE_DEBUG
  stk_cout << _T("FuncSpectraDesign<Size_>::FuncSpectraDesign constructor terminated\n");
  stk_cout << " allKernel_ =\n" << allKernel_.transpose();
  stk_cout << " builder_.coefficients() =\n" << builder_.coefficients();
#endif
}
/* destructor */
template<int Size_>
FuncSpectraDesign<Size_>::~FuncSpectraDesign(){}

/* compute the design matrix for an individual (a pixel) using the
 *  sampling times of the spectra, the inverse square roots of the
 *  covariance matrix and weight the resulting matrices by the inverse
 *  square root of the covariance matrix.
 *  @param tt the vector with the date of the observation
 *  @param xt,yt the weighted design matrix xt = \Sigma_i^{-1/2} B^i
 *  @return determinant of the covariance matrix
 **/
template<int Size_>
Real FuncSpectraDesign<Size_>::computeDesign(CPointX const& tt, ArrayXX& xt, ArraySpectra& yt) const
{
  // get design matrix xt
  xt.move(getDesignMatrix(tt));
  // build covariance matrix with kernel
  CSquareX c(tt.range());
  for(int t=tt.begin(); t< tt.end(); ++t)
  {
    for(int u=t; u< tt.end(); ++u)
    { c(t,u) = allKernel_[(int)std::abs(tt[t] - tt[u])];}
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
ArrayXX FuncSpectraDesign<Size_>::getDesignMatrix(CPointX const& tt) const
{
  //stk_cout << "tt=" << tt;
  ArrayXX xt(tt.range(), builder_.coefficients().cols());
  for (int t=tt.begin(); t<tt.end(); ++t)
  { xt.row(t) = builder_.coefficients().row( (int)tt[t] );}
  return xt;
}


} // namespace STK

#endif /* STK_FUNCSPECTRADESIGN_H */
