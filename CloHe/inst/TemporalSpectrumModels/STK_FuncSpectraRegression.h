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
 * Project:  stkpp::Regress
 * created on: 27 oct. 2010
 * Purpose:  .
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_FuncSpectraRegression.h
 *  @brief In this file we define the FuncSpectraRegression class.
 **/

#ifndef STK_FUNCSPECTRAREGRESSION_H
#define STK_FUNCSPECTRAREGRESSION_H

#include <Arrays/include/STK_Array2D.h>
#include <Algebra/include/STK_InvertSymMatrix.h>
#include <Regress/include/STK_IRegression.h>

namespace STK
{

/** @brief The @c FuncSpectraRegression class allows to regress a multidimensional
 *  output variable among a multivariate explanation variable.
 */
template<int Size_, class Weights_>
class FuncSpectraRegression: public IRegression< CArray<Real, UnknownSize, Size_>, CArrayXX, Weights_>
                           , public IStatModelBase
{
  public:
    typedef CArrayPoint<Real, Size_> Spectra;
    typedef CArray<Real, UnknownSize, Size_> YArray;
    typedef CArrayXX XArray;
    typedef IRegression<YArray, XArray, Weights_> Base;

    using Base::p_x_;
    using Base::p_y_;
    using Base::predicted_;
    using Base::residuals_;

    /** Default constructor. */
    FuncSpectraRegression();
    /** Constructor.
     * @param y,x Variates to explain and co-variates
     */
    FuncSpectraRegression( YArray const& y, XArray const& x);
    /** Destructor. */
    virtual ~FuncSpectraRegression() {}
    /** @return the coefficients */
    inline CArrayXX const& coefs() const { return coefs_;}
    /** @return the variance of the residuals */
    inline Spectra const& sigma2() const { return sigma2_;}
    /** @return the variance of the residuals */
    inline int nbFreeParameter() const { return IStatModelBase::nbFreeParameter();}
    /** @return the extrapolated values y from the value @c x.
     *  Given the data set @c x will compute the values \f$ y = x.\hat{\beta} \f$.
     *  The coefficients @c coefs_ have to be estimated previously.
     *  @param x the input data set
     */
    virtual YArray extrapolate(XArray const& x) const;

  protected:
    CArrayXX coefs_;
    Spectra sigma2_;

  private:
    /** compute the regression function. */
    virtual bool regressionStep();
    /** compute the weighted regression function.
     * @param weights the weights of the samples
     **/
    virtual bool regressionStep(Weights_ const& weights);
    /** Compute the predicted outputs by the regression function. */
    virtual bool predictionStep();
    /** Compute the residual variance and the likelihood of a Gaussian model. */
    virtual bool finalizeStep();
    /** Compute the number of parameter of the regression function.
     * @return the number of parameter of the regression function
     **/
    inline virtual int computeNbFreeParameter() const
    { return coefs_.sizeCols() * coefs_.sizeRows(); }
};

/* default constructor */
template<int Size_, class Weights_>
FuncSpectraRegression<Size_,Weights_>::FuncSpectraRegression(): Base(), IStatModelBase(), coefs_()
{}

template<int Size_, class Weights_>
FuncSpectraRegression<Size_,Weights_>::FuncSpectraRegression( YArray const& y, XArray const& x)
                                                            : Base(y, x)
                                                            , IStatModelBase(x.sizeRows(), x.sizeCols())
                                                            , coefs_()
{}

/* compute the regression function. */
template<int Size_, class Weights_>
bool FuncSpectraRegression<Size_,Weights_>::regressionStep()
{
  // compute a'a
  CSquareX prod = p_x_->transpose() * *p_x_;
  // compute (a'a)^{-1}
  lapack::SymEigen<CSquareX> decomp(prod);
  decomp.run();
  // compute (a'a)^{-1}b'a
  if (p_x_->sizeRows() < p_y_->sizeCols())
  { coefs_ = ((decomp.ginv(prod) * p_x_->transpose()) * *p_y_);}
  else
  { coefs_ = (decomp.ginv(prod) * (p_x_->transpose() * *p_y_));}
  return true;
//  Range arows = p_x_->rows(), acols = p_x_->cols();
//  Range brows = p_y_->rows(), bcols = p_y_->cols();
//  // auxiliary arrays
//  CArrayXX a = *p_x_, b = *p_y_;
//  CArrayXX b(brows.size()<acols.size() ? acols : brows, bcols);
//  b.sub(brows, bcols) = *p_y_;
//  int m   = a.sizeRows(), n= a.sizeCols(), nrhs = b.sizeCols();
//  int lda = a.sizeRows(), ldb = b.sizeRows();
//  CVectorX s(arows, 0.);
  // shift to 0
//  a.shift(0,0);
//  b.shift(0,0);
//  s.shift(0);
//  // get sizes
//  Real wkopt, rcond_ = -1;
//  int iwkopt, rank_,lWork = -1;
//  int info = lapack::gelsd( m, n, nrhs
//                          , a.p_data(), lda
//                          , b.p_data(), ldb
//                          , s.p_data()
//                          , &rcond_, &rank_
//                          , &wkopt, lWork, &iwkopt
//                          );
//  if( info != 0 )
//  {
//    if (info>0)
//    { this->msg_error_ = STKERROR_NO_ARG(FuncSpectraRegression::regressionStep,internal error);
//      return false;
//    }
//    this->msg_error_= STKERROR_1ARG(FuncSpectraRegression::regressionStep get,-info,error parameter);
//    return false;
//  }
//  // get working sizes
//  lWork = (int)wkopt;
//  Real* work = new Real[lWork];
//  int* iwork = new int[iwkopt];
//  // Solve the least square problem
//  info = lapack::gelsd( m, n, nrhs
//                      , a.p_data(), lda
//                      , b.p_data(), ldb
//                      , s.p_data()
//                      , &rcond_, &rank_
//                      , work, lWork, iwork
//                      );
//  delete[] iwork;
//  delete[] work;
//  if( info != 0 )
//  {
//    if (info>0)
//    { this->msg_error_ = STKERROR_NO_ARG(FuncSpectraRegression::regressionStep,internal error);
//      return false;
//    }
//    this->msg_error_= STKERROR_1ARG(FuncSpectraRegression::regressionStep get,-info,error parameter);
//    return false;
//  }
//  // shift back and get coefficients
//  b.shift(p_y_->beginRows(), p_y_->beginCols());
//  coefs_ = b.sub(acols,bcols);
//  return true;
}

/* compute the regression function. */
template<int Size_, class Weights_>
bool FuncSpectraRegression<Size_, Weights_>::regressionStep(Weights_ const& weights)
{
  // compute X'WX
  ArraySquareX prod;
  prod.move(weightedMultLeftTranspose(p_x_->asDerived(), weights));
  // compute (X'WX)^{-1}
  GInvertSymMatrix()(prod);
  // compute (X'WX)^{-1}X'WY
  coefs_.move(mult(prod, wmultLeftTranspose(p_x_->asDerived(), p_y_->asDerived(), weights)));
  return true;
}

/* Compute the predicted outputs by the regression function. */
template<int Size_, class Weights_>
bool FuncSpectraRegression<Size_,Weights_>::predictionStep()
{
  predicted_ = (*p_x_) * coefs_;
  return true;
}
/* Compute the predicted outputs by the regression function. */
template<int Size_, class Weights_>
bool FuncSpectraRegression<Size_,Weights_>::finalizeStep()
{
  sigma2_= Stat::varianceByCol(residuals_, false);
  this->setNbFreeParameter(sigma2_.size() + coefs_.sizeArray());
  this->setLnLikelihood( nbSample() * sigma2_.log().sum()  // \sum T_i \log(\sigma^2_k)
                       + (residuals_*sigma2_.asDiagonal().inverse()).norm2());
  return true;
}


/* @brief Extrapolate the the values @c y from the value @c x.
 *  Given the data set @c x will compute the values \f$ y = \hat{f}(x) \f$.
 *  The regression function @e f have to be estimated previously.
 */
template<int Size_, class Weights_>
typename FuncSpectraRegression<Size_,Weights_>::YArray
  FuncSpectraRegression<Size_,Weights_>::extrapolate( XArray const& x) const
  { return(x * coefs_);}

} // namespace STK

#endif /* STK_FUNCSPECTRAREGRESSION_H */
