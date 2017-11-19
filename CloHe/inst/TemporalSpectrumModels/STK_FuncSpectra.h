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

/** @file STK_FuncSpectrum.h
 *  @brief In this file we define the FuncSpectra class.
 **/

#ifndef STK_FUNCSPECTRUM_H
#define STK_FUNCSPECTRUM_H

#ifdef _OPENMP
#include <omp.h>
#endif


#include <Classification/include/STK_IClassification.h>
#include "STK_CloHe_Util.h"
#include "STK_FuncSpectraRegression.h"
#include "STK_FuncSpectraDataHandler.h"
#include "STK_FuncSpectraDesign.h"

namespace STK
{
/** @brief This class allows to estimate the Func(mu_t,Sigma_t) model.
 *  @tparam Size_ the number of spectrum
 *  In each class, the y values are stored in arrays of size
 *  \f$ d_k= \left[\sum_{Z_i=k} T_i \right] \f$ by @c Size_.
 *
 */
template<int Size_, class Weights_>
class FuncSpectra: public IClassification< Labels
                                         , typename hidden::ClassifTraits<Size_ >::ArraySeriesSpectra
                                         , Weights_>
{
  public:
    // spectrum container (CArrayPoint<double, Size_>)
    typedef typename hidden::ClassifTraits<Size_ >::Spectra Spectra;

    // array of spectra (CArray<Real, UnknownSize, Size_>)
    typedef typename hidden::ClassifTraits<Size_ >::ArraySpectra ArraySpectra;
    // array of values during time for a given pixel (CArrayPoint<Spectra>)
    typedef typename hidden::ClassifTraits<Size_ >::SeriesSpectra SeriesSpectra;
    // 3D array of all values (Array1D< CArrayPoint<Spectra>  >)
    typedef typename hidden::ClassifTraits<Size_ >::ArraySeriesSpectra ArraySeriesSpectra;

    // Base class
    typedef IClassification< Labels, ArraySeriesSpectra, Weights_> Base;
    // data handler class
    typedef FuncSpectraDataHandler<Size_> Handler;
    enum
    {
      size_ =  hidden::ClassifTraits<Size_ >::size_,
    };

    /** typedef for the y values. Array (of array) of size K */
    typedef Array1D< ArraySpectra >  YArrays;
    /** typedef for the parameters values. */
    typedef Array1D< Spectra >   PArrays;

    using Base::nbFreeParameter_;

    /** Default constructor
     *  @param handler data handler
     *  @param dim number of base functions to use
     **/
    FuncSpectra( Handler const& handler
               , String const& kernelName, Real const& h
               , int dim, int degree, String const& posKnots
               , String const& criterion);
    /** Default constructor
     *  @param handler data handler
     *  @param dim number of base functions to use
     **/
    FuncSpectra( Handler const& handler
               , String const& kernelName, Real const& h
               , CPointXi const& dim, int degree, String const& posKnots
               , String const& criterion);
    /** destructor */
    virtual ~FuncSpectra() {}
    // data handler
    /** @return the number of sample */
    inline int nbSample() const { return handler_.nbSample();}
    /** @return the number of spectrum */
    inline int nbSpectrum() const { return handler_.nbSpectrum();};
    /** @return the number of class */
    inline int nbClass() const { return handler_.nbClass();};
    /** number of pixels in each classes */
    inline CPointXi const& nk() const { return handler_.nk();}
    /** number of available sampling for each sample */
    inline CPointXi const& ni() const { return handler_.ni();}
    /** @return labels of all pixels */
    Labels const& labels() const { return handler_.labels();}
    /** @return labels of pixel i*/
    inline int labels(int i) const { return handler_.labels()[i];}
    /** @return sampling times for all pixels */
    Times const& times() const { return handler_.times();}
    /** @return sampling times of pixel i */
    CPointX const& times(int i) const { return handler_.times()[i];}
    /** spectrum of the pixel i at time t */
    Spectra const& data(int i, int t) const { return handler_.data()[i][t];}
    // designer
    /** @return a reference on the array of knots */
    VectorX const& knots() const { return designer_.knots();}
	/** @return a reference on the array of coefficients */
    ArrayXX const& coefficients() const { return designer_.coefficients();}
    // getters
    /** @return regression coefficients */
    YArrays const& alpha() const {return alpha_;}
    /** @return estimated variance */
    PArrays const& sigma2() const { return sigma2_;}
    /** @return lnLikelihhod of each class */
    CPointX const& lnLikelihood() const { return lnLikelihood_;}
    /** @return criterion of each class */
    CPointX const& criterion() const { return criterion_;}
    /** @return criterion of each class */
    CPointXi const& nbFreeParameters() const { return nbFreeParameters_;}
    /** @return the determinant of the covariance matrix of each pixel */
    CVectorX const& det() const { return det_;}
    /** @return the prior probabilities of each class */
    CPointX const& pk() const { return pk_;}
    /** @return the posterior probabilities of each class */
    CArrayXX const& tik() const { return tik_;}
    /** @return a reference on the array of dimensions */
    CPointXi const& dim() const { return dim_;}
    /** @return the predicted labels of each pixel */
    Labels const& zi() const { return zi_;}

    /** classify the label given a new observation (t,y) */
    int classify( CPointX const& ti, ArraySpectra const& yi) const;

    /** compute the mean curve of each class */
    YArrays mean(Real first=1, Real last=365, int nbSampling=100) const;

  protected:
    /** reference on the data handler */
    Handler const& handler_;
    /** utility class allowing to design the matrices xi */
    FuncSpectraDesign<Size_> designer_;
    /** vector of dimensions for each class */
    CPointXi dim_;
    Model::TypeCriterion type_;

    /** predicted class after estimation */
    Labels zi_;

    /** coefficients array for each class */
    YArrays alpha_;
    /** variances for each class */
    PArrays sigma2_;
    /** log-Likelihood of each class */
    CPointX lnLikelihood_;
    /** criterion of each class */
    CPointX criterion_;
    /** Size of the design matrix in each class */
    CPointXi nbFreeParameters_;
    /** determinant of the covariance matrix for each pixels (nbSample) */
    CVectorX det_;
    /** prior probabilities */
    CPointX pk_;
    /** posterior probabilities */
    CArrayXX tik_;

  private:
    /** matrix X^tX */
    Array1D<CSquareX> xtx_jj_;
    /** matrix X^ty */
    Array1D< ArraySpectra > xty_js_;

    /** Size of the design matrix in each class */
    CVectorXi sizek_;

    /** Compute the design matrices using BSpline*/
    virtual bool initializeStep();
    /** Estimate the parameters */
    virtual bool estimationStep();
    /** Estimate the parameters */
    virtual bool finalizeStep();
    /** predict the posterior probabilities tik and class partnership zi */
    virtual bool predictionStep();

    virtual bool estimationStep( Weights_ const& weights)
    { /* Nothing to do, no weights available */
      return true;
    }
    /** compute the number of free parameters */
    virtual int computeNbFreeParameter() const;
    /** compute the matrices X^tX and X^tY */
    void computeXtXandXtY();
    /** compute the matrices predicted values, the residuals and the variances */
    void computeVariances();
};
/* default constructor */
template<int Size_, class Weights_>
FuncSpectra<Size_, Weights_>::FuncSpectra( Handler const& handler
                                         , String const& kernelName, Real const& h
                                         , int dim, int degree, String const& posKnots
                                         , String const& criterion
                                         )
                                         : Base(handler.labels(), handler.data())
                                         , handler_(handler)
                                         , designer_(handler.times(), kernelName, h, dim, degree, posKnots)
                                         , dim_(handler_.nbClass(), dim)
                                         , type_(Model::stringToTypeCriterion(criterion))
                                         , zi_(handler_.nbSample())
                                         , alpha_(handler_.nbClass())
                                         , sigma2_(handler_.nbClass())
                                         , lnLikelihood_(handler_.nbClass(), -Arithmetic<Real>::infinity())
                                         , criterion_(handler_.nbClass(), Arithmetic<Real>::infinity())
                                         , nbFreeParameters_(handler_.nbClass(), 0)
                                         , det_(handler_.nbSample())
                                         , pk_(handler_.nbClass())
                                         , tik_(handler_.nbSample(),handler_.nbClass())
                                         , sizek_(handler_.nbClass(), 0)
                                         , xtx_jj_(handler_.nbClass())
                                         , xty_js_(handler_.nbClass())
{
#ifdef STK_CLOHE_DEBUG
  stk_cout << _T("FuncSpectra<Size_, Weights_> created\n");
#endif
}

/* default constructor */
template<int Size_, class Weights_>
FuncSpectra<Size_, Weights_>::FuncSpectra( Handler const& handler
                                         , String const& kernelName, Real const& h
                                         , CPointXi const& dim, int degree, String const& posKnots
                                         , String const& criterion
                                         )
                                         : Base(handler.labels(), handler.data())
                                         , handler_(handler)
                                         , designer_(handler.times(), kernelName, h, dim, degree, posKnots)
                                         , dim_(dim)
                                         , type_(Model::stringToTypeCriterion(criterion))
                                         , zi_(handler_.nbSample())
                                         , alpha_(handler_.nbClass())
                                         , sigma2_(handler_.nbClass())
                                         , lnLikelihood_(handler_.nbClass(), -Arithmetic<Real>::infinity())
                                         , criterion_(handler_.nbClass(), Arithmetic<Real>::infinity())
                                         , nbFreeParameters_(handler_.nbClass(), 0)
                                         , det_(handler_.nbSample())
                                         , pk_(handler_.nbClass())
                                         , tik_(handler_.nbSample(),handler_.nbClass())
                                         , sizek_(handler_.nbClass(), 0)
                                         , xtx_jj_(handler_.nbClass())
                                         , xty_js_(handler_.nbClass())
{
#ifdef STK_CLOHE_DEBUG
  stk_cout << _T("FuncSpectra<Size_, Weights_> created\n");
#endif
}

/* compute the mean curve of each class */
template<int Size_, class Weights_>
typename FuncSpectra<Size_, Weights_>::YArrays
              FuncSpectra<Size_, Weights_>::mean( Real first, Real last, int nbSampling) const
{
  CPointX ti(nbSampling);
  YArrays mu(nbClass());
  Real delta = (last-first)/nbSampling, val = first;
  for (int i=ti.begin(); i<ti.end(); ++i) { ti[i] = val; val += delta;}
  //  mu[k] = xi*alpha_[k];
  for (int k= alpha_.begin(); k < alpha_.end(); ++k)
  { mu[k] = designer_.getDesignMatrix(ti)*alpha_[k];}
  return mu;
}

/* Compute the design matrices using BSpline*/
template<int Size_, class Weights_>
bool FuncSpectra<Size_, Weights_>::initializeStep()
{
#ifdef STK_CLOHE_DEBUG
  stk_cout << _T("Entering FuncSpectrum<Size_, Weights_>::initializeStep()\n");
#endif

  // compute proportions pk
  for (int k= pk_.begin(); k < pk_.end(); k++)
  {
    pk_[k] = Real(nk()[k])/Real(nbSample());
    xtx_jj_[k].resize(dim_[k]) = 0.;
    xty_js_[k].resize(dim_[k], nbSpectrum()) = 0.;
    sigma2_[k].resize(nbSpectrum()) = 0.;
  }

#ifdef STK_CLOHE_DEBUG
  stk_cout << _T("FuncSpectrum<Size_, Weights_>::initializeStep() terminated\n");
#endif
  return true;
}

/* Compute the design matrices using BSpline*/
template<int Size_, class Weights_>
bool FuncSpectra<Size_, Weights_>::estimationStep()
{
#ifdef STK_CLOHE_DEBUG
  stk_cout << _T("Entering FuncSpectrum<Size_, Weights_>::estimationStep()\n");
#endif

  bool flag = true;
  computeXtXandXtY();

  // this is the main loop over the classes
  // can be parallelized
#ifdef _OPENMP
//#pragma omp for schedule(dynamic,1)
#endif
  for (int k=pk_.begin(); k < pk_.end(); k++)
  {
    lapack::SymEigen<CSquareX> decomp(xtx_jj_[k]);
    decomp.run();
    // compute (x'x)^{-1}x'j
    alpha_[k] = decomp.ginv(xtx_jj_[k]) * xty_js_[k];
  }
  // compute residal variances
  computeVariances();
  return flag;

#ifdef STK_CLOHE_DEBUG
  stk_cout << _T("Terminating FuncSpectrum<Size_, Weights_>::estimationStep()\n");
#endif
}

/* @brief Compute the lnLikelihood of each class */
template<int Size_, class Weights_>
bool FuncSpectra<Size_, Weights_>::predictionStep()
{
#ifdef STK_CLOHE_DEBUG
  stk_cout << _T("Entering FuncSpectrum<Size_, Weights_>::predictionStep()\n");
#endif


  for (int i=zi_.begin(); i<zi_.end(); ++i)
  {
	ArrayXX      xt;
	ArraySpectra yt ;
    yt.move(handler_.getArraySpectra(i));
	Real det = designer_.computeDesign(times(i), xt, yt);

    CPointX lnComp(sigma2_.range());
    for (int k=sigma2_.begin(); k< sigma2_.end(); k++)
    {
      lnComp[k] = std::log(pk_[k])
                - 0.5*( (Stat::sumByCol( (yt - xt * alpha_[k]).square() )/sigma2_[k]).sum()
                       + std::log(det) + sigma2_[k].log().sum() * times(i).size()
                      );
    }
    int kmax;
    Real max = lnComp.maxElt(kmax);
    //  sum_k p_k exp{lnCom_k - lnComp_kmax}
    Real sum =  (tik_.row(i) = (lnComp - max).exp()).sum();
    tik_.row(i) /= sum;
    zi_[i] = kmax;
  }

#ifdef STK_CLOHE_DEBUG
  stk_cout << _T("Terminating FuncSpectrum<Size_, Weights_>::predictionStep()\n");
#endif
  return true;
}

template<int Size_, class Weights_>
bool FuncSpectra<Size_, Weights_>::finalizeStep()
{
#ifdef STK_CLOHE_DEBUG
  stk_cout << _T("Entering FuncSpectrum<Size_, Weights_>::finalizeStep()\n");
  stk_cout << _T("compute nbFreeParameters_, lnLikelihood_ and criterion_\n");
#endif

 for(int k=criterion_.begin(); k<criterion_.end(); ++k)
 {
   nbFreeParameters_[k] = alpha_[k].sizeArray() + sigma2_[k].size();
   lnLikelihood_[k]     = nk()[k] * (sigma2_[k].log().sum() + nbSpectrum())/2.;
   switch (type_)
   {
     case Model::aic_:
       criterion_[k] = -2. * lnLikelihood_[k] + nbFreeParameters_[k];
       break;
     case Model::bic_:
       criterion_[k] = -2. * lnLikelihood_[k] + nbFreeParameters_[k] * std::log(nk()[k]);
       break;
     case Model::ml_:
       criterion_[k] = -2*lnLikelihood_[k];
       break;
     default: // use BIC
       criterion_[k] = -2. * lnLikelihood_[k] + nbFreeParameters_[k] * std::log(nk()[k]);
       break;
   }
 }

#ifdef STK_CLOHE_DEBUG
  stk_cout << _T("Terminating FuncSpectrum<Size_, Weights_>::finalizeStep()\n");
#endif
  return true;
}

/* compute the number of free parameters */
template<int Size_, class Weights_>
int FuncSpectra<Size_, Weights_>::computeNbFreeParameter() const
{
  int sum = 0;
  for (int k= nbFreeParameters_.begin(); k < nbFreeParameters_.end(); ++k)
  { sum+=nbFreeParameters_[k];}
  return sum;
}

/* classify the label given a new observation (t,y) */
template<int Size_, class Weights_>
int FuncSpectra<Size_, Weights_>::classify( CPointX const& ti, ArraySpectra const& yi) const
{
  ArrayXX      xt;
  ArraySpectra yt = yi;
  Real det = designer_.computeDesign(ti, xt, yt);

  CPointX lnComp(sigma2_.range());
  for (int k=sigma2_.begin(); k< sigma2_.end(); k++)
  {
    lnComp[k] = std::log(pk_[k])
              - 0.5*( (Stat::sumByCol( (yt - xt * alpha_[k]).square() )/sigma2_[k]).sum()
                     + std::log(det) + sigma2_[k].log().sum() * ti.size()
                    );
  }
  int kmax;
  Real max = lnComp.maxElt(kmax);
  return kmax;

}

/* compute the matrices X^tX and X^tY */
template<int Size_, class Weights_>
void FuncSpectra<Size_, Weights_>::computeXtXandXtY()
{
  for(int i = labels().begin(); i<labels().end(); ++i)
  {
    int k=labels(i);
    ArraySpectra y_ts;
    ArrayXX      x_tj;
    y_ts.move( handler_.getArraySpectra(i) );
    det_[i] = designer_.computeDesign(times(i), x_tj, y_ts);

    xtx_jj_[k] += x_tj.transpose() * x_tj;
    xty_js_[k] += x_tj.transpose() * y_ts;
    sizek_[k]  += times(i).size();
  }
}

/* compute the matrices X^tX and X^tY */
template<int Size_, class Weights_>
void FuncSpectra<Size_, Weights_>::computeVariances()
{
  for(int i = labels().begin(); i<labels().end(); ++i)
  {
    int k=labels(i);
    ArraySpectra y_ts;
    ArrayXX      x_tj;
    y_ts.move( handler_.getArraySpectra(i) );
    det_[i] = designer_.computeDesign(times(i), x_tj, y_ts);

    sigma2_[k]  += Stat::sumByCol( (y_ts - x_tj * alpha_[k]).square() );
  }
  //
  for (int k=sigma2_.begin(); k< sigma2_.end(); ++k)
  { sigma2_[k] /= sizek_[k];}
}


} // namespace STK

#endif /* STK_FUNCSPECTRA_H */
