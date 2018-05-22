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
 * Project:  CloHe::
 * created on: 8 août 2016
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Gaussian_mut_Sigmat.cpp
 *  @brief In this file we implement the Gaussian_mut_Sigmat class.
 **/

#include "../inst/TemporalSpectrumModels/STK_Gaussian_mut_Sigmat.h"

namespace STK
{

Gaussian_mut_Sigmat::Gaussian_mut_Sigmat( ArraySp4 const& data
                                        , CArrayXX const& clouds
                                        )
                                        : IStatModelBase(data.sizeRows(), data.sizeCols())
                                        , datait_(data)
                                        , residualsit_(data.rows(), data.cols(), 0.)
                                        , cloudsit_(clouds)
                                        , mut_(datait_.cols()), sigmat_(datait_.cols())
{ setNbFreeParameter(4*nbVariable()+16*nbVariable());}

Gaussian_mut_Sigmat::~Gaussian_mut_Sigmat() {}

/* compute the means and the variances */
bool Gaussian_mut_Sigmat::run()
{
  // put missing values when necessary
  for (int j= datait_ .beginCols(); j< datait_.endCols(); ++j)
  {
    for (int i=datait_.beginRows(); i< datait_.endRows(); ++i)
    { if(cloudsit_(i,j)) { datait_(i,j) = Arithmetic<double>::NA();}}
  }
  // compute means
  computeMeans();
  // compute covariances
  computeCovariances();
  // compute residuals
  computeResiduals();
  // compute likelihood
  computeLogLikelihood();
  return true;
}

void Gaussian_mut_Sigmat::computeMeans()
{
  for (int t= datait_ .beginCols(); t< datait_.endCols(); ++t)
  {
    mut_[t]  = 0 ; //* Const::Vector<double, 4>();
    int nobs = datait_.sizeRows();
    for (int i=datait_.beginRows(); i< datait_.endRows(); ++i)
    {
      if (!datait_(i,t).isNA().any()) { mut_[t] += datait_(i,t);}
      else --nobs;
    }
    mut_[t]/= nobs;
  }
}
void Gaussian_mut_Sigmat::computeCovariances()
{
  for (int t= datait_ .beginCols(); t< datait_.endCols(); ++t)
  {
    int nobs = datait_.sizeRows();
    CovSpectra4 var; var=0.0;
    Spectra4 sum = 0.0;
    for (int i=datait_.beginRows(); i< datait_.endRows(); ++i)
    {
      Spectra4 dev;
      if (!datait_(i,t).isNA().any())
      {
        sum += (dev = datait_(i,t) - mut_[t]);
        var += dev * dev.transpose();
      }
      else --nobs;
    }
    sigmat_[t]  = (var - ((sum * sum.transpose())/(double)nobs) )/(double)(nobs) ;
  }
}

void Gaussian_mut_Sigmat::computeResiduals()
{ residualsit_ = datait_- Const::VectorXd(datait_.rows()) * mut_;}

// compute log-likelihoods
void Gaussian_mut_Sigmat::computeLogLikelihood()
{
  double sum =0.0;
  for (int t= datait_ .beginCols(); t< datait_.endCols(); ++t)
  {
    InvertSymMatrix<CovSpectra4, 4> inverter(sigmat_[t]);
    if(!inverter.isInvertible())
    {
#ifdef STK_CLOHE_DEBUG
    std::cout << "sigma_["<<t<<"] is not invertible\n";
#endif
    }
    else
    {
      for (int i=datait_.beginRows(); i< datait_.endRows(); ++i)
      {
        if(cloudsit_(i,t) == 0) // check if value is available
        {
          sum -= 0.5 *( ( residualsit_(i,t).transpose()
                         * inverter.inv()
                         * residualsit_(i,t)
                         )
                         + std::log(inverter.det())
                      );
        }
//      else { sum -= 0.5 * std::log(inverter.det());}
      }
    }
  }
  setLnLikelihood(sum);
}


/* Compute the log-likelihood of the original data. If there is clouds,
 *  value is replaced by its mean.
 **/
void Gaussian_mut_Sigmat::computeLnLikelihood( RMatrix<double> const& xb
                                             , RMatrix<double> const& xg
                                             , RMatrix<double> const& xr
                                             , RMatrix<double> const& xi
                                             , RMatrix<double> const& clouds
                                             , CVectorXd& laux
                                            )
{
  laux = 0.;
  for (int t= clouds.beginCols(); t< clouds.endCols(); ++t)
  {
    InvertSymMatrix<CovSpectra4, 4> inverter(sigmat_[t]);
    if(!inverter.isInvertible())
    { std::cout << "sigma_["<<t<<"] is not invertible\n";}
    else
    {
      for (int i=clouds.beginRows(); i< clouds.endRows(); ++i)
      {
        if(clouds(i,t) == 0) // check if value is available
        {
          Spectra4 xit;
          xit << xb(i,t), xg(i,t), xr(i,t), xi(i,t);
          laux[i] -= 0.5 *( (  (xit - mut_[t]).transpose()
                             * inverter.inv()
                             * (xit - mut_[t])
                            )
                            + std::log(inverter.det())
                           );
        }
//      else { laux[i] -= 0.5 * std::log(inverter.det());}
      }
    }
  }
}

} /* namespace STK */
