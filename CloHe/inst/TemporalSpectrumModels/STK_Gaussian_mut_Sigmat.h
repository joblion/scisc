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

/** @file STK_Gaussian_mut_Sigmat.h
 *  @brief In this file we define the Gaussian_mut_Sigmat class.
 **/

#ifndef STK_GAUSSIAN_MUT_SIGMAT_H_
#define STK_GAUSSIAN_MUT_SIGMAT_H_

#include "STK_CloHe_Util.h"

namespace STK
{

/** @brief This class allows to estimate the Gaussian(mu_t,Sigma_t) model
 *
 */
class Gaussian_mut_Sigmat: public IStatModelBase, public IRunnerBase
{
  public:
    //typedef IClassification<ArraySp4, CVectorXi, CVectorX> Base;
    /** Default constructor
     *  @param data a matrix with for each individuals the spectrum
     *  @param clouds a matrix with the clouds, shadows,...
     * */
    Gaussian_mut_Sigmat( ArraySp4 const& data, CArrayXX const& clouds);
    /** destructor */
    virtual ~Gaussian_mut_Sigmat();

    /** @return mean at each times */
    inline PointSpectra4 const& mut() const { return mut_;}
    /** @return covarainces at each times */
    inline PointCov4 const& sigmat() const { return sigmat_;}
    /** @return residuals at each times */
    inline ArraySp4 const& residuals() const { return residualsit_;}

    /** compute the means and the variances */
    virtual bool run();
    /** Compute the log-likelihood of the original data. If there is clouds,
     *  value is replaced by its mean.
     **/
    void computeLnLikelihood( RMatrix<double> const& xb
                            , RMatrix<double> const& xg
                            , RMatrix<double> const& xr
                            , RMatrix<double> const& xi
                            , RMatrix<double> const& clouds
                            , CVectorXd& likelihoods
                            );

  protected:
    /** Array with the spectrums */
    ArraySp4 datait_;
    /** Array with the the residuals of the spectrums */
    ArraySp4 residualsit_;
    /** Array with the clouds, shadows,... */
    CArrayXXd cloudsit_;
    /** means of the spectrum at each instants */
    PointSpectra4 mut_;
    /** covariances of the spectrum at each instants */
    PointCov4 sigmat_;

  private:
    void computeMeans();
    void computeCovariances();
    void computeResiduals();
    void computeLogLikelihood();
};

} // namespace STK

#endif /* STK_GAUSSIAN_MUT_SIGMAT_H_ */
