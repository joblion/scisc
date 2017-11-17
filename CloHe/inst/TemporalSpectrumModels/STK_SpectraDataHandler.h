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
 * created on: 5 août 2016
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_SpectraDataHandler.h
 *  @brief In this file we define the data handler for functionnal models.
 **/

#ifndef STK_SPECTRADATAHANDLER_H
#define STK_SPECTRADATAHANDLER_H

#include "STK_CloHe_Util.h"

namespace STK
{
/** @brief This facade create the arrays of data of each class and
 *  launch the computation of the models parameters.
 */
class SpectraDataHandler: IRunnerBase
{
  public:
    // Factor class for labels
    typedef Stat::Factor< STK::RVector<int> > Factor;
    /** constructor */
    SpectraDataHandler( STK::RVector<int> const& labels,
                          STK::RVector<double> const& times,
                          STK::RMatrix<double> const& xb,
                          STK::RMatrix<double> const& xg,
                          STK::RMatrix<double> const& xr,
                          STK::RMatrix<double> const& xi,
                          STK::RMatrix<double> const& clouds
                        );
    /** destructor */
    virtual ~SpectraDataHandler();

    /** compute the data sets */
    virtual bool run();
    /** @return the factors of the labels */
    Factor const& fact() const { return fact_;}
    /** @return set of data with all spectrum */
    Data4DSet const& dataSet() const { return dataSet_;}
    /** @return set of data with all spectrum */
    Data1DSet const& maskSet() const { return maskSet_;}
    /** @return a vector of data with blue spectrum */
    Data1DSet const& xbSet() const { return xbSet_;}
    /** @return a vector of data with green spectrum */
    Data1DSet const& xgSet() const { return xgSet_;}
    /** @return a vector of data with red spectrum */
    Data1DSet const& xrSet() const { return xrSet_;}
    /** @return a vector of data with infrared spectrum */
    Data1DSet const& xiSet() const { return xiSet_;}

  protected:
    /** Factor of the labels */
    Factor fact_;
    /** sets with the spectrums */
    Data4DSet dataSet_;
    /** sets with the clouds */
    Data1DSet maskSet_;
    /** sets with the blue spectrum */
    Data1DSet xbSet_;
    /** sets with the green spectrum */
    Data1DSet xgSet_;
    /** sets with the red spectrum */
    Data1DSet xrSet_;
    /** sets with the infrared spectrum */
    Data1DSet xiSet_;

  private:
    STK::RVector<int>    labels_;
    STK::RVector<double> times_;
    STK::RMatrix<double> xb_;
    STK::RMatrix<double> xg_;
    STK::RMatrix<double> xr_;
    STK::RMatrix<double> xi_;
    STK::RMatrix<double> clouds_;
};

} // namespace STK

#endif /* STK_SPECTRADATAHANDLER_H */
