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

/** @file STK_GaussianModelsFacade.cpp
 *  @brief In this file xe implement the SpectraDataHandler.
 **/

#include "../inst/TemporalSpectrumModels/STK_SpectraDataHandler.h"

namespace STK
{

SpectraDataHandler::SpectraDataHandler( STK::RVector<int>    const& labels,
                                        STK::RVector<double> const& times,
                                        STK::RMatrix<double> const& xb,
                                        STK::RMatrix<double> const& xg,
                                        STK::RMatrix<double> const& xr,
                                        STK::RMatrix<double> const& xi,
                                        STK::RMatrix<double> const& clouds
                                      )
                                      : fact_(labels), labels_(labels), times_(times)
                                      , xb_(xb), xg_(xg), xr_(xr), xi_(xi)
                                      , clouds_(clouds)
{ fact_.run();}

SpectraDataHandler::~SpectraDataHandler() {}

bool SpectraDataHandler::run()
{
  typedef Stat::Factor< STK::RVector<int> >::EncodingMap EncodingMap;
  typedef Stat::Factor< STK::RVector<int> >::DecodingMap DecodingMap;
  //typedef EncodingMap::const_iterator mapIterator;

  // get reference on the the encoder and decoder maps
  EncodingMap const& encoder = fact_.encoder();//.front();
  DecodingMap const& decoder = fact_.decoder();//.front();

  // copy levels, counts and nbClass for Factor clas
  VectorXi const& levels = fact_.levels();//.front() ;
  VectorXi const& counts = fact_.counts();//.front();
  CVectorXi labelsAsInteger = fact_.asInteger();//.col(0);

  int nbClass     = fact_.nbLevels();//.front();

  // get number of sampling
  Range nbSamples = times_.range();

  // prepare data sets to receive data
  dataSet_.resize(levels.range());
  maskSet_.resize(levels.range());
  xbSet_.resize(levels.range());
  xgSet_.resize(levels.range());
  xrSet_.resize(levels.range());
  xiSet_.resize(levels.range());
  for (int k=dataSet_.begin(); k< dataSet_.end(); ++k)
  {
    dataSet_[k].resize(counts[k], nbSamples);
    maskSet_[k].resize(counts[k], nbSamples);
    xbSet_[k].resize(counts[k], nbSamples);
    xgSet_[k].resize(counts[k], nbSamples);
    xrSet_[k].resize(counts[k], nbSamples);
    xiSet_[k].resize(counts[k], nbSamples);
  }

  CVectorXi row(dataSet_.range(), 0);
  for (int i=labelsAsInteger.begin(); i<labelsAsInteger.end(); ++i)
  {
    int k = labelsAsInteger[i]; // get number of the class for label i
    int irow = row[k];
    for (int j = dataSet_[k].beginCols(); j<dataSet_[k].endCols(); ++j)
    {
      dataSet_[k](irow,j) << xb_(i,j), xg_(i,j), xr_(i,j), xi_(i,j);
      maskSet_[k](irow,j) = clouds_(i,j);
      xbSet_[k](irow,j) = xb_(i,j);
      xgSet_[k](irow,j) = xg_(i,j);
      xrSet_[k](irow,j) = xr_(i,j);
      xiSet_[k](irow,j) = xi_(i,j);
    }
    row[k]++;
  }
  return true;
}

} /* namespace STK */
