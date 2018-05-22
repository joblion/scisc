/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2017  Serge Iovleff, Université Lille 1, Inria

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
 * created on: Oct 26, 2017
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file BSplineWCoefficients.cpp
 *  @brief In this file .
 **/
#include "../inst/TemporalSpectrumModels/STK_BSplineWCoefficients.h"

namespace STK
{

BSplineWCoefficients::BSplineWCoefficients( Times const& times
                                          , Labels const& labels
                                          , int nbControlPoints
                                          , int degree
                                          , Regress::KnotsPosition position
                                          )
                                          : Base(0, nbControlPoints, degree, position)
                                          , times_(times)
                                          , labels_(labels)
                                          , allTimes_(366)
                                          , tmin_(0), tmax_(365)
                                          , classNumber_(1)
{}

/* copy constructor.
 *  @param coefs the coefficients to copy
 **/
BSplineWCoefficients::BSplineWCoefficients( BSplineWCoefficients const& coefs)
                                          : Base(coefs)
                                          , times_(coefs.times_)
                                          , labels_(coefs.labels_)
                                          , tmin_(coefs.tmin_), tmax_(coefs.tmax_)
                                          , classNumber_(coefs.classNumber_)
{}

bool BSplineWCoefficients::BSplineWCoefficients::run()
{
  tmin_ = 366; tmax_ = 0;
  // for each individual
  for (int i = labels_.begin(); i < labels_.end(); ++i)
  {
    if (labels_[i] == classNumber_)
    {
      for(int t=times_[i].begin(); t< times_[i].end(); ++t)
      {
        tmin_ = std::min(tmin_, times_[i][t]);
        tmax_ = std::max(tmax_, times_[i][t]);
      }
    }
  }

  // create data set with integer sampled dates in the range [tmin_, tmax_]
  allTimes_.resize( Range((int)tmin_, (int)tmax_, 0));
  for (int i = allTimes_.begin(); i < allTimes_.end(); ++i) { allTimes_[i] = i;}
  p_data_ = &allTimes_;

  // for density and periodic knots, use base class method
  if (position_ != Regress::densityKnotsPositions_) { return Base::run();}

  // initialize
  if (!initializeStep()) return false;
  // compute the knots and coefficients
  if (computeWeightedKnots()) { computeCoefficients();}
  else { return false;}
  this->hasRun_ = true;
  return true;
}

bool BSplineWCoefficients::computeWeightedKnots()
{
#ifdef STK_CLOHE_DEBUG
  try
  {
#endif
    // count the sampled days
    CPointXi count( allTimes_.range(), 0 ), cumcount( allTimes_.range(), 0);
    for (int i = times_.begin(); i < times_.end(); ++i)
    {
      if (labels_[i] == classNumber_)
      {
        for(int t = times_[i].begin(); t < times_[i].end(); ++t)
        { count[ (int)times_[i][t] ]++;} // Note: must be between 0 and 355
      }
    }
    // create cumulated counts
    cumcount.front() = count.front();
    for(int t = count.begin()+1; t < count.end(); ++t)
    { cumcount[t] = cumcount[t-1] + count[t];}

    // set density knots
    Real step = (Real)cumcount.back()/(Real)(lastControlPoint_-degree_+1);
    Real oldt = tmin_;
    int t = cumcount.begin();
    // set internal knots
    for (int k =  degree_ + 1, kcell =1; k < dim_; k++, kcell++)
    {
      Real prop = kcell * step;
      while (cumcount[t]<prop) { ++t;} // t is such that cumcount[t]>=prop

      Real freq = (count[t]==0) ? 0. : (cumcount[t]-prop)/count[t]; // count[t]==0 should be always false ?
      knots_[k] = (t-1) *freq + (1.-freq)*t;
    }
    if (knots_[degree_+1] < tmin_) { knots_[degree_+1] = tmin_;}
    // set external knots
    for ( int k=0, j = dim_; k <= degree_; j++, k++)
    {
      knots_[k] = (int)tmin_;
      knots_[j] = (int)tmax_;
    }
#ifdef STK_CLOHE_DEBUG
  }
  catch (Exception const& e)
  {
    msg_error_ = STKERROR_NO_ARG(BSplineWCoefficients::computeWeightedKnots, e.error());
  }
#endif
  return true;
}

} /* namespace STK */
