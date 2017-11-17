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
                                          , int nbControlPoints
                                          , int degree
                                          , Regress::KnotsPosition position
                                          )
                                          : Base(0, nbControlPoints, degree, position)
                                          , times_(times)
                                          , allTimes_(366)
                                          , tmin_(0)
                                          , tmax_(365)
{}

/* copy constructor.
 *  @param coefs the coefficients to copy
 **/
BSplineWCoefficients::BSplineWCoefficients( BSplineWCoefficients const& coefs)
                                          : times_(coefs.times_), Base(coefs)
                                          , tmin_(coefs.tmin_), tmax_(coefs.tmax_) {}

bool BSplineWCoefficients::BSplineWCoefficients::run()
{
  tmin_ = 366; tmax_ = 0;
  for (int i = times_.begin(); i < times_.end(); ++i)
  { for(int t=times_[i].begin(); t< times_[i].end(); ++t)
    {
      tmin_ = std::min(tmin_, times_[i][t]);
      tmax_ = std::max(tmax_, times_[i][t]);
    }
  }
  allTimes_.resize( Range((int)tmin_, (int)tmax_, 0));
  // initialize data
  for (int i = allTimes_.begin(); i < allTimes_.end(); ++i) { allTimes_[i] = i;}
  p_data_ = &allTimes_;

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
    Range timeRange =allTimes_.range();
    timeRange.incEnd();
    CPointXi count( allTimes_.range(), 0 ), cumcount( allTimes_.range(), 0);
    for (int i = times_.begin(); i < times_.end(); ++i)
    {
      for(int t = times_[i].begin(); t < times_[i].end(); ++t)
      {
        count[ (int)times_[i][t] ]++; // Note: must be between 0 and 355
      }
    }
    // create cumulated counts
    cumcount.front() = count.front();
    for(int t = count.begin()+1; t < count.end(); ++t)
    { cumcount[t] = cumcount[t-1] + count[t];}

    // set density knots
    Real step = (Real)cumcount.back()/(Real)(lastControlPoint_-degree_+1);
    int t = cumcount.begin();
    // set internal knots
    for (int k =  degree_ + 1, kcell =1; k <= lastControlPoint_; k++, kcell++)
    {
      Real prop = kcell * step;
      while (cumcount[t]<prop) { ++t;} // t is such that cumcount[t]>=prop
      int s = t;
      while( count[s]==0 ) { --s; if (s==0) break;} // s is such that count[s]!=0
      Real freq = (s==0) ? 0. : (cumcount[t]-prop)/count[s];
      knots_[k] = s *freq + (1.-freq)*t;
    }
    // set external knots
    for ( int k=0, j = lastControlPoint_+1; k <= degree_; j++, k++)
    {
      knots_[k] = (int)tmin_;
      knots_[j] = (int)tmax_;
    }
#ifdef STK_CLOHE_DEBUG
  }
  catch (Exception const& e)
  {
    msg_error_ = STKERROR_NO_ARG(BSplineWCoefficients::computeWeightedKnots,e.error());
    return false;
  }
#endif
  return true;
}

} /* namespace STK */
