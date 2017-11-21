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

/** @file STK_FuncSpectraDatahandler.cpp
 *  @brief In this file we implement the data handler run method
 **/

#include "../inst/TemporalSpectrumModels/STK_FuncSpectraDataHandler.h"

namespace STK
{

template<>
bool FuncSpectraDataHandler<4>::run()
{
  typedef Stat::Factor< STK::RVector<int> >::EncodingMap EncodingMap;
  typedef Stat::Factor< STK::RVector<int> >::DecodingMap DecodingMap;

  // get reference on the the encoder and decoder maps
  EncodingMap const& encoder = fact_.encoder();
  DecodingMap const& decoder = fact_.decoder();

  // get number of individuals in each class
  // Some checks
  if (nk_.size() < 1) { this->msg_error_ = "No class given"; return false;}
  if (nbSample_ < 3*nbClass_) { this->msg_error_ = "No enough samples"; return false;}

  // prepare data sets to receive data
  datait_.resize(nbSample_);
  timesit_.resize(nbSample_);
  nti_.resize(nbSample_);

  // create index and iterators, loop over
  int iSample = classLabels_.begin();
  tmin_=355; tmax_=0;

  Rcpp::List::iterator it=Rspectra_.begin();
  Rcpp::List Rx1_ = *it; ++it;
  Rcpp::List Rx2_ = *it; ++it;
  Rcpp::List Rx3_ = *it; ++it;
  Rcpp::List Rx4_ = *it;
  for( Rcpp::List::iterator it_times = Rtimes_.begin()
     ,                      it1    = Rx1_.begin()
     ,                      it2    = Rx2_.begin()
     ,                      it3    = Rx3_.begin()
     ,                      it4    = Rx4_.begin()
     ,                      it_clouds = Rclouds_.begin()
     ; it_times!=Rtimes_.end()
     ; ++it_times, ++it1, ++it2, ++it3, ++it4, ++it_clouds)
  {
#ifdef STK_CLOHE_DEBUG
  stk_cout << _T("create Spectra and clouds matrices\n");
#endif
    // All lists have the same number of records
    Rcpp::NumericVector Rtimes = *it_times;
    // matrices of size number of sample of this year by the number of times
    Rcpp::NumericMatrix Rx1 = *it1, Rx2 = *it2, Rx3 = *it3, Rx4 = *it4;
    // clouds/shadows/etc...
    Rcpp::NumericMatrix Rclouds = *it_clouds;
    // loop over the individuals of this passage
    for (int i=0; i< Rclouds.nrow(); ++i, iSample++)
    {
      // compute the number of available sampling for individuals i
      nti_[iSample] = 0;
      for (int j=0; j< Rclouds.ncol(); ++j) { if (Rclouds(i,j) == 0) ++nti_[iSample];}

      // fill data sets
      datait_[iSample].resize(nti_[iSample]);
      timesit_[iSample].resize(nti_[iSample]);

      for (int j=0, t=datait_[iSample].begin(); j< Rclouds.ncol(); ++j)
      {
        if (Rclouds(i,j) == 0) // use only no clouds
        {
          datait_[iSample][t] << Rx1(i,j), Rx2(i,j), Rx3(i,j), Rx4(i,j);
          timesit_[iSample][t] = Rtimes[j];
          ++t;
        }
      }
      tmin_ = std::min(tmin_, timesit_[iSample].minElt());
      tmax_ = std::max(tmax_, timesit_[iSample].maxElt());
    }
  }
  return true;
}

template<>
bool FuncSpectraDataHandler<10>::run()
{
  typedef Stat::Factor< STK::RVector<int> >::EncodingMap EncodingMap;
  typedef Stat::Factor< STK::RVector<int> >::DecodingMap DecodingMap;

  // get reference on the the encoder and decoder maps
  EncodingMap const& encoder = fact_.encoder();
  DecodingMap const& decoder = fact_.decoder();

  // get number of individuals in each class
  // Some checks
  if (nk_.size() < 1) { this->msg_error_ = "No class given"; return false;}
  if (nbSample_ < 3*nbClass_) { this->msg_error_ = "No enough samples"; return false;}

  // prepare data sets to receive data
  datait_.resize(nbSample_);
  timesit_.resize(nbSample_);
  nti_.resize(nbSample_);

  // create iterators over the ten spectra
  // For some reason an array of Rcpp::List does not work
  Rcpp::List::iterator it=Rspectra_.begin();
  Rcpp::List Rx1_ = *it; ++it;
  Rcpp::List Rx2_ = *it; ++it;
  Rcpp::List Rx3_ = *it; ++it;
  Rcpp::List Rx4_ = *it; ++it;
  Rcpp::List Rx5_ = *it; ++it;
  Rcpp::List Rx6_ = *it; ++it;
  Rcpp::List Rx7_ = *it; ++it;
  Rcpp::List Rx8_ = *it; ++it;
  Rcpp::List Rx9_ = *it; ++it;
  Rcpp::List Rx10_ = *it;

  // loop over the years
  int iSample = classLabels_.begin();
  tmin_=355; tmax_=0;
  for( Rcpp::List::iterator it_times = Rtimes_.begin()
     ,                      it1    = Rx1_.begin()
     ,                      it2    = Rx2_.begin()
     ,                      it3    = Rx3_.begin()
     ,                      it4    = Rx4_.begin()
     ,                      it5    = Rx5_.begin()
     ,                      it6    = Rx6_.begin()
     ,                      it7    = Rx7_.begin()
     ,                      it8    = Rx8_.begin()
     ,                      it9    = Rx9_.begin()
     ,                      it10   = Rx10_.begin()
     ,                      it_clouds = Rclouds_.begin()
     ; it_times!=Rtimes_.end()
     ; ++it_times, ++it1, ++it2, ++it3, ++it4, ++it5, ++it6, ++it7, ++it8, ++it9, ++it10, ++it_clouds)
  {
#ifdef STK_CLOHE_DEBUG
  stk_cout << _T("Current year. create times vector, Spectra and clouds matrices\n");
#endif
    // All lists have the same number of records
    Rcpp::NumericVector Rtimes = *it_times;
    // matrices of size number of sample of this year by the number of times
    Rcpp::NumericMatrix Rx1 = *it1
                      , Rx2 = *it2
                      , Rx3 = *it3
                      , Rx4 = *it4
                      , Rx5 = *it5
                      , Rx6 = *it6
                      , Rx7 = *it7
                      , Rx8 = *it8
                      , Rx9 = *it9
                      , Rx10 = *it10 ;
    // clouds/shadows/etc...
    Rcpp::NumericMatrix Rclouds = *it_clouds;
    // loop over the individuals of this passage
    for (int i=0; i< Rclouds.nrow(); ++i, iSample++)
    {
      // compute the number of available sampling for individuals i
      nti_[iSample] = 0;
      for (int j=0; j< Rclouds.ncol(); ++j) { if (Rclouds(i,j) == 0) ++nti_[iSample];}

      // fill data sets
      datait_[iSample].resize(nti_[iSample]);
      timesit_[iSample].resize(nti_[iSample]);

      for (int j=0, t=datait_[iSample].begin(); j< Rclouds.ncol(); ++j)
      {
        if (Rclouds(i,j) == 0) // use only no clouds
        {
          datait_[iSample][t] << Rx1(i,j), Rx2(i,j), Rx3(i,j), Rx4(i,j)
                               , Rx5(i,j), Rx6(i,j), Rx7(i,j), Rx8(i,j)
                               , Rx9(i,j), Rx10(i,j);
          timesit_[iSample][t] = Rtimes[j];
          ++t;
        }
      }
      tmin_ = std::min(tmin_, timesit_[iSample].minElt());
      tmax_ = std::max(tmax_, timesit_[iSample].maxElt());
    }
  }
  return true;
}


} // namespace STK

