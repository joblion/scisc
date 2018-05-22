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

/** @file STK_FuncSpectraDatahandler.h
 *  @brief In this file we define the data handler for the FuncSpectrum model.
 **/

#ifndef STK_FUNCSPECTRADATAHANDLER_H
#define STK_FUNCSPECTRADATAHANDLER_H

#include "STK_CloHe_Util.h"

namespace STK
{
/** @brief This handler creates the arrays of data needed.
 */
template<int Size_>
class FuncSpectraDataHandler: IRunnerBase
{
  public:
    // Factor class for labels
    typedef Stat::Factor< CVectorXi > Factor;
    // spectra container (CArrayPoint<double, size_>)
    typedef typename hidden::ClassifTraits<Size_ >::Spectra Spectra;

    // array of spectra (CArray<Real, UnknownSize, Size_>)
    typedef typename hidden::ClassifTraits<Size_ >::ArraySp ArraySp;
    // array of values during time for a given pixel (CArrayVector<Spectra>)
    typedef typename hidden::ClassifTraits<Size_ >::SeriesSpectra SeriesSpectra;
    // 3D array of all values (Array1D<SeriesSpectra>)
    typedef typename hidden::ClassifTraits<Size_ >::ArraySeriesSpectra ArraySeriesSpectra;
    enum
    {
      size_ =  hidden::ClassifTraits<Size_ >::size_,
    };

    /** constructor.
     *  All parameters are list of records by years.
     *  @param labels a list of labels for each year.
     * */
    FuncSpectraDataHandler( Rcpp::List const& labels,
                            Rcpp::List const& times,
                            Rcpp::List const& spectra,
                            Rcpp::List const& clouds
                          );
    /** destructor */
    virtual ~FuncSpectraDataHandler() {}

    /** @return the number of sample */
    inline int nbSample() const { return nbSample_;}
    /** @return the number of spectrum */
    inline int nbSpectrum() const { return nbSpectrum_;};
    /** @return the number of class */
    inline int nbClass() const { return nbClass_;};

    /** @return the factors of the labels */
    inline Factor const& fact() const { return fact_;}
    /** @return a vector with for each individuals its label computed by fact */
    inline Labels const& labels() const { return fact().asInteger();}
    /** @return a vector with for each individuals the original class labels */
    inline Labels const& classLabels() const { return classLabels_;}
    /** @return set of data with for each individuals the times series of the spectra */
    inline ArraySeriesSpectra const& data() const { return datait_;}
    /** @return a vector of data with for each individuals the times sampling */
    inline Times const& times() const { return timesit_;}
    /** @return first time sampling */
    inline Real const& tmin() const { return tmin_;}
    /** @return last time sampling */
    inline Real const& tmax() const { return tmax_;}
    /** @return number of individuals in each class */
    inline CPointXi const& nk() const { return nk_;}
    /** @return number of sampling times for each individuals */
    inline CPointXi const& ni() const { return nti_;}

    /** @param i,yi index and spectra of the individual i
     *  @return the spectra of the individual i */
    ArraySp getArraySp(int i) const;

    /** compute the data sets */
    virtual bool run();

    /** set the maximal value of the spectra */
    void setMaxValue(double maxValue) { maxValue_ = maxValue;}
  protected:
    /** number of sample */
    int nbSample_;
    /** number of spectra */
    int nbSpectrum_;
    /** number of class */
    int nbClass_;
    /** Factor of the labels */
    Factor fact_;

    /** array with for each individuals the class labels (original labelling) */
    Labels classLabels_;
    /** sets with the spectra */
    ArraySeriesSpectra datait_;
    /** Vector array with for each individual the sampling times */
    Times timesit_;
    /** number of individuals in each classes */
    CPointXi nk_;
    /** number of available sampling for each individual */
    CPointXi nti_;
    /** first time */
    Real tmin_;
    /** last time */
    Real tmax_;

  private:
    Rcpp::List Rlabels_;
    Rcpp::List Rtimes_;
    Rcpp::List Rspectra_;
    Rcpp::List Rclouds_;
    /** maximal value of the spectra */
    double maxValue_;
};

template<int Size_>
FuncSpectraDataHandler<Size_>::FuncSpectraDataHandler( Rcpp::List const& labels,
                                                       Rcpp::List const& times,
                                                       Rcpp::List const& spectra,
                                                       Rcpp::List const& clouds
                                                     )
                                                     : nbSample_(0)
                                                     , nbSpectrum_(size_)
                                                     , nbClass_(0)
                                                     , fact_()
                                                     , datait_(), timesit_()
                                                     , Rlabels_(labels), Rtimes_(times)
                                                     , Rspectra_(spectra)
                                                     , Rclouds_(clouds)
                                                     , tmin_(0)
                                                     , tmax_(355)
                                                     , maxValue_(STK::Arithmetic<Real>::max())
{
#ifdef STK_CLOHE_DEBUG
  stk_cout << _T("Entering FuncSpectraDataHandler<Size_> constructor\n") ;
#endif
  // compute nbSample
  for(Rcpp::List::iterator it = Rlabels_.begin(); it!=Rlabels_.end(); ++it)
  {
    Rcpp::IntegerVector l = *it;
    nbSample_ += l.length();
  }
  // get labels
#ifdef STK_CLOHE_DEBUG
  stk_cout << _T("FuncSpectraDataHandler<Size_>::nbSample_ =") << nbSample_ << _T("\n");
  stk_cout << _T("build classLabels_ vector\n");
#endif

  classLabels_.resize(nbSample_);
  int iLabel = classLabels_.begin();
  for(Rcpp::List::iterator it = Rlabels_.begin(); it!=Rlabels_.end(); ++it)
  {
    Rcpp::IntegerVector l = *it;
    for( int icur=0; icur<l.length(); ++icur)
    { classLabels_[iLabel] = (int)l[icur]; ++iLabel;}
  }

#ifdef STK_CLOHE_DEBUG
  stk_cout << _T("original labels copied\n");
  stk_cout << _T("classLabels_[1:10]=") << classLabels_.sub(STK::Range(0,10)).transpose();
  stk_cout << _T("create and launch Factor class\n");
#endif

  fact_.setData(classLabels_);
  fact_.run();
  nk_       = fact_.counts();
  nbClass_  = nk_.size();

#ifdef STK_CLOHE_DEBUG
  stk_cout << _T("FuncSpectraDataHandler<Size_>::nbClass_ =") << nbClass_ << _T("\n");
  stk_cout << _T("FuncSpectraDataHandler<Size_>::nbSpectrum_ =") << nbSpectrum_ << _T("\n");
  stk_cout << _T("FuncSpectraDataHandler<Size_>::nk_ =") << nk_;
  stk_cout << _T("fact_.asInteger()[1:10]=")             << fact_.asInteger().sub(STK::Range(0,10)).transpose();
  stk_cout << _T("Terminating FuncSpectraDataHandler<Size_> constructor\n") ;
#endif
}


template<>
bool FuncSpectraDataHandler<4>::run();
template<>
bool FuncSpectraDataHandler<10>::run();

/* @param i individual inumber
 *  @return the times sampling and the spectra */
template<int Size_>
typename FuncSpectraDataHandler<Size_>::ArraySp FuncSpectraDataHandler<Size_>::getArraySp(int i) const
{
  // get sampled  spectra
  ArraySp yt(datait_[i].range(), nbSpectrum());
  for(int t=yt.beginRows(); t< yt.endRows(); ++t)
  { yt.row(t) = datait_[i][t];}
  return yt;
}


} // namespace STK

#endif /* STK_FUNCSPECTRADATAHANDLER_H */
