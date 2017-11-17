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
 * Project:  stkpp::Classif::CloHe
 * created on: 6 juin 2017
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_CloHe_Util.h
 *  @brief In this file we define utilities methods and constants for the CloHe project.
 **/


#ifndef STK_CLOHE_UTIL_H
#define STK_CLOHE_UTIL_H

#include <RTKpp.h>

namespace STK
{
// array with for each individuals the times sampling
typedef CArrayVector<CPointX> Times;
// array with for each individuals the labels
typedef CVectorXi Labels;

// typedef for Gaussian_mut_Sigmat model
// values (b,g,r,ir)
typedef CArrayVector<double, 4> Spectra4;
// covariances values (b,g,r,ir)
typedef CArraySquare<double, 4> CovSpectra4;
// array of values
typedef CArray<Spectra4> ArraySpectra4;
// array of mean values
typedef CArrayPoint<Spectra4> PointSpectra4;
// array of values
typedef CArrayPoint<CovSpectra4> PointCov4;

// All classes data sets
typedef Array1D< ArraySpectra4 > Data4DSet;
// All class clouds data sets
typedef Array1D< CArrayXX > Data1DSet;



namespace hidden
{
/** @ingroup
 *  @brief ClassifTraits enclose all structures used by the project CloHe
 **/
template<int Size_> struct ClassifTraits;

/** Specialization for FuncSpectrum models */
template<int Size_ >
struct ClassifTraits
{
  enum
  {
    size_ = Size_
  };

  // if Size_ = 4, this will be the values (b,g,r,ir)
  typedef CArrayPoint<Real, size_> Spectra;
  // times series of the spectra of a given individual
  typedef CArrayVector<Spectra> SeriesSpectra;
  // Array of spectra
  typedef CArray<Real, UnknownSize, Size_> ArraySpectra;
  // array with for each individuals, the times series of the spectra
  typedef Array1D<SeriesSpectra> ArraySeriesSpectra;

};

} // namespace hidden



}

#endif /* STK_CLOHE_UTIL_H */
