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

/** @file STK_launchFuncSpectra10.cpp
 *  @brief In this file we launch the computation of the model.
 **/

#ifdef _OPENMP
#include <omp.h>
#endif

#include "../inst/TemporalSpectrumModels/STK_CloHe_Util.h"
#include "../inst/TemporalSpectrumModels/STK_FuncSpectra.h"
#include "../inst/TemporalSpectrumModels/STK_FuncSpectraDataHandler.h"

/** Main function.
 *  @param r_labels vector with the number of the class
 *  @param r_times vector with the a sample times of each pixel
 *  @param r_spectra spectra
 *  @param r_clouds clouds of each pixel at each time
 *  @param r_params model parameters
 *  @param r_models vector with S4 classes models
 **/
RcppExport SEXP launchFuncSpectra10( SEXP r_labels
                                   , SEXP r_times
                                   , SEXP r_spectra
                                   , SEXP r_clouds
                                   , SEXP r_params
                                   , SEXP r_models
                                   )
{
BEGIN_RCPP
#ifdef STK_CLOHE_DEBUG
  stk_cout << _T("Entering launchFuncSpectrum\n");
#endif
  // typedefs
  typedef STK::FuncSpectraDataHandler<10> DataHandler;
  typedef STK::FuncSpectra<10, STK::CVectorX> FuncModel;

  typedef FuncModel::ArraySeriesSpectra ArraySeriesSpectra;
  typedef FuncModel::YArrays YArrays;

  // convert input to Rcpp format
  Rcpp::List Rlabels  = r_labels;
  Rcpp::List Rtimes   = r_times;
  Rcpp::List Rspectra = r_spectra;
  Rcpp::List Rclouds  = r_clouds;
  Rcpp::List Rparams  = r_params;

#ifdef STK_CLOHE_DEBUG
  stk_cout << _T("In launchFuncSpectrum. Read parameters\n");
#endif

  // get parameter values
  std::string kernelName = Rcpp::as<std::string>(Rparams["kernelName"]);
  double width           = Rcpp::as<double>(Rparams["width"]);
  int dim                = Rcpp::as<int>(Rparams["dim"]);
  std::string posKnots   = Rcpp::as<std::string>(Rparams["posKnots"]);
  int degree             = Rcpp::as<int>(Rparams["degree"]);
  std::string criterion  = Rcpp::as<std::string>(Rparams["criterion"]);

#ifdef STK_CLOHE_DEBUG
  stk_cout << _T("kernelName =") << kernelName << _T("\n");
  stk_cout << _T("width =")      << width      << _T("\n");
  stk_cout << _T("dim =")        << dim        << _T("\n");
  stk_cout << _T("posKnots =")   << posKnots   << _T("\n");
  stk_cout << _T("degree =")     << degree     << _T("\n");
  stk_cout << _T("criterion =")  << criterion  << _T("\n");
#endif

  // create output
  Rcpp::List Rmodels = r_models;
  // create handler and launch data sets creation
  DataHandler handler(Rlabels, Rtimes, Rspectra, Rclouds);
  handler.run();

#ifdef STK_CLOHE_DEBUG
  stk_cout << _T("handler.nbSample()=")      << handler.nbSample() <<_T("\n");
  stk_cout << _T("handler.nbClass()=")       << handler.nbClass() <<_T("\n");
  stk_cout << _T("handler.nbSpectrum()=")    << handler.nbSpectrum() <<_T("\n");
  stk_cout << _T("handler.nk()=")            << handler.nk();
  stk_cout << _T("handler.tmin()=")          << handler.tmin() << _T("\n");
  stk_cout << _T("handler.tmax()=")          << handler.tmax() << _T("\n");
  stk_cout << _T("\ncreate and launch FuncModel\n");
#endif

#ifdef _OPENMP
//  omp_set_num_threads(omp_get_num_procs());
  omp_set_num_threads(2);
#endif

  FuncModel model(handler, kernelName, width, dim, degree, posKnots, criterion);
  model.run();
  // get mean value

  YArrays mu;
  mu.move(model.mean(handler.tmin(), handler.tmax(), handler.tmax() - handler.tmin() +1));

#ifdef STK_CLOHE_DEBUG
  stk_cout << "get results\n";
#endif

  for (int k=0; k < handler.nbClass(); ++k)
  {
    Rcpp::S4 modelk = Rmodels[k];
    modelk.slot("classNumber")     = k;
    modelk.slot("classLabel")      = handler.fact().decoder().find(k)->second;
    modelk.slot("nbSpectrum")      = handler.nbSpectrum();
    modelk.slot("lnLikelihood")    = model.lnLikelihood()[k];
    modelk.slot("criterion")       = model.criterion()[k];
    modelk.slot("nbFreeParameter") = model.nbFreeParameters()[k];
    modelk.slot("nk")              = model.nk()[k];
    modelk.slot("pk")              = model.pk()[k];
    modelk.slot("muk")             = STK::wrap(mu[k]);
    modelk.slot("sigma2k")         = STK::wrap(model.sigma2()[k]);
  }
  // get classification and class label
  STK::Labels z(model.zi());
  for (int i=z.begin(); i<z.end(); ++i)
  { z[i] = handler.fact().decoder().find(z[i])->second;}

  // return results
  return Rcpp::List::create( Rcpp::Named("nbSample")     = handler.nbSample()
                           , Rcpp::Named("nbClass")      = handler.nbClass()
                           , Rcpp::Named("nbSpectrum")   = handler.nbSpectrum()
                           , Rcpp::Named("labels")       = STK::wrap(handler.classLabels())
                           , Rcpp::Named("zi")           = STK::wrap(z)
                           , Rcpp::Named("tik")          = STK::wrap(model.tik())
                           , Rcpp::Named("nk")           = STK::wrap(model.nk())
                           , Rcpp::Named("pk")           = STK::wrap(model.pk())
                           , Rcpp::Named("lnLikelihood") = STK::wrap(model.lnLikelihood())
                           , Rcpp::Named("criterion")    = STK::wrap(model.criterion())
                           , Rcpp::Named("knots")        = STK::wrap(model.knots())
                           , Rcpp::Named("coefficients") = STK::wrap(model.coefficients())
                           , Rcpp::Named("tMin")         = handler.tmin()
                           , Rcpp::Named("tMax")         = handler.tmax()
                           , Rcpp::Named("params")       = Rparams
                           , Rcpp::Named("models")       = Rmodels
                           );

END_RCPP
}

