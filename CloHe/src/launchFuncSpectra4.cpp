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

/** @file STK_launchFuncSpectra4.cpp
 *  @brief In this file we launch the computation of the model.
 **/

#ifdef _OPENMP
#include <omp.h>
#endif

#include "../inst/TemporalSpectrumModels/STK_FuncSpectra.h"

/** Main function with 4 spectra.
 *  @param r_data list with the labels, times, spectra and clouds
 *  @param r_params model parameters
 *  @param r_models list of S4 FuncSpectra class
 *  @param r_test list with the data test
 **/
RcppExport SEXP launchFuncSpectra4( SEXP r_data, SEXP r_params, SEXP r_models, SEXP r_test)
{
BEGIN_RCPP
#ifdef STK_CLOHE_DEBUG
  stk_cout << _T("Entering launchFuncSpectrum\n");
#endif
  // typedefs
  typedef STK::FuncSpectraDataHandler<4> DataHandler;
  typedef STK::FuncSpectraDesigner<4> Designer;
  typedef STK::FuncSpectra<4, STK::CVectorX> FuncModel;

  typedef FuncModel::ArraySeriesSpectra ArraySeriesSpectra;
  typedef FuncModel::ArraySp ArraySp;
  typedef FuncModel::YArrays YArrays;

  // convert input to Rcpp format
  Rcpp::List Rdata    = r_data;
  Rcpp::List Rlabels  = Rdata["labels"];
  Rcpp::List Rtimes   = Rdata["times"];
  Rcpp::List Rspectra = Rdata["spectra"];
  Rcpp::List Rclouds  = Rdata["clouds"];
  Rcpp::List Rparams  = r_params;

  double maxValue = STK::Arithmetic<double>::max();
  if (!Rf_isNull(Rparams["maxValue"])) { maxValue = Rcpp::as<double>(Rparams["maxValue"]);};

#ifdef STK_CLOHE_DEBUG
  stk_cout << _T("Creating and launching DataHandler\n");
#endif
  // create handler and launch data sets creation
  DataHandler handler(Rlabels, Rtimes, Rspectra, Rclouds);
  handler.setMaxValue(maxValue);
  handler.run();

#ifdef STK_CLOHE_DEBUG
  stk_cout << _T("handler.nbSample()=")      << handler.nbSample() <<_T("\n");
  stk_cout << _T("handler.nbClass()=")       << handler.nbClass() <<_T("\n");
  stk_cout << _T("handler.nbSpectrum()=")    << handler.nbSpectrum() <<_T("\n");
  stk_cout << _T("handler.nk()=")            << handler.nk();
  stk_cout << _T("handler.tmin()=")          << handler.tmin() << _T("\n");
  stk_cout << _T("handler.tmax()=")          << handler.tmax() << _T("\n");
  stk_cout << _T("In launchFuncSpectrum. Read parameters\n");
#endif
  // get parameter values
  std::string kernelName = Rcpp::as<std::string>(Rparams["kernelName"]);
  double width           = Rcpp::as<double>(Rparams["width"]);
  int dim                = Rcpp::as<int>(Rparams["dim"]);
  STK::CPointXi dims(handler.nbClass(), dim);

  std::string basisName  = Rcpp::as<std::string>(Rparams["basisName"]);
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
  stk_cout << _T("\ncreate and launch FuncModel\n");
#endif

  /* utility class allowing to design the matrices xi */
  Designer designer(handler.times(), handler.labels(), kernelName, width, basisName, dims, degree, posKnots);

#ifdef _OPENMP
//  omp_set_num_threads(omp_get_num_procs());
  omp_set_num_threads(2);
#endif
  FuncModel model(handler, designer, criterion);
  model.run();

  // get mean values
  YArrays mu;
  mu.move(model.mean( handler.tmax() - handler.tmin() -1));

#ifdef STK_CLOHE_DEBUG
  stk_cout << "get results\n";
#endif

  // create output
  Rcpp::List Rmodels = r_models;
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
    modelk.slot("mu")              = STK::wrap(mu[k]);
    modelk.slot("sigma2")          = STK::wrap(model.sigma2()[k]);
    modelk.slot("alpha")           = STK::wrap(model.alpha()[k]);
    modelk.slot("eigenvalues")     = STK::wrap(model.eigenvalues()[k]);
    modelk.slot("knots")           = STK::wrap(model.knots()[k]);
    modelk.slot("coefficients")    = STK::wrap(model.coefficients()[k]);
    modelk.slot("tmin")            = model.tmin()[k];
    modelk.slot("tmax")            = model.tmax()[k];
  }
  // get classification and class label
  STK::Labels z(model.zi());
  for (int i=z.begin(); i<z.end(); ++i)
  { z[i] = handler.fact().decoder().find(z[i])->second;}
  // classify data test
  STK::Labels ztest, ltest;
  if (!Rf_isNull(r_test))
  {
    Rcpp::List Rtest         = r_test;
    Rcpp::List Rlabels_test  = Rtest["labels"];
    Rcpp::List Rtimes_test   = Rtest["times"];
    Rcpp::List Rspectra_test = Rtest["spectra"];
    Rcpp::List Rclouds_test  = Rtest["clouds"];

  	DataHandler handlerTest(Rlabels_test, Rtimes_test, Rspectra_test, Rclouds_test);
    handlerTest.setMaxValue(maxValue);
  	handlerTest.run();
  	ztest.resize(handlerTest.nbSample())= 0;
#ifdef STK_CLOHE_DEBUG
  stk_cout << "Launch classification for test data\n";
#endif
  	for( int i=ztest.begin(); i<ztest.end(); ++i)
  	{
  	  ArraySp yt;
      yt.move(handlerTest.getArraySp(i));
      ztest[i] = model.classify(handlerTest.times()[i], yt);
      ztest[i] = handler.fact().decoder().find(ztest[i])->second;
  	}
  	ltest = handlerTest.classLabels();
  }
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
                           , Rcpp::Named("tMin")         = handler.tmin()
                           , Rcpp::Named("tMax")         = handler.tmax()
                           , Rcpp::Named("params")       = Rparams
                           , Rcpp::Named("models")       = Rmodels
                           , Rcpp::Named("labelstest")   = STK::wrap(ltest)
                           , Rcpp::Named("zitest")       = STK::wrap(ztest)
                           );

END_RCPP
}

