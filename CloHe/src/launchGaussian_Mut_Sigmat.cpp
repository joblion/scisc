


#include "../inst/TemporalSpectrumModels/STK_Gaussian_mut_Sigmat.h"
#include "../inst/TemporalSpectrumModels/STK_SpectraDataHandler.h"

RcppExport SEXP launchGaussian_mut_sigmat( SEXP r_labels
                                         , SEXP r_times
                                         , SEXP r_xb, SEXP r_xg, SEXP r_xr, SEXP r_xi
                                         , SEXP r_clouds
                                         , SEXP r_res)
{
BEGIN_RCPP

  STK::RVector<int>    labels = r_labels;
  STK::RVector<double> times  = r_times;
  STK::RMatrix<double> xb = r_xb;
  STK::RMatrix<double> xg = r_xg;
  STK::RMatrix<double> xr = r_xr;
  STK::RMatrix<double> xi = r_xi;
  STK::RMatrix<double> clouds = r_clouds;

  Rcpp::List res = r_res;
  // compute data sets
  STK::SpectraDataHandler dataHandler(labels, times, xb, xg, xr, xi, clouds);
  dataHandler.run();

  // intermediate arrays for computing the predicted class
  STK::CVectorXi kmax(labels.range());
  STK::CVectorXd lmax(labels.range(), - STK::Arithmetic<double>::max());
  // compute models
  for (int k=0; k < res.size(); ++k)
  {
    // get model
    Rcpp::S4 modelk = res[k];
    modelk.slot("classNumber")     = k;
    modelk.slot("classLabel")      = dataHandler.fact().decoder().find(k)->second;

    // run computations
    STK::Gaussian_mut_Sigmat model(dataHandler.dataSet()[k], dataHandler.maskSet()[k]);
    model.run();

    // save current model
    Rcpp::List r_mut(model.mut().size());
    Rcpp::List r_sigmat(model.sigmat().size());
    for (int t = model.mut().begin(); t<model.mut().end(); ++t)
    {
      r_mut[t]    = STK::wrap(model.mut()[t]);
      r_sigmat[t] = STK::wrap(model.sigmat()[t]);
    }
    modelk.slot("mut")    = r_mut;
    modelk.slot("sigmat") = r_sigmat;
    modelk.slot("lnLikelihood") = model.lnLikelihood();
    modelk.slot("nbFreeParameter") = model.nbFreeParameter();

    // update predictions
    STK::CVectorXd laux(labels.range(), 0.);
    model.computeLnLikelihood(xb, xg, xr, xi, clouds, laux);
    for(int i = lmax.begin(); i< lmax.end(); i++)
    {
      if (lmax[i] < laux[i])
      {
        lmax[i] = laux[i];
        kmax[i] = k;
      }
    }
  }
  // return results
  return Rcpp::List::create( Rcpp::Named("predict") = STK::wrap(kmax)
                           , Rcpp::Named("models")  = r_res
                           );
  //return Rcpp::wrap(r_res);

  END_RCPP
}
