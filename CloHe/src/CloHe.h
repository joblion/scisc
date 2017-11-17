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
 * created on: 10 juil. 2017
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file rclohe.h
 *  @brief In this file .
 **/


#ifndef SRC_RCLOHE_H
#define SRC_RCLOHE_H

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#ifdef __cplusplus
extern "C"
{
#endif
SEXP launchFuncSpectra4( SEXP r_labels
                       , SEXP r_times
                       , SEXP r_spectra
                       , SEXP r_clouds
                       , SEXP r_params
                       , SEXP r_res
                       );
SEXP launchFuncSpectra10( SEXP r_labels
                        , SEXP r_times
                        , SEXP r_spectra
                        , SEXP r_clouds
                        , SEXP r_params
                        , SEXP r_res
                         );
SEXP launchGaussian_mut_sigmat( SEXP r_labels
                              , SEXP r_times
                              , SEXP r_xb, SEXP r_xg, SEXP r_xr, SEXP r_xi
                              , SEXP r_clouds
                              , SEXP r_res
                              );
#ifdef __cplusplus
extern "C"
}
#endif

#endif /* SRC_RCLOHE_H_ */
