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
 * Project:  CloHe
 * created on: 10 juil. 2017
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file init_clohe.c
 *  @brief In this file we .
 **/

// RegisteringDynamic Symbols

#include "CloHe.h"

// declare functions
static const R_CallMethodDef callMethods[]  =
{
  {"launchFuncSpectra4", (DL_FUNC) &launchFuncSpectra4, 6},
  {"launchFuncSpectra10", (DL_FUNC) &launchFuncSpectra10, 6},
  {"launchGaussian_mut_sigmat", (DL_FUNC) &launchGaussian_mut_sigmat, 8},
  {NULL}
};


void R_init_myRoutines(DllInfo *info)
{
  /* Register the .Call routines.
  No .C  .Fortran() or .External() routines,
  so pass those arrays as NULL.
  */
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}


