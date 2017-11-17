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

/** @file BSplineWCoefficients.h
 *  @brief In this file we define the BSplineWCoefficients class.
 **/

#ifndef BSPLINEWCOEFFICIENTS_H
#define BSPLINEWCOEFFICIENTS_H

#include "STK_CloHe_Util.h"
#include <Regress/include/STK_BSplineCoefficients.h>

namespace STK
{

/** @brief derivation of the BSplineCefficients class allowing to compute
 *  density positioned knots with weights.
 */
class BSplineWCoefficients: public BSplineCoefficients<CPointX>
{
  public:
    typedef BSplineCoefficients<CPointX> Base;
    /** @brief Default constructor : initialize the data members with default
     *  values.
     *  The number of knots is given by the formula
     *  nbKnots = nbControlPoints + degree +1.
     *  @param times the input times values
     *  @param nbControlPoints number of control points
     *  @param degree degree of the B-spline curves
     *  @param position method to use for positioning the knots
     **/
    BSplineWCoefficients( Times const& times
                        , int nbControlPoints =1
                        , int degree = 3
                        , Regress::KnotsPosition position = Regress::uniformKnotsPositions_
                        );
    /** copy constructor.
     *  @param coefs the coefficients to copy
     **/
    BSplineWCoefficients( BSplineWCoefficients const& coefs);
    /** Destructor. */
    inline virtual ~BSplineWCoefficients() {}
    /** clone pattern implementation */
    inline BSplineWCoefficients* clone() const { return new BSplineWCoefficients(*this);}
    /** run the computations. */
    virtual bool run();
    /** run weighted computations */
    template<class Weights_>
    bool run(Weights_ const& w);

  protected:
    /** constant reference on the sampled days */
    Times const& times_;
    /** first time */
    Real tmin_;
    /** last time */
    Real tmax_;

  private:
    /** */
    CPointX allTimes_;
    /** Compute the coefficients of the B-spline curves.*/
    bool computeWeightedKnots();
};

} // namespace STK

#endif /* BSPLINEWCOEFFICIENTS_H */
