/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/

#pragma once

/**
 * @file IntegralInvariantBarycenterEstimatorFromSurfelsFunctors.h
 * @author David Coeurjolly (\c david.coeurjolly@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systemes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 *
 * @date 2013/05/30
 *
 * Header file for module IntegralInvariantBarycenterEstimatorFromSurfelsFunctors.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(IntegralInvariantBarycenterEstimatorFromSurfelsFunctors_RECURSES)
#error Recursive header files inclusion detected in IntegralInvariantBarycenterEstimatorFromSurfelsFunctors.h
#else // defined(IntegralInvariantBarycenterEstimatorFromSurfelsFunctors_RECURSES)
/** Prevents recursive inclusion of headers. */
#define IntegralInvariantBarycenterEstimatorFromSurfelsFunctors_RECURSES

#if !defined IntegralInvariantBarycenterEstimatorFromSurfelsFunctors_h
/** Prevents repeated inclusion of headers. */
#define IntegralInvariantBarycenterEstimatorFromSurfelsFunctors_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/kernel/NumberTraits.h"
#include "DGtal/topology/CSCellEmbedder.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{
  namespace functors
  {
  template<typename TSurfel, typename TSCellEmbedder>
  struct IntegralInvariantBarycenterEstimatorFromSurfels
  {
  public:

    template<typename Space>
    struct BarycenterStruct
    {
      typename Space::RealVector barycenter;
      typename Space::RealVector center;

      inline bool operator==(const BarycenterStruct& a)
      {
       if (a.barycenter==barycenter && a.center== center)
          return true;
       else
          return false;
      }

      inline bool operator<(const BarycenterStruct& a)
      {
       if (a.center < center)
          return true;
       else
          return false;
      }

      inline bool operator>(const BarycenterStruct& a)
      {
       if (a.center > center)
          return true;
       else
          return false;
      }

      inline bool operator!=(const BarycenterStruct& a)
      {
        return !operator==( a );
      }
    };



    ///Surfel type
    typedef TSurfel Surfel;
    typedef TSCellEmbedder SCellEmbedder;
    BOOST_CONCEPT_ASSERT(( concepts::CSCellEmbedder<SCellEmbedder> ));
    typedef typename SCellEmbedder::Space Space;
    typedef typename Space::RealVector RealVector;
    typedef typename Space::RealPoint RealPoint;

    ///Embedder type
    ///Type of output values
    typedef BarycenterStruct<Space> Quantity;

    /**
     * Constructor.
     *
     * @param [in] anEmbedder any model of CSCellEmbedder.
     * @param [in] h a grid step
     */
    IntegralInvariantBarycenterEstimatorFromSurfels(ConstAlias<SCellEmbedder> anEmbedder ,
                              const double h):
      myEmbedder(&anEmbedder), myH(h)
    {
      myFirstSurfel = true;
      myN = 0;
      mySumX = RealVector::diagonal(0.0);
    }

    /**
     * Destructor
     */
    ~IntegralInvariantBarycenterEstimatorFromSurfels() {}

    /**
     * Push a surfel to the estimator. For this dummy estimator,
     * we just count the number of surfels.
     */
    void pushSurfel(const Surfel &aSurfel,
                    const double aDistance)
    {
      BOOST_VERIFY(aDistance == aDistance);
      BOOST_VERIFY(aSurfel == aSurfel);

      if (myFirstSurfel)
      {
        myReceiver = myEmbedder->operator()(aSurfel);
        myFirstSurfel = false;
      }
      //else
      {
        RealPoint p = myEmbedder->operator()(aSurfel);

        mySumX += p;

        ++myN;
      }
    }

    /**
     * @return the estimated quantity.
     */
    Quantity eval( ) const
    {
      Quantity res;

      res.barycenter = (mySumX*(1.0/myN))*myH;
      res.center = myReceiver*myH;

      return res;
    }

    /**
     * Reset the estimator.
     */
    void reset()
    {
      myFirstSurfel = true;
      myN = 0;
      mySumX = RealVector::diagonal(0.0);
    }

  private:

    /**
     * Private default constructor.
     */
    IntegralInvariantBarycenterEstimatorFromSurfels();

    ///ConstAlias of the Embedder
    const SCellEmbedder * myEmbedder;

    RealVector myReceiver;

    ///Surfel counter.
    bool myFirstSurfel;

    RealVector mySumX;
    unsigned int myN;

    ///Grid step
    double myH;
  };
  }

} // namespace DGtal



#endif // !defined IntegralInvariantBarycenterEstimatorFromSurfelsFunctors_h

#undef IntegralInvariantBarycenterEstimatorFromSurfelsFunctors_RECURSES
#endif // else defined(IntegralInvariantBarycenterEstimatorFromSurfelsFunctors_RECURSES)
