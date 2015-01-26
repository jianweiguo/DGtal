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
 * @file IntegralInvariantMeanCurvatureFromSurfaceAreaEstimator.h
 * @brief Computes the true quantity to each element of a range associated to a parametric shape.
 * @author David Coeurjolly (\c david.coeurjolly@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 *
 * @date 2015/01/24
 *
 * Header file for module IntegralInvariantMeanCurvatureFromSurfaceAreaEstimator.cpp
 *
 * This file is part of the DGtal library.
 *
 * @see testLengthEstimators.cpp, testTrueLocalEstimator.cpp
 */

#if defined(IntegralInvariantMeanCurvatureFromSurfaceAreaEstimator_RECURSES)
#error Recursive header files inclusion detected in IntegralInvariantMeanCurvatureFromSurfaceAreaEstimator.h
#else // defined(IntegralInvariantMeanCurvatureFromSurfaceAreaEstimator_RECURSES)
/** Prevents recursive inclusion of headers. */
#define IntegralInvariantMeanCurvatureFromSurfaceAreaEstimator_RECURSES

#if !defined IntegralInvariantMeanCurvatureFromSurfaceAreaEstimator_h
/** Prevents repeated inclusion of headers. */
#define IntegralInvariantMeanCurvatureFromSurfaceAreaEstimator_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <DGtal/base/Common.h>
#include "DGtal/geometry/surfaces/estimation/CSurfelLocalEstimator.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{
  namespace functors
  {
    /////////////////////////////////////////////////////////////////////////////
    // template class IntegralInvariantMeanCurvatureFromSurfaceAreaEstimator
    /**
     * Description of template class 'IntegralInvariantMeanCurvatureFromSurfaceAreaEstimator' <p>
     * \brief Aim: Estimate the mean curvature using integral
     * invariant estimator from local patch surface area. This
     * functor needs a normal vector estimator (or cache) to evaluate
     * the surface area of the given patch.
     *
     * model of CLocalEstimatorFromSurfelFunctor.
     *
     *
     * @tparam TSurfel type of surfels
     * @tparam TEmbedder type of functors which embed surfel to @f$ \mathbb{R}^3@f$
     * @tparam TNormalVectorEstimator the type of normal vector
     * estimator to consider (e.g. cached version with EstimatorCache class).
     *
     * @see testSphereFitting.cpp
     **/
    template <typename TSurfel,
              typename TEmbedder,
              typename TNormalVectorEstimator>
    class IntegralInvariantMeanCurvatureFromSurfaceAreaEstimator
    {
    public:

      ///Surfel type
      typedef TSurfel Surfel;
      ///Surfel embedder
      typedef TEmbedder SCellEmbedder;
      ///Real vector to represent normal vectors
      typedef typename SCellEmbedder::RealPoint RealPoint;
      ///Normal vector estimator type
      typedef TNormalVectorEstimator NormalVectorEstimator;

      //BOOST_CONCEPT_ASSERT(( concepts::CSurfelLocalEstimator<TNormalVectorEstimator> ));
      
      
      ///Quantity return type
      typedef double Quantity;
  

      /**
       * Constructor.
       *
       * @param [in] anEmbedder embedder to map surfel to R^n.
       * @param [in] h gridstep.
       * @param [in] radius integration ball radius (Euclidean)
       * @param [in] anEstimator a normal vector estimator
       */
      IntegralInvariantMeanCurvatureFromSurfaceAreaEstimator(ConstAlias<SCellEmbedder> anEmbedder,
                                                             const double h,
                                                             const double radius,
                                                             ConstAlias<NormalVectorEstimator> anEstimator):
        myEmbedder(&anEmbedder), myH2(h*h), myRadius(radius), myNormalEsitmatorCache(&anEstimator), myArea(0.0)
      { }


      /**
       * Destructor.
       */
      ~IntegralInvariantMeanCurvatureFromSurfaceAreaEstimator( )
      { }
      
      /**
       * Add the geometrical embedding of a surfel to the point list
       *
       * @param [in] aSurf a surfel to add
       * @param [in] aDistance of aSurf to the neighborhood boundary
       */
      void pushSurfel(const Surfel & aSurf,
                      const double aDistance)
      {
        BOOST_VERIFY(aDistance==aDistance);
        
        RealPoint elementary;
        Dimension i = myEmbedder->space().sOrthDir ( aSurf );
        elementary[ i ] = myEmbedder->space().sDirect ( aSurf, i ) ? 1 : -1;
        RealPoint estimatedNormal = myNormalEsitmatorCache->eval( &aSurf );          

        myArea += elementary.dot(estimatedNormal);
      }

      /**
       * Evaluate the surfaca area of the local patch.
       *
       * @return the surface area
       */
      Quantity eval() const
      {
        return (2.0/myRadius - myArea*myH2 / (M_PI*myRadius*myRadius*myRadius)) ;
      }
                             
      
      /**
       * Reset the point list.
       *
       */
      void reset()
      {
        myArea = 0.0;
      }


    private:

      ///Alias of the geometrical embedder
      const SCellEmbedder * myEmbedder;

      ///Square of the grid step
      const double myH2;

      ///Ball radius
      const double myRadius;

      ///NormalVectorCache
      const NormalVectorEstimator *myNormalEsitmatorCache;

      ///Surface area
      double myArea;
      
    }; // end of class IntegralInvariantMeanCurvatureFromSurfaceAreaEstimator
  }
} // namespace DGtal


//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined IntegralInvariantMeanCurvatureFromSurfaceAreaEstimator_h

#undef IntegralInvariantMeanCurvatureFromSurfaceAreaEstimator_RECURSES
#endif // else defined(IntegralInvariantMeanCurvatureFromSurfaceAreaEstimator_RECURSES)
