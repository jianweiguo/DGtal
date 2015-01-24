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
 * @file LocalSurfaceAreaEstimator.h
 * @brief Computes the true quantity to each element of a range associated to a parametric shape.
 * @author David Coeurjolly (\c david.coeurjolly@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 *
 * @date 2015/01/24
 *
 * Header file for module LocalSurfaceAreaEstimator.cpp
 *
 * This file is part of the DGtal library.
 *
 * @see testLengthEstimators.cpp, testTrueLocalEstimator.cpp
 */

#if defined(LocalSurfaceAreaEstimator_RECURSES)
#error Recursive header files inclusion detected in LocalSurfaceAreaEstimator.h
#else // defined(LocalSurfaceAreaEstimator_RECURSES)
/** Prevents recursive inclusion of headers. */
#define LocalSurfaceAreaEstimator_RECURSES

#if !defined LocalSurfaceAreaEstimator_h
/** Prevents repeated inclusion of headers. */
#define LocalSurfaceAreaEstimator_h

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
    // template class LocalSurfaceAreaEstimator
    /**
     * Description of template class 'LocalSurfaceAreaEstimator' <p>
     * \brief Aim: Estimate the surface area of a local patch. This
     * estimation is given by summing a dot product between estimated
     * normal vectors and surfel elementary normal vector.
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
    class LocalSurfaceAreaEstimator
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
       * @param [in] anEstimator a normal vector estimator
       */
      LocalSurfaceAreaEstimator(ConstAlias<SCellEmbedder> anEmbedder,
                                const double h,
                                ConstAlias<NormalVectorEstimator> anEstimator):
        myEmbedder(&anEmbedder), myH(h), myNormalEsitmatorCache(&anEstimator), myArea(0.0)
      { }


      /**
       * Destructor.
       */
      ~LocalSurfaceAreaEstimator( )
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
        RealPoint estimatedNormal = myNormalEsitmatorCache->eval( aSurf );          

        myArea += elementary.dot(estimatedNormal);
      }

      /**
       * Evaluate the surfaca area of the local patch.
       *
       * @return the surface area
       */
      Quantity eval() const
      {
        return myArea*myH*myH;
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

      ///Grid step
      double myH;

      ///NormalVectorCache
      const NormalVectorEstimator *myNormalEsitmatorCache;

      ///Surface area
      double myArea;
      
    }; // end of class LocalSurfaceAreaEstimator
  }
} // namespace DGtal


//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined LocalSurfaceAreaEstimator_h

#undef LocalSurfaceAreaEstimator_RECURSES
#endif // else defined(LocalSurfaceAreaEstimator_RECURSES)
