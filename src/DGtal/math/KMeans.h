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
 * @file IntegralInvariantGaussianCurvatureEstimator.h
 * @author Jeremy Levallois (\c jeremy.levallois@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), INSA-Lyon, France
 * LAboratoire de MAthématiques - LAMA (CNRS, UMR 5127), Université de Savoie, France
 *
 * @date 2014/03/03
 *
 * Header file for module KMeans.ih
 *
 * @brief 
 *
 * This file is part of the DGtal library.
 */

#if defined(KMeans_RECURSES)
#error Recursive header files inclusion detected in KMeans.h
#else // defined(KMeans_RECURSES)
/** Prevents recursive inclusion of headers. */
#define KMeans_RECURSES

#if !defined KMeans_h
/** Prevents repeated inclusion of headers. */
#define KMeans_h

#include <iostream>
#include "DGtal/base/Common.h"

namespace DGtal
{
   /**
   * @brief The KMeans class
   */
  class KMeans
  {
  public:
    /**
     * @brief compute K-Means algorithm on set of Cells
     * @details [long description]
     * 
     * @tparam Quantity a type of quantity (Point, Cell, double, etc.)
     * @tparam DistQuantityFunctor a type of distance functor between two quantities
     * 
     * @param[in] v_inputQuantities input vector of quantities. (size: n = v_inputQuantities.size())
     * @param[in] nbBins number of bins (each interval is called a bin).
     * @param[in] distFunctor functor of distance between two quantities.
     * @param[out] v_registration vector of mapping quantities - bin (each value of the vector is the position of the corresponding bin inside v_centroid) (size: n)
     * @param[out] v_centroid vector of bins (size: nbBins)
     */
    template< typename Quantity, typename DistQuantityFunctor >
    static void computeKMeans(
        const std::vector< Quantity > &v_inputQuantities,
        const Dimension nbBins,
        const DistQuantityFunctor &distFunctor,
        std::vector< Dimension > &v_registration,
        std::vector< Quantity > &v_centroid
        );
  };
} // namespace DGtal

///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "DGtal/math/KMeans.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined KMeans_h

#undef KMeans_RECURSES
#endif // else defined(KMeans_RECURSES)
