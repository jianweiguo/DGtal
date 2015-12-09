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

/**
 * @file curvature-estimation-II-2d.cpp
 * @ingroup Examples
 * @author Jérémy Levallois (\c jeremy.levallois@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), INSA-Lyon, France
 * LAboratoire de MAthématiques - LAMA (CNRS, UMR 5127), Université de Savoie, France
 *
 * @date 2012/12/17
 *
 * An example file named curvature-estimation-II-2d.
 * Documentation: \link http://liris.cnrs.fr/dgtal/doc/stable/moduleIntegralInvariant.html \endlink
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "DGtal/base/Common.h"
// Shape construction
#include "DGtal/shapes/parametric/Flower2D.h"
// Digitization
#include "DGtal/shapes/GaussDigitizer.h"
// Digital surface construction
#include "DGtal/topology/LightImplicitDigitalSurface.h"
#include "DGtal/topology/DigitalSurface.h"
// Digital surface exploration
#include "DGtal/graph/DepthFirstVisitor.h"
#include "DGtal/graph/GraphVisitorRange.h"
// Estimator
#include "DGtal/geometry/surfaces/estimation/IIGeometricFunctors.h"
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantVolumeEstimator.h"
// Drawing
#include "DGtal/io/boards/Board2D.h"
#include "DGtal/io/colormaps/GradientColorMap.h"

///////////////////////////////////////////////////////////////////////////////
using namespace DGtal;
///////////////////////////////////////////////////////////////////////////////

int main( int /* argc */, char** /* argv */ )
{
    trace.beginBlock ( "Example curvature-estimation-II-2d" );

    trace.beginBlock( "Shape construction..." );
    typedef Flower2D< Z2i::Space > MyShape;

    MyShape shape( 0, 0, 20, 10, 6, 3.0 ); // => x_0, y_0, radius, small_radius, nb_extremeties, phase
    trace.endBlock();

    trace.beginBlock( "Digitization..." );
    typedef GaussDigitizer< Z2i::Space, MyShape > MyDigitalShape;

    double h = 0.5;
    MyDigitalShape digShape;
    digShape.attach( shape );
    digShape.init( shape.getLowerBound() + Z2i::Point::diagonal(-1), shape.getUpperBound() + Z2i::Point::diagonal(1), h );
    Z2i::Domain domainShape = digShape.getDomain();
    Z2i::KSpace K;
    K.init( domainShape.lowerBound(), domainShape.upperBound(), true );
    trace.endBlock();

    trace.beginBlock( "Digital surface construction..." );
    typedef LightImplicitDigitalSurface< Z2i::KSpace, MyDigitalShape > MyBoundary;
    typedef DigitalSurface< MyBoundary > MyDigitalSurface;
    typedef typename MyDigitalSurface::ConstIterator MyConstIterator;

    Z2i::SCell bel = Surfaces< Z2i::KSpace >::findABel( K, digShape, 10000 );
    MyBoundary boundary( K, digShape, SurfelAdjacency< Z2i::KSpace::dimension >( true ), bel );
    MyDigitalSurface surf( boundary );
    trace.endBlock();

    trace.beginBlock( "Digital surface exploration..." );
    typedef DepthFirstVisitor< MyDigitalSurface > MyVisitor;
    typedef GraphVisitorRange< MyVisitor > MyVisitorRange;
    typedef typename MyVisitorRange::ConstIterator MyVisitorConstIterator;
    typedef Z2i::KSpace::SurfelSet MySet;

    MyVisitorRange range( new MyVisitor( surf, *surf.begin() ));
    MyVisitorConstIterator abegin = range.begin();
    MyVisitorConstIterator aend = range.end();
    MySet set;
    set.insert( abegin, aend );
    trace.endBlock();

    trace.beginBlock( "Integral Invariant computation..." );
    //! [IntegralInvariantUsage]
    double radius = 5.0; // (Euclidean) radius of the convolution kernel. This is the only parameter of this method.

    typedef functors::IICurvatureFunctor< Z2i::Space > MyIICurvatureFunctor;
    typedef IntegralInvariantVolumeEstimator< Z2i::KSpace, MyDigitalShape, MyIICurvatureFunctor > MyIICurvatureEstimator;
    typedef MyIICurvatureFunctor::Value Value;

    MyIICurvatureFunctor functor; // Functor used to convert volume -> curvature
    functor.init( h, radius ); // Initialisation for a grid step and a given Euclidean radius of convolution kernel

    MyIICurvatureEstimator estimator( functor );
    estimator.attach( K, digShape ); // Setting a KSpace and a predicate on the object to evaluate
    estimator.setParams( radius/h ); // Setting the digital radius of the convolution kernel. Used to precompute some kernels.
    estimator.init( h, set.begin(), set.end() ); // Initialisation for a given h

    std::vector< Value > results;
    std::back_insert_iterator< std::vector< Value > > resultsIt( results ); // output iterator for results of Integral Invariant curvature computation
    estimator.eval( set.begin(), set.end(), resultsIt ); // Computation for a range of surfels
    //! [IntegralInvariantUsage]
    trace.endBlock();

    trace.beginBlock( "Drawing results..." );
    auto minmax = std::minmax_element(std::begin(results), std::end(results));
    Board2D board;

    typedef GradientColorMap< Value > Gradient;
    Gradient cmap_grad( *minmax.first, *minmax.second );
    cmap_grad.addColor( Color( 50, 50, 255 ) );
    cmap_grad.addColor( Color( 255, 0, 0 ) );
    cmap_grad.addColor( Color( 255, 255, 10 ) );

    board << SetMode( (*set.begin()).className(), "Paving" );
    std::string specificStyle = (*set.begin()).className() + "/Paving";
    MySet::const_iterator it = set.begin();
    for( Value &val : results )
    {
        Z2i::KSpace::SCell currentCell = K.sIndirectIncident( *it, *K.sOrthDirs( *it ) ); // We apply the color to the inner spel (more visible than surfel)
        board << CustomStyle( specificStyle, new CustomColors( Color::Black, cmap_grad( val )))
              << currentCell;
        ++it;
    }
    board.saveSVG ( "curvature-estimation-II-2d.svg" );
    trace.endBlock();

    trace.endBlock();
    return 0;
}
