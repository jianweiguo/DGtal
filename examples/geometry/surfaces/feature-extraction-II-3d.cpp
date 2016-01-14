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
 * @file curvature-estimation-II-3d.cpp
 * @ingroup Examples
 * @author Jérémy Levallois (\c jeremy.levallois@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), INSA-Lyon, France
 * LAboratoire de MAthématiques - LAMA (CNRS, UMR 5127), Université de Savoie, France
 *
 * @date 2012/12/17
 *
 * An example file named curvature-estimation-II-3d.
 * Documentation: \link http://liris.cnrs.fr/dgtal/doc/stable/moduleIntegralInvariant.html \endlink
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "DGtal/base/Common.h"

// Shape construction
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "DGtal/images/IntervalForegroundPredicate.h"
#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/topology/LightImplicitDigitalSurface.h"
#include "DGtal/images/ImageHelper.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/graph/DepthFirstVisitor.h"
#include "DGtal/graph/GraphVisitorRange.h"
#include "ConfigExamples.h"

/// Estimator
#include "DGtal/geometry/surfaces/estimation/IIGeometricFunctors.h"
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantVolumeEstimator.h"

// Drawing
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/colormaps/GradientColorMap.h"

///////////////////////////////////////////////////////////////////////////////
using namespace DGtal;
///////////////////////////////////////////////////////////////////////////////

int main( int argc, char** argv )
{
    trace.beginBlock ( "Example IntegralInvariantCurvature3D" );

    std::string inputFilename = examplesPath + "samples/Al.100.vol";
    int thresholdMin = 0;
    int thresholdMax = 1;
    const double h = 1;

    trace.info() << "File             = " << inputFilename << std::endl;
    trace.info() << "Min image thres. = " << thresholdMin << std::endl;
    trace.info() << "Max image thres. = " << thresholdMax << std::endl;

    trace.beginBlock( "Loading image into memory." );
    typedef ImageSelector< Z3i::Domain, unsigned char >::curvature-estimation-II-3dcurvature-estimation-II-3dType MyImage;
    typedef functors::IntervalForegroundPredicate< MyImage > MyDigitalShape;

    MyImage image = GenericReader< MyImage >::import( inputFilename );
    MyDigitalShape digShape( image, thresholdMin, thresholdMax );
    Z3i::Domain domainShape = image.domain();
    Z3i::KSpace K;
    K.init( domainShape.lowerBound(), domainShape.upperBound(), true );
    trace.endBlock();

    trace.beginBlock( "Digital surface construction..." );
    typedef LightImplicitDigitalSurface< Z3i::KSpace, MyDigitalShape > MyBoundary;
    typedef DigitalSurface< MyBoundary > MyDigitalSurface;
    typedef typename MyDigitalSurface::ConstIterator MyConstIterator;

    Z3i::SCell bel = Surfaces< Z3i::KSpace >::findABel( K, digShape, 10000 );
    MyBoundary boundary( K, digShape, SurfelAdjacency< Z3i::KSpace::dimension >( true ), bel );
    MyDigitalSurface surf( boundary );
    trace.endBlock();

    trace.beginBlock( "Digital surface exploration..." );
    typedef DepthFirstVisitor< MyDigitalSurface > MyVisitor;
    typedef GraphVisitorRange< MyVisitor > MyVisitorRange;
    typedef typename MyVisitorRange::ConstIterator MyVisitorConstIterator;
    typedef Z3i::KSpace::SurfelSet MySet;

    MyVisitorRange range( new MyVisitor( surf, *surf.begin() ));
    MyVisitorConstIterator abegin = range.begin();
    MyVisitorConstIterator aend = range.end();
    MySet set;
    set.insert( abegin, aend );
    trace.endBlock();

    trace.beginBlock( "Integral Invariant computation..." );
    //! [IntegralInvariantUsage]
    const double radius = 5; // (Euclidean) radius of the convolution kernel. This is the only parameter of this method.

    typedef functors::IIMeanCurvature3DFunctor< Z3i::Space > MyIICurvatureFunctor; // Mean curvature functor. see @IIGeometricFunctors
    typedef IntegralInvariantVolumeEstimator< Z3i::KSpace, MyDigitalShape, MyIICurvatureFunctor > MyIICurvatureEstimator;
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

    QApplication application( argc, argv );
    typedef Viewer3D< Z3i::Space, Z3i::KSpace > Viewer;
    Viewer viewer( K );
    viewer.setWindowTitle("example Integral Invariant 3D");
    viewer.show();

    typedef GradientColorMap< Value > Gradient;
    Gradient cmap_grad( *minmax.first, *minmax.second );
    cmap_grad.addColor( Color( 50, 50, 255 ) );
    cmap_grad.addColor( Color( 255, 0, 0 ) );
    cmap_grad.addColor( Color( 255, 255, 10 ) );

    Z3i::KSpace::Cell dummy_cell;
    viewer << SetMode3D( dummy_cell.className(), "Basic" );
    MySet::const_iterator it = set.begin();

    for( Value &val : results )
    {
        viewer << CustomColors3D( Color::Black, cmap_grad( val ))
               << K.unsigns( *it );
        ++it;
    }
    trace.endBlock();

    viewer << Viewer3D<>::updateDisplay;

    trace.endBlock();
    return application.exec();
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
