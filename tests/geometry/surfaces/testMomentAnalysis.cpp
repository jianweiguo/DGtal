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
 * @file testIntegralInvariantCovarianceEstimator.cpp
 * @ingroup Tests
 * @author Jérémy Levallois (\c jeremy.levallois@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), INSA-Lyon, France
 * LAboratoire de MAthématiques - LAMA (CNRS, UMR 5127), Université de Savoie, France
 *
 * @date 2014/06/26
 *
 * Functions for testing class IntegralInvariantCovarianceEstimator and IIGeometricFunctor.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "DGtal/base/Common.h"

 /// Shape
#include "DGtal/shapes/implicit/ImplicitBall.h"
#include "DGtal/shapes/implicit/ImplicitHyperCube.h"

 /// Digitization
#include "DGtal/shapes/GaussDigitizer.h"
#include "DGtal/topology/LightImplicitDigitalSurface.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/graph/DepthFirstVisitor.h"
#include "DGtal/graph/GraphVisitorRange.h"

/// Estimator
#include "DGtal/geometry/surfaces/estimation/IIGeometricFunctors.h"
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantBarycenterEstimator.h"
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantCovarianceEstimator.h"

#include "DGtal/shapes/EuclideanShapesDecorator.h"


///////////////////////////////////////////////////////////////////////////////


using namespace DGtal;

///////////////////////////////////////////////////////////////////////////////
// Functions for testing class IntegralInvariantCovarianceEstimator and IIGeometricFunctor.
///////////////////////////////////////////////////////////////////////////////

bool testCubeSphere()
{
  typedef ImplicitHyperCube<Z3i::Space> Cube;
  typedef ImplicitBall<Z3i::Space> Sphere;
  typedef EuclideanShapesCSG< Cube, Sphere > CubeSphere;
  typedef GaussDigitizer<Z3i::Space, CubeSphere> DigitalShape;
  typedef LightImplicitDigitalSurface<Z3i::KSpace,DigitalShape> Boundary;
  typedef DigitalSurface< Boundary > MyDigitalSurface;
  typedef DepthFirstVisitor< MyDigitalSurface > Visitor;
  typedef GraphVisitorRange< Visitor > VisitorRange;
  typedef VisitorRange::ConstIterator VisitorConstIterator;

  typedef functors::IIBarycenterSpeedFunctor<Z3i::Space> MyIIBarycenterSpeedFunctor;
  typedef IntegralInvariantBarycenterEstimator< Z3i::KSpace, DigitalShape, MyIIBarycenterSpeedFunctor > MyIIBarycenterEstimator;
  typedef MyIIBarycenterSpeedFunctor::Value Value;

  typedef functors::IIEigenvalues3DFunctor<Z3i::Space> MyIIPrincipalCurvaturesFunctor;
  typedef IntegralInvariantCovarianceEstimator< Z3i::KSpace, DigitalShape, MyIIPrincipalCurvaturesFunctor > MyIIPrincipalCurvatureEstimator;
  typedef MyIIPrincipalCurvaturesFunctor::Value Value_k;

  const double h  = 1;//0.5;

  trace.beginBlock( "Shape initialisation ..." );

  Z3i::Point p1( -100, -100, -100 );
  Z3i::Point p2( 100, 100, 100 );
  Cube cube( Z3i::RealPoint( 0, 0, 0 ), 30 );
  Sphere sphere( Z3i::RealPoint( 30, 0, 0 ), 15 );
  CubeSphere cubesphere( cube );
  cubesphere.op_union( sphere );
  DigitalShape dshape;
  dshape.attach( cubesphere );
  dshape.init( p1, p2, h );

  Z3i::KSpace K;
  if ( !K.init( p1, p2, true ) )
  {
    trace.error() << "Problem with Khalimsky space" << std::endl;
    return false;
  }

  Z3i::KSpace::Surfel bel = Surfaces<Z3i::KSpace>::findABel( K, dshape, 1000000 );
  Boundary boundary( K, dshape, SurfelAdjacency<Z3i::KSpace::dimension>( true ), bel );
  MyDigitalSurface surf ( boundary );
  trace.info() << "MyDigitalSurface " << surf.size() << std::endl;

  trace.endBlock();

  trace.beginBlock("Saving border");
  std::vector<Z3i::SCell> v_border;

  {
    VisitorRange range( new Visitor( surf, *surf.begin() ));
    VisitorConstIterator ibegin = range.begin();
    VisitorConstIterator iend = range.end();

    while( ibegin != iend )
    {
      v_border.push_back( *ibegin );
      ++ibegin;
    }
  }

  trace.endBlock();

  double re = 5;

  {
    trace.beginBlock( "Curvature estimator initialisation ...");

    MyIIPrincipalCurvaturesFunctor kFunctor;
    kFunctor.init( h, re );

    MyIIPrincipalCurvatureEstimator kEstimator( kFunctor );
    kEstimator.attach( K, dshape );
    kEstimator.setParams( re/h );
    kEstimator.init( h, v_border.begin(), v_border.end() );

    trace.endBlock();

    trace.beginBlock( "Curvature evaluation ...");

    std::vector< Value_k > results_k;
    std::back_insert_iterator< std::vector< Value_k > > resultsIt( results_k );
    kEstimator.eval( v_border.begin(), v_border.end(), resultsIt );

    trace.endBlock();
  }

  {
    trace.beginBlock( "Barycenter speed estimator initialisation ...");

    MyIIBarycenterSpeedFunctor barycenterSpeedFunctor;
    barycenterSpeedFunctor.init( h, re );

    MyIIBarycenterEstimator barycenterEstimator( barycenterSpeedFunctor );
    barycenterEstimator.attach( K, dshape );
    barycenterEstimator.setParams( re/h );
    barycenterEstimator.init( h, v_border.begin(), v_border.end() );

    trace.endBlock();

    trace.beginBlock( "Barycenter speed evaluation ...");

    std::vector< Value > results_bar;
    std::back_insert_iterator< std::vector< Value > > resultsIt( results_bar );
    barycenterEstimator.eval( v_border.begin(), v_border.end(), resultsIt );

    trace.endBlock();
  }


  return true;
}

///////////////////////////////////////////////////////////////////////////////
// Standard services - public :

int main( int argc, char** argv )
{
  trace.beginBlock ( "Testing class IntegralInvariantBarycenterEstimator and 3d functors (TestMomentAnalysisViewer)" );
  trace.info() << "Args:";
  for ( int i = 0; i < argc; ++i )
  {
    trace.info() << " " << argv[ i ];
  }
  trace.info() << std::endl;

  bool res = testCubeSphere();
  trace.emphase() << ( res ? "Passed." : "Error." ) << std::endl;
  trace.endBlock();
  return res ? 0 : 1;
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
