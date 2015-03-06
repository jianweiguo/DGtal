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
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/images/IntervalForegroundPredicate.h"

 /// Digitization
#include "DGtal/shapes/GaussDigitizer.h"
#include "DGtal/topology/LightImplicitDigitalSurface.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/graph/DepthFirstVisitor.h"
#include "DGtal/graph/GraphVisitorRange.h"

/// Estimator
#include "DGtal/geometry/surfaces/estimation/LocalEstimatorFromSurfelFunctorAdapter.h"
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantBarycenterEstimator.h"

#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/boards/Board3D.h"
#include "DGtal/io/Color.h"
#include "DGtal/io/colormaps/GradientColorMap.h"

#include "DGtal/shapes/EuclideanShapesDecorator.h"


///////////////////////////////////////////////////////////////////////////////


using namespace DGtal;

///////////////////////////////////////////////////////////////////////////////
// Functions for testing class IntegralInvariantCovarianceEstimator and IIGeometricFunctor.
///////////////////////////////////////////////////////////////////////////////


bool testT(const double radius, const double radius_kernel, const double alpha, const double beta, int argc, char** argv)
{

  typedef ImplicitHyperCube<Z3i::Space> Cube;
  typedef ImplicitBall<Z3i::Space> Sphere;
  typedef EuclideanShapesCSG< Cube, Sphere > Shape;
  typedef GaussDigitizer<Z3i::Space,Shape> Gauss;

  typedef LightImplicitDigitalSurface<Z3i::KSpace,Gauss> SurfaceContainer;
  typedef DigitalSurface<SurfaceContainer> Surface;
  //typedef Surface::SurfelConstIterator ConstIterator;
  //typedef Surface::Tracker Tracker;
  typedef typename Surface::Surfel Surfel;

  const double h = 1.0;

  trace.beginBlock("Creating Surface");
  Z3i::Point p1( -(1.5*radius + 20 ), -(1.5*radius + 20 ), -(1.5*radius + 20 ) );
  Z3i::Point p2( (1.5*radius + 20 ), (1.5*radius + 20 ), (1.5*radius + 20 ) );
  Z3i::KSpace K;
  if( !K.init( p1, p2, true ))
    return false;

  //Shape
  Cube cube( Z3i::RealPoint( 0, 0, 0 ), radius );
  Sphere sphere( Z3i::RealPoint( radius, 0, 0 ), radius/2.0 );
  Shape shape( cube );
  shape.op_union( sphere );
  Gauss gauss;
  gauss.attach(shape);
  gauss.init(p1,p2,h);

  //Surface
  Surfel bel = Surfaces<Z3i::KSpace>::findABel( K, gauss, 10000 );
  SurfaceContainer* surfaceContainer = new SurfaceContainer
    ( K, gauss, SurfelAdjacency<Z3i::KSpace::dimension>( true ), bel );
  Surface surface( surfaceContainer ); // acquired
  trace.endBlock();

  trace.beginBlock("Creating  adapters");
  typedef functors::IntegralInvariantBarycenterEstimatorFromSurfels<Surfel, CanonicSCellEmbedder<Z3i::KSpace> > Functor;
  typedef Functor::Quantity Quantity;

  typedef functors::ConstValue< double > ConvFunctor;
  typedef LocalEstimatorFromSurfelFunctorAdapter<SurfaceContainer, Z3i::L2Metric, Functor, ConvFunctor> Reporter;

  CanonicSCellEmbedder<Z3i::KSpace> embedder(surface.container().space());
  Functor estimator(embedder,h);

  ConvFunctor convFunc(1.0);
  Reporter reporter;
  reporter.attach(surface);
  reporter.setParams(Z3i::l2Metric, estimator , convFunc, radius_kernel/h);

  reporter.init(h, surface.begin() , surface.end());

  std::vector<Quantity> values;
  reporter.eval( surface.begin(), surface.end(), std::back_insert_iterator<std::vector<Quantity> >(values));

  std::vector<double> distance;
  for(uint i = 0; i < values.size(); ++i)
  {
    double distanceBary = (values[i].barycenter*h - values[i].center*h).norm();

    // double a = ( distanceBary * values[i].eigenvalues[0] ) / ( radius_kernel * values[i].eigenvalues[2] );
    // double dist = 1.0 / ( alpha + beta * a * a);

    // double dist = distanceBary;
    
    double dist = 1.0 / ( alpha + beta * ( distanceBary / radius_kernel ) * ( distanceBary / radius_kernel ));

    // if( a != a )
    // {
    //   trace.error() << "--------------------" << std::endl;
    //   // trace.error() << distanceBary << std::endl;
    //   trace.error() << values[i].eigenvalues[0] << " " << values[i].eigenvalues[2] << std::endl;
    //   // trace.error() << " a != a " << a << std::endl;
    //   dist = 0;
    // }

    distance.push_back( dist );

    // trace.info() << values[i].eigenvalues[0] << " " << values[i].eigenvalues[1] << " " << values[i].eigenvalues[2] << std::endl;
  }

  //distance[0] = 0;

  double maxval = *std::max_element(distance.begin(), distance.end());
  double minval = *std::min_element(distance.begin(), distance.end());
  trace.info() << "Min/max= "<< minval<<"/"<<maxval<<std::endl;
  QApplication application( argc, argv );
  typedef Viewer3D<Z3i::Space, Z3i::KSpace> Viewer;
  Viewer viewer( K );
  viewer.setWindowTitle("Features from Moment Analysis");
  viewer.show();

  typedef GradientColorMap< double > Gradient;
  Gradient cmap_grad( minval, maxval );
  cmap_grad.addColor( Color( 50, 50, 255 ) );
  cmap_grad.addColor( Color( 50, 255, 50 ) );
  // cmap_grad.addColor( Color( 255, 255, 10 ) );
  cmap_grad.addColor( Color( 255, 0, 0 ) );

  // cmap_grad.addColor( Color( 255, 50, 50 ) );
  // cmap_grad.addColor( Color( 50, 255, 50 ) );


  viewer << SetMode3D((*(surface.begin())).className(), "Basic" );
  
  unsigned int i=0;
  for(typename Surface::ConstIterator it = surface.begin(), itend=surface.end();
      it!= itend;
      ++it, ++i)
    {
      viewer << CustomColors3D( Color::Black, cmap_grad( distance[ i ] ))
             <<  *it ;    
    }
  
  
  viewer << Viewer3D<>::updateDisplay;
  
  trace.endBlock();
  application.exec();

  return true;
}

bool testT2(const double radius, const double radius_kernel, const double alpha, const double beta, int argc, char** argv)
{

  typedef ImplicitBall<Z3i::Space> Cube;
  typedef ImplicitBall<Z3i::Space> Sphere;
  typedef EuclideanShapesCSG< Cube, Sphere > Shape;
  typedef GaussDigitizer<Z3i::Space,Shape> Gauss;

  typedef LightImplicitDigitalSurface<Z3i::KSpace,Gauss> SurfaceContainer;
  typedef DigitalSurface<SurfaceContainer> Surface;
  //typedef Surface::SurfelConstIterator ConstIterator;
  //typedef Surface::Tracker Tracker;
  typedef typename Surface::Surfel Surfel;

  const double h = 1.0;

  trace.beginBlock("Creating Surface");
  Z3i::Point p1( -(3*radius + 20 ), -(3*radius + 20 ), -(3*radius + 20 ) );
  Z3i::Point p2( (3*radius + 20 ), (3*radius + 20 ), (3*radius + 20 ) );
  Z3i::KSpace K;
  if( !K.init( p1, p2, true ))
    return false;

  //Shape
  Cube cube( Z3i::RealPoint( 0, 0, 0 ), radius );
  Sphere sphere( Z3i::RealPoint( radius, 0, 0 ), radius/2.0 );
  Sphere sphere2( Z3i::RealPoint( radius + radius/2.0, 0, 0 ), radius/4.0 );
  Sphere sphere3( Z3i::RealPoint( radius + radius/2.0 + radius / 4.0, 0, 0 ), radius/8.0 );
  Shape shape( cube );
  shape.op_union( sphere );
  shape.op_union( sphere2 );
  shape.op_union( sphere3 );
  Gauss gauss;
  gauss.attach(shape);
  gauss.init(p1,p2,h);

  //Surface
  Surfel bel = Surfaces<Z3i::KSpace>::findABel( K, gauss, 10000 );
  SurfaceContainer* surfaceContainer = new SurfaceContainer
    ( K, gauss, SurfelAdjacency<Z3i::KSpace::dimension>( true ), bel );
  Surface surface( surfaceContainer ); // acquired
  trace.endBlock();

  trace.beginBlock("Creating  adapters");
  typedef functors::IntegralInvariantBarycenterEstimatorFromSurfels<Surfel, CanonicSCellEmbedder<Z3i::KSpace> > Functor;
  typedef Functor::Quantity Quantity;

  typedef functors::ConstValue< double > ConvFunctor;
  typedef LocalEstimatorFromSurfelFunctorAdapter<SurfaceContainer, Z3i::L2Metric, Functor, ConvFunctor> Reporter;

  CanonicSCellEmbedder<Z3i::KSpace> embedder(surface.container().space());
  Functor estimator(embedder,h);

  ConvFunctor convFunc(1.0);
  Reporter reporter;
  reporter.attach(surface);
  reporter.setParams(Z3i::l2Metric, estimator , convFunc, radius_kernel/h);

  reporter.init(h, surface.begin() , surface.end());

  std::vector<Quantity> values;
  reporter.eval( surface.begin(), surface.end(), std::back_insert_iterator<std::vector<Quantity> >(values));

  std::vector<double> distance;
  for(uint i = 0; i < values.size(); ++i)
  {
    double distanceBary = (values[i].barycenter*h - values[i].center*h).norm();

    // double a = ( distanceBary * values[i].eigenvalues[0] ) / ( radius_kernel * values[i].eigenvalues[2] );
    // double dist = 1.0 / ( alpha + beta * a * a);

    // double dist = distanceBary;
    
    double dist = 1.0 / ( alpha + beta * ( distanceBary / radius_kernel ) * ( distanceBary / radius_kernel ));

    // if( a != a )
    // {
    //   trace.error() << "--------------------" << std::endl;
    //   // trace.error() << distanceBary << std::endl;
    //   trace.error() << values[i].eigenvalues[0] << " " << values[i].eigenvalues[2] << std::endl;
    //   // trace.error() << " a != a " << a << std::endl;
    //   dist = 0;
    // }

    distance.push_back( dist );

    // trace.info() << values[i].eigenvalues[0] << " " << values[i].eigenvalues[1] << " " << values[i].eigenvalues[2] << std::endl;
  }

  //distance[0] = 0;

  double maxval = *std::max_element(distance.begin(), distance.end());
  double minval = *std::min_element(distance.begin(), distance.end());
  trace.info() << "Min/max= "<< minval<<"/"<<maxval<<std::endl;
  QApplication application( argc, argv );
  typedef Viewer3D<Z3i::Space, Z3i::KSpace> Viewer;
  Viewer viewer( K );
  viewer.setWindowTitle("Features from Moment Analysis");
  viewer.show();

  typedef GradientColorMap< double > Gradient;
  Gradient cmap_grad( minval, maxval );
  cmap_grad.addColor( Color( 50, 50, 255 ) );
  cmap_grad.addColor( Color( 50, 255, 50 ) );
  // cmap_grad.addColor( Color( 255, 255, 10 ) );
  cmap_grad.addColor( Color( 255, 0, 0 ) );

  // cmap_grad.addColor( Color( 255, 50, 50 ) );
  // cmap_grad.addColor( Color( 50, 255, 50 ) );


  viewer << SetMode3D((*(surface.begin())).className(), "Basic" );
  
  unsigned int i=0;
  for(typename Surface::ConstIterator it = surface.begin(), itend=surface.end();
      it!= itend;
      ++it, ++i)
    {
      viewer << CustomColors3D( Color::Black, cmap_grad( distance[ i ] ))
             <<  *it ;    
    }
  
  
  viewer << Viewer3D<>::updateDisplay;
  
  trace.endBlock();
  application.exec();

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

  double radius = 100;
  double radius_kernel = 15;
  double alpha = 1;
  double beta = 5;
  if( argc == 5 )
  {
    radius = atof(argv[1]);
    radius_kernel = atof(argv[2]);
    alpha = atof(argv[3]);
    beta = atof(argv[4]);
  }

  bool res = testT2( radius, radius_kernel, alpha, beta, argc, argv );
  trace.emphase() << ( res ? "Passed." : "Error." ) << std::endl;
  trace.endBlock();
  return res ? 0 : 1;
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
