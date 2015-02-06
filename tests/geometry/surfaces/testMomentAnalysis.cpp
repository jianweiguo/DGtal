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
#include "DGtal/geometry/surfaces/estimation/IIGeometricFunctors.h"
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantBarycenterEstimator.h"
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantCovarianceEstimator.h"

#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/Color.h"
#include "DGtal/io/colormaps/GradientColorMap.h"

#include "DGtal/shapes/EuclideanShapesDecorator.h"


///////////////////////////////////////////////////////////////////////////////


using namespace DGtal;

///////////////////////////////////////////////////////////////////////////////
// Functions for testing class IntegralInvariantCovarianceEstimator and IIGeometricFunctor.
///////////////////////////////////////////////////////////////////////////////

bool testCube( double radius, double alpha, double beta, int argc, char** argv )
{
  typedef ImageContainerBySTLVector<Z3i::Domain,unsigned char> ImageVol;
  typedef functors::IntervalForegroundPredicate<ImageVol> DigitalShape;

  trace.info() << "InputFile" << std::endl;
  ImageVol image = GenericReader<ImageVol>::import ( "bunny.vol" );//VolReader<ImageVol>::importVol(vm["input-file"].as< std::string >());
  DigitalShape dshape( image, 0, 255 );

  Z3i::Domain domain = image.domain();
  trace.info() << domain << std::endl;

  // typedef ImplicitHyperCube<Z3i::Space> Cube;
  // typedef ImplicitBall<Z3i::Space> Sphere;
  // typedef EuclideanShapesUnion< Cube, Sphere > CubeSphere;
  // typedef GaussDigitizer<Z3i::Space, CubeSphere> DigitalShape;
  // typedef LightImplicitDigitalSurface<Z3i::KSpace,DigitalShape> Boundary;
  
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

  // Cube cube( Z3i::RealPoint( 0, 0, 0 ), (3.0/4.0)*radius );
  // Sphere sphere( Z3i::RealPoint( 0, 0, 0 ), radius );
  // CubeSphere cubesphere( cube, sphere );
  // DigitalShape dshape;
  // dshape.attach( cubesphere );
  // dshape.init( Z3i::RealPoint( -15.0, -15.0, -15.0 ), Z3i::RealPoint( 15.0, 15.0, 15.0 ), h );

  Z3i::KSpace K;
  // if ( !K.init( dshape.getLowerBound(), dshape.getUpperBound(), true ) )
  // {
  //   trace.error() << "Problem with Khalimsky space" << std::endl;
  //   return false;
  // }
  if ( !K.init( domain.lowerBound(), domain.upperBound(), true ) )
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

  unsigned int nbTest = 1;
  // std::vector<std::vector< double > > results(nbTest);
  std::vector< double > results;

  for( unsigned int i = 0; i < nbTest; ++i )
  {
    double re = 5;//(6.0 + i) * h;

    std::vector< Value_k > results_k;

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
  
      std::back_insert_iterator< std::vector< Value_k > > resultsIt( results_k );
      kEstimator.eval( v_border.begin(), v_border.end(), resultsIt );
  
      trace.endBlock();
    }

    {
      trace.beginBlock( "Barycenter speed estimator initialisation ...");
      
      // VisitorRange range( new Visitor( surf, *surf.begin() ));
      // VisitorConstIterator ibegin = range.begin();
      // VisitorConstIterator iend = range.end();

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
      // while(ibegin != iend)
      // {
      //   *resultsIt = barycenterEstimator.eval(ibegin);
      //   ++resultsIt;
      //   ++ibegin;
      // }
      barycenterEstimator.eval( v_border.begin(), v_border.end(), resultsIt );

      for(unsigned int j = 0; j < results_bar.size(); ++j )
      {
        if( results_k[j][1] > results_k[j][2] )
        {
          trace.error() << "Wrong sens lambda1 lambda2 " << results_k[j][0] << " " << results_k[j][1] << " " << results_k[j][2] << std::endl;
        }

        double value = results_bar[j].norm();
        // double value = 1.0 / ( alpha + beta * ( results_bar[j].norm() / re ) * ( results_bar[j].norm() / re ));
        // double value = 1.0 / ( alpha + beta * (( results_bar[j].norm() * results_k[j][0] ) / ( re * results_k[j][2] )) * (( results_bar[j].norm() * results_k[j][0] ) / ( re * results_k[j][2] )));
        // double value = 1.0 / ( alpha + beta * (( results_bar[j].norm() * results_k[j][2] ) / ( re * results_k[j][0] )) * (( results_bar[j].norm() * results_k[j][2] ) / ( re * results_k[j][0] )));
        results.push_back(value);
      }

      trace.endBlock();
    }
  }


  {
    trace.beginBlock ( "Comparing results ..." );
    double maxval = *std::max_element(results.begin(), results.end());
    double minval = *std::min_element(results.begin(), results.end());
    trace.info() << "Min/max= "<< minval<<"/"<<maxval<<std::endl;
    QApplication application( argc, argv );
    typedef Viewer3D<Z3i::Space, Z3i::KSpace> Viewer;
    Viewer viewer( K );
    viewer.setWindowTitle("Features from Tensor Voting");
    viewer.show();

    typedef GradientColorMap< double > Gradient;
    Gradient cmap_grad( minval, maxval );
    cmap_grad.addColor( Color( 255, 50, 50 ) );
    cmap_grad.addColor( Color( 50, 255, 50 ) );
    // cmap_grad.addColor( Color( 255, 255, 10 ) );

    // VisitorRange range( new Visitor( surf, *surf.begin() ));
    // VisitorConstIterator ibegin = range.begin();
    // VisitorConstIterator iend = range.end();


    viewer << SetMode3D((*(v_border.begin())).className(), "Basic" );
    unsigned int i=0;
    auto it_value = results.begin();
    
    for( auto it = v_border.begin(), itend=v_border.end();
      it!= itend;
      ++it, ++it_value, ++i)
    {
      viewer << CustomColors3D( Color::Black, cmap_grad( *it_value ))
             <<  *it ;

      trace.progressBar( i, v_border.size() );
    }
    
    
    viewer << Viewer3D<>::updateDisplay;
    
    trace.endBlock();
    application.exec();
  }
}

///////////////////////////////////////////////////////////////////////////////
// Standard services - public :

int main( int argc, char** argv )
{
  trace.beginBlock ( "Testing class IntegralInvariantBarycenterEstimator and 3d functors" );
    bool res = testCube( 10, 1, 5, argc, argv );// && testSphere( 10 );
    trace.emphase() << ( res ? "Passed." : "Error." ) << std::endl;
  trace.endBlock();
  return res ? 0 : 1;
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
