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
 * @file exampleIntegralInvariantCurvature2D.cpp
 * @ingroup Examples
 * @author Jérémy Levallois (\c jeremy.levallois@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), INSA-Lyon, France
 * LAboratoire de MAthématiques - LAMA (CNRS, UMR 5127), Université de Savoie, France
 *
 * @date 2012/12/17
 *
 * An example file named exampleIntegralInvariantCurvature2D.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "DGtal/base/Common.h"

// Shape construction
#include "DGtal/shapes/parametric/Flower2D.h"
#include "DGtal/shapes/GaussDigitizer.h"
#include "DGtal/topology/LightImplicitDigitalSurface.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/graph/DepthFirstVisitor.h"
#include "DGtal/graph/GraphVisitorRange.h"

// Integral Invariant includes
#include "DGtal/images/ImageHelper.h"
#include "DGtal/geometry/surfaces/FunctorOnCells.h"
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantMeanCurvatureEstimator.h"

#include "DGtal/math/Statistic.h"
#include "DGtal/geometry/curves/ArithmeticalDSS.h"
#include "DGtal/geometry/curves/SaturatedSegmentation.h"
#include "DGtal/math/KMeans.h"

// Drawing
#include "DGtal/io/boards/Board2D.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/colormaps/RandomColorMap.h"
#include "DGtal/io/colormaps/HueShadeColorMap.h"



///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;

void checkSizeRadius( double & re,
                      const double h,
                      const double minRadiusAABB )
{
  if(( re / h ) < 5.0 ) /// 	ridiculously small radius check
  {
    re = 5.0 * h;
  }
  if( re > ( 0.75 * minRadiusAABB ))
  {
    re = 0.75 * minRadiusAABB;
  }
}

struct Euclidean
{
  double distance(const double &a, const double &b) const
  {
    return (double)std::sqrt((b-a)*(b-a));
  }

  template< typename Point >
  double distance(const Point &a, const Point &b) const
  {
    return (double)(b-a).norm();
  }
};


template <typename KSpace, typename Iterator>
void analyseAllLengthMS( std::vector< Statistic<double> > & statE,
                         Iterator itb, Iterator ite )
{
  typedef typename KSpace::Space Space;
  typedef typename Space::Point Point;
  typedef typename Space::Vector Vector;
  typedef ArithmeticalDSSComputer< Iterator, int, 4 > SegmentComputer;
  typedef SaturatedSegmentation< SegmentComputer > Decomposition;
  typedef typename Decomposition::SegmentComputerIterator SegmentComputerIterator;
  typedef std::vector< SegmentComputerIterator > VectorOfSegmentComputerIterator;
  typedef std::map< Point, VectorOfSegmentComputerIterator > Pmap;

  // Computes the tangential cover
  SegmentComputer algo;
  Decomposition theDecomposition( itb, ite, algo);

  Pmap map;

  Iterator itc = itb;
  //for( Iterator itc = itb; itc != ite; ++itc )
  do
  {
    map.insert( std::pair< Point, VectorOfSegmentComputerIterator >( *itc, VectorOfSegmentComputerIterator() ) );
    ++itc;
  } while ( itc != ite );


  for ( SegmentComputerIterator scIt = theDecomposition.begin(), scItEnd = theDecomposition.end();
        scIt != scItEnd; ++scIt )
  {
    const SegmentComputer & sc = *scIt;
    for ( Iterator ptIt = sc.begin(), ptItEnd = sc.end(); ptIt != ptItEnd; ++ptIt )
    {
      typename Pmap::iterator mloc = map.find( *ptIt );
      if( mloc != map.end() )
      {
        mloc->second.push_back( scIt );
      }
    }
  }

  Dimension ii = 0;

  itc = itb;
  //for( Iterator itc = itb; itc != ite; ++itc )
  do
  {
    //statD[ii].clear();
    statE[ii].clear();
    typename Pmap::iterator mloc = map.find( *itc );
    ASSERT(( mloc != map.end() ));

    /////////////
    for( typename VectorOfSegmentComputerIterator::iterator scIt = mloc->second.begin(), scItEnd = mloc->second.end(); scIt != scItEnd; ++scIt )
    {
      const SegmentComputer & sc = *(*scIt);
      /*int64_t l = 0;
          for ( Iterator ptIt = sc.begin(), ptItEnd = sc.end(); ptIt != ptItEnd; ++ptIt )
            ++l;
          statD[ii].addValue( (double) l );*/

      double v = (sc.back( ) - sc.front()).norm1();
      statE[ii].addValue( v );
    }
    /////////////

    ++ii;
    ++itc;
  } while( itc != ite );
}

template <typename KSpace, typename Iterator>
void analyseLengthMS( /*Statistic<double> & statD,*/ Statistic<double> & statE,
                      Iterator itb, Iterator ite )
{
  typedef typename KSpace::Space Space;
  typedef typename Space::Point Point;
  typedef typename Space::Vector Vector;
  typedef ArithmeticalDSSComputer< Iterator, int, 4 > SegmentComputer;
  typedef SaturatedSegmentation< SegmentComputer > Decomposition;
  typedef typename Decomposition::SegmentComputerIterator SegmentComputerIterator;
  // Computes the tangential cover
  SegmentComputer algo;
  Decomposition theDecomposition( itb, ite, algo);
  //statD.clear();
  statE.clear();
  for ( SegmentComputerIterator scIt = theDecomposition.begin(), scItEnd = theDecomposition.end();
        scIt != scItEnd; ++scIt )
  {
    const SegmentComputer & sc = *scIt;
    /*int64_t l = 0;
      for ( Iterator ptIt = sc.begin(), ptItEnd = sc.end(); ptIt != ptItEnd; ++ptIt )
        ++l;
      statD.addValue( (double) l );*/
    double v = (sc.back( ) - sc.front()).norm1();
    statE.addValue( v );
  }
}

void suggestedRadiusForIntegralInvariantEstimators( const std::vector< double > & radius,
                                                    std::vector< Dimension > & registration,
                                                    std::vector< double > & chosenRadius,
                                                    const Dimension nbRadius )
{
  KMeans<double, Euclidean>( radius, nbRadius, Euclidean(), registration, chosenRadius);
}


///////////////////////////////////////////////////////////////////////////////

int main( int argc, char** argv )
{
    trace.beginBlock ( "Example IntegralInvariantCurvature2D" );
    trace.info() << "Args:";
    for ( int i = 0; i < argc; ++i )
        trace.info() << " " << argv[ i ];
    trace.info() << endl;

    /// Construction of the shape + digitalization
    double h = 0.1;

    unsigned int nbKernels = 5;

    typedef Flower2D< Z2i::Space > MyShape;
    typedef GaussDigitizer< Z2i::Space, MyShape > MyGaussDigitizer;
    typedef Z2i::KSpace::Surfel Surfel;
    typedef Z2i::KSpace::SCell SCell;
    typedef LightImplicitDigitalSurface< Z2i::KSpace, MyGaussDigitizer > LightImplicitDigSurface;
    typedef DigitalSurface< LightImplicitDigSurface > MyDigitalSurface;

    MyShape shape( 0, 0, 20, 7, 6, 0.2 );

    MyGaussDigitizer digShape;
    digShape.attach( shape );
    digShape.init( shape.getLowerBound(), shape.getUpperBound(), h );
    Z2i::Domain domainShape = digShape.getDomain();
    Z2i::KSpace KSpaceShape;
    bool space_ok = KSpaceShape.init( domainShape.lowerBound(), domainShape.upperBound(), true );
    if ( !space_ok )
    {
        trace.error() << "Error in the Khamisky space construction." << std::endl;
        return 2;
    }

    typedef ImageSelector< Z2i::Domain, unsigned int >::Type Image;
    Image image( domainShape );
    DGtal::imageFromRangeAndValue( domainShape.begin(), domainShape.end(), image );

    SurfelAdjacency<Z2i::KSpace::dimension> SAdj( true );
    Surfel bel = Surfaces<Z2i::KSpace>::findABel( KSpaceShape, digShape, 100000 );

    ///////////////////////
    std::vector< SCell > points;
    Surfaces< Z2i::KSpace >::track2DBoundary( points, KSpaceShape, SAdj, digShape, bel );
    GridCurve< Z2i::KSpace > gridcurve;
    gridcurve.initFromSCellsVector( points );
    typedef typename GridCurve< Z2i::KSpace >::PointsRange PointsRange;
    PointsRange pointsRange2 = gridcurve.getPointsRange();

    const Dimension pr2size = (pointsRange2.size());
    std::vector< Statistic< double > > v_statMSEL(pr2size);
    for(Dimension ii = 0; ii < pr2size; ++ii )
    {
      v_statMSEL[ii] = Statistic<double>(true);
    }

    trace.beginBlock("Analyse segments and Mapping segments <-> Surfels...");

    analyseAllLengthMS<Z2i::KSpace>( v_statMSEL, pointsRange2.c(), pointsRange2.c() );

    trace.endBlock();
    ///////////////////////


    LightImplicitDigSurface LightImplDigSurf( KSpaceShape, digShape, SAdj, bel );
    MyDigitalSurface digSurf( LightImplDigSurf );

    typedef DepthFirstVisitor< MyDigitalSurface > Visitor;
    typedef GraphVisitorRange< Visitor > VisitorRange;
    typedef VisitorRange::ConstIterator SurfelConstIterator;

    VisitorRange range( new Visitor( digSurf, *digSurf.begin() ) );
    SurfelConstIterator abegin = range.begin();
    SurfelConstIterator aend = range.end();

    std::vector< SCell > contour;
    for( ; abegin != aend; ++abegin )
    {
      contour.push_back( *abegin );
    }

    std::vector< double > v_curvatures( contour.size() );
    std::vector< double > v_estimated_radius( contour.size() );

    trace.beginBlock("Computation of radius...");
    {
      std::map< double, unsigned int > nbKernelRadius;
      for( Dimension ii = 0; ii < pr2size; ++ii )
      {
        Dimension current_pos = pr2size - 1 - ii;
        double mean = v_statMSEL[ current_pos ].mean();
        double re = (0.1 * (mean * mean)) * h;
//        checkSizeRadius( re, h, minRadiusAABB );
        v_estimated_radius[ii] = re;

        if( nbKernels != 0 )
        {
          if( nbKernelRadius.find(re) == nbKernelRadius.end() )
          {
            nbKernelRadius[ re ] = 1;
          }
          else
          {
            nbKernelRadius[ re ] += 1;
          }
        }
      }
      if( nbKernelRadius.size() < nbKernels )
      {
        nbKernels = nbKernelRadius.size();
      }
    }
    trace.endBlock();

    trace.beginBlock("Sorting radius & pre-computing estimators...");
    std::vector< double > v_radius;
    std::vector< Dimension > v_registration;

    if( nbKernels > 0 )
    {
      suggestedRadiusForIntegralInvariantEstimators( v_estimated_radius, v_registration, v_radius, nbKernels );
    }
    trace.endBlock();

    /// Drawing results
    typedef double Quantity;
    Quantity min = numeric_limits < Quantity >::max();
    Quantity max = numeric_limits < Quantity >::min();
    for ( unsigned int i = 0; i < v_estimated_radius.size(); ++i )
    {
      double current = v_registration[i];
        if ( current < min )
        {
            min = current;
        }
        else if ( current > max )
        {
            max = current;
        }
    }
    Board2D board;
    VisitorRange range2( new Visitor( digSurf, *digSurf.begin() ) );
    abegin = range2.begin();

//    typedef  RandomColorMap Gradient;
//    Gradient cmap_grad( 0, 20 );
//    cmap_grad.addColor( Color::Blue );
//    cmap_grad.addColor( Color::Green );
//    cmap_grad.addColor( Color::Black );
//    cmap_grad.addColor( Color::Lime );
//    cmap_grad.addColor( Color::Purple );

//    typedef GradientColorMap< Quantity > Gradient;
//    Gradient cmap_grad( min, max, 10 );
//    cmap_grad.addColor( Color( 50, 50, 255 ) );
//    cmap_grad.addColor( Color( 255, 0, 0 ) );
//    cmap_grad.addColor( Color( 255, 255, 10 ) );

    typedef HueShadeColorMap<double> Gradient;
    Gradient cmap_grad( min, max, 1 );

    board << SetMode( (*abegin).className(), "Paving" );
    string specificStyle = (*abegin).className() + "/Paving";
    for ( unsigned int i = 0; i < v_estimated_radius.size(); ++i )
    {
        SCell currentCell = KSpaceShape.sIndirectIncident( *abegin, *KSpaceShape.sOrthDirs( *abegin ) ); // We apply the color to the inner spel (more visible than surfel)
        board << CustomStyle( specificStyle, new CustomColors( Color::Black, cmap_grad( v_registration[i] )))
              << currentCell;
        ++abegin;
    }
    board.saveSVG ( "example-integralinvariant2D.svg" );
    trace.endBlock();
    return 0;
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
