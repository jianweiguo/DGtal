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
 * @file testKMeans.cpp
 * @ingroup Tests
 * @author Jeremy Levallois (\c jeremy.levallois@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), INSA-Lyon, France
 * LAboratoire de MAthématiques - LAMA (CNRS, UMR 5127), Université de Savoie, France
 *
 * @date 2014/03/03
 *
 * Functions for testing class KMeans.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <vector>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/math/KMeans.h"
#include "DGtal/io/boards/Board2D.h"

///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;

///////////////////////////////////////////////////////////////////////////////
// Functions for testing class KMeans.
///////////////////////////////////////////////////////////////////////////////

template< typename TPoint >
struct EuclideanDistancePointFunctor
{
  typedef TPoint Point;

  double operator()(const Point &a, const Point &b) const
  {
    return (double)(b-a).norm();
  }
};

struct MyDrawStyleCustomGreen : public DrawableWithBoard2D
{
  virtual void setStyle( DGtal::Board2D & aBoard ) const
   {
     aBoard.setFillColorRGBi(0,160,0);
     aBoard.setPenColorRGBi(80,0,0);
   }
};

bool testKMeans()
{
  unsigned int nbok = 0;
  unsigned int nb = 0;
  
  typedef Z2i::Point Point2D;
  std::vector< Point2D > v_points;
  for( int ii = -10; ii <= 20; ++ii )
  {
    v_points.push_back( Point2D(ii, ii));
  }

  EuclideanDistancePointFunctor< Point2D > functor;

  std::vector< Dimension > v_registration;
  std::vector< Point2D > v_centroids;
  Dimension k = 3;
  KMeans::computeKMeans< Point2D, EuclideanDistancePointFunctor< Point2D > >( v_points, k, functor, v_registration, v_centroids );
  std::cout << v_centroids[0] << " " << v_centroids[1] << " " << v_centroids[2] << std::endl;

  /////////////
  Board2D board;
  Z2i::Domain domain( Point2D(-11,-11), Point2D(21,21));
  board << SetMode( domain.className(), "Grid" ) << domain;

  std::vector<Color> v_colorMap(k);
  for( Dimension ii = 0; ii < k; ++ii )
  {
    v_colorMap[ii] = Color( std::rand()%255, std::rand()%255, std::rand()%255 );
  }

  std::string specificStyle = v_points[0].className();
  for( Dimension ii = 0; ii < v_points.size(); ++ii )
  {
    board << CustomStyle( specificStyle, new CustomColors( v_colorMap[v_registration[ii]], v_colorMap[v_registration[ii]] ) )
      << v_points[ii];
  }
  board.saveSVG("testKMeans.svg");
  /////////////


  ++nb, nbok += true ? 1 : 0;
  return nbok == nb;
}


///////////////////////////////////////////////////////////////////////////////
// Standard services - public :

int main( int, char** )
{
  trace.beginBlock ( "Testing class Histogram" );

  bool res = testKMeans();
  trace.emphase() << ( res ? "Passed." : "Error." ) << endl;

  trace.endBlock();
  return res ? 0 : 1;
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
