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
 * @file SCellsFunctors.h
 *
 * @author Tristan Roussillon (\c tristan.roussillon@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 * @date 2012/02/02
 *
 * This files contains several basic classes representing Functors
 * on points.
 *
 * This file is part of the DGtal library.
 */

#if defined(SCellsFunctors_RECURSES)
#error Recursive header files inclusion detected in SCellsFunctors.h
#else // defined(SCellsFunctors_RECURSES)
/** Prevents recursive inclusion of headers. */
#define SCellsFunctors_RECURSES

#if !defined SCellsFunctors_h
/** Prevents repeated inclusion of headers. */
#define SCellsFunctors_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <iterator>
#include "DGtal/base/Common.h"
#include "DGtal/kernel/SpaceND.h"
#include "DGtal/base/BasicBoolFunctions.h"

#include "DGtal/topology/KhalimskySpaceND.h"
#include "DGtal/base/IteratorAdapter.h"
#include "DGtal/base/ConstIteratorAdapter.h"
#include "DGtal/kernel/BasicPointFunctors.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

/////////////////////////////////////////////////////////////////////////////
// template class SCellToPoint
/**
   * Description of template class 'SCellToPoint' <p>
   * \brief Aim: transforms a scell into a point
   * @tparam KSpace the Khalimsky space
   * @code
  KSpace aKSpace;
  typename KSpace::SCell aSCell;
  typename KSpace::Point aPoint;
  SCellToPoint<KSpace> f(aKSpace);
  ...
  aPoint = f(aSCell);
   * @endcode
   * @see ConstIteratorAdapter KhalimskySpaceND PointVector
   */
template <typename KSpace>
class SCellToPoint
{
public:
  typedef typename KSpace::Point Output;
  typedef typename KSpace::SCell Input;

private:
  /**
       * Aliasing pointer on the Khalimsky space.
      */
  const KSpace* myK;

public:

  /**
       * Default constructor.
      */
  SCellToPoint() : myK(NULL) { }
  /**
       * Constructor.
       * @param aK a Khalimsky space
      */
  SCellToPoint(const KSpace& aK) : myK(&aK) { }

  /**
     * Copy constructor.
     * @param other any SCellToPoint functor
     */
  SCellToPoint(const SCellToPoint& other)
    : myK(other.myK) { }

  /**
     * Assignment.
     *
     * @param other the object to copy.
     * @return a reference on 'this'.
     */
  SCellToPoint& operator= ( const SCellToPoint & other )
  {
    if (this != &other)
    {
      myK = other.myK;
    }
    return *this;
  }

  /**
     * Returns a point (with integer coordinates)
     * from a scell (with khalimsky coordinates)
     * @param aSCell a scell
     * @return the corresponding point.
     */
  Output operator()(const Input& aSCell) const
  {
    ASSERT( myK );
    Input s = aSCell;
    while ( myK->sDim(s) > 0 )
    {
      Input tmp( myK->sIndirectIncident( s, *myK->sDirs( s ) ) );
      ASSERT( myK->sDim(tmp) < myK->sDim(s) );
      s = tmp;
    }
    return Output( myK->sCoords(s) );
  }

}; // end of class SCellToPoint

/**
   * SCellToMidPoint is now deprecated. Please use CanonicSCellEmbedder instead.
   */
namespace deprecated
{
/////////////////////////////////////////////////////////////////////////////
// template class SCellToMidPoint
/**
   * Description of template class 'SCellToMidPoint' <p>
   * \brief Aim: transforms a scell into a real point
   * (the coordinates are divided by 2)
   * @tparam KSpace the Khalimsky space
   * @code
  KSpace aKSpace;
  typename KSpace::SCell aSCell;
  typename KSpace::Space::RealPoint aPoint;
  SCellToMidPoint<KSpace> f(aKSpace);
  ...
  aPoint = f(aSCell);
   * @endcode
   * @see ConstIteratorAdapter KhalimskySpaceND PointVector
   */
template <typename KSpace>
class SCellToMidPoint
{

public:

  typedef typename KSpace::Space::RealPoint Output;
  typedef typename KSpace::SCell Input;

private:
  /**
       *  Aliasing pointer on the Khalimsky space.
      */
  const KSpace* myK;

public:

  /**
       * Default constructor.
      */
  SCellToMidPoint() : myK(NULL) { }
  /**
       *  Constructor.
       * @param aK a Khalimsky space
      */
  SCellToMidPoint(const KSpace& aK) : myK(&aK) { }

  /**
     *  Copy constructor.
     * @param other any SCellToMidPoint functor
     */
  SCellToMidPoint(const SCellToMidPoint& other)
    : myK(other.myK) { }

  /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     */
  SCellToMidPoint & operator= ( const SCellToMidPoint & other )
  {
    if (this != &other)
    {
      myK = other.myK;
    }
    return *this;
  }

  /**
     * Return a real point (double coordinates) from a scell (khalimsky coordinates)
     * @param s a scell
     * @return the corresponding point.
     */
  Output operator()(const Input& s) const
  {
    ASSERT( myK );
    Output o( myK->sKCoords(s) );
    o /= 2;
    for( unsigned int i = 0; i < o.dimension; ++i )
      o[i] -= 0.5;
    return o;
  }

}; // end of class SCellToMidPoint
} // end of namespace deprecated
/////////////////////////////////////////////////////////////////////////////
// template class SCellToArrow
/**
   * Description of template class 'SCellToArrow' <p>
   * \brief Aim: transforms a signed cell into an arrow,
   * ie. a pair point-vector
   * @tparam KSpace the Khalimsky space
   * @see SCellToPoint ConstIteratorAdapter KhalimskySpaceND PointVector
   */
template <typename KSpace>
class SCellToArrow
{

public:

  typedef typename KSpace::Point Point;
  typedef typename KSpace::Vector Vector;
  typedef std::pair<Point,Vector> Output;
  typedef typename KSpace::SCell Input;

private:
  /**
       *  Aliasing pointer on the Khalimsky space.
      */
  const KSpace* myK;

public:

  /**
       * Default constructor.
      */
  SCellToArrow() : myK(NULL) { }
  /**
       *  Constructor.
       * @param aK a Khalimsky space
      */
  SCellToArrow(const KSpace& aK) : myK(&aK) { }

  /**
     *  Copy constructor.
     * @param other any SCellToArrow modifier
     */
  SCellToArrow(const SCellToArrow& other)
    : myK(other.myK) { }

  /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     */
  SCellToArrow & operator= ( const SCellToArrow & other )
  {
    if (this != &other)
    {
      myK = other.myK;
    }
    return *this;
  }

  /**
     * Get an arrow, ie a pair point-vector with integer coordinates
     * from a scell in khalimsky coordinates
     * @param s a scell
     * @return the corresponding point.
     */
  Output operator()(const Input& s) const
  {
    ASSERT( myK );
    //starting point of the arrow
    Input pointel( myK->sIndirectIncident( s, *myK->sDirs( s ) ) );
    Point p( myK->sCoords( pointel ) );   //integer coordinates
    //displacement vector
    Vector v( myK->sKCoords( s ) - myK->sKCoords( pointel ) );
    return Output(p,v);
  }

}; // end of class SCellToArrow

/////////////////////////////////////////////////////////////////////////////
// template class SCellToInnerPoint
/**
   * Description of template class 'SCellToInnerPoint' <p>
   * \brief Aim: transforms a signed cell c into a point
   * corresponding to the signed cell of greater dimension
   * that is indirectly incident to c.
   *
   * For instance, a linel is mapped into the indirect incident pixel center
   * and a surfel is mapped into the indirect incident voxel center.
   *
   * @tparam KSpace the Khalimsky space
   * @see SCellToPoint SCellToOuterPoint ConstIteratorAdapter KhalimskySpaceND PointVector
   */
template <typename KSpace>
class SCellToInnerPoint
{

public:

  typedef typename KSpace::Point Output;
  typedef typename KSpace::SCell Input;

private:
  /**
       *  Aliasing pointer on the Khalimsky space.
      */
  const KSpace* myK;

public:

  /**
       * Default constructor.
      */
  SCellToInnerPoint() : myK(NULL) { }
  /**
       *  Constructor.
       * @param aK a Khalimsky space
      */
  SCellToInnerPoint(const KSpace& aK) : myK(&aK) { }

  /**
     *  Copy constructor.
     * @param other any SCellToInnerPoint functor
     */
  SCellToInnerPoint(const SCellToInnerPoint& other)
    : myK(other.myK) { }

  /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     */
  SCellToInnerPoint & operator= ( const SCellToInnerPoint & other )
  {
    if (this != &other)
    {
      myK = other.myK;
    }
    return *this;
  }

  /**
     * Return a point (integer coordinates) from a scell (khalimsky coordinates)
     * @param s a linel
     * @return the inner pixel center
     */
  Output operator()(const Input& s) const
  {
    ASSERT( myK );
    Input pixel( myK->sDirectIncident( s, *myK->sOrthDirs( s ) ) );
    return Output( myK->sCoords( pixel ) ); //integer coordinates
  }

}; // end of class SCellToInnerPoint

/////////////////////////////////////////////////////////////////////////////
// template class SCellToOuterPoint
/**
   * Description of template class 'SCellToOuterPoint' <p>
   * \brief Aim: transforms a signed cell c into a point
   * corresponding to the signed cell of greater dimension
   * that is directly incident to c.
   *
   * For instance, a linel is mapped into the direct incident pixel center
   * and a surfel is mapped into the direct incident voxel center.
   *
   * @tparam KSpace the Khalimsky space
   * @see SCellToPoint SCellToInnerPoint ConstIteratorAdapter KhalimskySpaceND PointVector
   */
template <typename KSpace>
class SCellToOuterPoint
{
public:

  typedef typename KSpace::Point Output;
  typedef typename KSpace::SCell Input;

private:
  /**
       *  Aliasing pointer on the Khalimsky space.
      */
  const KSpace* myK;

public:

  /**
       * Default constructor.
      */
  SCellToOuterPoint() : myK(NULL) { }
  /**
       *  Constructor.
       * @param aK a Khalimsky space
      */
  SCellToOuterPoint(const KSpace& aK) : myK(&aK) { }

  /**
     *  Copy constructor.
     * @param other any SCellToOuterPoint modifier
     */
  SCellToOuterPoint(const SCellToOuterPoint& other)
    : myK(other.myK) { }

  /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     */
  SCellToOuterPoint & operator= ( const SCellToOuterPoint & other )
  {
    if (this != &other)
    {
      myK = other.myK;
    }
    return *this;
  }

  /**
     * Return a point (integer coordinates) from a scell (khalimsky coordinates)
     * @param s a linel
     * @return the outer pixel center
     */
  Output operator()(const Input& s) const
  {
    ASSERT( myK );
    Input pixel( myK->sIndirectIncident( s, *myK->sOrthDirs( s ) ) );
    return Output( myK->sCoords( pixel ) ); //integer coordinates
  }

}; // end of class SCellToOuterPoint

/////////////////////////////////////////////////////////////////////////////
// template class SCellToIncidentPoints
/**
   * Description of template class 'SCellToIncidentPoints' <p>
   * \brief Aim: transforms a signed cell c into a pair of points
   * corresponding to the signed cells of greater dimension
   * that are indirectly and directly incident to c.
   *
   * For instance, a linel is mapped into the pair of incident pixel centers
   * and a surfel is mapped into the pair of incident voxel centers.
   *
   * @tparam KSpace the Khalimsky space
   * @see SCellToInnerPoint SCellToOuterPoint ConstIteratorAdapter KhalimskySpaceND PointVector
   */
template <typename KSpace>
class SCellToIncidentPoints
{

public:

  typedef typename KSpace::Point Point;
  typedef std::pair<Point,Point> Output;
  typedef typename KSpace::SCell Input;

private:
  /**
       *  Aliasing pointer on the Khalimsky space.
      */
  const KSpace* myK;

public:

  /**
       * Default constructor.
      */
  SCellToIncidentPoints() : myK(NULL) { }
  /**
       *  Constructor.
       * @param aK a Khalimsky space
      */
  SCellToIncidentPoints(const KSpace& aK) : myK(&aK) { }

  /**
     *  Copy constructor.
     * @param other any SCellToIncidentPoints functor
     */
  SCellToIncidentPoints(const SCellToIncidentPoints& other)
    : myK(other.myK) { }

  /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     */
  SCellToIncidentPoints & operator= ( const SCellToIncidentPoints & other )
  {
    if (this != &other)
    {
      myK = other.myK;
    }
    return *this;
  }

  /**
     * Get a pair of point (integer coordinates) from a scell (khalimsky coordinates)
     * @param s a linel
     * @return the pair of points
     */
  Output operator()(const Input& s) const
  {
    ASSERT( myK );
    //inner point
    Input innerPixel( myK->sDirectIncident( s, *myK->sOrthDirs( s ) ) );
    //outer point
    Input outerPixel( myK->sIndirectIncident( s, *myK->sOrthDirs( s ) ) );

    return Output(myK->sCoords( innerPixel ),myK->sCoords( outerPixel ));
  }

}; // end of class SCellToIncidentPoints

/////////////////////////////////////////////////////////////////////////////
// template class SCellToCode
/**
   * Description of template class 'SCellToCode' <p>
   * \brief Aim: transforms a 2d signed cell, basically a linel,
    * into a code (0,1,2 or 3),
   * @tparam KSpace the 2d Khalimsky space
   * @see ConstIteratorAdapter KhalimskySpaceND
   */
template <typename KSpace>
class SCellToCode
{

  BOOST_STATIC_ASSERT( KSpace::dimension == 2 );

public:

  typedef typename KSpace::Point Point;
  typedef typename KSpace::Vector Vector;
  typedef char Output;

  typedef typename KSpace::SCell Input;

private:
  /**
       *  Aliasing pointer on the Khalimsky space.
      */
  const KSpace* myK;

public:

  /**
       * Default constructor.
      */
  SCellToCode() : myK(NULL) { }
  /**
       *  Constructor.
       * @param aK a Khalimsky space
      */
  SCellToCode(const KSpace& aK) : myK(&aK) { }

  /**
     *  Copy constructor.
     * @param other any SCellToCode modifier
     */
  SCellToCode(const SCellToCode& other)
    : myK(other.myK) { }

  /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     */
  SCellToCode & operator= ( const SCellToCode & other )
  {
    if (this != &other)
    {
      myK = other.myK;
    }
    return *this;
  }

  /**
     * Return a code from a linel
     * @param s a linel
     * @return the corresponding code
     */
  Output operator()(const Input& s) const
  {
    ASSERT( myK );
    //starting point of the arrow
    Input pointel( myK->sIndirectIncident( s, *myK->sDirs( s ) ) );
    Point p( myK->sCoords( pointel ) );   //integer coordinates
    //displacement vector
    Vector v( myK->sKCoords( s ) - myK->sKCoords( pointel ) );
    if (v == Vector(1,0)) return '0';
    else if (v == Vector(0,1)) return '1';
    else if (v == Vector(-1,0)) return '2';
    else if (v == Vector(0,-1)) return '3';
    else return 'e'; //e for error!
  }

}; // end of class SCellToCode

/////////////////////////////////////////////////////////////////////////////
// template class SCellProjector
template < typename K = KhalimskySpaceND< 2, DGtal::int32_t > >
struct SCellProjector
{
  typedef K KSpace;
  typedef typename KSpace::Space::Dimension Dimension;
  BOOST_STATIC_CONSTANT( Dimension, dimension = KSpace::dimension );
  typedef typename KSpace::Integer Integer;
  typedef typename KSpace::SCell SCell;

  /**
     * Default constructor
     */
  SCellProjector(const Integer & aDefaultInteger = NumberTraits<Integer>::zero())
    : myDims(), myDefaultInteger(aDefaultInteger)
  { //default projection
    Dimension k = 0;
    for ( ; k < dimension; ++k)
    {
      myDims[k] = k;
    }
  }

  /**
     * Initialization of the array of relevant dimensions
     * @param itb begin iterator on dimensions.
     * @param ite end iterator on dimensions.
     */
  template<typename TIterator>
  void init ( const TIterator & itb, const TIterator & ite )
  {
    BOOST_STATIC_ASSERT ((boost::is_same< Dimension,
                          typename std::iterator_traits<TIterator>::value_type >::value));

    TIterator it = itb;
    Dimension k = 0;
    for ( ; ( (k < dimension)&&(it != ite) ); ++it, ++k)
    {
      myDims[k] = *it;
    }
    for ( ; k < dimension; ++k)
    {
      myDims[k] = dimension;
    }
  }

  /**
     *  Initialisation by removing a given dimension.
     *  @param dimRemoved the removed dimension.
     */
  void initRemoveOneDim ( const Dimension & dimRemoved )
  {
    std::vector<Dimension> vectDims;
    Dimension aCurrentPos = 0;
    Dimension k=0;
    for ( ; k < dimension; ++k)
    {
      if(k!=dimRemoved){
        vectDims.push_back(aCurrentPos);
      }
      aCurrentPos++;
    }
    init(vectDims.begin(), vectDims.end());
  }

  /**
     *  Initialisation by adding a given dimension.
     *  @param newDim the new dimension.
     */
  void initAddOneDim ( const Dimension & newDim )
  {
    std::vector<Dimension> vectDims;
    Dimension maxIndex = dimension;
    Dimension aCurrentPos = 0;
    Dimension k=0;
    for ( ; k < dimension; ++k)
    {
      if(k==newDim){
        vectDims.push_back(maxIndex);
      }else{
        vectDims.push_back(aCurrentPos);
        aCurrentPos++;
      }
    }
    init(vectDims.begin(), vectDims.end());
  }

  /**
     *  Main operator
     * @param aPoint any point.
     * @return the projected point.
     */
  template<typename TInputSCell>
  SCell operator()( const TInputSCell & aSCell ) const
  {
    BOOST_STATIC_ASSERT ((boost::is_same< typename TInputSCell::Integer, Integer >::value));
    SCell res;

#ifdef CPP11_ARRAY
    typename std::array<Dimension,dimension>::const_iterator it = myDims.begin();
    typename std::array<Dimension,dimension>::const_iterator itEnd = myDims.end();
#else
    typename boost::array<Dimension,dimension>::const_iterator it = myDims.begin();
    typename boost::array<Dimension,dimension>::const_iterator itEnd = myDims.end();
#endif

    Dimension k = 0;
    for ( ; it != itEnd; ++it, ++k)
    {
      Dimension l = *it;
      if (l < TInputSCell::Point::dimension)
        res.myCoordinates[k] = aSCell.myCoordinates[l];
      else
        res.myCoordinates[k] = myDefaultInteger;
    }

    return res;
  }

private:
  /**
     * Array storing the coordinates that are copied from
     * the input point to its projection (order matters)
     */
#ifdef CPP11_ARRAY
  std::array<Dimension, dimension> myDims;
#else
  boost::array<Dimension, dimension> myDims;
#endif
  /**
     * Default integer set to coordinates of the projected point
     * not in the input point
     */
  Integer myDefaultInteger;

}; // end of class SCellProjector

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
//#include "DGtal/kernel/SCellsFunctors.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined SCellsFunctors_h

#undef SCellsFunctors_RECURSES
#endif // else defined(SCellsFunctors_RECURSES)
