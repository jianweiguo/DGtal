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
 * @file EuclideanShapesDecorator.h
 * @author Jeremy Levallois (\c jeremy.levallois@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), INSA-Lyon, France
 * LAboratoire de MAthématiques - LAMA (CNRS, UMR 5127), Université de Savoie, France
 *
 * @date 2012/08/28
 *
 * This file is part of the DGtal library.
 */

#if defined(EuclideanShapesDecorator_RECURSES)
#error Recursive header files inclusion detected in EuclideanShapesDecorator.h
#else // defined(EuclideanShapesDecorator_RECURSES)
/** Prevents recursive inclusion of headers. */
#define EuclideanShapesDecorator_RECURSES

#if !defined EuclideanShapesDecorator_h
/** Prevents repeated inclusion of headers. */
#define EuclideanShapesDecorator_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/base/ConstAlias.h"

#include "DGtal/shapes/CEuclideanBoundedShape.h"
#include "DGtal/shapes/CEuclideanOrientedShape.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  template <typename ShapeA, typename ShapeB>
  class EuclideanShapesCSG
  {
  protected:
    enum e_operator
    {
      e_union,
      e_intersection,
      e_minus
    };

  public:
    BOOST_CONCEPT_ASSERT (( concepts::CEuclideanBoundedShape< ShapeA > ));
    BOOST_CONCEPT_ASSERT (( concepts::CEuclideanOrientedShape< ShapeA > ));

    typedef typename ShapeA::Space Space;
    typedef typename ShapeA::RealPoint RealPoint;

    EuclideanShapesCSG( )
    {
    }

    EuclideanShapesCSG ( const EuclideanShapesCSG & other )
    {
      myShapeA = other.myShapeA;
      v_shapes = other.v_shapes;

      myLowerBound = other.myLowerBound;
      myUpperBound = other.myUpperBound;

    }

    EuclideanShapesCSG & operator= ( const EuclideanShapesCSG & other )
    {
      myShapeA = other.myShapeA;
      v_shapes = other.v_shapes;

      myLowerBound = other.myLowerBound;
      myUpperBound = other.myUpperBound;
      return *this;
    }

    EuclideanShapesCSG( ShapeA* a )
      : myShapeA( a )
    {
      myLowerBound = myShapeA->getLowerBound();
      myUpperBound = myShapeA->getUpperBound();
    }

    void op_union( ShapeB* b )
    {
      BOOST_CONCEPT_ASSERT (( concepts::CEuclideanBoundedShape< ShapeB > ));
      BOOST_CONCEPT_ASSERT (( concepts::CEuclideanOrientedShape< ShapeB > ));
      std::pair<e_operator, ShapeB*> shape( e_union, b );

      for(uint i =0; i < Space::dimension; ++i)
      {
        myLowerBound[i] = std::min(myLowerBound[i], b->getLowerBound()[i]);
        myUpperBound[i] = std::max(myUpperBound[i], b->getUpperBound()[i]);
      }

      v_shapes.push_back(shape); 
    }

    void op_intersection( ShapeB* b )
    {
      BOOST_CONCEPT_ASSERT (( concepts::CEuclideanBoundedShape< ShapeB > ));
      BOOST_CONCEPT_ASSERT (( concepts::CEuclideanOrientedShape< ShapeB > ));
      std::pair<e_operator, ShapeB*> shape( e_intersection, b );

      for(uint i =0; i < Space::dimension; ++i)
      {
        myLowerBound[i] = std::max(myLowerBound[i], b->getLowerBound()[i]);
        myUpperBound[i] = std::min(myUpperBound[i], b->getUpperBound()[i]);
      }

      v_shapes.push_back(shape); 
    }

    void op_minus( ShapeB* b )
    {
      BOOST_CONCEPT_ASSERT (( concepts::CEuclideanBoundedShape< ShapeB > ));
      BOOST_CONCEPT_ASSERT (( concepts::CEuclideanOrientedShape< ShapeB > ));
      std::pair<e_operator, ShapeB*> shape( e_minus, b );
      v_shapes.push_back(shape); 

    }

    RealPoint getLowerBound() const
    {
      return myLowerBound;
    }

    RealPoint getUpperBound() const
    {
      return myUpperBound;
    }

    Orientation orientation( const RealPoint & p ) const
    {
      Orientation orient = myShapeA->orientation( p );

      for(unsigned int i = 0; i < v_shapes.size(); ++i)
      {
        if( v_shapes[i].first == e_minus )
        {
          if (( v_shapes[i].second->orientation( p ) == INSIDE ) || ( v_shapes[i].second->orientation( p ) == ON ))
          {
            orient = OUTSIDE;
          }
        }
        else if( v_shapes[i].first == e_intersection )
        {
          if (( orient == ON ) && ( v_shapes[i].second->orientation( p ) != OUTSIDE ))
          {
            orient = ON;
          }
          else if (( v_shapes[i].second->orientation( p ) == ON ) && ( orient != OUTSIDE ))
          {
            orient = ON;
          }
          else if (( orient == INSIDE ) && ( v_shapes[i].second->orientation( p ) == INSIDE ))
          {
            orient = INSIDE;
          }

          orient = OUTSIDE;
        }
        else /// e_union
        {
          if (( orient == OUTSIDE ) && ( v_shapes[i].second->orientation( p ) == OUTSIDE ))
          {
            orient = OUTSIDE;
          }
          else if(( orient == OUTSIDE ) && ( v_shapes[i].second->orientation( p ) != OUTSIDE ))
          {
            orient = v_shapes[i].second->orientation( p );
          }
          else if(( orient != OUTSIDE ) && ( v_shapes[i].second->orientation( p ) == OUTSIDE ))
          {
            orient = orient;
          }
          else if (( orient == INSIDE ) || ( v_shapes[i].second->orientation( p ) == INSIDE ))
          {
              orient = INSIDE;
          }
          else if (( orient == ON ) || ( v_shapes[i].second->orientation( p ) == ON ))
          {
              orient = ON;
          }
        }
      }

      return orient;
    }

  public:

    /**
     * Writes/Displays the object on an output stream.
     * @param out the output stream where the object is written.
     */
    void selfDisplay ( std::ostream & out ) const;

    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    bool isValid() const;

    // ------------------------- Hidden services ------------------------------
  protected:

    /**
     * Constructor.
     * Forbidden by default (protected to avoid g++ warnings).
     */
    // EuclideanShapesCSG();

  private:

    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    // EuclideanShapesCSG ( const EuclideanShapesCSG & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    // EuclideanShapesCSG & operator= ( const EuclideanShapesCSG & other );

    // ------------------------- Internals ------------------------------------
  private:
    ShapeA * myShapeA;
    std::vector< std::pair<e_operator, ShapeB*> > v_shapes;

    RealPoint myLowerBound;
    RealPoint myUpperBound;

  };

 //namespace deprecated
// {
/////////////////////////////////////////////////////////////////////////////
// template class EuclideanShapesDecorator
/**
 * Description of template class 'EuclideanShapesDecorator' <p>
 * \brief Aim: Union between two models of CEuclideanBoundedShape and CEuclideanOrientedShape
 *
 * @tparam ShapeA type of the first shape. Must be a model of CEuclideanBoundedShape and CEuclideanOrientedShape
 * @tparam ShapeB type of the second shape. Must be a model of CEuclideanBoundedShape and CEuclideanOrientedShape
 */

  template <typename ShapeA, typename ShapeB>
  class EuclideanShapesUnion
  {
    // ----------------------- Standard services ------------------------------
  public:
    BOOST_CONCEPT_ASSERT (( concepts::CEuclideanBoundedShape< ShapeA > ));
    BOOST_CONCEPT_ASSERT (( concepts::CEuclideanOrientedShape< ShapeA > ));
    BOOST_CONCEPT_ASSERT (( concepts::CEuclideanBoundedShape< ShapeB > ));
    BOOST_CONCEPT_ASSERT (( concepts::CEuclideanOrientedShape< ShapeB > ));

    typedef typename ShapeA::Space Space;
    typedef typename ShapeA::RealPoint RealPoint;

    /**
      * Constructor.
      *
      * @param[in] a a model of CEuclideanBoundedShape and CEuclideanOrientedShape
      * @param[in] b a model of CEuclideanBoundedShape and CEuclideanOrientedShape
      */
    EuclideanShapesUnion( ConstAlias< ShapeA > a, ConstAlias< ShapeB > b )
      : myShapeA( a ),
        myShapeB( b )
    {
      RealPoint shapeALowerBoundary = myShapeA.getLowerBound();
      RealPoint shapeBLowerBoundary = myShapeB.getLowerBound();
      RealPoint shapeAUpperBoundary = myShapeA.getUpperBound();
      RealPoint shapeBUpperBoundary = myShapeB.getUpperBound();
      for ( unsigned int i = 0; i < myLowerBound.size(); ++i )
      {
        myLowerBound[ i ] = std::min( shapeALowerBoundary[ i ], shapeBLowerBoundary[ i ] );
        myUpperBound[ i ] = std::max( shapeAUpperBoundary[ i ], shapeBUpperBoundary[ i ] );
      }
    }

    /**
     * @return the lower bound of the shape bounding box.
     *
     */
    RealPoint getLowerBound() const
    {
      return myLowerBound;
    }

    /**
     * @return the upper bound of the shape bounding box.
     *
     */
    RealPoint getUpperBound() const
    {
      return myUpperBound;
    }

    /**
     * Return the orientation of a point with respect to a shape.
     *
     * @param[in] p input point
     *
     * @return the orientation of the point (0 = INSIDE, 1 = ON, 2 = OUTSIDE)
     */
    Orientation orientation( const RealPoint & p ) const
    {
      if (( myShapeA.orientation( p ) == INSIDE ) || ( myShapeB.orientation( p ) == INSIDE ))
        {
            return INSIDE;
        }
      else if (( myShapeA.orientation( p ) == ON ) || ( myShapeB.orientation( p ) == ON ))
        {
            return ON;
        }
        return OUTSIDE;
    }

    /**
     * Destructor.
     */
    ~EuclideanShapesUnion(){}

    // ----------------------- Interface --------------------------------------
  public:

    /**
     * Writes/Displays the object on an output stream.
     * @param out the output stream where the object is written.
     */
    void selfDisplay ( std::ostream & out ) const;

    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    bool isValid() const;

    // ------------------------- Hidden services ------------------------------
  protected:

    /**
     * Constructor.
     * Forbidden by default (protected to avoid g++ warnings).
     */
    EuclideanShapesUnion();

  private:

    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    EuclideanShapesUnion ( const EuclideanShapesUnion & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    EuclideanShapesUnion & operator= ( const EuclideanShapesUnion & other );

    // ------------------------- Internals ------------------------------------
  private:
    const ShapeA & myShapeA;
    const ShapeB & myShapeB;

    RealPoint myLowerBound;
    RealPoint myUpperBound;

  }; // end of class EuclideanShapesUnion

  /////////////////////////////////////////////////////////////////////////////
  // template class EuclideanShapesIntersection
  /**
   * Description of template class 'EuclideanShapesIntersection' <p>
   * \brief Aim: Intersection between two models of CEuclideanBoundedShape and CEuclideanOrientedShape
   *
   * @tparam ShapeA type of the first shape. Must be a model of CEuclideanBoundedShape and CEuclideanOrientedShape
   * @tparam ShapeB type of the second shape. Must be a model of CEuclideanBoundedShape and CEuclideanOrientedShape
   */
  template <typename ShapeA, typename ShapeB>
  class EuclideanShapesIntersection
  {
    // ----------------------- Standard services ------------------------------
  public:
    BOOST_CONCEPT_ASSERT (( concepts::CEuclideanBoundedShape< ShapeA > ));
    BOOST_CONCEPT_ASSERT (( concepts::CEuclideanOrientedShape< ShapeA > ));
    BOOST_CONCEPT_ASSERT (( concepts::CEuclideanBoundedShape< ShapeB > ));
    BOOST_CONCEPT_ASSERT (( concepts::CEuclideanOrientedShape< ShapeB > ));

    typedef typename ShapeA::Space Space;
    typedef typename ShapeA::RealPoint RealPoint;

    /**
      * Constructor.
      *
      * @param[in] a a model of CEuclideanBoundedShape and CEuclideanOrientedShape
      * @param[in] b a model of CEuclideanBoundedShape and CEuclideanOrientedShape
      */
    EuclideanShapesIntersection( ConstAlias< ShapeA > a, ConstAlias< ShapeB > b )
      : myShapeA( a ),
        myShapeB( b )
    {
      RealPoint shapeALowerBoundary = myShapeA.getLowerBound();
      RealPoint shapeBLowerBoundary = myShapeB.getLowerBound();
      RealPoint shapeAUpperBoundary = myShapeA.getUpperBound();
      RealPoint shapeBUpperBoundary = myShapeB.getUpperBound();
      for ( unsigned int i = 0; i < myLowerBound.size(); ++i )
      {
        myLowerBound[ i ] = std::min( shapeALowerBoundary[ i ], shapeBLowerBoundary[ i ] );
        myUpperBound[ i ] = std::max( shapeAUpperBoundary[ i ], shapeBUpperBoundary[ i ] );
      }
    }


    /**
     * @return the lower bound of the shape bounding box.
     *
     */
    RealPoint getLowerBound() const
    {
      return myLowerBound;
    }

    /**
     * @return the upper bound of the shape bounding box.
     *
     */
    RealPoint getUpperBound() const
    {
      return myUpperBound;
    }

    /**
     * Return the orientation of a point with respect to a shape.
     *
     * @param[in] p input point
     *
     * @return the orientation of the point (0 = INSIDE, 1 = ON, 2 = OUTSIDE)
     */
    Orientation orientation( const RealPoint & p ) const
    {
      if (( myShapeA.orientation( p ) == ON ) && ( myShapeB.orientation( p ) != OUTSIDE ))
      {
        return ON;
      }
      else if (( myShapeB.orientation( p ) == ON ) && ( myShapeA.orientation( p ) != OUTSIDE ))
      {
        return ON;
      }
      else if (( myShapeA.orientation( p ) == INSIDE ) && ( myShapeB.orientation( p ) == INSIDE ))
      {
        return INSIDE;
      }

      return OUTSIDE;
    }


    /**
     * Destructor.
     */
    ~EuclideanShapesIntersection(){}

    // ----------------------- Interface --------------------------------------
  public:

    /**
     * Writes/Displays the object on an output stream.
     * @param out the output stream where the object is written.
     */
    void selfDisplay ( std::ostream & out ) const;

    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    bool isValid() const;

    // ------------------------- Hidden services ------------------------------
  protected:

    /**
     * Constructor.
     * Forbidden by default (protected to avoid g++ warnings).
     */
    EuclideanShapesIntersection();

  private:

    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    EuclideanShapesIntersection ( const EuclideanShapesIntersection & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    EuclideanShapesIntersection & operator= ( const EuclideanShapesIntersection & other );

    // ------------------------- Internals ------------------------------------
  private:
    const ShapeA & myShapeA;
    const ShapeB & myShapeB;

    RealPoint myLowerBound;
    RealPoint myUpperBound;

  }; // end of class EuclideanShapesIntersection

  /////////////////////////////////////////////////////////////////////////////
  // template class EuclideanShapesMinus
  /**
   * Description of template class 'EuclideanShapesMinus' <p>
   * \brief Aim: Minus between two models of CEuclideanBoundedShape and CEuclideanOrientedShape
   *
   * @tparam ShapeA type of the first shape. Must be a model of CEuclideanBoundedShape and CEuclideanOrientedShape
   * @tparam ShapeB type of the second shape. Must be a model of CEuclideanBoundedShape and CEuclideanOrientedShape
   */
  template <typename ShapeA, typename ShapeB>
  class EuclideanShapesMinus
  {
    // ----------------------- Standard services ------------------------------
  public:
    BOOST_CONCEPT_ASSERT (( concepts::CEuclideanBoundedShape< ShapeA > ));
    BOOST_CONCEPT_ASSERT (( concepts::CEuclideanOrientedShape< ShapeA > ));
    BOOST_CONCEPT_ASSERT (( concepts::CEuclideanBoundedShape< ShapeB > ));
    BOOST_CONCEPT_ASSERT (( concepts::CEuclideanOrientedShape< ShapeB > ));

    typedef typename ShapeA::Space Space;
    typedef typename ShapeA::RealPoint RealPoint;

    /**
      * Constructor.
      *
      * @param[in] a a model of CEuclideanBoundedShape and CEuclideanOrientedShape
      * @param[in] b a model of CEuclideanBoundedShape and CEuclideanOrientedShape
      */
    EuclideanShapesMinus( ConstAlias< ShapeA > a, ConstAlias< ShapeB > b )
      : myShapeA( a ),
        myShapeB( b )
    {
      RealPoint shapeALowerBoundary = myShapeA.getLowerBound();
      RealPoint shapeBLowerBoundary = myShapeB.getLowerBound();
      RealPoint shapeAUpperBoundary = myShapeA.getUpperBound();
      RealPoint shapeBUpperBoundary = myShapeB.getUpperBound();
      for ( unsigned int i = 0; i < myLowerBound.size(); ++i )
      {
        myLowerBound[ i ] = std::min( shapeALowerBoundary[ i ], shapeBLowerBoundary[ i ] );
        myUpperBound[ i ] = std::max( shapeAUpperBoundary[ i ], shapeBUpperBoundary[ i ] );
      }
    }

    /**
     * @return the lower bound of the shape bounding box.
     *
     */
    RealPoint getLowerBound() const
    {
      return myLowerBound;
    }

    /**
     * @return the upper bound of the shape bounding box.
     *
     */
    RealPoint getUpperBound() const
    {
      return myUpperBound;
    }

    /**
     * Return the orientation of a point with respect to a shape.
     *
     * @param[in] p input point
     *
     * @return the orientation of the point (0 = INSIDE, 1 = ON, 2 = OUTSIDE)
     */
    Orientation orientation( const RealPoint & p ) const
    {
      if (( myShapeB.orientation( p ) == INSIDE ) || ( myShapeB.orientation( p ) == ON ))
      {
        return OUTSIDE;
      }
      return myShapeA.orientation( p );
    }


    /**
     * Destructor.
     */
    ~EuclideanShapesMinus(){}

    // ----------------------- Interface --------------------------------------
  public:

    /**
     * Writes/Displays the object on an output stream.
     * @param out the output stream where the object is written.
     */
    void selfDisplay ( std::ostream & out ) const;

    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    bool isValid() const;

    // ------------------------- Hidden services ------------------------------
  protected:

    /**
     * Constructor.
     * Forbidden by default (protected to avoid g++ warnings).
     */
    EuclideanShapesMinus();

  private:

    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    EuclideanShapesMinus ( const EuclideanShapesMinus & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    EuclideanShapesMinus & operator= ( const EuclideanShapesMinus & other );

    // ------------------------- Internals ------------------------------------
  private:
    const ShapeA & myShapeA;
    const ShapeB & myShapeB;

    RealPoint myLowerBound;
    RealPoint myUpperBound;

  }; // end of class EuclideanShapesMinus


 //}


  /**
   * Overloads 'operator<<' for displaying objects of class 'EuclideanShapesDecorator'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'EuclideanShapesDecorator' to write.
   * @return the output stream after the writing.
   */
  template <typename ShapeA, typename ShapeB>
  std::ostream&
  operator<< ( std::ostream & out, const /*deprecated::*/EuclideanShapesUnion<ShapeA, ShapeB> & object );

  template <typename ShapeA, typename ShapeB>
  std::ostream&
  operator<< ( std::ostream & out, const /*deprecated::*/EuclideanShapesIntersection<ShapeA, ShapeB> & object );

  template <typename ShapeA, typename ShapeB>
  std::ostream&
  operator<< ( std::ostream & out, const /*deprecated::*/EuclideanShapesMinus<ShapeA, ShapeB> & object );

} // namespace DGtal


//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined EuclideanShapesDecorator_h

#undef EuclideanShapesDecorator_RECURSES
#endif // else defined(EuclideanShapesDecorator_RECURSES)
