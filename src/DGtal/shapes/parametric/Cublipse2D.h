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
 * @file Cublipse2D.h
 * @author Jérémy Levallois (\c jeremy.levallois@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), INSA-Lyon, France
 * LAboratoire de MAthématiques - LAMA (CNRS, UMR 5127), Université de Savoie, France
 *
 * @date 2013/12/05
 *
 * Header file for module Cublipse2D.cpp
 *
 * Half cube, half ellipse => Cublipse !
 *
 * This file is part of the DGtal library.
 */

#if defined(Cublipse2D_RECURSES)
#error Recursive header files inclusion detected in Cublipse2D.h
#else // defined(Cublipse2D_RECURSES)
/** Prevents recursive inclusion of headers. */
#define Cublipse2D_RECURSES

#if !defined Cublipse2D_h
/** Prevents repeated inclusion of headers. */
#define Cublipse2D_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/shapes/parametric/StarShaped2D.h"
#include <cmath>
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // template class Cublipse2D
  /**
   * Description of template class 'Cublipse2D' <p>
   * \brief Aim: Model of the concept StarShaped
   * represents any ellipse in the plane.
   *
   * NB: A backport from [ImaGene](https://gforge.liris.cnrs.fr/projects/imagene).
   */
  template <typename TSpace>
  class Cublipse2D:  public StarShaped2D<TSpace>
  {
    // ----------------------- Standard services ------------------------------
  public:

    typedef TSpace Space;
    typedef typename Space::Point Point;
    typedef typename Space::RealPoint RealPoint2D;
    typedef typename Space::RealVector RealVector2D;
      
    /**
     * Destructor.
     */
    ~Cublipse2D();
    
    /**
     * Constructor. 
     * @param x0 the x-coordinate of the circle center.
     * @param y0 the y-coordinate of the circle center.
     * @param a0 the half big axis of the ellipse.
     * @param a1 the half small axis of the ellipse.
     * @param theta the orientation of the ellipse.
     */
    Cublipse2D( const double x0, const double y0,
         const double a0, const double a1, const double theta);

    /**
     * Constructor. 
     * @param aPoint the circle center.
     * @param a0 the half big axis of the ellipse.
     * @param a1 the half small axis of the ellipse.
     * @param theta the orientation of the ellipse.
     */
    Cublipse2D(const RealPoint2D &aPoint,
        const double a0, const double a1, const double theta);

    /**
     * Constructor. 
     * @param aPoint the circle center.
     * @param a0 the half big axis of the ellipse.
     * @param a1 the half small axis of the ellipse.
     * @param theta the orientation of the ellipse.
     */
    Cublipse2D(const Point &aPoint,
        const double a0, const double a1, const double theta);

    
  // ------------- Implementation of 'StarShaped' services ------------------
  public:

    RealPoint interiorPoint() const
    {
      return center();
    }

    /**
     * @return the lower bound of the shape bounding box.
     *
     */
    RealPoint2D getLowerBound() const
    {
      return RealPoint2D(myCenter[0] - myAxis1, myCenter[1] - myAxis1);
    }

    /**
     * @return the upper bound of the shape bounding box.
     *
     */
    RealPoint2D getUpperBound() const
    {
      return RealPoint2D(myCenter[0] + myAxis1, myCenter[1] + myAxis1);
    }

    /**
     * @return the center of the star-shaped object.
     */
    RealPoint2D center() const
    {
      return myCenter;
    }
   
    /**
     * @param p any point in the plane.
     *
     * @return the angle parameter between 0 and 2*Pi corresponding to
     * this point for the shape.
     */
    double parameter( const RealPoint2D & p ) const;


    /**
     * @param t any angle between 0 and 2*Pi.
     *
     * @return the vector (x(t),y(t)) which is the position on the
     * shape boundary.
     */
    RealPoint2D x( const double t ) const;

    /**
     * @param t any angle between 0 and 2*Pi.
     *
     * @return the vector (x'(t),y'(t)) which is the tangent to the
     * shape boundary.
     */
    RealVector2D xp( const double t ) const;

    /**
     * @param t any angle between 0 and 2*Pi.
     *
     * @return the vector (x''(t),y''(t)).
     */
    RealVector2D xpp( const double t ) const;
    

    // ------------------------- data ----------------------------
  private:

    /**
     * Center of the circle.
     */
    RealPoint2D myCenter;
    
    /**
     * First axis.
     */
    double myAxis1;


    /**
     * Second axis.
     */
    double myAxis2;
    
    /**
     * Orientation (radian).
     */
    double myTheta;

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
    Cublipse2D();

  private:

    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    //  Cublipse2D ( const Cublipse2D & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    Cublipse2D & operator= ( const Cublipse2D & other );

    // ------------------------- Internals ------------------------------------
  private:

  }; // end of class Cublipse2D


  /**
   * Overloads 'operator<<' for displaying objects of class 'Cublipse2D'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'Cublipse2D' to write.
   * @return the output stream after the writing.
   */
  template <typename T>
  std::ostream&
  operator<< ( std::ostream & out, const Cublipse2D<T> & object );

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "DGtal/shapes/parametric/Cublipse2D.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined Cublipse2D_h

#undef Cublipse2D_RECURSES
#endif // else defined(Cublipse2D_RECURSES)

