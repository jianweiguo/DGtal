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
 * @file ArithmeticalDSS.h @author Tristan Roussillon (\c
 * tristan.roussillon@liris.cnrs.fr ) Laboratoire d'InfoRmatique en
 * Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS,
 * France
 *
 * @date 2010/07/01
 *
 * Header file for module ArithmeticalDSS.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(ArithmeticalDSS_RECURSES)
#error Recursive header files inclusion detected in ArithmeticalDSS.h
#else // defined(ArithmeticalDSS_RECURSES)
/** Prevents recursive inclusion of headers. */
#define ArithmeticalDSS_RECURSES

#if !defined ArithmeticalDSS_h
/** Prevents repeated inclusion of headers. */
#define ArithmeticalDSS_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include "DGtal/base/Exceptions.h"
#include "DGtal/base/Common.h"
#include "DGtal/kernel/PointVector.h"
#include "DGtal/kernel/CInteger.h"
#include "DGtal/io/DGtalBoard.h"
//////////////////////////////////////////////////////////////////////////////


namespace DGtal
{


/////////////////////////////////////////////////////////////////////////////
// class ArithmeticalDSS
/**
 * Description of class 'ArithmeticalDSS' <p>
 * \brief Aim:
 * Dynamic recognition of a digital straight segment (DSS)
 * defined as the sequence of connected points (x,y)
 * such that mu <= ax - by < mu + omega.  
 * (see Debled and Reveilles [1995]).
 *
 * This class is templated by the typename 'TInteger',
 * the type of the coordinates of the points (satisfying CInteger) 
 * and by the integer 'connectivity', which may be equal to 
 * 4 for standard (4-connected) DSS or 8 for naive (8-connected) DSS. 
 * Any other integers act as 8. The following shortcuts are provided in 
 * the namespace Z2i:
 * @code 
 typedef ArithmeticalDSS<Coordinate,4> > DSS4;
 typedef ArithmeticalDSS<Coordinate,8> > DSS8;
 * @endcode
 *
 * Here is a short example of how to use this class:
 * @code 

	// Input points
	std::vector<Point> contour;
	contour.push_back(Point(0,0));
	contour.push_back(Point(1,0));
	contour.push_back(Point(2,0));
	contour.push_back(Point(3,0));
	contour.push_back(Point(3,1));
	contour.push_back(Point(4,1));
	contour.push_back(Point(5,1));
	contour.push_back(Point(5,2));

  
  // Add points while it is possible
	DSS4 theDSS4(contour.at(0));		
	unsigned int i = 1;
	while ( (i<contour.size())
				&&(theDSS4.extend(contour.at(i))) ) {
		i++;
	}

	// Output parameters
	cout << theDSS4 << endl;
 //You must see:
 //[ArithmeticalDSS]
 //Parameters (a,b,mu,omega)=(2, 5, 0, 7)
 //First point (0,0) Last point (5,2)
 //Leaning points:
 //   Uf (0,0)
 //   Ul (5,2)
 //   Lf (3,0)
 //   Ll (3,0)
 //Steps:
 //   (1,0)
 //   (0,1)
 //[End ArithmeticalDSS]

 * @endcode
 */
template <typename TInteger, int connectivity>
class ArithmeticalDSS
{


    // ----------------------- Types ------------------------------
public:

  //2D point and 2D vector
  BOOST_CONCEPT_ASSERT(( CInteger<TInteger> ) );
  //does not compile
  typedef TInteger Integer;
  typedef DGtal::PointVector<2,Integer> Point;
  typedef DGtal::PointVector<2,Integer> Vector;
  
  typedef DGtal::PointVector<2,double> PointD;
  

//////////////////////////////////////////////////////////////////////////////
// generic class for computing the norm of (a,b)
// max(a,b) for 8-connectivity 
// |a|+|b| for 4-connectivity

		//default (for 8-connectivity) 
		template <typename TInt, int c>
		struct Tools
		{
			static TInt norm(const TInt& a, const TInt& b) 
			{
				DGtal::PointVector<2,TInt> v(a,b);
				return v.norm(DGtal::PointVector<2,TInt>::L_infty);
			}
			static TInt dualNorm(const TInt& a, const TInt& b) 
			{
				//does not compile with GMP type mpz_class
				//DGtal::PointVector<2,TInt> v(a,b);
				//return v.norm(DGtal::PointVector<2,TInt>::L_1);
				if (a > 0) {
					if (b > 0) {
						return (a+b);
					} else {
						return (a-b);		
					}
				} else {
					if (b > 0) {
						return (-a+b);						
					} else {
						return (-a-b);		
					}
				}
			}
		};

		//specialisation for 4-connectivity
		template <typename TInt>
		struct Tools<TInt,4> 
		{
			static TInt norm(const TInt& a, const TInt& b) 
			{
				//does not compile with GMP type mpz_class
				//DGtal::PointVector<2,TInt> v(a,b);
				//return v.norm(DGtal::PointVector<2,TInt>::L_1);
				if (a > 0) {
					if (b > 0) {
						return (a+b);
					} else {
						return (a-b);		
					}
				} else {
					if (b > 0) {
						return (-a+b);						
					} else {
						return (-a-b);		
					}
				}
			}

			static TInt dualNorm(const TInt& a, const TInt& b) 
			{
				//does not compile with GMP type mpz_class
				//DGtal::PointVector<2,TInt> v(a,b);
				//return v.norm(DGtal::PointVector<2,TInt>::L_infty);
				//why???
				TInt x;
				if (a>=0) x = a;
				else x = -a;
				TInt y;
				if (b>=0) y = b;
				else y = -b;
				return (x>=y)?x:y;
			}

		};

    // ----------------------- Standard services ------------------------------
public:


    /**
     * Default constructor.
     */
    ArithmeticalDSS();

    /**
     * Constructor.
		 * @param aPoint  
		 * a point
     */
    ArithmeticalDSS(const Point& aPoint);


    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    ArithmeticalDSS ( const ArithmeticalDSS & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    ArithmeticalDSS & operator= ( const ArithmeticalDSS & other );

    /**
     * Equality operator.
     * @param other the object to compare with.
     * @return 'true' either if the leaning points perfectly match
	 	 * or if the first leaning points match to the last ones
     * (same DSS scanned in the conversed way) 
     * and 'false' otherwise
     */
		bool operator==( const ArithmeticalDSS & other ) const;

    /**
     * Difference operator.
     * @param other the object to compare with.
     * @return 'false' if equal
     * 'true' otherwise
     */
		bool operator!=( const ArithmeticalDSS & other ) const;

    /**
     * Destructor.
     */
    ~ArithmeticalDSS(){};

    // ----------------------- Interface --------------------------------------
public:
     


    /**
		 * Tests whether the union between a point 
     * (adding to the front of the DSS 
     * with respect to the scan orientaion) 
		 * and a DSS is a DSS. 
     * Computes the parameters of the new DSS 
     * with the adding point if true.
     * @param aPoint the new pixel (connected to the DSS) 
     * @return 'true' if the union is a DSS, 'false' otherwise.
     */
    bool extend(const Point & aPoint);

    /**
		 * Removes the first point of a DSS
     * (located at the back with respect to 
     * the scan orientaion)
	   * if the DSS has more than two points
     * @return 'true' if the first point is removed, 'false' otherwise.
     */
    bool retract();

    /**
		 * Computes the sequence of (connected) points
		 * belonging to the DSL(a,b,mu,omega)
     * between the first and last point of the DSS
		 * Nb: in O(n)
     * @return the computed sequence of points.
     */
    std::vector<Point> retrieve() const;


    /**
		 * Computes the remainder of a point
     * (that does not necessarily belong to the DSS)
     * @param aPoint the point whose remainder is returned 
     * @return the remainder.
     */
    Integer getRemainder( const Point& aPoint ) const;

    /**
		 * Checks whether a point is in the DSL
     * of parameters (myA,myB,myMu,myOmega)
     * @return 'true' if yes, 'false' otherwise
     */
    bool isInDSL( const Point& aPoint ) const;

    /**
		 * Checks whether a point belongs to the DSS or not
     * @return 'true' if yes, 'false' otherwise
     */
    bool isInDSS( const Point& aPoint ) const;

    // ------------------------- Accessors ------------------------------
    /**
		 * myA accessor
     * @return an Integer of value myA.
     */
    Integer getA() const;
    /**
		 * myB accessor
     * @return an Integer of value myB.
     */
    Integer getB() const;
    /**
		 * myMu accessor
     * @return an Integer of value myMu.
     */
    Integer getMu() const;
    /**
		 * myOmega accessor
     * @return an Integer of value myOmega.
     */
    Integer getOmega() const;
    /**
		 * Accessor to the first upper leaning point
     * @return first upper leaning point.
     */
    Point getUf() const;
    /**
		 * Accessor to the last upper leaning point
     * @return last upper leaning point.
     */
    Point getUl() const;
    /**
		 * Accessor to the first lower leaning point
     * @return first lower leaning point.
     */
    Point getLf() const;
    /**
		 * Accessor to the last lower leaning point
     * @return last lower leaning point.
     */
    Point getLl() const;
    /**
		 * Accessor to the first added point to the DSS
     * @return point.
     */
    Point getBackPoint() const;
    /**
		 * Accessor to the last added point to the DSS
     * @return point.
     */
    Point getFrontPoint() const;

    /**
     * Checks the validity/consistency of the object.
     * NB: O(log(max(|myA|,|myB|)))
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    bool isValid() const;

// ------------------ Display ------------------------------------------

	public:
    /**
   * Projects the point [m] onto the average straight line (ie (mu+nu)/2).
   * @param m any point expressed in the local reference frame (may not be part of the segment).
   * @return the projected point.
   */
    PointD project( const Point & m ) const;

    /**
     * Projects the point [m] onto the straight line whose points have
     * remainder [r].
     *
     * @param m any point expressed in the local reference frame (may not
     * be part of the segment).
     *
     * @param r the remainder (may not be an integer).
     * @return the projected point.
     */
    PointD project( const Point & m, double r ) const;
    
    /**
     * Projects the point [m] onto the straight line going through point [p].
     *
     * @param m any point expressed in the local reference frame (may not
     * be part of the segment).
     *
     * @param p any point expressed in the local reference frame (may not
     * be part of the segment).
     *
     * @return the projected point.
     */
    PointD project( const Point & m, const Point & p ) const;
  

    /**
     * Defined as: norm( project(cp_n) - project(c_n) )
     * @return the projected length of the segment.
     * @see projectRegularly
     */
    double projectedSegmentLength() const;
    

    /**
     * Writes/Displays the object on an output stream.
     * @param out the output stream where the object is written.
     */
    void selfDisplay ( std::ostream & out ) ;

    
    /**
     * Draw the retrieved digital points of the DSS linked into a 
     * polygonal line on a LiBoard board
     * @param board the output board where the object is drawn.
     * @tparam Functor a Functor to specialize the Board style
     */
    
    void selfDrawAsDigitalPoints( DGtalBoard & board ) const;
    
    
    /**
     * Draw the bounding box of the DSS on a LiBoard board
     * @param board the output board where the object is drawn.
     * @tparam Functor a Functor to specialize the Board style
     */
    void selfDrawAsBoundingBox( DGtalBoard & board ) const;
    
    
    
    // ------------------------- Private Datas --------------------------------
 private:

    /**
     * Default style.
     */
    struct DefaultDrawStyleBB : public DrawableWithDGtalBoard
    {
      virtual void selfDraw(DGtalBoard & aBoard) const
      {
				// Set board style
				aBoard.setLineStyle(DGtalBoard::Shape::SolidStyle);
				aBoard.setPenColor(DGtalBoard::Color::Red);
				aBoard.setLineWidth(1);
				aBoard.setFillColor(DGtalBoard::Color::None);
      }
    };
    
    /**
       * Default style.
       */
    struct DefaultDrawStylePoints : public DrawableWithDGtalBoard
    {
        virtual void selfDraw(DGtalBoard & aBoard) const
        {
				// Set board style
				aBoard.setLineStyle(DGtalBoard::Shape::SolidStyle);
				aBoard.setPenColor(DGtalBoard::Color::Black);
				aBoard.setLineWidth(2);
				aBoard.setFillColor(DGtalBoard::Color::None);
			  }
    };

    // --------------- CDrawableWithDGtalBoard realization --------------------
  public:
    
    /**
     * Default drawing style object.
     * @return the dyn. alloc. default style for this object.
     */
    DrawableWithDGtalBoard* defaultStyle( std::string mode = "" ) const;
    
    /**
     * @return the style name used for drawing this object.
     */
    std::string styleName() const;





    /**
     * Draw the DSS on a LiBoard board as its bounding box and the
     * polyline of its points 
     * @see BoundingBoxDraw
     * @see DigitalPointsDraw
     * @param board the output board where the object is drawn.
     *
     */
    void selfDraw(DGtalBoard & board ) const;
    
    




    // ------------------------- Protected Datas ------------------------------
protected:

	//parameters of the DSS
	Integer myA, myB, myMu, myOmega;

	//leaning points
	Point myUf, myUl, myLf, myLl;

	//first and last added points
	Point myF, myL;

	//steps of the DSS 
  //e.g. right and up in the first octant
	std::vector<Vector> step;

    // ------------------------- Private Datas --------------------------------
	
private:

    // ------------------------- Hidden services ------------------------------
private:

		/**
		 * Checks whether the DSS has less or more
		 * than two displacement vectors (steps)
		 * between two consecutive points
		 * (must be called only in the main stage)
		 * @param aStep, the last displacement vector. 
		 * @return 'true' if less or equal, 'false' otherwise.
		 */
    bool hasLessThanTwoSteps(const Vector& aStep) const;


    /**
		 * Returns the 2D vector 
		 * starting at a point of remainder 0
		 * and pointing at the closer point of
     * remainder omega
     * @return the 2D vector.
     */
    Vector vectorFrom0ToOmega() const;

    /**
		 * Returns the point that follows a given 
     * point of a DSL
		 * @param a, b, mu, omega, the parameters 
     * of the DSL and aPoint, a given point of 
     * the DSL. 
     * @return the next point.
     */
    Point next(const Point& aPoint) const;




}; // end of class ArithmeticalDSS


/**
 * Modifier class in a DGtalBoard stream. Realizes the concept
 * CDrawableWithDGtalBoard.
 */
 struct DrawDSSBoundingBox : public DrawWithBoardModifier {
   void selfDraw( DGtalBoard & board ) const
   {
     board.myModes[ "ArithmeticalDSS" ] = "BoundingBox";
   }
 };
 
 /**
  * Modifier class in a DGtalBoard stream. Realizes the concept
  * CDrawableWithDGtalBoard.
  */
 struct DrawDSSPoints : public DrawWithBoardModifier {
   void selfDraw( DGtalBoard & board ) const
   {
     board.myModes[ "ArithmeticalDSS" ] = "Points";
   }
 };

/**
 * Overloads 'operator<<' for displaying objects of class 'ArithmeticalDSS'.
 * @param out the output stream where the object is written.
 * @param object the object of class 'ArithmeticalDSS' to write.
 * @return the output stream after the writing.
 */
template <typename TInteger, int connectivity>
std::ostream&
operator<< ( std::ostream & out,  ArithmeticalDSS<TInteger,connectivity> & object )
  {
      object.selfDisplay( out);
      return out;
    }


} // namespace DGtal



///////////////////////////////////////////////////////////////////////////////
// Includes inline functions/methods if necessary.
#if defined(INLINE)
#include "DGtal/geometry/2d/ArithmeticalDSS.ih"
#endif

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined ArithmeticalDSS_h

#undef ArithmeticalDSS_RECURSES
#endif // else defined(ArithmeticalDSS_RECURSES)
