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
 * @file IntegralInvariantGaussianCurvatureEstimator_0memory.h
 * @author Jeremy Levallois (\c jeremy.levallois@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), INSA-Lyon, France
 * LAboratoire de MAthématiques - LAMA (CNRS, UMR 5127), Université de Savoie, France
 *
 * @date 2012/07/10
 *
 * Header file for module IntegralInvariantGaussianCurvatureEstimator_0memory.cpp
 *
 * @brief Compute Gaussian curvature on border of shapes of n-dimension, based on integral invariant.
 *
 * @see related article:
 *       Coeurjolly, D.; Lachaud, J.O; Levallois, J., (2013). Integral based Curvature
 *       Estimators in Digital Geometry. DGCI 2013. Retrieved from
 *       https://liris.cnrs.fr/publis/?id=5866
 *
 * This file is part of the DGtal library.
 */

#if defined(IntegralInvariantGaussianCurvatureEstimator_0memory_RECURSES)
#error Recursive header files inclusion detected in IntegralInvariantGaussianCurvatureEstimator_0memory.h
#else // defined(IntegralInvariantGaussianCurvatureEstimator_0memory_RECURSES)
/** Prevents recursive inclusion of headers. */
#define IntegralInvariantGaussianCurvatureEstimator_0memory_RECURSES

#if !defined IntegralInvariantGaussianCurvatureEstimator_0memory_h
/** Prevents repeated inclusion of headers. */
#define IntegralInvariantGaussianCurvatureEstimator_0memory_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include "DGtal/base/Common.h"

#include "DGtal/shapes/GaussDigitizer.h"
#include "DGtal/shapes/Shapes.h"

#include "DGtal/geometry/surfaces/DigitalSurfaceConvolver.h"
#include "DGtal/geometry/surfaces/DigitalSurfaceConvolver_0memory.h"
#include "DGtal/shapes/EuclideanShapesDecorator.h"

#include "DGtal/shapes/implicit/ImplicitBall.h"

#include "DGtal/math/EigenValues3D.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

template< typename Matrix3x3 >
class GaussianCurvatureFunctor3 : std::unary_function <double,double>
{
public:
    typedef double Value;
    typedef EigenValues3D< Value >::Vector3 Vector3;

    GaussianCurvatureFunctor3(){}

    void init( const double & h, const double & r )
    {
        double r3 = r * r * r;
        double r6 = r3 * r3;
        d6_PIr6 = 6.0 / ( M_PI * r6 );
        d8_5r = 8.0 / ( 5.0 * r );
        double h2 = h * h;
        dh5 = h2 * h2 * h;
    }

    Value operator()(const Matrix3x3 & aInput)
    {
        Matrix3x3 cp_matrix = aInput;
//        std::cout << cp_matrix << std::endl;
        cp_matrix *= dh5;
        Value k1,k2;
        Matrix3x3 eigenVectors;
        Vector3 eigenValues;
        evalk1k2( cp_matrix, eigenVectors, eigenValues, k1, k2 );
        return k1 * k2;
    }

protected:
    void evalk1k2(
            Matrix3x3 & matrix,
            Matrix3x3 & eigenVectors,
            Vector3 & eigenValues,
            Value & k1,
            Value & k2 )
    {
        EigenValues3D< Value >::getEigenDecomposition( matrix, eigenVectors, eigenValues );

        ASSERT ( eigenValues[ 0 ] == eigenValues[ 0 ] ); // NaN
        ASSERT ( (eigenValues[ 0 ] <= eigenValues[ 2 ]) && (eigenValues[ 0 ] <= eigenValues[ 1 ]) && (eigenValues[ 1 ] <= eigenValues[ 2 ]) );

        k1 = d6_PIr6 * ( eigenValues[ 1 ] - ( 3.0 * eigenValues[ 2 ] )) + d8_5r;
        k2 = d6_PIr6 * ( eigenValues[ 2 ] - ( 3.0 * eigenValues[ 1 ] )) + d8_5r;
    }

private:
    Value dh5;
    Value d6_PIr6;
    Value d8_5r;

};

template< typename Matrix2x2 >
class GaussianCurvatureFunctor2 : std::unary_function <double,double>
{
public:
    typedef double Value;

    GaussianCurvatureFunctor2(){}

    void init( const double & h, const double & r )
    {
        dh2 = h * h;
        d3_r = 3.0 / r;
        dPI_2 = M_PI / 2.0;
        d1_r2 = 1.0 / ( r * r );
    }

    Value operator()(const Matrix2x2 & aInput)
    {
        Matrix2x2 cp_matrix = aInput;
//        std::cout << cp_matrix << std::endl;
        cp_matrix *= dh2;
        Value k = d3_r * ( dPI_2 - d1_r2 * cp_matrix[ 0 ] );

        return k;
    }

private:
    Value dh2; /// h*h
    Value d3_r; /// 3/r
    Value dPI_2; /// PI/2
    Value d1_r2; /// 1/r^2

};


/////////////////////////////////////////////////////////////////////////////
// template class IntegralInvariantGaussianCurvatureEstimator_0memory
/**
   * Description of template class 'IntegralInvariantMeanCurvatureEstimator' <p>
   * \brief Aim: This class implement a Integral Invariant Gaussian curvature estimation.
   *
   * @see related article:
   *       Coeurjolly, D.; Lachaud, J.O; Levallois, J., (2013). Integral based Curvature
   *       Estimators in Digital Geometry. DGCI 2013. Retrieved from
   *       https://liris.cnrs.fr/publis/?id=5866
   *
   * The algorithm we propose uses covariance matrix to approximate Gaussian curvature.
   * To compute the covariance matrix on each cell of the surface, we convolve a kernel (2D: Ball2D, 3D: Ball3D) around
   * the surface.
   * The result covariance matrix give us principal curvatures and principal directions by computing
   * respectively eigenvectors and eigenvalues on it.
   * Experimental results showed a multigrid convergence.
   *
   * Some optimization are available when we set a range of 0-adjacent surfels to the estimator.
   *
   * @tparam TKSpace space in which the shape is defined.
   * @tparam TShapeFunctor TFunctor a model of a functor for the shape ( f(x) ).
   * @tparam dimension dimension of the shape. Let default value to use the correct specialization.
   *
   * @see exampleIntegralInvariantGaussianCurvature3D.cpp testIntegralInvariantGaussianCurvature3D.cpp
   */
template <typename TKSpace, typename TShapeFunctor, Dimension dimension = TKSpace::dimension>
class IntegralInvariantGaussianCurvatureEstimator_0memory
{
public:
    typedef TKSpace KSpace;
    typedef typename Z2i::Domain Domain;
    typedef typename KSpace::Space::RealPoint RealPoint;
    typedef typename Z2i::DigitalSet DigitalSet;
    typedef typename KSpace::SCell Spel;
    typedef typename KSpace::SurfelSet SurfelSet;
    typedef typename SurfelSet::const_iterator ConstIteratorKernel;

    typedef double Quantity;
    typedef int Value;

    typedef ImplicitBall<Z2i::Space> KernelSupport;
    typedef EuclideanShapesMinus< KernelSupport, KernelSupport > EuclideanMinus;

    typedef TShapeFunctor ShapeCellFunctor;
    typedef ConstValueFunctor<Value> KernelCellFunctor;
    typedef DigitalSurfaceConvolver_0memory<ShapeCellFunctor, KernelCellFunctor, KSpace, ConstIteratorKernel, EuclideanMinus> Convolver;
    typedef typename Convolver::PairIterators PairIterators;

    //  BOOST_CONCEPT_ASSERT (( CCellFunctor< ShapeCellFunctor > ));

    // ----------------------- Standard services ------------------------------
public:
    /**
     * Constructor.
     *
     * @param space space in which the shape is defined.
     * @param f functor on cell of the shape.
     */
    IntegralInvariantGaussianCurvatureEstimator_0memory ( const KSpace & space, const ShapeCellFunctor & f );

    /**
     * Destructor.
     */
    ~IntegralInvariantGaussianCurvatureEstimator_0memory()
    {}

    // ----------------------- Interface --------------------------------------
public:

    /**
      * Initialise the IntegralInvariantGaussianCurvatureEstimator_0memory with a specific Euclidean kernel re, and grid step h.
      *
      * @param _h precision of the grid
      * @param re Euclidean radius of the kernel support
      */
    void init ( const double _h, const double re );

    /**
      * Compute the integral invariant Gaussian curvature to cell *it of a shape.
      *
      * @tparam ConstIteratorOnCells iterator on a Cell
      *
      * @param it iterator of a cell (from a shape) we want compute the integral invariant curvature.
      *
      * @return quantity of the result of Integral Invariant estimator at position *it
      */
    template<typename ConstIteratorOnCells>
    Quantity eval ( const ConstIteratorOnCells & it );

    /**
      * Compute the integral invariant Gaussian curvature from two cells (from *itb to *ite (exclude) ) of a shape.
      * Return the result on an OutputIterator (param).
      *
      * @tparam ConstIteratorOnCells iterator on a Cell
      * @tparam OutputIterator iterator of a list of Quantity
      *
      * @param ite iterator of the begin position on the shape where we compute the integral invariant curvature.
      * @param itb iterator of the end position (excluded) on the shape where we compute the integral invariant curvature.
      * @param result iterator of results of the computation.
      */
    template< typename ConstIteratorOnCells, typename OutputIterator >
    void eval ( const ConstIteratorOnCells & itb,
                const ConstIteratorOnCells & ite,
                OutputIterator & result );

    /**
      * @return iterator of the begin spel of the kernel support
      */
    const ConstIteratorKernel & beginKernel() const;

    /**
      * @return iterator of the end spel of the kernel support
      */
    const ConstIteratorKernel & endKernel() const;

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

    // ------------------------- Private Datas --------------------------------
private:
    /// array of shifting masks.
    std::vector< SurfelSet > kernels;
    /// array of begin/end iterator of shifting masks.
    std::vector< PairIterators > kernelsIterators;

    /// origin spel of the kernel support
    Spel myOrigin;

    /// kernel functor
    const KernelCellFunctor myKernelFunctor;

    /// convolver
    Convolver myConvolver;

    /// precision of the grid
    double h;

    /// Euclidean radius of the kernel
    double radius;

private:

    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    IntegralInvariantGaussianCurvatureEstimator_0memory ( const IntegralInvariantGaussianCurvatureEstimator_0memory & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    IntegralInvariantGaussianCurvatureEstimator_0memory & operator= ( const IntegralInvariantGaussianCurvatureEstimator_0memory & other );

}; // end of class IntegralInvariantGaussianCurvatureEstimator_0memory

/**
      * Specialization for dimension = 2
      */
template <typename TKSpace, typename TShapeFunctor>
class IntegralInvariantGaussianCurvatureEstimator_0memory<TKSpace, TShapeFunctor, 2>
{
public:
    typedef TKSpace KSpace;
    typedef typename Z2i::Domain Domain;
    typedef typename KSpace::Space::RealPoint RealPoint;
    typedef typename Z2i::DigitalSet DigitalSet;
    typedef typename KSpace::SCell Spel;
    typedef typename KSpace::SurfelSet SurfelSet;
    typedef typename SurfelSet::const_iterator ConstIteratorKernel;

    typedef double Quantity;
    typedef int Value;

    typedef ImplicitBall<Z2i::Space> KernelSupport;
    typedef EuclideanShapesMinus< KernelSupport, KernelSupport > EuclideanMinus;
    typedef GaussDigitizer< Z3i::Space, KernelSupport > DigitalShapeKernel;
    typedef GaussDigitizer< Z2i::Space, EuclideanMinus > DigitalShape;

    typedef TShapeFunctor ShapeCellFunctor;
    typedef ConstValueFunctor<Value> KernelCellFunctor;
    typedef DigitalSurfaceConvolver_0memory<ShapeCellFunctor, KernelCellFunctor, KSpace, ConstIteratorKernel, DigitalShapeKernel> Convolver;
    typedef typename Convolver::PairIterators PairIterators;

    typedef typename Convolver::CovarianceMatrix Matrix2x2;



    //  BOOST_CONCEPT_ASSERT (( CCellFunctor< ShapeCellFunctor > ));

    // ----------------------- Standard services ------------------------------
public:
    /**
     * Constructor.
     *
     * @param space space in which the shape is defined.
     * @param f functor on cell of the shape.
     */
    IntegralInvariantGaussianCurvatureEstimator_0memory ( const KSpace & space, const ShapeCellFunctor & f );

    /**
     * Destructor.
     */
    ~IntegralInvariantGaussianCurvatureEstimator_0memory()
    {
        for( unsigned int i = 0; i < kernelsSet.size(); ++i )
        {
            delete kernelsSet[ i ];
        }
        kernelsSet.clear();
        delete kernel;
        delete digKernel;
    }

    // ----------------------- Interface --------------------------------------
public:

    /**
      * Initialise the IntegralInvariantGaussianCurvatureEstimator_0memory with a specific Euclidean kernel radius re, and grid step h.
      *
      * @param _h precision of the grid
      * @param re Euclidean radius of the kernel support
      *
      * @bug known bug with radius of kernel. Small hack for the moment.
      */
    void init ( const double _h, const double re );

    /**
      * Compute the integral invariant Gaussian curvature to cell *it of a shape.
      *
      * @tparam ConstIteratorOnCells iterator on a Cell
      *
      * @param it iterator of a cell (from a shape) we want compute the integral invariant curvature.
      *
      * @return quantity of the result of Integral Invariant estimator at position *it
      */
    template<typename ConstIteratorOnCells>
    Quantity eval ( const ConstIteratorOnCells & it );

    /**
      * Compute the integral invariant Gaussian curvature to cell *it of a shape.
      *
      * @tparam ConstIteratorOnCells iterator on a Cell
      *
      * @param it iterator of a cell (from a shape) we want compute the integral invariant curvature.
      *
      * @return quantity of the result of Integral Invariant estimator at position *it
      */
    template<typename ConstIteratorOnCells, typename Shape>
    Quantity eval ( const ConstIteratorOnCells & it,
                    const Shape & shape );

    /**
      * Compute the integral invariant Gaussian curvature from two cells (from *itb to *ite (exclude) ) of a shape.
      * Return the result on an OutputIterator (param).
      *
      * @tparam ConstIteratorOnCells iterator on a Cell
      * @tparam OutputIterator iterator of a list of Quantity
      *
      * @param ite iterator of the begin position on the shape where we compute the integral invariant curvature.
      * @param itb iterator of the end position (excluded) on the shape where we compute the integral invariant curvature.
      * @param result iterator of results of the computation.
      */
    template< typename ConstIteratorOnCells, typename OutputIterator >
    void eval ( const ConstIteratorOnCells & itb,
                const ConstIteratorOnCells & ite,
                OutputIterator & result );

    /**
      * Compute the integral invariant Gaussian curvature from two cells (from *itb to *ite (exclude) ) of a shape.
      * Return the result on an OutputIterator (param).
      *
      * @tparam ConstIteratorOnCells iterator on a Cell
      * @tparam OutputIterator iterator of a list of Quantity
      *
      * @param ite iterator of the begin position on the shape where we compute the integral invariant curvature.
      * @param itb iterator of the end position (excluded) on the shape where we compute the integral invariant curvature.
      * @param result iterator of results of the computation.
      */
    template< typename ConstIteratorOnCells, typename OutputIterator, typename Shape >
    void eval ( const ConstIteratorOnCells & itb,
                const ConstIteratorOnCells & ite,
                OutputIterator & result,
                const Shape & shape );


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

    // ------------------------- Private Datas --------------------------------
private:

    /// array of shifting masks. Size = 9 for each shiftings (0-adjacent and full kernel included)
    std::vector< PairIterators > kernels;
    std::vector< DigitalSet* > kernelsSet;
//    std::vector< EuclideanMinus* > masks;
//    std::vector< KernelSupport* > kernelsShifted;
    KernelSupport * kernel;
    DigitalShapeKernel * digKernel;

    /// origin spel of the kernel support.
    Spel myOrigin;

    /// kernel functor
    const KernelCellFunctor myKernelFunctor;

    /// convolver
    Convolver myConvolver;

    /// precision of the grid
    double h;

    /// Euclidean radius of the kernel
    double radius;

    /// kernel's radius-dependant variable. Used to compute IntegralInvariant.
    /*Quantity dh2; /// h*h
    Quantity d3_r; /// 3/r
    Quantity dPI_2; /// PI/2
    Quantity d1_r2; /// 1/r^2*/

    GaussianCurvatureFunctor2< Matrix2x2 > gaussFunctor;

private:

    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    IntegralInvariantGaussianCurvatureEstimator_0memory ( const IntegralInvariantGaussianCurvatureEstimator_0memory & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    IntegralInvariantGaussianCurvatureEstimator_0memory & operator= ( const IntegralInvariantGaussianCurvatureEstimator_0memory & other );
}; // end of class IntegralInvariantGaussianCurvatureEstimator_0memory for dimension = 2

/**
    * Specialization for dimension = 3
    */
template <typename TKSpace, typename TShapeFunctor>
class IntegralInvariantGaussianCurvatureEstimator_0memory<TKSpace, TShapeFunctor, 3>
{
public:
    typedef TKSpace KSpace;
    typedef typename Z3i::Domain Domain;
    typedef typename KSpace::Space::RealPoint RealPoint;
    typedef typename Z3i::DigitalSet DigitalSet;
    typedef typename KSpace::SCell Spel;
    typedef typename KSpace::SurfelSet SurfelSet;
    typedef typename SurfelSet::const_iterator ConstIteratorKernel;

    typedef double Quantity;
    typedef int Value;

    typedef ImplicitBall<Z3i::Space> KernelSupport;
    typedef EuclideanShapesMinus< KernelSupport, KernelSupport > EuclideanMinus;
    typedef GaussDigitizer< Z3i::Space, KernelSupport > DigitalShapeKernel;
    typedef GaussDigitizer< Z3i::Space, EuclideanMinus > DigitalShape;

    typedef TShapeFunctor ShapeCellFunctor;
    typedef ConstValueFunctor<Value> KernelCellFunctor;
    typedef DigitalSurfaceConvolver_0memory<ShapeCellFunctor, KernelCellFunctor, KSpace, ConstIteratorKernel, DigitalShapeKernel> Convolver;
    typedef typename Convolver::PairIterators PairIterators;

    typedef typename Convolver::CovarianceMatrix Matrix3x3;
    typedef EigenValues3D< Quantity >::Vector3 Vector3;

    //  BOOST_CONCEPT_ASSERT (( CCellFunctor< ShapeCellFunctor > ));

    // ----------------------- Standard services ------------------------------
public:
    /**
     * Constructor.
     *
     * @param space space in which the shape is defined.
     * @param f functor on cell of the shape.
     */
    IntegralInvariantGaussianCurvatureEstimator_0memory ( const KSpace & space, const ShapeCellFunctor & f );

    /**
     * Destructor.
     */
    ~IntegralInvariantGaussianCurvatureEstimator_0memory()
    {
        /*for( unsigned int i = 0; i < masks.size(); ++i )
        {
            delete masks[ i ];
        }
        for( unsigned int i = 0; i < kernelsShifted.size(); ++i )
        {
            delete kernelsShifted[ i ];
        }
        kernelsShifted.clear();
        masks.clear();*/
        for( unsigned int i = 0; i < kernelsSet.size(); ++i )
        {
            delete kernelsSet[ i ];
        }
        kernelsSet.clear();
        delete kernel;
        delete digKernel;
    }

    // ----------------------- Interface --------------------------------------
public:

    /**
      * Initialise the IntegralInvariantGaussianCurvatureEstimator_0memory with a specific Euclidean kernel radius re, and grid step k.
      *
      * @param _h precision of the grid
      * @param re Euclidean radius of the kernel support
      *
      * @bug known bug with radius of kernel. Small hack for the moment.
      */
    void init ( const double _h, const double re );

    /**
      * Compute the integral invariant Gaussian curvature to cell *it of a shape.
      *
      * @tparam ConstIteratorOnCells iterator on a Cell
      *
      * @param it iterator of a cell (from a shape) we want compute the integral invariant curvature.
      *
      * @return quantity of the result of Integral Invariant estimator at position *it
      */
    template<typename ConstIteratorOnCells>
    Quantity eval ( const ConstIteratorOnCells & it );

    /**
      * Compute the integral invariant Gaussian curvature to cell *it of a shape.
      *
      * @tparam ConstIteratorOnCells iterator on a Cell
      *
      * @param it iterator of a cell (from a shape) we want compute the integral invariant curvature.
      *
      * @return quantity of the result of Integral Invariant estimator at position *it
      */
    template<typename ConstIteratorOnCells, typename Shape>
    Quantity eval ( const ConstIteratorOnCells & it,
                    const Shape & shape );

    /**
      * Compute the integral invariant Gaussian curvature from two cells (from *itb to *ite (exclude) ) of a shape.
      * Return the result on an OutputIterator (param).
      *
      * @tparam ConstIteratorOnCells iterator on a Cell
      * @tparam OutputIterator iterator of a list of Quantity
      *
      * @param ite iterator of the begin position on the shape where we compute the integral invariant curvature.
      * @param itb iterator of the end position (excluded) on the shape where we compute the integral invariant curvature.
      * @param result iterator of results of the computation.
      */
    template< typename ConstIteratorOnCells, typename OutputIterator >
    void eval ( const ConstIteratorOnCells & itb,
                const ConstIteratorOnCells & ite,
                OutputIterator & result );

    /**
      * Compute the integral invariant Gaussian curvature from two cells (from *itb to *ite (exclude) ) of a shape.
      * Return the result on an OutputIterator (param).
      *
      * @tparam ConstIteratorOnCells iterator on a Cell
      * @tparam OutputIterator iterator of a list of Quantity
      *
      * @param ite iterator of the begin position on the shape where we compute the integral invariant curvature.
      * @param itb iterator of the end position (excluded) on the shape where we compute the integral invariant curvature.
      * @param result iterator of results of the computation.
      */
    template< typename ConstIteratorOnCells, typename OutputIterator, typename Shape >
    void eval ( const ConstIteratorOnCells & itb,
                const ConstIteratorOnCells & ite,
                OutputIterator & result,
                const Shape & shape );

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

    // ------------------------- Private Datas --------------------------------
private:

    /// array of shifting masks. Size = 27 for each shiftings (0-adjacent and full kernel included)
    std::vector< PairIterators > kernels;
    std::vector< DigitalSet* > kernelsSet;
//    std::vector< EuclideanMinus* > masks;
//    std::vector< KernelSupport* > kernelsShifted;
    KernelSupport * kernel;
    DigitalShapeKernel * digKernel;

    /// origin spel of the kernel support.
    Spel myOrigin;

    /// kernel functor
    const KernelCellFunctor myKernelFunctor;

    /// convolver
    Convolver myConvolver;

    /// precision of the grid
    double h;

    /// Euclidean radius of the kernel
    double radius;

    /*Quantity d6_PIr6; /// 6/PI*r^6
    Quantity d8_5r; /// 8/5r
    Quantity dh5; /// h^5*/

    GaussianCurvatureFunctor3< Matrix3x3 > gaussFunctor;

private:

    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    IntegralInvariantGaussianCurvatureEstimator_0memory ( const IntegralInvariantGaussianCurvatureEstimator_0memory & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    IntegralInvariantGaussianCurvatureEstimator_0memory & operator= ( const IntegralInvariantGaussianCurvatureEstimator_0memory & other );

}; // end of specialization for dimension = 3









/**
   * Overloads 'operator<<' for displaying objects of class 'IntegralInvariantGaussianCurvatureEstimator_0memory'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'IntegralInvariantGaussianCurvatureEstimator_0memory' to write.
   * @return the output stream after the writing.
   */
template <typename TKS, typename TSF, Dimension dimension>
std::ostream&
operator<< ( std::ostream & out, const IntegralInvariantGaussianCurvatureEstimator_0memory<TKS, TSF, dimension> & object );

template <typename TKS, typename TSF>
std::ostream&
operator<< ( std::ostream & out, const IntegralInvariantGaussianCurvatureEstimator_0memory<TKS, TSF, 2> & object );

template <typename TKS, typename TSF>
std::ostream&
operator<< ( std::ostream & out, const IntegralInvariantGaussianCurvatureEstimator_0memory<TKS, TSF, 3> & object );

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantGaussianCurvatureEstimator_0memory.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined IntegralInvariantGaussianCurvatureEstimator_0memory_h

#undef IntegralInvariantGaussianCurvatureEstimator_0memory_RECURSES
#endif // else defined(IntegralInvariantGaussianCurvatureEstimator_0memory_RECURSES)
