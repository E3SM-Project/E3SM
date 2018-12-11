/*
 * LEGAL NOTICE
 * This computer software was prepared by Battelle Memorial Institute,
 * hereinafter the Contractor, under Contract No. DE-AC05-76RL0 1830
 * with the Department of Energy (DOE). NEITHER THE GOVERNMENT NOR THE
 * CONTRACTOR MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY
 * LIABILITY FOR THE USE OF THIS SOFTWARE. This notice including this
 * sentence must appear on any copies of this computer software.
 * 
 * EXPORT CONTROL
 * User agrees that the Software will not be shipped, transferred or
 * exported into any country or used in any manner prohibited by the
 * United States Export Administration Act or any other applicable
 * export laws, restrictions or regulations (collectively the "Export Laws").
 * Export of the Software may require some form of license or other
 * authority from the U.S. Government, and failure to obtain such
 * export control license may result in criminal liability under
 * U.S. laws. In addition, if the Software is identified as export controlled
 * items under the Export Laws, User represents and warrants that User
 * is not a citizen, or otherwise located within, an embargoed nation
 * (including without limitation Iran, Syria, Sudan, Cuba, and North Korea)
 *     and that User is not otherwise prohibited
 * under the Export Laws from receiving the Software.
 * 
 * Copyright 2011 Battelle Memorial Institute.  All Rights Reserved.
 * Distributed as open-source under the terms of the Educational Community 
 * License version 2.0 (ECL 2.0). http://www.opensource.org/licenses/ecl2.php
 * 
 * For further details, see: http://www.globalchange.umd.edu/models/gcam/
 * 
 */

/* cost_curve.h
 * Created: 02/06/2007
 * Version: 02/21/2007
 *
 * This software, which is provided in confidence, was prepared by employees
 * of Pacific Northwest National Laboratory operated by Battelle Memorial
 * Institute. Battelle has certain unperfected rights in the software
 * which should not be copied or otherwise disseminated outside your
 * organization without the express written authorization from Battelle.
 * All rights to the software are reserved by Battelle.   Battelle makes no
 * warranty, express or implied, and assumes no liability or responsibility
 * for the use of this software.
 */

#if !defined( __COST_CURVE_H )
#define __COST_CURVE_H     // prevent multiple includes

// include files ***********************************************************

#include <cmath>

// namespaces **************************************************************

namespace ObjECTS {

// class: TCostCurve<T> ****************************************************

/*! \ingroup Objects
 *  \brief TCostCurve<T> is a class template that is be used to compute
 *         cost curves.
 *  \details The class template TCostCurve<T> is used to compute cost
 *           curves based on a price, mid price and curve exponent.
 *           In the case of multiple calculations using a variable price and
 *           constant mid price and curve exponent, part of the calculation
 *           is cached, thus reducing the number of function calls to the
 *           standard C++ math library.
 * 
 *  \author Kevin Walker
 */
template <class T = double>
class TCostCurve
{
public :

   typedef T value_type;

   // Constructors
   TCostCurve( void )
      : mMidprice( 0 ),
        mCurveExponent( 0 ),
        mCachedValue( 0 ),
        mDirty( true ) {}
   TCostCurve(
      const value_type& midprice,
      const value_type& curveExponent )
      : mMidprice( midprice ),
        mCurveExponent( curveExponent ),
        mCachedValue( 0 ),
        mDirty( true ) {}
   TCostCurve(const TCostCurve<T>& other)
      : mMidprice( other.mMidprice ),
        mCurveExponent( other.mCurveExponent ),
        mCachedValue( other.mCachedValue ),
        mDirty( other.mDirty ) {}

   // Destructor
   virtual ~TCostCurve(void) {}

   // Assignment operator
   TCostCurve<T>& operator = (const TCostCurve<T>& other)
   {
      if ( &other != this )
      {
         mMidprice      = other.mMidprice;
         mCurveExponent = other.mCurveExponent;
         mCachedValue   = other.mCachedValue;
         mDirty         = other.mDirty;
      }
      return *this;
   }

   /*! Calculate the cost curve value for the specified price
    *  \param aPrice the current price
    *  \return the cost curve value for the specified price
    */
   value_type operator () ( const value_type& aPrice ) const
   { return calculate( aPrice ); }
   virtual value_type calculate( const value_type& aPrice ) const;

   /*! Calculate the cost curve value for the specified price, mid price
    *  and curve exponent
    *  \param aPrice the current price
    *  \param aMidprice the price where 50% of the maximum supply is brought
    *         to market
    *  \param aCurveExponent the steepness parameter of the cost curve
    *  \return the cost curve value for the specified price, mid price
    *          and curve exponent
    */
   static value_type calculate(
      const value_type& aPrice,
      const value_type& aMidprice,
      const value_type& aCurveExponent );

   //! Get the steepness parameter of the cost curve
   virtual const value_type& getCurveExponent( void ) const
   { return mCurveExponent; }

   //! Get the price where 50% of the maximum supply is brought to market
   virtual const value_type& getMidprice( void ) const { return mMidprice; }

   //! Set the steepness parameter of the cost curve
   virtual void setCurveExponent( const value_type& aCurveExponent )
   {
      mCurveExponent = aCurveExponent;
      mDirty         = true;
   }

   //! Set the price where 50% of the maximum supply is brought to market
   virtual void setMidprice( const value_type& aMidprice )
   {
      mMidprice = aMidprice;
      mDirty    = true;
   }

private :

   //! The price where 50% of the maximum supply is brought to market
   value_type           mMidprice;

   //! The steepness parameter of the cost curve
   value_type           mCurveExponent;

   //! pow( midprice, curve-exponent ) for optimisation
   mutable value_type   mCachedValue;

   //! Flag if mCachedValue needs to be computed
   mutable bool         mDirty;
};

// TCostCurve<T>::calculate ************************************************

template <class T>
inline typename TCostCurve<T>::value_type TCostCurve<T>::calculate(
   const value_type& aPrice ) const
{
   // Validate the price
   if ( aPrice <= value_type( 0 ) )
   {
      return value_type( 0 );
   }

   // Check if to update the cached value
   if ( mDirty )
   {
      mCachedValue = std::pow( mMidprice, mCurveExponent );
      mDirty       = false;
   }

   value_type  temp = std::pow( aPrice, mCurveExponent );
   return temp / ( mCachedValue + temp );
}

template <class T>
inline typename TCostCurve<T>::value_type TCostCurve<T>::calculate(
   const value_type& aPrice,
   const value_type& aMidprice,
   const value_type& aCurveExponent )
{
   return TCostCurve<T>( aMidprice, aCurveExponent )( aPrice );
}

} // namespace ObjECTS

#endif   // __COST_CURVE_H

// end of cost_curve.h *****************************************************

