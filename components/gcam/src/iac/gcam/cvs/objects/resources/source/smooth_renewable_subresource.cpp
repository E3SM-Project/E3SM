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

/*!
 * smooth_renewable_subresource.cpp
 * Created: 02/02/2007
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

// include files ***********************************************************

#include "util/base/include/definitions.h"

#include "resources/include/smooth_renewable_subresource.h"

#include <cassert>
#include <cmath>

// SmoothRenewableSubresource::sXMLName ************************************

const std::string SmoothRenewableSubresource::sXMLName
   = "smooth-renewable-subresource";

// Constructor: SmoothRenewableSubresource: ********************************

SmoothRenewableSubresource::SmoothRenewableSubresource(void)
   : parent(),
     mCostCurve(),
     mPriceExponent( 0.01 ),
     mMidPrice( 0 )
{
}

// Destructor: SmoothRenewableSubresource **********************************

SmoothRenewableSubresource::~SmoothRenewableSubresource(void)
{
}

// SmoothRenewableSubresource::annualsupply ********************************

void SmoothRenewableSubresource::annualsupply(
   int        aPeriod,
   const GDP* aGDP,
   double     aPrice,
   double     aPrevPrice )
{
   // Compute the fraction of the total possible supply that is
   // available at a given price
   double fractionAvailable = mCostCurve( aPrice );

   // Make supply increase continuously with price to improve convergence.
   // Default mPriceExponent value is very small so as not to significantly change resource base 
   // The factor of 5 below is arbitary, but was chosen so as to not change results signifiantly.
   // The equation below changes max resource value (using default  mPriceExponent) by 1% at 2 * mid-price.
   fractionAvailable *=
      std::pow( ( 1 + ( aPrice / ( 5.0 * mCostCurve.getMidprice() ) ) ), mPriceExponent );
    
   // Calculate expansion in supply due to GDP increase
   double gpdSupplyExpansion = std::pow(
      aGDP->getApproxGDP( aPeriod ) / aGDP->getApproxGDP( 0 ),
      gdpSupplyElasticity );

   // now convert to absolute value of production
   mAnnualProd[ aPeriod ] =
      fractionAvailable * maxSubResource * gpdSupplyExpansion;
}

// SmoothRenewableSubresource::completeInit ********************************

void SmoothRenewableSubresource::completeInit( const IInfo* aSectorInfo )
{
#if !defined( _MSC_VER )
   using std::exit;
#endif

   parent::completeInit( aSectorInfo );

   if ( !( mCostCurve.getMidprice() > 0 && mCostCurve.getCurveExponent() > 0 ) )
   // Invalid input parameter
   {
      ILogger& mainLog = ILogger::getLogger( "main_log" );
      mainLog.setLevel( ILogger::ERROR );
      mainLog << "Invalid input parameter(s) to "
         << getXMLNameStatic() << std::endl;
      exit( -1 );
   }
}

/*! \brief Perform any initializations needed for each period.
* \details Any initializations or calculations that only need to be done once per
*          period(instead of every iteration) should be placed in this function.
* \author Kate Calvin
* \param aRegionName Region name.
* \param aResourceName Resource name.
* \param aPeriod Model aPeriod
*/
void SmoothRenewableSubresource::initCalc(const std::string& aRegionName, const std::string& aResourceName, 
                                          const int aPeriod)
{
    SubResource::initCalc( aRegionName, aResourceName, aPeriod);

    // Reset the mid price to account for technical change
    // The mid price is the price at which 50% of potential supply is utilized
    // This implementation assumes techChange shifts the supply curve
    if( aPeriod == 0 ){
        mCumulativeTechChange[ aPeriod ] = 1.0;
    }
    else {
        const Modeltime* modeltime = scenario->getModeltime();
        mCumulativeTechChange[ aPeriod ] = mCumulativeTechChange[ aPeriod - 1 ] * 
                    pow( ( 1.0 + mTechChange[ aPeriod ] ), modeltime->gettimestep( aPeriod ) );
    }

    mCostCurve.setMidprice( mMidPrice / mCumulativeTechChange[ aPeriod ] );

}

// SmoothRenewableSubresource::getXMLName **********************************

const std::string& SmoothRenewableSubresource::getXMLName( void ) const
// Pre:
// Modifies:
// Post: Return the XML tag name
{
   return getXMLNameStatic();
}

// SmoothRenewableSubresource::getXMLNameStatic ****************************

const std::string& SmoothRenewableSubresource::getXMLNameStatic( void )
// Pre:
// Modifies:
// Post: Return the XML tag name
{
   return sXMLName;
}

// SmoothRenewableSubresource::toXMLforDerivedClass ************************

void SmoothRenewableSubresource::toXMLforDerivedClass(
   std::ostream& out,
   Tabs*         tabs ) const
{
   parent::toXMLforDerivedClass( out, tabs );

   XMLWriteElementCheckDefault(
      mMidPrice, "mid-price", out, tabs, 0.0 );
   XMLWriteElementCheckDefault(
      mCostCurve.getCurveExponent(), "curve-exponent", out, tabs, 0.0 );
   XMLWriteElementCheckDefault(
      mPriceExponent, "price-exponent", out, tabs, 0.01 );
}

// SmoothRenewableSubresource::XMLDerivedClassParse ************************

bool SmoothRenewableSubresource::XMLDerivedClassParse(
   const std::string&      nodeName,
   const xercesc::DOMNode* node )
{
   bool didParse = parent::XMLDerivedClassParse( nodeName, node );
   if ( !didParse )
   {
      if( nodeName == "mid-price" )
      {
          mMidPrice = XMLHelper<double>::getValue( node );
          mCostCurve.setMidprice( mMidPrice );
          didParse  = true;
       }
       else if( nodeName == "curve-exponent" )
      {
          mCostCurve.setCurveExponent( XMLHelper<double>::getValue( node ) );
          didParse       = true;
       }
      else if ( nodeName == "price-exponent" )
      {
         mPriceExponent = XMLHelper<double>::getValue( node );
         didParse       = true;
      }
   }

   return didParse;
}

// end of smooth_renewable_subresource.cpp *********************************


