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
 * \file building_final_demand.cpp
 * \ingroup Objects
 * \brief BuildingFinalDemand class source file.
 * \author Steve Smith, Josh Lurz
 */

#include "util/base/include/definitions.h"
#include <string>
#include <cmath>
#include <algorithm>

#include <xercesc/dom/DOMNode.hpp>
#include "sectors/include/building_final_demand.h"
#include "util/base/include/xml_helper.h"
#include "util/base/include/model_time.h"
#include "marketplace/include/marketplace.h"
#include "containers/include/scenario.h"
#include "containers/include/gdp.h"
#include "util/base/include/model_time.h"
#include "demographics/include/demographic.h"
#include "sectors/include/sector_utils.h"
#include "containers/include/iinfo.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;

/*!
 * \brief Constructor.
 * \author Steve Smith, Josh Lurz
 */
BuildingFinalDemand::BuildingFinalDemand()
{
}

//! Destructor
BuildingFinalDemand::~BuildingFinalDemand() {
}

bool BuildingFinalDemand::XMLDerivedClassParse( const string& nodeName, const DOMNode* curr ) {
    if( nodeName == "saturation-elasticity" ){
        mSaturationElasticity = XMLHelper<double>::getValue( curr );
    }
    else if( nodeName == "saturation-point" ){
        mSaturationPoint = XMLHelper<double>::getValue( curr );
    }
    else if( EnergyFinalDemand::XMLDerivedClassParse( nodeName, curr ) ){
        return false;
    }
    // If was true somewhere above then noce was parsed
    return true;
}

void BuildingFinalDemand::toInputXMLDerived( ostream& out, Tabs* tabs ) const {  
    EnergyFinalDemand::toInputXMLDerived( out, tabs );
    XMLWriteElement( mSaturationElasticity, "saturation-elasticity", out, tabs );
    XMLWriteElement( mSaturationPoint, "saturation-point", out, tabs );
}   

void BuildingFinalDemand::toDebugXMLDerived( const int period, ostream& out, Tabs* tabs ) const {
    EnergyFinalDemand::toDebugXMLDerived( period, out, tabs );
    XMLWriteElement( mSaturationElasticity, "saturation-elasticity", out, tabs );
    XMLWriteElement( mSaturationPoint, "saturation-point", out, tabs );
}

const string& BuildingFinalDemand::getXMLName() const {
    return getXMLNameStatic();
}

/*! \brief Get the XML node name in static form for comparison when parsing XML.
*
* This public function accesses the private constant string, XML_NAME.
* This way the tag is always consistent for both read-in and output and can be easily changed.
* The "==" operator that is used when parsing, required this second function to return static.
* \note A function cannot be static and virtual.
* \author Josh Lurz, James Blackwood
* \return The constant XML_NAME as a static.
*/
const string& BuildingFinalDemand::getXMLNameStatic() {
    static const string XML_NAME = "building-final-demand";
    return XML_NAME;
}

void BuildingFinalDemand::completeInit( const string& aRegionName,
                                         const IInfo* aRegionInfo )
{
    // Perhaps this should be in the base class completeInit, but not sure if there could be
    // demands without a FinalEnergyConsumer?
    if ( !mFinalEnergyConsumer.get() ) {
        mFinalEnergyConsumer.reset( new FinalEnergyConsumer( mName ) );
    }

    EnergyFinalDemand::completeInit( aRegionName, aRegionInfo );
}

void BuildingFinalDemand::initCalc( const string& aRegionName,
                                    const GDP* aGDP,
                                    const Demographic* aDemographics,
                                    const int aPeriod )
{
    EnergyFinalDemand::initCalc( aRegionName, aGDP, aDemographics, aPeriod );
}

/*!
* \brief Calculate the macro-economic scaler for the service demand.
* \details
*
* \return The macro-economic scaler.
* \todo This method is not currently different from the aggregate per capita 
*       demand function.  Reserved for potential change.
*/
double BuildingFinalDemand::calcMacroScaler( const string& aRegionName,
                                             const Demographic* aDemographics,
                                             const GDP* aGDP,
                                             const int aPeriod ) const
{
    if( aPeriod == 0 ){
        // No changes in price, income and population scales.
        return 1;
    }

    int previousPeriod = 0;
    if( aPeriod > 0 ){
        previousPeriod = aPeriod - 1;
    }

    double priceRatio = SectorUtils::calcPriceRatio( aRegionName, mName,
                                                     previousPeriod, aPeriod );

    double GDPperCapRatio = aGDP->getGDPperCap( aPeriod )
                          / aGDP->getGDPperCap( aPeriod - 1);

    double populationRatio = aDemographics->getTotal( aPeriod )
                           / aDemographics->getTotal( aPeriod - 1);

    // TODO: A null saturation elasticity results in the reduction of
    // the unscaled service demand by 1/2 and must be accounted for 
    // in the base scaler.  For non-calibrated periods, the base scaler is
    // 1 and the saturation formulation with null elasticity will reduce the
    // demand by a half.  Reformulate the saturation expression so that it is
    // a multiple of 1 instead.
    // SHK 7/11/07
    double macroScaler = pow( priceRatio, mPriceElasticity[ aPeriod ] )
                         * pow( GDPperCapRatio, mIncomeElasticity[ aPeriod ] )
                       //  / ( 1 + pow( aGDP->getGDPperCap( aPeriod ) / mSaturationPoint, mSaturationElasticity.get() ) )
                         * populationRatio;

    return macroScaler;
}
