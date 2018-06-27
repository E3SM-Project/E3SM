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
 * All rights to use the Software are granted on condition that such
 * rights are forfeited if User fails to comply with the terms of
 * this Agreement.
 * 
 * User agrees to identify, defend and hold harmless BATTELLE,
 * its officers, agents and employees from all liability involving
 * the violation of such Export Laws, either directly or indirectly,
 * by User.
 */

/*! 
 * \file rcp_emissions_visitor.h
 * \ingroup Objects
 * \brief The RCPEmissionsVisitor class source file.
 *
 * \author Pralit Patel
 */

#include "util/base/include/definitions.h"

#include "containers/include/scenario.h"
#include "util/base/include/model_time.h"
#include "util/logger/include/ilogger.h"

#include "reporting/include/rcp_emissions_visitor.h"
#include "containers/include/region.h"
#include "resources/include/aresource.h"
#include "sectors/include/sector.h"
#include "sectors/include/subsector.h"
#include "technologies/include/technology.h"
#include "emissions/include/aghg.h"

extern Scenario* scenario; // for modeltime

using namespace std;
using namespace objects;

//! Default Constructor
RCPEmissionsVisitor::RCPEmissionsVisitor()
{
    // Create a mapping to deal with the fact that GCAM may use multiple names
    // to refer to the same gas.
    mGHGMap[ "CH4" ] = "CH4";
    mGHGMap[ "CH4_AWB" ] = "CH4";
    mGHGMap[ "CH4_AGR" ] = "CH4";
    mGHGMap[ "N20" ] = "N20";
    mGHGMap[ "N20_AWB" ] = "N20";
    mGHGMap[ "N20_AGR" ] = "N20";
    mGHGMap[ "BC" ] = "BC";
    mGHGMap[ "BC_AWB" ] = "BC";
    mGHGMap[ "BC_AGR" ] = "BC";
    mGHGMap[ "OC" ] = "OC";
    mGHGMap[ "OC_AWB" ] = "OC";
    mGHGMap[ "OC_AGR" ] = "OC";
    mGHGMap[ "CO" ] = "CO";
    mGHGMap[ "CO_AWB" ] = "CO";
    mGHGMap[ "CO_AGR" ] = "CO";
    mGHGMap[ "NOx" ] = "NOx";
    mGHGMap[ "NOx_AWB" ] = "NOx";
    mGHGMap[ "NOx_AGR" ] = "NOx";
    mGHGMap[ "NMVOC" ] = "NMVOC";
    mGHGMap[ "NMVOC_AWB" ] = "NMVOC";
    mGHGMap[ "NMVOC_AGR" ] = "NMVOC";
    mGHGMap[ "NH3" ] = "NH3";
    mGHGMap[ "NH3_AWB" ] = "NH3";
    mGHGMap[ "NH3_AGR" ] = "NH3";
    mGHGMap[ "SO2_1" ] = "SO2";
    mGHGMap[ "SO2_1_AWB" ] = "SO2";
    mGHGMap[ "SO2_2" ] = "SO2";
    mGHGMap[ "SO2_2_AWB" ] = "SO2";
    mGHGMap[ "SO2_3" ] = "SO2";
    mGHGMap[ "SO2_3_AWB" ] = "SO2";
    mGHGMap[ "SO2_4" ] = "SO2";
    mGHGMap[ "SO2_4_AWB" ] = "SO2";
    
    // Define the GCAM to RCP category mappings which may need to cut sectors,
    // subsectors, or technologies in different ways.
    // TODO: a better way of doing this?
    // Note:  Commented out mappings are relevent in later versions of GCAM and are
    // left in for reference/future use.
    RCPMapping tempMapping ( "ENE" );
    tempMapping.mSectorName = "regional crude oil";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mSectorName = "regional natural gas";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mSectorName = "regional unconventional oil";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mSectorName = "regional coal";
    mRCPMappingsList.push_back( tempMapping );
    /*tempMapping.mSectorName = "regional oil";
    mRCPMappingsList.push_back( tempMapping );*/
    tempMapping.mSectorName = "electricity";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mSectorName = "elect_td_bld";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mSectorName = "elect_td_ind";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mSectorName = "elect_td_trn";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mSectorName = "base load generation";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mSectorName = "intermediate generation";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mSectorName = "subpeak generation";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mSectorName = "peak generation";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mSectorName = "H2 Central Production";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mSectorName = "H2 Forecourt Production";
    mRCPMappingsList.push_back( tempMapping );
    /*tempMapping.mSectorName = "district heat";
    mRCPMappingsList.push_back( tempMapping );*/
    tempMapping.mSectorName = "refined liquids enduse";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mSectorName = "refined liquids industrial";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mSectorName = "refined liquids electricity"; // Removed in future versions
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mSectorName = "crude oil";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mSectorName = "coal";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mSectorName = "natural gas";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mSectorName = "unconventional oil";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mSectorName = "gas processing";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mSectorName = "regional biomass";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mSectorName = "regional sugar for ethanol";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mSectorName = "regional biomassOil";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mSectorName = "regional corn for ethanol";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mSectorName = "regional sugarbeet for ethanol";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mSectorName = "backup_electricity";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mSectorName = "csp_backup";
    mRCPMappingsList.push_back( tempMapping );
    /*tempMapping.mSectorName = "unconventional oil production";
    mRCPMappingsList.push_back( tempMapping );*/
    
    tempMapping = RCPMapping( "IND" );
    tempMapping.mSectorName = "industrial energy use";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mSectorName = "industrial processes";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mSectorName = "process heat cement";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mSectorName = "industrial feedstocks";
    mRCPMappingsList.push_back( tempMapping );
    
    tempMapping = RCPMapping( "DOM" );
    tempMapping.mSectorName = "building";
    mRCPMappingsList.push_back( tempMapping );
    /*tempMapping.mSectorName = "comm heating";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mSectorName = "comm cooling";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mSectorName = "comm others";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mSectorName = "resid heating";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mSectorName = "resid cooling";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mSectorName = "resid others";
    mRCPMappingsList.push_back( tempMapping );*/
    
    tempMapping = RCPMapping( "SHIP" );
    tempMapping.mSectorName = "trn_shipping_intl";
    mRCPMappingsList.push_back( tempMapping );
    
    tempMapping = RCPMapping( "AIR" );
    tempMapping.mSectorName = "trn_passenger";
    tempMapping.mSubsectorName = "air";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mSectorName = "trn_freight";
    mRCPMappingsList.push_back( tempMapping );
    
    tempMapping = RCPMapping( "TRA" );
    tempMapping.mSectorName = "trn_passenger";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mSectorName = "trn_freight";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mSectorName = "trn_pass_road";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mSectorName = "gas pipeline";
    mRCPMappingsList.push_back( tempMapping );
    // TODO: gas pipeline?
    
    tempMapping = RCPMapping( "WST" );
    tempMapping.mSectorName = "urban processes";
    mRCPMappingsList.push_back( tempMapping );
    
    tempMapping = RCPMapping( "SLV" );
    tempMapping.mSectorName = "industrial processes";
    tempMapping.mSubsectorName = "solvents";
    mRCPMappingsList.push_back( tempMapping );
    
    tempMapping = RCPMapping( "AGR" );
    tempMapping.mGHGName = "CH4_AGR";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mGHGName = "N20_AGR";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mGHGName = "BC_AGR";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mGHGName = "OC_AGR";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mGHGName = "CO_AGR";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mGHGName = "NOx_AGR";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mGHGName = "NMVOC_AGR";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mGHGName = "CH4_AGR";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mGHGName = "NH3_AGR";
    mRCPMappingsList.push_back( tempMapping );
    
    tempMapping = RCPMapping( "AWB" );
    tempMapping.mGHGName = "CH4_AWB";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mGHGName = "N20_AWB";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mGHGName = "BC_AWB";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mGHGName = "OC_AWB";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mGHGName = "CO_AWB";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mGHGName = "NOx_AWB";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mGHGName = "NMVOC_AWB";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mGHGName = "CH4_AWB";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mGHGName = "NH3_AWB";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mGHGName = "SO2_1_AWB";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mGHGName = "SO2_2_AWB";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mGHGName = "SO2_3_AWB";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mGHGName = "SO2_4_AWB";
    mRCPMappingsList.push_back( tempMapping );
    
    // Note in the future we will have a 0 prefixed for aezs < 10.
    tempMapping = RCPMapping( "LCF" );
    tempMapping.mTechnologyName = "Forest FiresAEZ1";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "Forest FiresAEZ2";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "Forest FiresAEZ3";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "Forest FiresAEZ4";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "Forest FiresAEZ5";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "Forest FiresAEZ6";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "Forest FiresAEZ7";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "Forest FiresAEZ8";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "Forest FiresAEZ9";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "Forest FiresAEZ10";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "Forest FiresAEZ11";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "Forest FiresAEZ12";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "Forest FiresAEZ13";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "Forest FiresAEZ14";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "Forest FiresAEZ15";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "Forest FiresAEZ16";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "Forest FiresAEZ17";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "Forest FiresAEZ18";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "DeforestationAEZ1";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "DeforestationAEZ2";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "DeforestationAEZ3";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "DeforestationAEZ4";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "DeforestationAEZ5";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "DeforestationAEZ6";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "DeforestationAEZ7";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "DeforestationAEZ8";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "DeforestationAEZ9";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "DeforestationAEZ10";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "DeforestationAEZ11";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "DeforestationAEZ12";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "DeforestationAEZ13";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "DeforestationAEZ14";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "DeforestationAEZ15";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "DeforestationAEZ16";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "DeforestationAEZ17";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "DeforestationAEZ18";
    mRCPMappingsList.push_back( tempMapping );
    
    tempMapping = RCPMapping( "SAV" );
    tempMapping.mTechnologyName = "Savannah BurningAEZ1";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "Savannah BurningAEZ2";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "Savannah BurningAEZ3";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "Savannah BurningAEZ4";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "Savannah BurningAEZ5";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "Savannah BurningAEZ6";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "Savannah BurningAEZ7";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "Savannah BurningAEZ8";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "Savannah BurningAEZ9";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "Savannah BurningAEZ10";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "Savannah BurningAEZ11";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "Savannah BurningAEZ12";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "Savannah BurningAEZ13";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "Savannah BurningAEZ14";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "Savannah BurningAEZ15";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "Savannah BurningAEZ16";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "Savannah BurningAEZ17";
    mRCPMappingsList.push_back( tempMapping );
    tempMapping.mTechnologyName = "Savannah BurningAEZ18";
    mRCPMappingsList.push_back( tempMapping );
    
    
}

/*!
 * \brief Get RCP emissions by the given region, RCP category, gas, and year.
 * \param aRegionName The GCAM region to get emissions for.
 * \param aSectorCategory The RCP category to get emissions for.
 * \param aGasName The RCP GHG gas name to get emissions for.
 * \param aYear The year in which to get emissions
 * \return The emissions for the requested category or 0 if unavailable.
 */
double RCPEmissionsVisitor::getEmissions( const string& aRegionName,
                                          const string& aSectorCategory,
                                          const string& aGasName,
                                          const int aYear ) const
{
    const Modeltime* modeltime = scenario->getModeltime();
    const int period = modeltime->getyr_to_per( aYear );
    // TODO: interpolate if not a model year?
    if( period != 0 ) {
        // Do a sequence of finds on the nest of maps in mRCPEmissions to get to
        // the requested emissions value.
        map<string, map<string, map<string, PeriodVector<Value> > > >::const_iterator regionIt =
            mRCPEmissions.find( aRegionName );
        if( regionIt != mRCPEmissions.end() ) {
            typedef map<string, map<string, PeriodVector<Value> > > CategoryMap;
            const CategoryMap& categoryMap = (*regionIt).second;
            CategoryMap::const_iterator categoryIt = categoryMap.find( aSectorCategory );
            if( categoryIt != categoryMap.end() ) {
                typedef map<string, PeriodVector<Value> > GasMap;
                const GasMap& gasMap = (*categoryIt).second;
                GasMap::const_iterator gasIt = gasMap.find( aGasName );
                if( gasIt != gasMap.end() ) {
                    return (*gasIt).second[ period ];
                }
            }
        }
    }
    return 0;
}

void RCPEmissionsVisitor::startVisitRegion( const Region* aRegion, const int aPeriod ) {
    mCurrentRegion = aRegion->getName();
}

void RCPEmissionsVisitor::startVisitResource( const AResource* aResource, const int aPeriod ) {
    vector<const RCPMapping*> resourceFilter;
    for( vector<RCPMapping>::const_iterator mapIt = mRCPMappingsList.begin(); mapIt != mRCPMappingsList.end(); ++mapIt ) {
        if( (*mapIt).mSectorName == aResource->getName() ) {
            resourceFilter.insert( resourceFilter.begin(), &(*mapIt) );
        }
        else if( (*mapIt).mSectorName.empty() ) {
            resourceFilter.push_back( &(*mapIt) );
        }
    }
    mCurrentMappings.push( resourceFilter );
}

void RCPEmissionsVisitor::endVisitResource( const AResource* aResource, const int aPeriod ) {
    mCurrentMappings.pop();
}

void RCPEmissionsVisitor::startVisitSector( const Sector* aSector, const int aPeriod ) {
    vector<const RCPMapping*> sectorFilter;
    for( vector<RCPMapping>::const_iterator mapIt = mRCPMappingsList.begin(); mapIt != mRCPMappingsList.end(); ++mapIt ) {
        if( (*mapIt).mSectorName == aSector->getName() ) {
            sectorFilter.insert( sectorFilter.begin(), &(*mapIt) );
        }
        else if( (*mapIt).mSectorName.empty() ) {
            sectorFilter.push_back( &(*mapIt) );
        }
    }
    mCurrentMappings.push( sectorFilter );
}

void RCPEmissionsVisitor::endVisitSector( const Sector* aSector, const int aPeriod ) {
    mCurrentMappings.pop();
}

void RCPEmissionsVisitor::startVisitSubsector( const Subsector* aSubsector, const int aPeriod ) {
    vector<const RCPMapping*> subsectorFilter;
    const vector<const RCPMapping*>& currFilter = mCurrentMappings.top();
    for( vector<const RCPMapping*>::const_iterator mapIt = currFilter.begin(); mapIt != currFilter.end(); ++mapIt ) {
        if( (*mapIt)->mSubsectorName == aSubsector->getName() ) {
            subsectorFilter.insert( subsectorFilter.begin(), *mapIt );
        }
        else if( (*mapIt)->mSubsectorName.empty() ) {
            subsectorFilter.push_back( *mapIt );
        }
    }
    mCurrentMappings.push( subsectorFilter );
}

void RCPEmissionsVisitor::endVisitSubsector( const Subsector* aSubsector, const int aPeriod ) {
    mCurrentMappings.pop();
}

void RCPEmissionsVisitor::startVisitTechnology( const Technology* aTechnology, const int aPeriod ) {
    vector<const RCPMapping*> technologyFilter;
    const vector<const RCPMapping*>& currFilter = mCurrentMappings.top();
    for( vector<const RCPMapping*>::const_iterator mapIt = currFilter.begin(); mapIt != currFilter.end(); ++mapIt ) {
        if( (*mapIt)->mTechnologyName == aTechnology->getName() ) {
            technologyFilter.insert( technologyFilter.begin(), *mapIt );
        }
        else if( (*mapIt)->mTechnologyName.empty() ) {
            technologyFilter.push_back( *mapIt );
        }
    }
    mCurrentMappings.push( technologyFilter );
}

void RCPEmissionsVisitor::endVisitTechnology( const Technology* aTechnology, const int aPeriod ) {
    mCurrentMappings.pop();
}

void RCPEmissionsVisitor::startVisitGHG( const AGHG* aGHG, const int aPeriod ) {
    // GCAM uses multiple names to refer to the same GHG so we may need to translate
    // it.  Note the mappings assume the GCAM GHG name.
    const string rcpGHGName = getRCPGHGName( aGHG->getName() );

    const vector<const RCPMapping*>& currFilter = mCurrentMappings.top();    
    const RCPMapping* ghgFilter = 0;
    // TODO: do some error checking
    for( vector<const RCPMapping*>::const_iterator mapIt = currFilter.begin(); mapIt != currFilter.end(); ++mapIt ) {
        if( (*mapIt)->mGHGName == aGHG->getName() ) {
            ghgFilter = *mapIt;
        }
        else if( (*mapIt)->mGHGName.empty() && !ghgFilter) {
            ghgFilter = *mapIt;
        }
    }
    
    if( ghgFilter ) {
        mRCPEmissions[ mCurrentRegion ][ ghgFilter->mRCPName ][ rcpGHGName ][ aPeriod ] += aGHG->getEmission( aPeriod );
    }
}

/*!
 * \brief Retrieve the RCP GHG name given the GCAM GHG name.
 * \details This will do a check in the GHG map to see if GCAM has an alternate name
 *          for the gas.  If so that alternate name is returned, otherwise the GCAM
 *          name will be used.
 * \param aGasName The GCAM name for a gas which may need to be translated.
 * \return The appropriate GHG name to use.
 */
const string& RCPEmissionsVisitor::getRCPGHGName( const string& aGasName ) const {
    map<string, string>::const_iterator nameIt = mGHGMap.find( aGasName );
    return nameIt == mGHGMap.end() ? aGasName : (*nameIt).second;
}
