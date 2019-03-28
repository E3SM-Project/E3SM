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
* \file get_glm_data.h
* \ingroup Objects
* \brief The GetGLMData class source file.
*
* \author Pralit Patel
* \author Sonny Kim
*/

#include <boost/lexical_cast.hpp>

#include "util/base/include/definitions.h"

#include "reporting/include/get_glm_data.h"
#include "containers/include/region.h"
#include "land_allocator/include/land_leaf.h"
#include "land_allocator/include/aland_allocator_item.h"
#include "technologies/include/ag_production_technology.h"
#include "ccarbon_model/include/icarbon_calc.h"
#include "containers/include/scenario.h"
#include "util/base/include/model_time.h"

using namespace std;
using namespace objects;

extern Scenario* scenario;

//! Default Constructor
GetGLMData::GetGLMData()
{
     //TODO: Hard coding renames here for the moment.  These should be read-in
     // through data or otherwise set exogenously to allow felxability in changes
     // to regional or land cover categorization.
     // mRegionMap[ GCAM Region ] = GLM Region
     mRegionMap[ "USA" ] = "USA";
     mRegionMap[ "Canada" ] = "Canada";
     mRegionMap[ "Western Europe" ] = "Western Europe";
     mRegionMap[ "Japan" ] = "Japan";
     mRegionMap[ "Australia_NZ" ] = "Australia_NZ";
     mRegionMap[ "Former Soviet Union" ] = "Former Soviet Union";
     mRegionMap[ "China" ] = "China";
     mRegionMap[ "Middle East" ] = "Middle East";
     mRegionMap[ "Africa" ] = "Africa";
     mRegionMap[ "Latin America" ] = "Latin America";
     mRegionMap[ "Southeast Asia" ] = "Southeast Asia";
     mRegionMap[ "Eastern Europe" ] = "Eastern Europe";
     mRegionMap[ "Korea" ] = "Korea";
     mRegionMap[ "India" ] = "India";

     mLandMap[ "biomass" ] = "Cropland";
     mLandMap[ "eucalyptus" ] = "Cropland";
     mLandMap[ "miscanthus" ] = "Cropland";
     mLandMap[ "willow" ] = "Cropland";
     mLandMap[ "Jatropha" ] = "Cropland";

     mLandMap[ "Corn" ] = "Cropland";
     mLandMap[ "FiberCrop" ] = "Cropland";
     mLandMap[ "OtherArableLand" ] = "Cropland";
     mLandMap[ "OtherGrain" ] = "Cropland";
     mLandMap[ "Rice" ] = "Cropland";
     mLandMap[ "SugarCrop" ] = "Cropland";
     mLandMap[ "Wheat" ] = "Cropland";

     mLandMap[ "OilCrop" ] = "Cropland";
     mLandMap[ "PalmFruit" ] = "Cropland";

     mLandMap[ "MiscCrop" ] = "Cropland";
     mLandMap[ "Root_Tuber" ] = "Cropland";

     mLandMap[ "FodderHerb" ] = "Cropland";
     mLandMap[ "FodderGrass" ] = "Cropland";

     mLandMap[ "Forest" ] = "Forest";
     mLandMap[ "UnmanagedForest" ] = "Forest";

     mLandMap[ "Grassland" ] = "Grassland";
     mLandMap[ "Shrubland" ] = "Grassland";
     mLandMap[ "Tundra" ] = "Grassland";

     mLandMap[ "Pasture" ] = "Pasture";
     mLandMap[ "UnmanagedPasture" ] = "Pasture";

     mLandMap[ "UrbanLand" ] = "Build-up";

     mLandMap[ "RockIceDesert" ] = "Other Land";

     // Crop name mapping for crop production
     mCropMap[ "Forest" ] = "Forest";

     // Land name mapping for carbon densities
     // TODO: create these mappings
}

/*!
 * \brief Get the land cover for a given GLM region and land category in the given year.
 * \details If data was requested for and unknown region/AEZ/category/year -1 will be returned.
 * \param aGLMRegionName A region name.
 * \param aAEZ The AEZ number.
 * \param aGLMLandCategory A land cover category.
 * \param aYear A GCAM model year in which to get data.  TODO: do we want to allow interpolating?
 * \return The requested land cover or -1 for error.
 */
double GetGLMData::getLandCover( const string& aGLMRegionName,
                                 const int aAEZ,
                                 const string& aGLMLandCategory,
                                 const int aYear ) const
{
    const Modeltime* modeltime = scenario->getModeltime();
    const int modelPeriod = modeltime->getyr_to_per( aYear );
    if( modelPeriod == 0 ) {
        return -1;
    }
    RegAEZLandMap::const_iterator regionIt = mLandCoverData.find( aGLMRegionName );
    if( regionIt != mLandCoverData.end() ) {
        map<int, map<string, PeriodVector<double> > >::const_iterator aezIt = 
            (*regionIt).second.find( aAEZ );
        if( aezIt != (*regionIt).second.end() ) {
            map<string, PeriodVector<double> >::const_iterator categoryIt = 
                (*aezIt).second.find( aGLMLandCategory );
            if( categoryIt != (*aezIt).second.end() ) {
                return (*categoryIt).second[ modelPeriod ];
            }
        }
    }

    // Did not find the region/category
    return -1;
}

/*!
 * \brief Get the crop production in carbon for a given GLM region and category in the given year.
 * \details If data was requested for and unknown region/AEZ/category/year -1 will be returned.
 * \param aGLMRegionName A region name.
 * \param aAEZ The AEZ number.
 * \param aGLMLandCategory A crop category.
 * \param aYear A GCAM model year in which to get data.  TODO: do we want to allow interpolating?
 * \return The requested crop production or -1 for error.
 * \warning Only Forest is currently available and it's units are in tC.
 */
double GetGLMData::getProductionInCarbon( const string& aGLMRegionName,
                                          const int aAEZ,
                                          const string& aGLMCrop,
                                          const int aYear ) const
{
    const Modeltime* modeltime = scenario->getModeltime();
    const int modelPeriod = modeltime->getyr_to_per( aYear );
    if( modelPeriod == 0 ) {
        return -1;
    }
    RegAEZLandMap::const_iterator regionIt = mProductionData.find( aGLMRegionName );
    if( regionIt != mProductionData.end() ) {
        map<int, map<string, PeriodVector<double> > >::const_iterator aezIt = 
            (*regionIt).second.find( aAEZ );
        if( aezIt != (*regionIt).second.end() ) {
            map<string, PeriodVector<double> >::const_iterator categoryIt = 
                (*aezIt).second.find( aGLMCrop );
            if( categoryIt != (*aezIt).second.end() ) {
                return (*categoryIt).second[ modelPeriod ];
            }
        }
    }

    // Did not find the region/category
    return -1;
}

/*!
 * \brief Get the above and below carbon densities for a given GLM region and land category in the given year.
 * \details If data was requested for and unknown region/AEZ/category/year -1 will be returned.
 * \param aGLMRegionName A region name.
 * \param aAEZ The AEZ number.
 * \param aGLMLandCategory A land cover category.
 * \param aYear A GCAM model year in which to get data.  TODO: do we want to allow interpolating?
 * \return A pair of (above, below) carbon densities.  In case of error both values with be -1.
 */
pair<double, double> GetGLMData::getCarbonDensity( const string& aGLMRegionName,
                                                   const int aAEZ,
                                                   const string& aGLMLandCategory,
                                                   const int aYear ) const
{
    CarbonDensityMap::const_iterator regionIt = mCarbonDensityData.find( aGLMRegionName );
    if( regionIt != mCarbonDensityData.end() ) {
        map<int, map<string, map<int, pair<double, double> > > >::const_iterator aezIt = 
            (*regionIt).second.find( aAEZ );
        if( aezIt != (*regionIt).second.end() ) {
            map<string, map<int, pair<double, double> > >::const_iterator categoryIt = 
                (*aezIt).second.find( aGLMLandCategory );
            if( categoryIt != (*aezIt).second.end() ) {
                map<int, pair<double, double> >::const_iterator yearIt = 
                    (*categoryIt).second.find( aYear );
                if( yearIt != (*categoryIt).second.end() ) {
                    return (*yearIt).second;
                }
            }
        }
    }

    // Did not find the region/category
    return make_pair<double, double>( -1, -1 );
}

void GetGLMData::startVisitRegion( const Region* aRegion, const int aPeriod ) {
    // Store GLM region for later use when we have data to store
    mCurrentGLMRegionName = doMapLookup( mRegionMap, aRegion->getName() );
}

void GetGLMData::startVisitLandLeaf( const LandLeaf* aLandLeaf, const int aPeriod ) {
    const Modeltime* modeltime = scenario->getModeltime();
    // Allow summing of individual periods or all periods (the -1 flag)
    const int startPeriod = aPeriod == -1 ? 0 : aPeriod;
    const int endPeriod = aPeriod == -1 ? modeltime->getmaxper() - 1 : aPeriod;
    pair<string, int> nameAEZ = parseTypeAndAEZFromName( aLandLeaf->getName() );
    PeriodVector<double>& data = mLandCoverData[ mCurrentGLMRegionName ][ nameAEZ.second ][ doMapLookup( mLandMap, nameAEZ.first ) ];
    for( int period = startPeriod; period <= endPeriod; ++period ) {
        data[ period ] += aLandLeaf->getLandAllocation( aLandLeaf->getName(), period );
    }

    // TODO: behavior if different GCAM land types that are mapped to the same GLM land type
    // have different carbon densities this will likely not give the expeceted result.
    map<int, pair<double, double> >& carbonDensities = mCarbonDensityData[ mCurrentGLMRegionName ][ nameAEZ.second ][ doMapLookup( mLandMapForDensities, nameAEZ.first ) ];
    const int startYear = modeltime->getper_to_yr( startPeriod );
    const int endYear = modeltime->getper_to_yr( endPeriod );
    for( int year = startYear; year <= endYear; ++year ) {
        pair<double, double>& density = carbonDensities[ year ];
        density.first = aLandLeaf->getCarbonContentCalc()->getActualAboveGroundCarbonDensity( year );
        density.second = aLandLeaf->getCarbonContentCalc()->getActualBelowGroundCarbonDensity( year );
    }
}

void GetGLMData::startVisitAgProductionTechnology( const AgProductionTechnology* aAgTech, const int aPeriod ) {
    pair<string, int> nameAEZ = parseTypeAndAEZFromName( aAgTech->getName() );
    if( nameAEZ.first == "Forest" ) {
        // Note the hard coding of the forest to carbon conversion factor
        // The output of the forest technology is in billion m^3 and we want
        // the production in tC so this conversion factor is tC / billion m^3
        const double B_CUBIC_METER_FOREST_TO_CARBON = 0.288 * 1e9;

        const Modeltime* modeltime = scenario->getModeltime();
        // Technologies have special behavior when visiting such that aPeriod is really
        // the end period and we must loop from 0.
        PeriodVector<double>& data = mProductionData[ mCurrentGLMRegionName ][ nameAEZ.second ][ doMapLookup( mCropMap, nameAEZ.first ) ];
        for( int period = 0; period <= aPeriod; ++period ) {
            data[ period ] +=
                aAgTech->getOutput( period ) * B_CUBIC_METER_FOREST_TO_CARBON;
        }
    }
}

/*!
 * \brief A helper method to do a rename lookup from a GCAM name to a GLM name.
 * \details If an explicit rename was not provided the GCAM name will be used.
 * \param The name lookup map.
 * \param The name to translate.
 * \return The GLM name if a lookup exists otherwise just the GCAM name.
 */
const string& GetGLMData::doMapLookup( const map<string, string>& aLookupMap, const string& aKey ) {
    map<string, string>::const_iterator lookupIt = aLookupMap.find( aKey );
    if( lookupIt != aLookupMap.end() ) {
        return (*lookupIt).second;
    }
    else {
        // TODO: warn?
        return aKey;
    }
}

pair<string, int> GetGLMData::parseTypeAndAEZFromName( const string& aGCAMName ) {
    size_t index = aGCAMName.find( "AEZ" );
    if( index == string::npos ) {
        // TODO: warn?
        return make_pair<string, int>( aGCAMName, 0 );
    }
    else {
        string typeName( aGCAMName, 0, index);
        string aezStr( aGCAMName, index + 3, string::npos );
        int aez = boost::lexical_cast<int>( aezStr );
        return make_pair<string, int>( typeName, aez );
    }
}

