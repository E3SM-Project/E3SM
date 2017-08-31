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
* \file set_carbon_density.h
* \ingroup Objects
* \brief The SetCarbonDensity class source file.
*
* \author Pralit Patel
* \author Sonny Kim
*/

#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>

#include <boost/lexical_cast.hpp>

#include "util/base/include/definitions.h"

#include "util/base/include/xml_helper.h"
#include "util/logger/include/ilogger.h"
#include "util/base/include/util.h"

#include "land_allocator/include/set_carbon_density.h"
#include "containers/include/region.h"
#include "land_allocator/include/land_leaf.h"
#include "ccarbon_model/include/icarbon_calc.h"
#include "technologies/include/ag_production_technology.h"

using namespace std;
using namespace xercesc;

//! Default Constructor
SetCarbonDensity::SetCarbonDensity()
{
     //TODO: Hard coding renames here for the moment.  These should be read-in
     // through data or otherwise set exogenously to allow felxability in changes
     // to regional or land cover categorization.
     // mRegionMap[ GLM Region ] = GCAM Region
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

     /*
     vector<string> gcamLandTypes;
     gcamLandTypes.push_back( "Forest" );
     gcamLandTypes.push_back( "UnmanagedForest" );
     mLandMap[ "Forest" ] = gcamLandTypes;
     gcamLandTypes.clear();

     gcamLandTypes.push_back( "Corn" );
     mLandMap[ "Corn" ] = gcamLandTypes;
     */
}

bool SetCarbonDensity::XMLParse( const DOMNode* aNode ) {
    // assume we were passed a valid node.
    assert( aNode );

    // get the first child the node.
    DOMNode* curr = aNode->getFirstChild();
 
    // loop through the children
    while( curr ) {
        string nodeName = XMLHelper<string>::safeTranscode( curr->getNodeName() );

        if( nodeName == "#text" ) {
            curr = curr->getNextSibling();
            continue;
        }
        else if( nodeName == "data" ) {
            string regionName = XMLHelper<string>::getAttr( curr, "region" );
            assert( !regionName.empty() );

            int aez = XMLHelper<int>::getAttr( curr, "AEZ" );
            assert( aez != 0 );

            string landCategory = XMLHelper<string>::getAttr( curr, "crop" );
            assert( !regionName.empty() );

            double aboveGroundAdj = XMLHelper<double>::getAttr( curr, "above-adjust" );
            assert( aboveGroundAdj != 0 );

            double belowGroundAdj = XMLHelper<double>::getAttr( curr, "below-adjust" );
            assert( belowGroundAdj != 0 );

            int year = XMLHelper<int>::getAttr( curr, "year" );
            assert( year != 0 );

            // push values
            setCarbonDensityToPush( regionName, aez, landCategory, aboveGroundAdj, belowGroundAdj, year );
        }
        else {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            mainLog << "Unrecognized text string: " << nodeName << " found while parsing "
                << "set-carbon-density" << "." << endl;
        }
        
        // Get the next child of aNode to process.
        curr = curr->getNextSibling();
    }
    return true;
}

void SetCarbonDensity::toInputXML( ostream& aOut, Tabs* aTabs ) const {
    XMLWriteOpeningTag( "set-carbon-density", aOut, aTabs );

    for( RegionCarbonMap::const_iterator regionIter = mCarbonDensityAdjToSet.begin(); regionIter != mCarbonDensityAdjToSet.end(); ++regionIter ) {
        const string currRegion = (*regionIter).first;
        const LandTypeCarbonMap& landCarbonMap = (*regionIter).second;
        for( LandTypeCarbonMap::const_iterator landIter = landCarbonMap.begin(); landIter != landCarbonMap.end(); ++landIter ) {
            const string fullName = (*landIter).first;
            int aezIndex = fullName.find( "AEZ" );
            const string currCrop = fullName.substr( 0, aezIndex );
            const string currAEZ = fullName.substr( aezIndex + 3 );
            const YearCarbonMap& yearMap = (*landIter).second;
            for( YearCarbonMap::const_iterator yearIter = yearMap.begin(); yearIter != yearMap.end(); ++yearIter ) {
                aTabs->writeTabs( aOut );
                aOut << "<data region=\"" << currRegion << "\" crop=\"" << currCrop << "\" AEZ=\"" << currAEZ
                     << "\" year=\"" << (*yearIter).first << "\" above-adjust=\"" << (*yearIter).second.first
                     << "\" below-adjust=\"" << (*yearIter).second.second << "\" />" << endl;
            }
        }
    }

    XMLWriteClosingTag( "set-carbon-density", aOut, aTabs );
}

/*!
 * \brief Interpolate yearly between any read in carbon density adjustments that
 *        are set to push.
 * \details This allows us to only set data in 15 year time-steps.  This is only
 *          necessary for experiment 1 of iESM.
 */
void SetCarbonDensity::doInterpolations() {
    for( RegionCarbonMap::iterator regionIter = mCarbonDensityAdjToSet.begin(); regionIter != mCarbonDensityAdjToSet.end(); ++regionIter ) {
        LandTypeCarbonMap& landCarbonMap = (*regionIter).second;
        for( LandTypeCarbonMap::iterator landIter = landCarbonMap.begin(); landIter != landCarbonMap.end(); ++landIter ) {
            YearCarbonMap& yearMap = (*landIter).second;
            for( YearCarbonMap::iterator yearIter = yearMap.begin(); yearIter != yearMap.end(); ) {
                const int prevYear = (*yearIter).first;
                const pair<double, double> prevYearCarbon = (*yearIter).second;
                ++yearIter;
                if( yearIter != yearMap.end() ) {
                    const int nextYear = (*yearIter).first;
                    const pair<double, double> nextYearCarbon = (*yearIter).second;
                    for( int year = prevYear + 1; year < nextYear; ++year ) {
                        pair<double, double> currCarbon;
                        currCarbon.first = util::linearInterpolateY( year, prevYear, nextYear,
                                prevYearCarbon.first, nextYearCarbon.first );
                        currCarbon.second = util::linearInterpolateY( year, prevYear, nextYear,
                                prevYearCarbon.second, nextYearCarbon.second );
                        yearIter = yearMap.insert( make_pair( year, currCarbon ) ).first;
                    }
                }
            }
        }
    }
}

/*!
 * \brief Stores new carbon densities adjustments to push into GCAM when accept is called
 *        this visitor.
 * \details We asssume that we are given GLM region and land cover names and will
 *          map to the appropriate GCAM names.  Note that for land cover type
 *          there may be a one to many relationship from GLM names to GCAM in which
 *          case the carbon densities adjustment will be applied to all the related GCAM names.
 * \param aRegionName A region name.
 * \param aAEZ The AEZ number.
 * \param aLandCategory A land cover category.
 * \param aCarbonDensityAboveAdj The above ground carbon density adjustment to set.
 * \param aCarbonDensityBelowAdj The below ground carbon density adjustment to set.
 * \param aYear A GCAM model year in which to set data.
 */
void SetCarbonDensity::setCarbonDensityToPush( const string& aRegionName,
                                               const int aAEZ,
                                               const string& aLandCategory,
                                               const double aCarbonDensityAboveAdj,
                                               const double aCarbonDensityBelowAdj,
                                               const int aYear )
{
    map<string, string>::const_iterator it = mRegionMap.find( aRegionName );
    string regionName = it == mRegionMap.end() ? aRegionName : (*it).second;

    vector<string>& gcamLandTypes = mLandMap[ aLandCategory ];
    if( gcamLandTypes.empty() ) {
        // If there was no mapping set for the given land category just use the given
        // category directly.
        mCarbonDensityAdjToSet[ regionName ][ convertToGCAMAEZName( aLandCategory, aAEZ ) ][ aYear ] =
            make_pair( aCarbonDensityAboveAdj, aCarbonDensityBelowAdj );
    }
    else {
        // Set the data for each mapped land type.
        for( vector<string>::const_iterator landIter = gcamLandTypes.begin(); landIter != gcamLandTypes.end(); ++landIter ) {
            mCarbonDensityAdjToSet[ regionName ][ convertToGCAMAEZName( *landIter, aAEZ ) ][ aYear ] =
                make_pair( aCarbonDensityAboveAdj, aCarbonDensityBelowAdj );
        }
    }
}

void SetCarbonDensity::startVisitRegion( const Region* aRegion, const int aPeriod ) {
    // Store region for later use
    mCurrentRegionIter = mCarbonDensityAdjToSet.find( aRegion->getName()) ;
}

void SetCarbonDensity::startVisitLandLeaf( const LandLeaf* aLandLeaf, const int aPeriod ) {
    if( mCurrentRegionIter != mCarbonDensityAdjToSet.end() ) {
	    mCurrentLandTypeIter = (*mCurrentRegionIter).second.find( aLandLeaf->getName() );
    }
}

void SetCarbonDensity::startVisitAgProductionTechnology( const AgProductionTechnology* aAgTech,
                                                         const int aPeriod )
{
    if( mCurrentRegionIter != mCarbonDensityAdjToSet.end() ) {
        LandTypeCarbonMap::const_iterator landIter = (*mCurrentRegionIter).second.find( aAgTech->getName() );
        if( landIter != (*mCurrentRegionIter).second.end() ) {
            YearCarbonMap::const_iterator yearIt = (*landIter).second.find( aAgTech->year );
            if( yearIt != (*landIter).second.end() ) {
                // we adjust yields the same as above ground carbon densities
                const_cast<AgProductionTechnology*>( aAgTech )->mYieldAdjForDensity =
                    (*yearIt).second.first;
            }
        }
    }
}

void SetCarbonDensity::startVisitCarbonCalc( const ICarbonCalc* aCarbonCalc, const int aPeriod ) {
    if( mCurrentRegionIter != mCarbonDensityAdjToSet.end()
        && mCurrentLandTypeIter != (*mCurrentRegionIter).second.end() )
    {
        for( YearCarbonMap::const_iterator yearIt = (*mCurrentLandTypeIter).second.begin();
             yearIt != (*mCurrentLandTypeIter).second.end(); ++yearIt )
        {
            const double currCarbonDensityAbove = 
                aCarbonCalc->getActualAboveGroundCarbonDensity( (*yearIt).first );
            const double currCarbonDensityBelow =
                aCarbonCalc->getActualBelowGroundCarbonDensity( (*yearIt).first );
            // We really need this as non-const however a non-const visitor does not currently exist
            const_cast<ICarbonCalc*>( aCarbonCalc )->setActualAboveGroundCarbonDensity(
                    currCarbonDensityAbove * (*yearIt).second.first, (*yearIt).first );
            const_cast<ICarbonCalc*>( aCarbonCalc )->setActualBelowGroundCarbonDensity(
                    currCarbonDensityBelow * (*yearIt).second.second, (*yearIt).first );
        }
    }
}

/*!
 * \brief Convert a land type and aez number to the GCAM style naming.
 * \param aLandType The land type name.
 * \param aAEZ The aez number to combine with.
 * \return Combined name which can be used to refernce GCAM objects.
 */
string SetCarbonDensity::convertToGCAMAEZName( const string& aLandType, const int aAEZ ) {
    return aLandType + "AEZ" + boost::lexical_cast<string>( aAEZ );
}

