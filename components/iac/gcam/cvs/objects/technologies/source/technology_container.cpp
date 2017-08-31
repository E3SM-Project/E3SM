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
 * \file technology_container.cpp
 * \ingroup Objects
 * \brief TechnologyContainer class source file.
 * \author Pralit Patel
 */
#include "util/base/include/definitions.h"
#include <string>
#include <cassert>
#include <xercesc/dom/DOMNodeList.hpp>

#include "util/base/include/util.h"
#include "technologies/include/technology_container.h"
#include "util/base/include/xml_helper.h"
#include "util/base/include/interpolation_rule.h"
#include "technologies/include/stub_technology_container.h"

// technology includes
#include "technologies/include/technology.h"
#include "technologies/include/itechnology.h"
#include "technologies/include/default_technology.h"
#include "technologies/include/building_generic_dmd_technology.h"
#include "technologies/include/intermittent_technology.h"
#include "technologies/include/wind_technology.h"
#include "technologies/include/solar_technology.h"
#include "technologies/include/nuke_fuel_technology.h"
#include "technologies/include/tran_technology.h"
#include "technologies/include/ag_production_technology.h"
#include "technologies/include/unmanaged_land_technology.h"
#include "technologies/include/empty_technology.h"

extern Scenario* scenario;

using namespace std;
using namespace xercesc;

//! Constructor
TechnologyContainer::TechnologyContainer()
:mInitialAvailableYear( -1 ),
mFinalAvailableYear( -1 )
{
}

//! Destructor
TechnologyContainer::~TechnologyContainer() {
    for( CVintageIterator vintageIt = mVintages.begin(); vintageIt != mVintages.end(); ++vintageIt ) {
        delete ( *vintageIt ).second;
    }
    mVintages.clear();
    
    // just in case null out the period vector as well
    for( int period = 0; period < mVintagesByPeriod.size(); ++period ) {
        mVintagesByPeriod[ period ] = 0;
    }
    
    clearInterpolationRules();
}


/*!
 * \brief Create a deep clone of this technology container.
 * \note that initializations such as in completeInit may still need
 *          to be done.
 * \return A pointer to the cloned technology container.
 */
ITechnologyContainer* TechnologyContainer::clone() const {
    TechnologyContainer* clonedTechContainer = new TechnologyContainer();
    
    clonedTechContainer->mName = mName;
    clonedTechContainer->mInitialAvailableYear = mInitialAvailableYear;
    clonedTechContainer->mFinalAvailableYear = mFinalAvailableYear;
    
    for( CInterpRuleIterator ruleIter = mShareWeightInterpRules.begin(); ruleIter != mShareWeightInterpRules.end(); ++ruleIter ) {
        clonedTechContainer->mShareWeightInterpRules.push_back( new InterpolationRule( *( *ruleIter ) ) );
    }
    
    for( CVintageIterator vintageIter = mVintages.begin(); vintageIter != mVintages.end(); ++vintageIter ) {
        clonedTechContainer->mVintages[ ( *vintageIter ).first ] = ( *vintageIter ).second->clone();
    }
    
    return clonedTechContainer;
}

/*!
 * \brief Factory method to determine which technologies this container knows about.
 * \param aTechNodeName The xml name to check to see if it is a known technology
 *          type.
 * \return True if the given type is a known technology, false otherwise.
 * \note The list of known technologies here needs to be kept in sync with
 *       the ones found in TechnologyContainer::createAndParseVintage.
 */
bool TechnologyContainer::hasTechnologyType( const string& aTechNodeName ) {
    return ( aTechNodeName == DefaultTechnology::getXMLNameStatic() ||
             aTechNodeName == BuildingGenericDmdTechnology::getXMLNameStatic() ||
             aTechNodeName == IntermittentTechnology::getXMLNameStatic() ||
             aTechNodeName == WindTechnology::getXMLNameStatic() ||
             aTechNodeName == SolarTechnology::getXMLNameStatic() ||
             aTechNodeName == NukeFuelTechnology::getXMLNameStatic() ||
             aTechNodeName == TranTechnology::getXMLNameStatic() ||
             aTechNodeName == AgProductionTechnology::getXMLNameStatic() ||
             aTechNodeName == UnmanagedLandTechnology::getXMLNameStatic() );
}

/*!
 * \brief Create and parse a vintage of a technology from the given XML and of the given
 *        type.
 * \details Creates a new vintage in the year determined by the year attribute of
 *          aNode if one does not already exist for that year.  That vintage will
 *          then have XMLParse called on it.
 * \param aNode The XML which defines the vintage including year.
 * \param aTechType The type of technology which would need to be created.
 * \return Whether the creation and parsing of the vintage was successful.
 * \note The types of technologies that may be created should be kept in sync with
 *       TechnologyContainer::hasTechnologyType
 */
bool TechnologyContainer::createAndParseVintage( const DOMNode* aNode, const string& aTechType ) {
    /*! \pre Tech type should be known. */
    assert( hasTechnologyType( aTechType ) );
    
    ITechnology* newVintage;
    
    const int techYear = XMLHelper<int>::getAttr( aNode, "year" );
    if( techYear == 0 ) {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::ERROR );
        mainLog << "Could not determine year for technology " << mName << " while parsing " << aTechType << endl;
        return false;
    }
    
    // check if this vintage has already been created
    VintageIterator vintageIt = mVintages.find( techYear );
    if( vintageIt == mVintages.end() ) {
        // has not been added yet so we must create a new one
        if( aTechType == DefaultTechnology::getXMLNameStatic() ) {
            newVintage = new DefaultTechnology( mName, techYear );
        }
        else if( aTechType == BuildingGenericDmdTechnology::getXMLNameStatic() ) {
            newVintage = new BuildingGenericDmdTechnology( mName, techYear );
        }
        else if( aTechType == IntermittentTechnology::getXMLNameStatic() ) {
            newVintage = new IntermittentTechnology( mName, techYear );
        }
        else if( aTechType == WindTechnology::getXMLNameStatic() ) {
            newVintage = new WindTechnology( mName, techYear );
        }
        else if( aTechType == SolarTechnology::getXMLNameStatic() ) {
            newVintage = new SolarTechnology( mName, techYear );
        }
        else if( aTechType == NukeFuelTechnology::getXMLNameStatic() ) {
            newVintage = new NukeFuelTechnology( mName, techYear );
        }
        else if( aTechType == TranTechnology::getXMLNameStatic() ) {
            newVintage = new TranTechnology( mName, techYear );
        }
        else if( aTechType == AgProductionTechnology::getXMLNameStatic() ) {
            newVintage = new AgProductionTechnology( mName, techYear );
        }
        else if( aTechType == UnmanagedLandTechnology::getXMLNameStatic() ) {
            newVintage = new UnmanagedLandTechnology( mName, techYear );
        }
        else {
            // Getting an error message here implies that the known technologies in this method are
            // out of sync with hasTechnologyType.
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::ERROR );
            mainLog << "Could not determine a technology of type" << aTechType << endl;
            return false;
        }
        
        /*!
         * \note That we adding the vintage even if there are errors while parsing it.
         */
        mVintages[ techYear ] = newVintage;
    }
    else {
        // just fetch the previously created vintage
        newVintage = ( *vintageIt ).second;
    }
    
    return newVintage->XMLParse( aNode );
}

bool TechnologyContainer::XMLParse( const DOMNode* aNode ) {
    /*! \pre Make sure we were passed a valid node. */
    assert( aNode );
    
    // get the technology type
    string techType = XMLHelper<string>::safeTranscode( aNode->getNodeName() );
    
    // special case to accomodate overriding parameters from global technologies
    if( techType == StubTechnologyContainer::getXMLNameStatic() ) {
        /*!
         * \pre There must be atleast on vintage from the global technology.
         */
        assert( mVintages.begin() != mVintages.end() );
        techType = (*mVintages.begin()).second->getXMLName();
    }
    
    // get the name attribute.
    mName = XMLHelper<string>::getAttr( aNode, XMLHelper<void>::name() );
    
    if( mName.empty() ) {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::ERROR );
        mainLog << "Ignoring technology set because it does not have a name." << endl;
        return false;
    }
    
    // flag for return if parsing was successful
    // TODO: wait do we carry this through?
    bool parsingSuccessful = true;
    
    // get all child nodes.
    DOMNodeList* nodeList = aNode->getChildNodes();
    
    // loop through the child nodes.
    for( unsigned int i = 0; i < nodeList->getLength(); i++ ){
        DOMNode* curr = nodeList->item( i );
        string nodeName = XMLHelper<string>::safeTranscode( curr->getNodeName() );
        
        if( nodeName == XMLHelper<void>::text() ) {
            continue;
        }
        else if( nodeName == Technology::getXMLVintageNameStatic() ) {
            parsingSuccessful &= createAndParseVintage( curr, techType );
        }
        else if( nodeName == InterpolationRule::getXMLNameStatic() && XMLHelper<string>::getAttr( curr, "apply-to" ) == "share-weight" ) {
            // if the delete flag is set then for interpolation rules that means to clear
            // out any previously parsed rules
            if( XMLHelper<bool>::getAttr( curr, "delete" ) ) {
                clearInterpolationRules();
            }
            
            InterpolationRule* tempRule = new InterpolationRule();
            tempRule->XMLParse( curr );
            mShareWeightInterpRules.push_back( tempRule );
        }
        else if( nodeName == "initial-available-year" ) {
            mInitialAvailableYear = XMLHelper<int>::getValue( curr );
        }
        else if( nodeName == "final-available-year" ) {
            mFinalAvailableYear = XMLHelper<int>::getValue( curr );
        }
        else {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::ERROR );
            mainLog << "Unknown element " << nodeName << " encountered while parsing " << techType << endl;
            parsingSuccessful = false;
        }
    }
    
    return parsingSuccessful;
}

void TechnologyContainer::toInputXML( ostream& aOut, Tabs* aTabs ) const {
    // Note that we are assuming there is atleast 1 vintage and that the tech type
    // is the same for all vintages.
    const string techType = ( *mVintages.begin() ).second->getXMLName();
    XMLWriteOpeningTag( techType, aOut, aTabs, mName );
    
    XMLWriteElementCheckDefault( mInitialAvailableYear, "initial-available-year", aOut, aTabs, -1 );
    XMLWriteElementCheckDefault( mFinalAvailableYear, "final-available-year", aOut, aTabs, -1 );
    
    for( CInterpRuleIterator ruleIt = mShareWeightInterpRules.begin(); ruleIt != mShareWeightInterpRules.end(); ++ruleIt ) {
        ( *ruleIt )->toInputXML( aOut, aTabs );
    }
    
    for( CVintageIterator vintageIt = mVintages.begin(); vintageIt != mVintages.end(); ++vintageIt ) {
        // only write technologies that were not interpolated in toInputXML
        if( find( mInterpolatedTechYears.begin(), mInterpolatedTechYears.end(), ( *vintageIt ).first )
            == mInterpolatedTechYears.end() || techType == IntermittentTechnology::getXMLNameStatic() )
        {
            ( *vintageIt ).second->toInputXML( aOut, aTabs );
        }
    }
    
    XMLWriteClosingTag( techType, aOut, aTabs );
}

void TechnologyContainer::toDebugXML( const int aPeriod, ostream& aOut, Tabs* aTabs ) const {
    for( CVintageIterator vintageIt = mVintages.begin(); vintageIt != mVintages.end(); ++vintageIt ) {
        ( *vintageIt ).second->toDebugXML( aPeriod, aOut, aTabs );
    }
}

const string& TechnologyContainer::getName() const {
    return mName;
}

void TechnologyContainer::completeInit( const string& aRegionName,
                                        const string& aSectorName,
                                        const string& aSubsectorName,
                                        DependencyFinder* aDependencyFinder,
                                        const IInfo* aSubsecInfo,
                                        ILandAllocator* aLandAllocator )
{
    // Setup the period to vintage vector with the parsed technologies.  Also check
    // for technologies that are not in the initial year to final year range if
    // specified.
    const Modeltime* modeltime = scenario->getModeltime();
    for( VintageIterator vintageIt = mVintages.begin(); vintageIt != mVintages.end(); ) {
        if( mInitialAvailableYear != -1 && ( *vintageIt ).first < mInitialAvailableYear ) {
            // warn a technology exists before the intial year
            /*ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::DEBUG );
            mainLog << "Removing technology " << mName << " vintage " << ( *vintageIt ).first
                << " since it is before the intial year " << mInitialAvailableYear << endl;
            delete ( *vintageIt ).second;
            
            // The erase will invalidate the iterator so we must use the post-increment
            // operator to get the next value before erasing the iterator.
            mVintages.erase( vintageIt++ );*/
            // We can get away with not deleting previous technologies and they will
            // be necessary if the mInitialAvailableYear needs to be interpolated.
            vintageIt++;
        }
        else if( mFinalAvailableYear != -1 && ( *vintageIt ).first > mFinalAvailableYear ) {
            // warn a technology exists after the final investment year
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::DEBUG );
            mainLog << "Removing technology " << mName << " vintage " << ( *vintageIt ).first
                << " since it is after the final investment year " << mFinalAvailableYear << endl;
            delete ( *vintageIt ).second;
            
            // The erase will invalidate the iterator so we must use the post-increment
            // operator to get the next value before erasing the iterator.
            mVintages.erase( vintageIt++ );
        }
        else {
            // Add technologies that are on model years to the vintages by period vector
            // for easy lookup.
            if( modeltime->isModelYear( ( *vintageIt ).first ) ) {
                mVintagesByPeriod[ modeltime->getyr_to_per( ( *vintageIt ).first ) ] = ( *vintageIt ).second;
            }
            
            // did not erase so increment the iterator
            ++vintageIt;
        }
    }
    
    // Finish filling the vintages by period vector with empty technologies or
    // interpolate for missing years.
    for( int period = 0; period < mVintagesByPeriod.size(); ++period ) {
        const int year = modeltime->getper_to_yr( period );
        if( ( mInitialAvailableYear != -1 && year < mInitialAvailableYear ) ||
            ( mFinalAvailableYear != -1 && year > mFinalAvailableYear ) ) {
            mVintagesByPeriod[ period ] = EmptyTechnology::getInstance();
        }
        else if( !mVintagesByPeriod[ period ] ) {
            if( period == 0 ) {
                ILogger& mainLog = ILogger::getLogger( "main_log" );
                mainLog.setLevel( ILogger::ERROR );
                mainLog << "Missing a necessary technology in the base year." << endl;
                exit( 1 );
            }
            
            // We can interpolate a technology to fill this year.  Note that this
            // interpolated technology will not be written out in toInputXML
            CVintageIterator tempPrevTech, prevTech;
            tempPrevTech = prevTech = mVintages./*find*/lower_bound( 
                                      modeltime->getper_to_yr( period - 1 ) );
            // Use temporary iterator to get next tech to avoid changing prev tech.
            CVintageIterator nextTech = ++tempPrevTech;
            interpolateVintage( year, prevTech, nextTech );
        }
    }
    
    // Now that all interpolated technologies have been created we can call
    // completeInit.
    for( VintageIterator vintageIt = mVintages.begin(); vintageIt != mVintages.end(); ++vintageIt ) {
        // call complete init for all vintages, even those that are not a model year
        ( *vintageIt ).second->completeInit( aRegionName, aSectorName, aSubsectorName, aDependencyFinder,
                                             aSubsecInfo, aLandAllocator );
    }
}

void TechnologyContainer::initCalc( const string& aRegionName, const string& aSectorName,
                                    const IInfo* aSubsecInfo, const Demographic* aDemographic,
                                    const int aPeriod )
{
    
    // Pass forward any emissions information
    if( aPeriod > 1 && mVintagesByPeriod[ aPeriod - 1 ] != EmptyTechnology::getInstance() 
        && mVintagesByPeriod[ aPeriod ] != EmptyTechnology::getInstance() )
    {
        vector<string> ghgNames = mVintagesByPeriod[ aPeriod -1]->getGHGNames();
        int numberOfGHGs = mVintagesByPeriod[ aPeriod -1]->getNumbGHGs();
        for ( int j=0 ; j<numberOfGHGs; j++ ) {
            mVintagesByPeriod[ aPeriod ]->copyGHGParameters(
                mVintagesByPeriod[ aPeriod - 1 ]->getGHGPointer( ghgNames[j] ) );
        }
    }
    
    // Initialize the previous period info as having no input set and
    // cumulative Hicks neutral, energy of 1 and it is the first tech.
    PreviousPeriodInfo prevPeriodInfo = { 0, 1, true };
    // Warning: aPeriod is the current model period and not the technology vintage.
    // Currently calls initCalc on all vintages past and future.
    // TODO: Should not call initialization for all future technology vintages beyond the
    // current period but correction causing error (SHK).
    for( VintageIterator vintageIt = mVintages.begin(); vintageIt != mVintages.end(); ++vintageIt ) {
        ( *vintageIt ).second->initCalc( aRegionName, aSectorName, aSubsecInfo, aDemographic,
                                         prevPeriodInfo, aPeriod );
        prevPeriodInfo.mIsFirstTech = false;
    }
    
    // interpolate technology share weights
    interpolateShareWeights( aPeriod );
    
}

void TechnologyContainer::postCalc( const string& aRegionName, const int aPeriod ) {
    for( VintageIterator vintageIt = mVintages.begin(); vintageIt != mVintages.end(); ++vintageIt ) {
        ( *vintageIt ).second->postCalc( aRegionName, aPeriod );
    }
}

ITechnology* TechnologyContainer::getNewVintageTechnology( const int aPeriod ) {
    return mVintagesByPeriod[ aPeriod ];
}

const ITechnology* TechnologyContainer::getNewVintageTechnology( const int aPeriod ) const {
    return mVintagesByPeriod[ aPeriod ];
}

ITechnologyContainer::TechRangeIterator TechnologyContainer::getVintageBegin( const int aPeriod ) {
    const int year = scenario->getModeltime()->getper_to_yr( aPeriod );
    
    // Lower bound will give us the first technology which is not < year so we must
    // then make sure that it is not > year and is not beyond the final investment
    // year.  In those cases decrease the iterator to make sure we don't go include
    // a technology beyond aPeriod.
    VintageIterator vintageIter = mVintages.lower_bound( year );
    if( vintageIter == mVintages.end() || ( *vintageIter ).first > year ) {
        --vintageIter;
    }
    
    // Converting a forward iterator to a reverse in not completely inuitive.  We
    // need to decrease one from the converted forward iterator to be in the same place.
    return --TechRangeIterator( vintageIter );
}

ITechnologyContainer::CTechRangeIterator TechnologyContainer::getVintageBegin( const int aPeriod ) const {
    const int year = scenario->getModeltime()->getper_to_yr( aPeriod );
    
    // Lower bound will give us the first technology which is not < year so we must
    // then make sure that it is not > year and is not beyond the final investment
    // year.  In those cases decrease the iterator to make sure we don't go include
    // a technology beyond aPeriod.
    CVintageIterator vintageIter = mVintages.lower_bound( year );
    if( vintageIter == mVintages.end() || ( *vintageIter ).first > year ) {
        --vintageIter;
    }
    
    // Converting a forward iterator to a reverse in not completely inuitive.  We
    // need to decrease one from the converted forward iterator to be in the same place.
    return --CTechRangeIterator( vintageIter );
}

ITechnologyContainer::TechRangeIterator TechnologyContainer::getVintageEnd() {
    return mVintages.rend();
}

ITechnologyContainer::CTechRangeIterator TechnologyContainer::getVintageEnd() const {
    return mVintages.rend();
}

/*!
 * \brief A helper method to delete and clear the vector of interpolation rules.
 */
void TechnologyContainer::clearInterpolationRules() {
    for( CInterpRuleIterator ruleIter = mShareWeightInterpRules.begin(); ruleIter != mShareWeightInterpRules.end(); ++ruleIter ) {
        delete *ruleIter;
    }
    mShareWeightInterpRules.clear();
}

/*!
 * \brief Apply interpolation rules to fill share-weights for technologies.
 * \details Rules will only apply in the first period after calibration.  For
 *          technologies that were interpolated for time-steps that are left
 *          with uninitialized values after rules will have them linearly
 *          interpolated.  Any other uninitialized share-weights will be an error.
 * \param aPeriod The current model period.
 */
void TechnologyContainer::interpolateShareWeights( const int aPeriod ) {
    const Modeltime* modeltime = scenario->getModeltime();
    if( aPeriod != ( modeltime->getFinalCalibrationPeriod() + 1 ) ) {
        return;
    }
    
    objects::PeriodVector<Value> techShareWeights;
    objects::PeriodVector<Value> techParsedShareWeights;
    for( int period = 0; period < mVintagesByPeriod.size(); ++period ) {
        techParsedShareWeights[ period ] = mVintagesByPeriod[ period ]->getParsedShareWeight();
        if( period <= modeltime->getFinalCalibrationPeriod() ) {
            // be sure to get the calibrated share-weight
            techShareWeights[ period ] = mVintagesByPeriod[ period ]->getShareWeight();
        }
        else {
            // initialize the share-weights from the parsed ones
            techShareWeights[ period ] = techParsedShareWeights[ period ];
        }
    }
    for( CInterpRuleIterator ruleIter = mShareWeightInterpRules.begin(); ruleIter != mShareWeightInterpRules.end(); ++ruleIter ) {
        ( *ruleIter )->applyInterpolations( techShareWeights, techParsedShareWeights );
    }
    for( int period = 0; period < mVintagesByPeriod.size(); ++period ) {
        // Interpolated technologies may still have uninitialized share-weights,
        // particularly those between the last calibration year and the next non
        // interpolated technology.  We make the default behavior here to linearly
        // interpolate those values.
        const int year = modeltime->getper_to_yr( period );
        vector<int>::const_iterator interpYearsIt = find( mInterpolatedTechYears.begin(), mInterpolatedTechYears.end(), year );
        vector<int>::const_iterator interpYearsEnd = mInterpolatedTechYears.end();
        if( !techShareWeights[ period ].isInited() && interpYearsIt != interpYearsEnd ) {
            const int prevYear = modeltime->getper_to_yr( period - 1 );
            
            // Find the next non-interpolated technology to interpolate to.
            int nextPeriod = period;
            int nextYear = year;
            while( interpYearsIt != mInterpolatedTechYears.end() && nextPeriod < modeltime->getmaxper() ) {
                ++nextPeriod;
                // guard against going past the last model year
                if( nextPeriod < modeltime->getmaxper() ) {
                    nextYear = modeltime->getper_to_yr( nextPeriod );
                    interpYearsIt = find( interpYearsIt, interpYearsEnd, nextYear );
                }
            }
            if( nextPeriod == modeltime->getmaxper() ) {
                // If there was no next technology just copy the previous share-weight
                // forward.
                techShareWeights[ period ].set( techShareWeights[ period - 1] );
            }
            else {
                // linearly interpolate to fill this missing value
                techShareWeights[ period ].set( util::linearInterpolateY( year,
                    prevYear, nextYear, techShareWeights[ period - 1], techShareWeights[ nextPeriod ] ) );
            }
        }
        // All periods must have set a share weight value at this point, not having one is an error.
        if( period > modeltime->getFinalCalibrationPeriod() && !techShareWeights[ period ].isInited() ) {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::ERROR );
            mainLog << "Found uninitialized share weight in tech: " << mName
                << " in period " << aPeriod << endl;
            exit( 1 );
        }
        mVintagesByPeriod[ period ]->setShareWeight( techShareWeights[ period ] );
    }
}

/*!
 * \brief Interpolate a technology vintage given the previous and next vintage.
 * \details The vintage will be interpolated by cloning aPrevTech.  If aNextTech
 *          exists then doInterpolations will be called on the new vintage.  All
 *          datastructures will be updated with the new vintage.
 * \param aYear The year which needs to be interpolated.
 * \param aPrevTech The previous vintage to interpolate from.
 * \param aNextTech The next vintage to interpolate to.
 */
void TechnologyContainer::interpolateVintage( const int aYear, CVintageIterator aPrevTech,
                                              CVintageIterator aNextTech )
{
    /*!
     * \pre The given year does not already exist.
     */
    assert( mVintages.find( aYear ) == mVintages.end() );
    
    /*!
     * \pre The previous vintage must exist.
     */
    assert( aPrevTech != mVintages.end() );
    
    ITechnology* newTech = ( *aPrevTech ).second->clone();
    newTech->setYear( aYear );
    
    if( aNextTech != mVintages.end() ) {
        newTech->doInterpolations( static_cast<Technology*>( ( *aPrevTech ).second ),
                                   static_cast<Technology*>( ( *aNextTech ).second ) );
    }
    mVintages[ aYear ] = newTech;
    mInterpolatedTechYears.push_back( aYear );
    
    const Modeltime* modeltime = scenario->getModeltime();
    if( modeltime->isModelYear( aYear ) ) {
        mVintagesByPeriod[ modeltime->getyr_to_per( aYear ) ] = newTech;
    }
}

void TechnologyContainer::accept( IVisitor* aVisitor, const int aPeriod ) const {
    const Modeltime* modeltime = scenario->getModeltime();
    CVintageIterator end;
    int periodToUse;
    // If the period is -1 this means to update output containers for all periods.
    if( aPeriod == -1 ){
        end = mVintages.end();
        periodToUse = modeltime->getmaxper() - 1;
    }
    else {
        int lastTechYearToVisit = modeltime->getper_to_yr( aPeriod );
        // If the technology should not exist yet for this period then there is
        // nothing to visit.
        if( lastTechYearToVisit < mInitialAvailableYear ) {
            return;
        }
        
        // If the technology was no longer available for investment in the current
        // period then just visit to the last available period.
        if( mFinalAvailableYear != -1 && lastTechYearToVisit > mFinalAvailableYear ) {
            end = mVintages.end();
        }
        else {
            end = mVintages.find( modeltime->getper_to_yr( aPeriod ) );
            
            // aPeriod must correspond to a valid technology
            assert( end != mVintages.end() );
            
            // Increase end by one so that it is pointing to the first technology
            // that we don't want to visit.
            ++end;
        }
        periodToUse = aPeriod;
    }
    
    for( CVintageIterator vintageIt = mVintages.begin(); vintageIt != end; ++vintageIt ) {
        ( *vintageIt ).second->accept( aVisitor, periodToUse );
    }
}

void TechnologyContainer::interpolateAndParse( const DOMNode* aNode ) {
    // Get all child nodes.
    DOMNodeList* nodeList = aNode->getChildNodes();
    
    // Loop through the child nodes and interpolate any vintages which do not exist.
    for( unsigned int i = 0; i < nodeList->getLength(); i++ ){
        DOMNode* curr = nodeList->item( i );
        string nodeName = XMLHelper<string>::safeTranscode( curr->getNodeName() );
        
        if( nodeName == XMLHelper<void>::text() ) {
            continue;
        }
        
        const int year = XMLHelper<int>::getAttr( curr, "year" );
        if( year != 0 && mVintages.find( year ) == mVintages.end() ) {
            // Find the previous and next vintages so that we can interpolate.
            CVintageIterator prevTech = --mVintages.lower_bound( year );
            if( prevTech == mVintages.end() ) {
                ILogger& mainLog = ILogger::getLogger( "main_log" );
                mainLog.setLevel( ILogger::ERROR );
                mainLog << "Could not find a vintage before " << year << " to interplate from." << endl;
                exit( 1 );
            }
            CVintageIterator tempPrevTech = prevTech;
            CVintageIterator nextTech = ++tempPrevTech;
            interpolateVintage( year, prevTech, nextTech );
        }
    }
    
    // We can now parse values.
    XMLParse( aNode );
}
