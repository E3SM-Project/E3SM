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
* \file model_time.cpp
* \ingroup Objects
* \brief modeltime class source file.
* \author Sonny Kim
*/

#include "util/base/include/definitions.h"
#include <iostream>
#include <fstream>

#include <string>
#include <cassert>

// xml headers
#include <xercesc/dom/DOMNodeList.hpp>
#include "util/base/include/xml_helper.h"
#include "util/base/include/model_time.h"
#include "util/logger/include/ilogger.h"

using namespace std;
using namespace xercesc;

const Modeltime* Modeltime::getInstance() {
    const static Modeltime modeltime;
    return &modeltime;
}

//! Default constructor.
Modeltime::Modeltime()
:mIsInitialized( false ),
mStartYear( -1 ),
mEndYear( -1 ),
mFinalCalibrationYear( 2005 )
{
}

//! Set the data members from the XML input.
bool Modeltime::XMLParse( const DOMNode* aNode ) {
    
    // It is important that Modeltime only be parse once as the very first
    // object since the rest of the model will rely on it to initialize some
    // datastructures.
    if( mIsInitialized ) {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::ERROR );
        mainLog << "Modeltime can only be parsed once." << endl;
        return false;
    }
    
    // Temparary sorted map of years to their time step to keep track of
    // them while we are parsing.  After parsing we will use them to initialize
    // the Modeltime with the initMembers method.
    map<int, int> yearToTimeStep;
        
    // assume node is valid.
    assert( aNode );

    // get all children of the node.
    DOMNodeList* nodeList = aNode->getChildNodes();

    // loop through the children
    for ( unsigned int i = 0; i < nodeList->getLength(); ++i ){
        DOMNode* curr = nodeList->item( i );
        const string nodeName = XMLHelper<string>::safeTranscode( curr->getNodeName() );

        // select the type of node.
        if( nodeName == "#text" ) {
            continue;
        }

        else if ( nodeName == "start-year" ){
            mStartYear = XMLHelper<int>::getValue( curr );
            yearToTimeStep[ mStartYear ]  = XMLHelper<int>::getAttr( curr, "time-step" );
        } 
        else if ( nodeName == "inter-year" ){
            int interYear = XMLHelper<int>::getValue( curr );
            yearToTimeStep[ interYear ]  = XMLHelper<int>::getAttr( curr, "time-step" );
        } 
        else if ( nodeName == "end-year" ){
            mEndYear = XMLHelper<int>::getValue( curr );
            // the end year does not have a time step
            yearToTimeStep[ mEndYear ]  = -1;
        }
        else if ( nodeName == "final-calibration-year" ){
            int tempCalibrationYear = XMLHelper<int>::getValue( curr ); 
            // mFinalCalibrationYear is initialized to 2005
            if ( tempCalibrationYear != mFinalCalibrationYear ){
                ILogger& mainLog = ILogger::getLogger( "main_log" );
                mainLog.setLevel( ILogger::WARNING );
                mainLog << "Using read in final-calibration-year (" << tempCalibrationYear
                    << ") and not the last historical year (" << mFinalCalibrationYear << ")." << endl;
                mFinalCalibrationYear = tempCalibrationYear;
            }
        } 
        else {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            mainLog << "Unrecognized text string: " << nodeName << " found while parsing modeltime." << endl;
        }
    }
    
    // initialize the data members and finalize the Modeltime
    initMembers( yearToTimeStep );
    return true;
}

//! Write data members to datastream in XML format.
void Modeltime::toInputXML( ostream& aOut, Tabs* aTabs ) const {
    
    XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs );

    // note that all redundant inter-year will not be written back out even if they were read in
    map<string, int> attrs;
    for( int period = 0; period < mMaxPeriod; ++period ) {
        // skip writting out unnecessary inter-year elements
        if( period != 0 && period < mMaxPeriod - 1 && mPeriodToTimeStep[ period + 1 ] == mPeriodToTimeStep[ period ] ) {
            continue;
        }
        attrs[ "time-step" ] = period != mMaxPeriod - 1 ? mPeriodToTimeStep[ period + 1 ] : 0;
        if( period == 0 ) {
            XMLWriteElementWithAttributes( mStartYear, "start-year", aOut, aTabs, attrs );
        }
        else if( period == mMaxPeriod - 1 ) {
            XMLWriteElement( mEndYear, "end-year", aOut, aTabs );
        }
        else {
            XMLWriteElementWithAttributes( mPeriodToYear[ period ], "inter-year", aOut, aTabs, attrs );
        }
    }
    
    XMLWriteElement( mFinalCalibrationYear, "final-calibration-year", aOut, aTabs );
    
    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

//! Write out object to output stream for debugging.
void Modeltime::toDebugXML( const int aPeriod, ostream& aOut, Tabs* aTabs ) const {
    
    XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs );

    // write all info for the given period to debug
    map<string, int> attrs;
    attrs[ "time-step" ] = mPeriodToTimeStep[ aPeriod ];
    if( aPeriod == 0 ) {
        XMLWriteElementWithAttributes( mStartYear, "start-year", aOut, aTabs, attrs );
    }
    else if( aPeriod == mMaxPeriod - 1 ) {
        XMLWriteElement( mEndYear, "end-year", aOut, aTabs );
    }
    else {
        XMLWriteElementWithAttributes( mPeriodToYear[ aPeriod ], "inter-year", aOut, aTabs, attrs );
    }
    XMLWriteElement( mFinalCalibrationYear, "final-calibration-year", aOut, aTabs );

    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
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
const std::string& Modeltime::getXMLNameStatic() {
    const static string XML_NAME = "modeltime";
    return XML_NAME;
}

/*!
 * \brief Initialize all modeltime parameters based on years and the year increment
 *        from that year to break up the span of years in between.
 * \details This method can only be called one time and at the very beginning of the
 *          model (such as right after parsing modeltime data but before parsing
 *          data for anything else).  Error checking will take place and is not forgiving
 *          such that any error will leave the modeltime uninitialized.
 * \param aYearToTimeStep A sorted by year map that relates a year to a year interval
 *                        which is used to break up each span of years in the map.
 */
void Modeltime::initMembers( const map<int, int>& aYearToTimeStep ) {
    // start with some error checking
    ILogger& mainLog = ILogger::getLogger( "main_log" );
    mainLog.setLevel( ILogger::ERROR );
    if( mIsInitialized ) {
        mainLog << getXMLNameStatic() << " can only be initialized once." << endl;
        exit( 1 );
    }
    if( aYearToTimeStep.size() == 0 ) {
        mainLog << "No year information has been parsed in " << getXMLNameStatic() << endl;
        exit( 1 );
    }
    if( mStartYear == -1 ) {
        mainLog << "No start year information has been parsed in " << getXMLNameStatic() << endl;
        exit( 1 );
    }
    if( mEndYear == -1 ) {
        mainLog << "No end year information has been parsed in " << getXMLNameStatic() << endl;
        exit( 1 );
    }
    if( aYearToTimeStep.begin()->first != mStartYear ) {
        mainLog << "Parsed invalid year: " << aYearToTimeStep.begin()->first << " which is before the start year: "
            << mStartYear << endl;
        exit( 1 );
    }
    if( ( --aYearToTimeStep.end() )->first != mEndYear ) {
        mainLog << "Parsed invalid year: " << ( --aYearToTimeStep.end() )->first << " which is after the end year: "
            << mEndYear << endl;
        exit( 1 );
    }
    if( mFinalCalibrationYear < mStartYear || mFinalCalibrationYear > mEndYear ) {
        mainLog << "Parsed invalid final calibration year: " << mFinalCalibrationYear << endl;
        exit( 1 );
    }
    
    // start processing not to say no more errors are possible however
    mMaxPeriod = 0;
    // the timesteps are shifted by one and so the time step in period 0 does not make sense
    // note that the 15 is arbitrary here but a valid timestep is necessary for period 0
    mPeriodToTimeStep.push_back( 15 );
    for( map<int, int>::const_iterator yearIt = aYearToTimeStep.begin(); yearIt != --aYearToTimeStep.end(); ) {
        int currYear = yearIt->first;
        int currTimeStep = yearIt->second;
        int nextYear = ( ++yearIt )->first;
        if( ( ( nextYear - currYear ) % currTimeStep ) != 0 ) {
            mainLog << "Specified time step of " << currTimeStep << " does not evenly divide years from "
                << currYear << " to " << nextYear << endl;
            exit( 1 );
        }
        
        int numInbetweenPeriods = ( nextYear - currYear ) / currTimeStep;
        for( int periodOffset = 0; periodOffset < numInbetweenPeriods; ++periodOffset ) {
            int offsetYear = currYear + periodOffset * currTimeStep;
            mPeriodToYear.push_back( offsetYear );
            mPeriodToTimeStep.push_back( currTimeStep );
            mYearToPeriod[ offsetYear ] = mMaxPeriod++;
        }
    }
    
    // add info for the end year
    mPeriodToYear.push_back( mEndYear );
    mYearToPeriod[ mEndYear ] = mMaxPeriod++;
    
    // Fill non-model years in 1 year timesteps into the year to period map with the period of the next
    // model year. Required for the carbon box model.
    for( int currPeriod = 1; currPeriod < mMaxPeriod; ++currPeriod ) {
        for( int inBetweenYear = mPeriodToYear[ currPeriod - 1 ] + 1; inBetweenYear < mPeriodToYear[ currPeriod ]; ++inBetweenYear ) {
            mYearToPeriod[ inBetweenYear ] = currPeriod;
        }
    }
    
    // finally set the flag that we are initialized
    mIsInitialized = true;
}

//! Get the base period
int Modeltime::getBasePeriod() const {
    assert( mIsInitialized );
    
    return getyr_to_per( mStartYear );
}

//! Get the start year.
int Modeltime::getStartYear() const {
    assert( mIsInitialized );
    
    return mStartYear;
}

//! Get the end year.
int Modeltime::getEndYear() const {
    assert( mIsInitialized );
    
    return mEndYear;
}

/*!
* \brief Convert a period into a year.
* \details Converts the period into a year if it is valid. If it is not a valid
*          year the function will print a warning and return 0.
* \param aPeriod Model period.
* \return The first year of the period, 0 if the year is invalid.
*/
int Modeltime::getper_to_yr( const int aPeriod ) const {
    assert( mIsInitialized );
    
    if( aPeriod >= 0 && aPeriod < static_cast<int>( mPeriodToYear.size() ) ){
        return mPeriodToYear[ aPeriod ];
    }

    ILogger& mainLog = ILogger::getLogger( "main_log" );
    mainLog.setLevel( ILogger::ERROR );
    mainLog << "Invalid period " << aPeriod << " passed to Modeltime::getper_to_yr." << endl;
    return 0;
}

//! Convert a year to a period.
int Modeltime::getyr_to_per( const int aYear ) const {
    assert( mIsInitialized );
    
    map<int,int>::const_iterator iter = mYearToPeriod.find( aYear );
    // Check for an invalid time period.
    if( iter == mYearToPeriod.end() ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::ERROR );
        mainLog << "Invalid year: " << aYear << " passed to Modeltime::getyr_to_per. " << endl;
        return 0;
    }
    return iter->second; 
}

/*!
 * \brief Determine if the given year is a model year.
 * \param aYear The year to check.
 * \return True if aYear is a model year, false otherwise.
 */
bool Modeltime::isModelYear( const int aYear ) const {
    assert( mIsInitialized );
    
    // Note that simply checking the year to period map will not work as intended.
    // This is due to values being filled into the mYearToPeriod map for in between
    // years.  A work around is to convert the year to a period and back again to
    // make sure that year is the same as aYear.
    map<int, int>::const_iterator lookupYearIter = mYearToPeriod.find( aYear );
    return lookupYearIter != mYearToPeriod.end()
        && mPeriodToYear[ (*lookupYearIter).second ] == aYear;
}

/*!
 * \brief Get the final period in which base year calibration will occur.
 * \return Final period in which base year calibration will occur.
 */
int Modeltime::getFinalCalibrationPeriod() const {
    assert( mIsInitialized );
    
    return getyr_to_per( mFinalCalibrationYear );
}

