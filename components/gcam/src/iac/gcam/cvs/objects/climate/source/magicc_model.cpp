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
* \file magicc_model.cpp
* \ingroup Objects
* \brief This file contains the source for the MagiccModel class.
* \author Sonny Kim, Josh Lurz
*/
#include "util/base/include/definitions.h"
#include <iomanip>
#include <fstream>
#include <string>
#include <cassert>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>

#include "climate/include/magicc_model.h"
#include "containers/include/scenario.h"
#include "util/base/include/model_time.h"
#include "util/base/include/configuration.h"
#include "util/logger/include/ilogger.h"
#include "util/base/include/util.h"
#include "util/base/include/xml_helper.h"
#include "util/base/include/ivisitor.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;

#include "climate/include/ObjECTS_MAGICC.h"

// MAGICC 5.3 expects a 2000 year line in gas.emk
const int MagiccModel::GAS_EMK_CRIT_YEAR = 2000;

// Setup the gas name vector.
const string MagiccModel::sInputGasNames[] = { "CO2",
                                               "CO2NetLandUse",
                                               "CH4",
                                               "N2O",
                                               "SOXreg1",
                                               "SOXreg2",
                                               "SOXreg3",
                                               "CF4",
                                               "C2F6",
                                               "HFC125",
                                               "HFC134a",
                                               "HFC143a",
                                               "HFC227ea",
                                               "HFC245ca",
                                               "SF6",
                                               "NOx",
                                               "NMVOCs",
                                               "CO",
                                               "BC",
                                               "OC" };

// Setup the gas units vector.
const string MagiccModel::sInputGasUnits[] = { "(Pg C)",
                                               "(Pg C)",
                                               "(Tg)",
                                               "(Tg N)",
                                               "(Tg S)",
                                               "(Tg S)",
                                               "(Tg S)",
                                               "(kton)",
                                               "(kton)",
                                               "(kton)",
                                               "(kton)",
                                               "(kton)",
                                               "(kton)",
                                               "(kton)",
                                               "(kton)",
                                               "(Mt N)",
                                               "(Mt)",
                                               "(Mt)",
                                               "(Mt)",
                                               "(Mt)" };

/*! \brief Constructor
* \param aModeltime Pointer to the global modeltime object.
*/
MagiccModel::MagiccModel( const Modeltime* aModeltime ):
mModeltime( aModeltime ),
mIsValid( false ),
mGHGInputFileName(""),
mCarbonModelStartYear( 1975 ), // Need to have first model year here, but should be 1990 for MAGICC. FIX.
mSoilTempFeedback( -1 ),
mHumusTempFeedback( -1 ),
mGPPTempFeedback( -1 ),
mClimateSensitivity( -1 ),
mNetDeforestCarbFlux80s( -1 ),
mOceanCarbFlux80s( -1 )
{

}

/*! \brief Complete the initialization of the MagiccModel.
* \details This function first resizes the internal vectors which store
*          emissions by gas and period. It then reads in the default set of data
*          for emissions by gas from an input file. These emissions can be
*          overridden by calls to setEmission.
* \param aScenarioName Name of the current model scenario.
*/
void MagiccModel::completeInit( const string& aScenarioName ){
    // The climate model should not be reinitialized.
    assert( !mIsValid );

    mScenarioName = aScenarioName;
    // Resize to the number of gases based on the gas name vector.
    mEmissionsByGas.resize( getNumInputGases() );

    // Resize all the vectors to the number of data points for MAGICC.
    // Add comment here.
    const int numGasPoints = mModeltime->getmaxper() + getNumAdditionalGasPoints();
    for( unsigned int i = 0; i < mEmissionsByGas.size(); ++i ){
        mEmissionsByGas[ i ].resize( numGasPoints );
    }
	
	// Resize the LUC CO2 vector to the number of years plus the number of additional points
	const int numYears = mModeltime->getEndYear() - mModeltime->getStartYear()
							+ getNumAdditionalGasPoints();
	mLUCEmissionsByYear.resize( numYears );
	
    // Read in the data from the file. These values may later be overridden.
    readFile();
    mIsValid = false;
    
    // Set up correspondence for output gas names
    mOutputGasNameMap[ "total" ] = 0;
    mOutputGasNameMap[ "CO2" ] = 1;
    mOutputGasNameMap[ "CH4" ] = 2;  // Direct CH4 only
    mOutputGasNameMap[ "N2O" ] = 3;
    mOutputGasNameMap[ "C2F6" ] = 4;
    mOutputGasNameMap[ "HCFC125" ] = 5;
    mOutputGasNameMap[ "HCFC134A" ] = 6;
    mOutputGasNameMap[ "HCFC143A" ] = 7;
    mOutputGasNameMap[ "HCFC245fa" ] = 8;
    mOutputGasNameMap[ "SF6" ] = 9;
    mOutputGasNameMap[ "CF4" ] = 10;
    mOutputGasNameMap[ "HFC227ea" ] = 11;
    mOutputGasNameMap[ "OtherHC" ] = 12; // Other halocarbon
    mOutputGasNameMap[ "SO2" ] = 13; // Total SO2
    mOutputGasNameMap[ "DirSO2" ] = 14; // Direct SO2
    mOutputGasNameMap[ "TropO3" ] = 15; 
    mOutputGasNameMap[ "O3_CH4" ] = 16; // Ozone from CH4
    mOutputGasNameMap[ "H2O_CH4" ] = 17; // Water vapor from CH4
    mOutputGasNameMap[ "Montreal" ] = 18; // Montreal Gases
    mOutputGasNameMap[ "O3_CFC" ] = 19; // Ozone from CFCs
    mOutputGasNameMap[ "BioBurn" ] = 20;
    mOutputGasNameMap[ "FOC" ] = 21; // Fossil BC/OC
    mOutputGasNameMap[ "Land" ] = 22; // Land surface albedo
    mOutputGasNameMap[ "MinNOx" ] = 23; // Mineral Dust and Nitrate
    mOutputGasNameMap[ "BC" ] = 24; // Total BC (land and combustion)
    mOutputGasNameMap[ "OC" ] = 25; // Total OC (land and combustion)
    mOutputGasNameMap[ "EXTRA" ] = 26; // Extra forcing
    mOutputGasNameMap[ "RCP" ] = 27; // RCP radiative forcing (total - nitrate, albedo, mineral dust)
    
    overwriteMAGICCParameters( );
}

/*! \brief Overwrite MAGICC default parameters with new values.
*/
void MagiccModel::overwriteMAGICCParameters( ){
    // Override parameters in MAGICC if necessary
    int varIndex = 1;
    SETPARAMETERVALUES( varIndex, mClimateSensitivity );
    varIndex = 2;
    SETPARAMETERVALUES( varIndex, mSoilTempFeedback );
    varIndex = 3;
    SETPARAMETERVALUES( varIndex, mHumusTempFeedback );
    varIndex = 4;
    SETPARAMETERVALUES( varIndex, mGPPTempFeedback );
    varIndex = 5;
    SETPARAMETERVALUES( varIndex, mNetDeforestCarbFlux80s );
    varIndex = 6;
    SETPARAMETERVALUES( varIndex, mOceanCarbFlux80s );
}

//! parse MAGICC xml object
void MagiccModel::XMLParse( const DOMNode* node ){
    // make sure we were passed a valid node.
    assert( node );

    DOMNodeList* nodeList = node->getChildNodes();

    for( unsigned int i = 0; i < nodeList->getLength(); ++i ){
        DOMNode* curr = nodeList->item( i );

        // get the name of the node.
        string nodeName = XMLHelper<string>::safeTranscode( curr->getNodeName() );

        if( nodeName == "#text" ) {
            continue;
        }
        // GHG input file
        else if ( nodeName == "ghgInputFileName" ){
            mGHGInputFileName = XMLHelper<string>::getValue( curr );
        }
        // Climate Sensitivity
        else if ( nodeName == "climateSensitivity" ){
            mClimateSensitivity = XMLHelper<double>::getValue( curr );
        }
        // Soil Feedback Factor. 
        else if ( nodeName == "soilTempFeedback" ){
            mSoilTempFeedback = XMLHelper<double>::getValue( curr );
        }
        // Humus Feedback Factor. 
        else if ( nodeName == "humusTempFeedback" ){
            mHumusTempFeedback = XMLHelper<double>::getValue( curr );
        }
        // GPP Feedback Factor. 
        else if ( nodeName == "GPPTempFeedback" ){
            mGPPTempFeedback = XMLHelper<double>::getValue( curr );
        }
        // 1980s Ocean Uptake 
        else if ( nodeName == "oceanFlux80s" ){
            mOceanCarbFlux80s = XMLHelper<double>::getValue( curr );
        }
        // 1980s net terrestrial Deforestation 
        else if ( nodeName == "deforestFlux80s" ){
            mNetDeforestCarbFlux80s = XMLHelper<double>::getValue( curr );
        }
        else if( nodeName == "carbon-model-start-year" ){
            mCarbonModelStartYear = XMLHelper<int>::getValue( curr );
        }
        else {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            mainLog << "Unrecognized text string: " << nodeName << " found while parsing " << getXMLNameStatic() << endl;
        }
    }
}

//! Writes data members to data stream in XML format.
void MagiccModel::toInputXML( ostream& out, Tabs* tabs ) const {
    XMLWriteOpeningTag( getXMLNameStatic(), out, tabs );

    // Write out mClimateSensitivity
    XMLWriteElementCheckDefault( mClimateSensitivity, "climateSensitivity", out, tabs, -1.0 );

    // Write out Soil Feedback Factor.
    XMLWriteElementCheckDefault( mSoilTempFeedback, "soilTempFeedback", out, tabs, -1.0 );

    // Write out Humus Feedback Factor.
    XMLWriteElementCheckDefault( mHumusTempFeedback, "humusTempFeedback", out, tabs, -1.0 );

    // Write out GPP Feedback Factor.
    XMLWriteElementCheckDefault( mGPPTempFeedback, "GPPTempFeedback", out, tabs, -1.0 );

    // Write out 1980s Ocean Uptake.
    XMLWriteElementCheckDefault( mOceanCarbFlux80s, "oceanFlux80s", out, tabs, -1.0 );

    // Write out 1980s net terrestrial Deforestation.
    XMLWriteElementCheckDefault( mNetDeforestCarbFlux80s, "deforestFlux80s", out, tabs, -1.0 );

    // Write out Carbon model start year.
    XMLWriteElementCheckDefault( mCarbonModelStartYear, "carbon-model-start-year", out, tabs, 1975 );

    XMLWriteClosingTag( getXMLNameStatic(), out, tabs );
}

//! Writes data members to debugging data stream in XML format.
void MagiccModel::toDebugXML( const int period, ostream& out, Tabs* tabs ) const {

    XMLWriteOpeningTag( getXMLNameStatic(), out, tabs );

    // Write out mClimateSensitivity
    XMLWriteElement( mClimateSensitivity, "climateSensitivity", out, tabs );

    // Write out Soil Feedback Factor.
    XMLWriteElement( mSoilTempFeedback, "soilTempFeedback", out, tabs );

    // Write out Humus Feedback Factor.
    XMLWriteElement( mHumusTempFeedback, "humusTempFeedback", out, tabs );

    // Write out GPP Feedback Factor.
    XMLWriteElement( mGPPTempFeedback, "GPPTempFeedback", out, tabs );

    // Write out 1980s Ocean Uptake.
    XMLWriteElement( mOceanCarbFlux80s, "oceanFlux80s", out, tabs );

    // Write out 1980s net terrestrial Deforestation.
    XMLWriteElement( mNetDeforestCarbFlux80s, "deforestFlux80s", out, tabs );

    // Write out Carbon model start year.
    XMLWriteElement( mCarbonModelStartYear, "carbon-model-start-year", out, tabs );
    
    XMLWriteClosingTag( getXMLNameStatic(), out, tabs );
}

/*! \brief Checks that the year is a valid value for this climate model
* \details For the MAGICC model valid years range from 1765 to 2500
* \param aYear year
*/
bool MagiccModel::isValidClimateModelYear( const int aYear ) const {
    const int MIN_YEAR = 1765;
    const int MAX_YEAR = 2500;
    
    if ( aYear < MIN_YEAR || aYear > MAX_YEAR ) {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::ERROR );
        mainLog << "Invalid year " << aYear << " for MAGICC climate model." << endl;
        return false;    
    }
    
    return true;
}

/*! \brief Set the level of emissions for a particular for a period.
* \details This function overrides the read in emissions level for a given gas
*          with a specified value. Emissions are not added but reset each time.
* \param aGasName Name of the gas for which to set emissions.
* \param aPeriod Period in which to set emissions.
* \param aEmission Amount of emissions to set.
* \return Whether the emissions were successfully added.
*/
bool MagiccModel::setEmissions( const string& aGasName, const int aPeriod, const double aEmission ){
    // Lookup the gas index.
    const int gasIndex = getInputGasIndex( aGasName );
       
    // Check for error.
    if( gasIndex == INVALID_GAS_NAME ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::ERROR );
        mainLog << "Invalid gas name " << aGasName << " passed to the MAGICC model wrapper." << endl;
        return false;
    }

    // Check that the period is valid.
    if( aPeriod < 1 || aPeriod >= mModeltime->getmaxper() ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::ERROR );
        mainLog << "Invalid period " << aPeriod << " passed to MAGICC model wrapper." << endl;
        return false;
    }

    // Set the gas level. Base year values are not passed to MAGICC, so remove
    // one from the model period to determine the index.
    mEmissionsByGas[ gasIndex ][ aPeriod - 1 ] = aEmission;

    // Invalidate the output because an input has changed.
    mIsValid = false;

    return true;
}

/*! \brief Set the level of emissions for a particular for a year.
 * \details This function overrides the read in emissions level for a given gas
 *          with a specified value. Emissions are not added but reset each time.
 *          This is only relevant for LUC CO2 emissions.
 * \param aGasName Name of the gas for which to set emissions.
 * \param aYear Year in which to set emissions.
 * \param aEmission Amount of emissions to set.
 * \return Whether the emissions were successfully added.
 */
bool MagiccModel::setLUCEmissions( const string& aGasName, const int aYear, const double aEmission ){
	// Check for error.
    if( aGasName != "CO2NetLandUse" ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::ERROR );
        mainLog << "Invalid gas name " << aGasName << " passed to the MAGICC model." << endl;
        return false;
    }
	
    // Check that the period is valid.
    if( aYear < mModeltime->getStartYear() || aYear > mModeltime->getEndYear() ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::ERROR );
        mainLog << "Invalid year " << aYear << " passed to MAGICC model wrapper." << endl;
        return false;
    }
	
    // Set the gas level. Base year values are not passed to MAGICC, so remove
    // one from the model period to determine the index.
    mLUCEmissionsByYear[ aYear - mModeltime->getStartYear() - 1 ] = aEmission;
	
    // Invalidate the output because an input has changed.
    mIsValid = false;
	
    return true;
}

/*!
 * \brief Get the emissions that are currently set to be sent to MAGICC.
 * \details This method only retrieves the values that are currently set and
 *          provides no indication if they are the latest emissions from GCAM nor
 *          if they were the actual emissions MAGICC ran with.  This DOES NOT
 *          attempt to retrieve emissions from within MAGICC.
 * \param aGasName The name of the gas to get emissions for.
 * \param aYear Any valid year.  Emissions may be linearly interpolated if they
 *              are not stored in the given year.
 * \return The emissions for the given gas and year or -1 if an unknown gas or
 *         invalid year was given.
 */
double MagiccModel::getEmissions( const string& aGasName,
                                  const int aYear ) const
{
    // Check if a valid year is passed in.
    if( !isValidClimateModelYear( aYear ) ){
        return -1;
    }

    // Lookup the gas index.
    const int gasIndex = getInputGasIndex( aGasName );
    
    // Check for error.
    if( gasIndex == INVALID_GAS_NAME ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::ERROR );
        mainLog << "Invalid gas name " << aGasName << " passed to the MAGICC model wrapper." << endl;
        return -1;
    }

    const Modeltime* modeltime = scenario->getModeltime();

    int currPeriod = modeltime->getyr_to_per( aYear );

    // If invalid model period is passed, return -1 to indicate no valid emission data
    if( currPeriod <= 0 ) {
            return -1;
    }

    // Check if we don't need to interpolate. This will also handle the last
    // year of period 0.
    if( aYear == modeltime->getper_to_yr( currPeriod ) ){
        // NOTE: MAGICC indices are shifted by 1 position because
        //       the array does not have data for MiniCAM period 0.
        return mEmissionsByGas[ gasIndex ][ currPeriod - 1 ];
    }

    // Find the year which is the last year of the previous period.
    int x1 = modeltime->getper_to_yr( currPeriod - 1 );

    // Find the year which is the last year of the current period.
    int x2 = modeltime->getper_to_yr( currPeriod );

    // Emissions level in the last year of the current period.
    // NOTE: MAGICC indices are shifted by 1 position because
    //       the array does not have data for MiniCAM period 0.
    double y2 = mEmissionsByGas[ gasIndex ][ currPeriod - 1 ];

    // Emissions level in the last year of the previous period.
    double y1 = mEmissionsByGas[ gasIndex ][ currPeriod - 2 ];

    // Linearly interpolate the emissions for the current year based on the
    // emissions levels in the surrounding periods.
    return static_cast<double>( ( aYear - x1 ) * ( y2 - y1 ) ) / static_cast<double>( ( x2 - x1 ) ) + y1;
}

/*! \brief Run the MAGICC emissions model.
* \details This function will run the MAGICC model with the currently stored
*          emissions levels. It will first extrapolate future points for each
*          gas, set equal to the last period. It then prints the gases to the
*          MAGICC input file and calls MAGICC.
* \return Whether the model ran successfully.
* \todo Fix up the calculation of years past the end of the model.
*/
bool MagiccModel::runModel(){
    // Add on extra periods MAGICC needs. 
    // Loop through the gases and copy forward emissions.
    for( unsigned int gasNumber = 0; gasNumber < getNumInputGases(); ++gasNumber ){
        // Note the final period index is shifted by 1.
        const int finalPeriod = mModeltime->getmaxper() - 2;
        // Fill in data for extra periods
        fill( mEmissionsByGas[ gasNumber ].begin() + finalPeriod,
              mEmissionsByGas[ gasNumber ].end(),
              mEmissionsByGas[ gasNumber ][ finalPeriod ] );
    }
    
    // Open the gas file to write emissions into so that MAGICC can get them by
    // reading this file.
    ofstream gasFile;
    gasFile.open( Configuration::getInstance()->getFile( "climatFileName" ).c_str(), ios::out );
    // The first 5 lines of the gas file is skipped before the data is read by MAGICC.

	const int startYear = mModeltime->getper_to_yr( 1 );
	const int endYear = mModeltime->getEndYear();
	const int finalCalYear = mModeltime->getper_to_yr( mModeltime->getFinalCalibrationPeriod() ); 
	const int numCalPer = mModeltime->getFinalCalibrationPeriod();

    // Calculate the number of points which will be passed to magicc.
	// We are sending MAGICC annual emissions after the final calibration year 
	// (endYear - finalCalYear), emissions for each calibration period (numCalPer),
	// and some years after the final model period (getNumAdditionalGasPoints() )
	int numPoints = endYear - finalCalYear + getNumAdditionalGasPoints() + numCalPer;
	// In addition, MAGICC must have year 2000 as input.  If this isn't a model year,
	// we need to increment the number of points
	if ( !mModeltime->isModelYear(GAS_EMK_CRIT_YEAR) ) {
		numPoints += 1;
	}
	gasFile << numPoints << endl;
	
    /*line 2: Name of the scenario*/
    gasFile << " Scenario " << mScenarioName << endl;
    /*line 3: Blank line*/
    gasFile << endl;
    /*line 4: Gas names includes one extra space for commas*/
    gasFile << setw(5) << " ,"; // Year
    for( unsigned int gasNumber = 0; gasNumber < getNumInputGases(); ++gasNumber ) {
        gasFile << setw(9) << sInputGasNames[ gasNumber ] << ',';
    }
    gasFile << endl;
    /*line 5: Gas units*/
    gasFile << setw(5) << "Year,";
    for( unsigned int gasNumber = 0; gasNumber < getNumInputGases(); ++gasNumber ) {
        gasFile << setw(9) << sInputGasUnits[ gasNumber ] << ',';
    }
    gasFile << endl;
    
    // Setup the output format.
    gasFile.setf( ios::right, ios::adjustfield );
    gasFile.setf( ios::fixed, ios::floatfield );
    gasFile.setf( ios::showpoint );
    
    for( unsigned int year = startYear; year <= endYear; year++ ) {
		int period =  mModeltime->getyr_to_per( year );
		bool needToInterpolate = !mModeltime->isModelYear(year);
		
        // We want to write emissions to MAGICC annually to eliminate errors
		// with the LUC CO2 emissions
        if( needToInterpolate && ( year > finalCalYear || year == GAS_EMK_CRIT_YEAR ) ) {
            // note period is offset by 1
			int nextYear = mModeltime->getper_to_yr( period );
            int prevYear = mModeltime->getper_to_yr( period - 1);
            gasFile << setw( 4 ) << year << ",";
            // Write out all the gases.
            for( unsigned int gasNumber = 0; gasNumber < getNumInputGases(); ++gasNumber ){
				// We are passing LUC emissions to MAGICC annually.
				if ( sInputGasNames[ gasNumber ] == "CO2NetLandUse" ) {
					gasFile << setw( 8 ) << setprecision( 2 )
					<< mLUCEmissionsByYear[ year - mModeltime->getStartYear() - 1 ];
				
				}
				// For all other emissions, we need to linearly interpolate
				else {
					gasFile << setw( 8 ) << setprecision( 2 )
						<< util::linearInterpolateY( year, prevYear, nextYear,
											mEmissionsByGas[ gasNumber ][ period - 2 ],
											mEmissionsByGas[ gasNumber ][ period - 1 ] );
					
				}
                // Write a comma as long as this is not the last gas.
                if( gasNumber != getNumInputGases() - 1 ){
                    gasFile << ",";
                }
                else {
                    gasFile << endl;
                }
            }
        }
		else if ( !needToInterpolate ) {
			gasFile << setw( 4 ) << year << ",";
			// Write out all the gases.
			for( unsigned int gasNumber = 0; gasNumber < getNumInputGases(); ++gasNumber ){
				// We are passing LUC emissions to MAGICC annually.
				if ( sInputGasNames[ gasNumber ] == "CO2NetLandUse" ) {
					gasFile << setw( 8 ) << setprecision( 2 )
					<< mLUCEmissionsByYear[ year - mModeltime->getStartYear() - 1 ];
					
				}
				else {
					gasFile << setw( 8 ) << setprecision( 2 ) 
						<< mEmissionsByGas[ gasNumber ][ period - 1 ];
				}
				// Write a comma as long as this is not the last gas.
				if( gasNumber != getNumInputGases() - 1 ){
					gasFile << ",";
				}
				else {
					gasFile << endl;
				}
			}
		}
    }
	int year = endYear;
	int period = mModeltime->getmaxper();
	for ( unsigned int extra = 0; extra < getNumAdditionalGasPoints(); extra++ ) {
		year = year + 100;
		gasFile << setw( 4 ) << year << ",";
		// Write out all the gases.
		for( unsigned int gasNumber = 0; gasNumber < getNumInputGases(); ++gasNumber ){
			// We are passing LUC emissions to MAGICC annually.
			if ( sInputGasNames[ gasNumber ] == "CO2NetLandUse" ) {
				int index = mModeltime->getEndYear() - mModeltime->getStartYear() + extra - 1;
				gasFile << setw( 8 ) << setprecision( 2 )
				<< mLUCEmissionsByYear[ index ];
				
			}
			else {
				gasFile << setw( 8 ) << setprecision( 2 ) << mEmissionsByGas[ gasNumber ][ period + extra - 1 ];
			}
			// Write a comma as long as this is not the last gas.
			if( gasNumber != getNumInputGases() - 1 ){
				gasFile << ",";
			}
			else {
				gasFile << endl;
			}
		}
	}
	
    gasFile.close();

    // First overwrite parameters
    overwriteMAGICCParameters( );

    // now actually run the model.
    ILogger& mainLog = ILogger::getLogger( "main_log" );
    mainLog.setLevel( ILogger::NOTICE );
    mainLog << "Calling the climate model..."<< endl;
    CLIMAT();
    mainLog.setLevel( ILogger::DEBUG );
    mainLog << "Finished with CLIMAT()" << endl;
    mIsValid = true;
    return true;
}

/* \brief Get the number of input gases to MAGICC.
* \return The fixed number of gases.
*/
unsigned int MagiccModel::getNumInputGases(){
    // This is a C++ trick to find the number of elements in
    // an array by dividing its size by the size of the first element.
    return sizeof( sInputGasNames ) / sizeof( sInputGasNames[ 0 ] );
}

double MagiccModel::getConcentration( const string& aGasName,
                                      const int aYear ) const
{
    /*! \pre A gas was requested. */
    assert( !aGasName.empty() );

    // Check if the climate model has been run and that valid year is passed in.
    if( !mIsValid || !isValidClimateModelYear( aYear ) ){
        return -1;
    }

    int year = aYear;
    int gasNumber = util::searchForValue( mOutputGasNameMap, aGasName );
    if ( gasNumber != 0 ) {
        return GETGHGCONC( gasNumber, year );
    }
    return -1;
}

double MagiccModel::getTemperature( const int aYear ) const {
    // Check if the climate model has been run and that valid year is passed in.
    if( !mIsValid || !isValidClimateModelYear( aYear ) ){
        return -1;
    }

    // Need to store the year locally so it can be passed by reference.
    int year = aYear;
    return GETGMTEMP( year );
}

double MagiccModel::getForcing( const string& aGasName, const int aYear ) const {
    // Check if the climate model has been run and that valid year is passed in.
    if( !mIsValid || !isValidClimateModelYear( aYear ) ){
        return -1;
    }

    int year = aYear;
    int gasNumber = util::searchForValue( mOutputGasNameMap, aGasName );
    if ( gasNumber != 0 ) {
        return GETFORCING( gasNumber, year );
    }
    return -1;
}

double MagiccModel::getNetTerrestrialUptake( const int aYear ) const {
    // Check if the climate model has been run and that valid year is passed in.
    if( !mIsValid || !isValidClimateModelYear( aYear ) ){
        return -1;
    }
    
    // MAGICC array is not sized to include 1990 for net terrestrial uptake
    // so just return 0.
    if( aYear < 1990 ) {
        return 0;
    }

    int year = aYear;
    int itemNumber = 10;
    return GETCARBONRESULTS( itemNumber, year );
}

double MagiccModel::getNetOceanUptake( const int aYear ) const {
    // Check if the climate model has been run and that valid year is passed in.
    if( !mIsValid || !isValidClimateModelYear( aYear ) ){
        return -1;
    }
    
    // MAGICC array is not sized to include 1990 for net ocean uptake
    // so just return 0.
    if( aYear < 1990 ) {
        return 0;
    }

    int year = aYear;
    int itemNumber = 4;
    return GETCARBONRESULTS( itemNumber, year );
}

double MagiccModel::getNetLandUseChangeEmission( const int aYear ) const {
    // Check if the climate model has been run and that valid year is passed in.
    if( !mIsValid || !isValidClimateModelYear( aYear ) ){
        return -1;
    }
    
    // MAGICC array is not sized to include 1990 for net land use change emissions
    // so just return 0.
    if( aYear < 1990 ) {
        return 0;
    }

    int itemNumber = 2;
    int year = aYear;
    return GETCARBONRESULTS( itemNumber, year );
}

double MagiccModel::getTotalForcing( const int aYear ) const {
    // Check if the climate model has been run and that valid year is passed in.
    if( !mIsValid || !isValidClimateModelYear( aYear ) ){
        return -1;
    }

    // Need to store the year and gas number locally so it can be passed by reference.
    int year = aYear;
    int gasNumber = 0; // global forcing
    return GETFORCING( gasNumber, year );
}


int MagiccModel::getCarbonModelStartYear() const {
    return mCarbonModelStartYear;
}

/*! \brief Write out data from the emissions model to the main csv output file.*/
void MagiccModel::printFileOutput() const {
    void fileoutput3(string var1name,string var2name,string var3name,
        string var4name,string var5name,string uname,vector<double> dout);

    // Fill up a vector of CO2 concentrations.
    vector<double> data( mModeltime->getmaxper() );
    for( int period = 0; period < mModeltime->getmaxper(); ++period ){
        data[ period ] = getConcentration( "CO2", mModeltime->getper_to_yr( period ) );
    }
    // Write out the data to the text output function protocol
    fileoutput3( "global","MAGICC"," "," ","Concentration","PPM", data );
    
    // TEMP?
    for( int period = 0; period < mModeltime->getmaxper(); ++period ){
        data[ period ] = getEmissions( "CO2", mModeltime->getper_to_yr( period ) );
    }

    // Write out the data to the text output function protocol
    fileoutput3( "global","MAGICC"," "," ","CO2-Emissions","TgC", data );

    // Fill up a vector of total forcing.
    for( int period = 0; period < mModeltime->getmaxper(); ++period ){
        data[ period ] = getTotalForcing( mModeltime->getper_to_yr( period ) );
    }
    // Write out the data to the text output function protocol
    fileoutput3( "global","MAGICC"," "," ","Forcing","W/m^2", data );

    // Fill up a vector of CO2 forcing.
    for( int period = 0; period < mModeltime->getmaxper(); ++period ){
        data[ period ] = getForcing( "CO2", mModeltime->getper_to_yr( period ) );
    }
    // Write out the data to the text output function protocol
    fileoutput3( "global","MAGICC"," "," ","CO2-Forcing","W/m^2", data );

    // Fill up a vector of Global Mean Temperature.
    for( int period = 0; period < mModeltime->getmaxper(); ++period ){
        data[ period ] = getTemperature( mModeltime->getper_to_yr( period ) );
    }
    // Write out the data to the text output function protocol
    fileoutput3( "global","MAGICC"," "," ","Temperature","degC", data );

    // Fill up a vector of Net Terrestrial Uptake.
    for( int period = 0; period < mModeltime->getmaxper(); ++period ){
        data[ period ] = getNetTerrestrialUptake( mModeltime->getper_to_yr( period ) );
    }
    // Write out the data to the text output function protocol
    fileoutput3( "global","MAGICC"," "," ","NetTerUptake","GtC", data );

    // Fill up a vector of Ocean Uptake.
    for( int period = 0; period < mModeltime->getmaxper(); ++period ){
        data[ period ] = getNetOceanUptake( mModeltime->getper_to_yr( period ) );
    }
    // Write out the data to the text output function protocol
    fileoutput3( "global","MAGICC"," "," ","OceanUptake","GtC", data );
    
    // Fill up a vector of Net Land-Use Emissions.
    for( int period = 0; period < mModeltime->getmaxper(); ++period ){
        data[ period ] = getNetLandUseChangeEmission( mModeltime->getper_to_yr( period ) );
    }
    // Write out the data to the text output function protocol
    fileoutput3( "global","MAGICC"," "," ","netLUEm","GtC", data );
}

/*! \brief Write out data from the emissions model to the database. */
void MagiccModel::printDBOutput() const {
   void dboutput4(string var1name,string var2name,string var3name,string var4name,
        string uname,vector<double> dout);

    // Fill up a vector of concentrations.
    vector<double> data( mModeltime->getmaxper() );
    for( int period = 0; period < mModeltime->getmaxper(); ++period ){
        data[ period ] = getConcentration( "CO2", mModeltime->getper_to_yr( period ) );
    }
    // Write out the data to the database.
     dboutput4( "global", "General", "Concentrations", "Period", "PPM", data );
 
    // TEMP?
    for( int period = 0; period < mModeltime->getmaxper(); ++period ){
        data[ period ] =
            getEmissions( "CO2", mModeltime->getper_to_yr( period ) );
    }
    // Write out the data to the database.
     dboutput4( "global", "General", "CO2 Emissions", "Period", "TgC", data );

    // Fill up a vector of total forcing.
    for( int period = 0; period < mModeltime->getmaxper(); ++period ){
        data[ period ] = getTotalForcing( mModeltime->getper_to_yr( period ) );
    }
    dboutput4( "global", "General", "Forcing", "Period","W/m^2", data );

    // Fill up a vector of CO2 forcing.
    for( int period = 0; period < mModeltime->getmaxper(); ++period ){
        data[ period ] = getForcing( "CO2", mModeltime->getper_to_yr( period ) );
    }
     dboutput4( "global", "General", "CO2-Forcing", "Period","W/m^2", data );

     // Fill up a vector of Global Mean Temperature.
    for( int period = 0; period < mModeltime->getmaxper(); ++period ){
        data[ period ] = getTemperature( mModeltime->getper_to_yr( period ) );
    }

    // Write out the data to the database.
    dboutput4( "global", "General", "Temperature", "Period", "degC", data );

    // Fill up a vector of Net Terrestrial Uptake.
    for( int period = 0; period < mModeltime->getmaxper(); ++period ){
        data[ period ] = getNetTerrestrialUptake( mModeltime->getper_to_yr( period ) );
    }

    // Write out the data to the database.
    dboutput4( "global", "General", "NetTerUptake", "Period", "GtC", data );

    // Fill up a vector of Ocean Uptake.
    for( int period = 0; period < mModeltime->getmaxper(); ++period ){
        data[ period ] = getNetOceanUptake( mModeltime->getper_to_yr( period ) );
    }

    // Write out the data to the database.
    dboutput4( "global", "General", "OceanUptake", "Period", "GtC", data );

    // Fill up a vector of Net Land-Use Emissions.
    for( int period = 0; period < mModeltime->getmaxper(); ++period ){
        data[ period ] = getNetLandUseChangeEmission( mModeltime->getper_to_yr( period ) );
    }

    // Write out the data to the database.
    dboutput4( "global", "General", "netLUEm", "Period", "GtC", data );
}

/*! \brief Updates a visitor with information from the the climate model.
* \param aVisitor Visitor to update.
* \param aPeriod Period for which to perform the update, -1 means all periods.
* \todo Add a function to the magicc library to get the temperature.
*/
void MagiccModel::accept( IVisitor* aVisitor, const int aPeriod ) const {
    aVisitor->startVisitClimateModel( this, aPeriod );
    aVisitor->endVisitClimateModel( this, aPeriod );
}

/*! \brief Return the name of the XML node for this object.
* \return The XML element name for the object.
*/
const string& MagiccModel::getXMLNameStatic(){
    const static string XML_NAME = "MagiccModel";
    return XML_NAME;
}

/*! \brief Return the index of a gas name within the internal matrix or
*          INVALID_GAS_NAME if the gas name is not found.
* \param aGasName Name of the gas to search for.
* \return The index of the gas or INVALID_GAS_NAME if it is not found.
*/
int MagiccModel::getInputGasIndex( const string& aGasName ) const {
    // Search the gas name vector.
    for( unsigned int i = 0; i < getNumInputGases(); ++i ){
        if( sInputGasNames[ i ] == aGasName ){
            return i;
        }
    }
    // The gas does not exist.
    return INVALID_GAS_NAME;
}

/*! \brief Read initial emissions data from a file.
* \details This function reads in the initial values for each gas from an input
*          file specified by the configuration file.
* \warning input file data lines must end in a ","
* \warning The model currently assumes that all input files are in terms of model periods.
* \todo make routine more general to interpolate between read-in periods if necessary to get model data
*/
void MagiccModel::readFile(){
    // Open up the file with default emissions data and read it in.
    // Default to value given in config file, but use read-in value if that exists
    string gasFileName = Configuration::getInstance()->getFile( "GHGInputFileName" );
    if ( mGHGInputFileName != "" ) {
        gasFileName = mGHGInputFileName;
    }
    ifstream inputGasFile;
    inputGasFile.open( gasFileName.c_str(), ios::in ); // open input file for reading
    util::checkIsOpen( inputGasFile, gasFileName );

    // Read in all other gases from fossil fuels. The first line will tell us how many
    // lines to parse.  We are assuming there are a maximum of getNumAdditionalGasPoints()
    // years beyond the end of the model years.
    int remainingExtraGasPoints = getNumAdditionalGasPoints();
    int numInputGasPoints;
    inputGasFile >> numInputGasPoints;  //reads the # of data points in the first line
    
    // Ignore initial header lines of file.
    const int SKIP_LINES = 5;
    for ( unsigned int i = 0; i < SKIP_LINES; ++i ){
        inputGasFile.ignore( 200, '\n' );
    }
    
    // Now read in all the gases. 
    for ( int gasInputIndex = 0; gasInputIndex <= numInputGasPoints; ++gasInputIndex ) {
        // Parse the year as a double so that the comma after the year does
        // not cause parsing to fail.
        double year = 0;
        inputGasFile >> year;
        inputGasFile.ignore( 80, ',' ); // skip comma
        int yearAsInteger = static_cast<int>( floor( year ) );
        int gasYearIndex = -1;
        
        if( yearAsInteger > mModeltime->getEndYear() ) {
            // parsing a year beyond the last model year
            if( remainingExtraGasPoints == 0 ) {
                ILogger& mainLog = ILogger::getLogger( "main_log" );
                mainLog.setLevel( ILogger::WARNING );
                mainLog << "Invalid Year: " << year << " found while parsing MAGICC input file " << gasFileName.c_str() << endl;
                mainLog << "The maximum: " << getNumAdditionalGasPoints()
                    << " points beyond the last model year have already been parsed" << endl;
            }
            else {
                gasYearIndex = mModeltime->getmaxper() + getNumAdditionalGasPoints()
                    - remainingExtraGasPoints - 1;
                // Decrement the remaining allowed lines beyond the last model year.
                --remainingExtraGasPoints;
            }
        }
        else {
            // A potential model year.  If an invalid model year was parsed modeltime
            // will print the warning and return 0.  Note that 1975 is ignored and the
            // indices are then shift over by 1 so we must subtract 1 here.
            gasYearIndex = mModeltime->getyr_to_per( yearAsInteger ) - 1;
        }
        
        // Loop through all the gases.
        for( unsigned int gasNumber = 0; gasNumber < getNumInputGases(); ++gasNumber ){
            double value;
            inputGasFile >> value;
            // Only set the parsed value if we have a valid year, no warnings are necessary
            // as they should have already been given.
            if( gasYearIndex != -1 ) {
                mEmissionsByGas[ gasNumber ][ gasYearIndex ] = value;
            }
            inputGasFile.ignore( 80, ','); // skip comma
        }

        inputGasFile.ignore( 80, '\n' ); // next line
    }
    inputGasFile.close();
}

/*! \brief Return the name of the XML node for this object.
* \return The XML element name for the object.
*/
const string& MagiccModel::getnetDefor80sName(){
    const static string NET_DEF_80s_NAME = "netDef80s";
    return NET_DEF_80s_NAME;
}

/*!
 * \brief Get the maximum number of data points that we are expecting beyond the
 *        last model year.
 * \return The constant number of additional data points that will be read-in and
 *         written-out for MAGICC input files.
 */
int MagiccModel::getNumAdditionalGasPoints() {
    static const int ADDITIONAL_POINTS = 2;
    return ADDITIONAL_POINTS;
}
