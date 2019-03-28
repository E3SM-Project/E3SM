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
* \file energy_balance_table.cpp
* \ingroup Objects
* \brief The EnergyBalanceTable class source file.
*
* \author Pralit Patel
*/

#include "util/base/include/definitions.h"
#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "reporting/include/energy_balance_table.h"

#include "sectors/include/sector.h"
#include "sectors/include/subsector.h"
#include "resources/include/resource.h"
#include "resources/include/subresource.h"
#include "containers/include/scenario.h"
#include "util/base/include/model_time.h"
#include "technologies/include/technology.h"
#include "functions/include/minicam_input.h"
#include "reporting/include/storage_table.h"
#include "technologies/include/ioutput.h"
#include "sectors/include/energy_final_demand.h"
#include "technologies/include/iproduction_state.h"

using namespace std;
extern Scenario* scenario;

/*!
 * \brief Constructor which needs the region we will visit, a stream to write results to, whether to print a condensed,
 *		  and if we want to include non-calibrated values in the table.
 * \param aRegionName The name of the region we are collecting data for
 * \param aFile The output stream to write results to (not necessarily a file)
 * \param aPrintCondensed If we should print a condensed version of this table
 * \param aIncludeNonCalValues If we should include non-calibrated values in the table
 */
EnergyBalanceTable::EnergyBalanceTable( const string& aRegionName, ostream& aFile, const bool aPrintCondensed,
									    const bool aIncludeNonCalValues ):
mFile( aFile ), mRegionName( aRegionName ), mTable( new StorageTable ), mIsTechOperating( false ),
mCalOutput( -1 ), mPrintCondensed( aPrintCondensed ), mIncludeNonCalValues( aIncludeNonCalValues )
{
}

EnergyBalanceTable::~EnergyBalanceTable() {
}

/*!
 * \brief When finished gathering data this will write results to mFile.
 */
void EnergyBalanceTable::finish() const {
    if( mPrintCondensed ) {
        writeCondensedTable();
    }
    else {
        writeFullTable();
    }
}

/*!
 * \brief Write the complete version of the table.
 * \details This will dissaggregate a sector to differentiate between subsectors and technologies.
 */
void EnergyBalanceTable::writeFullTable() const {
    /*!
     * \pre Assumes that mPrintCondensed is false.
     */
    assert( !mPrintCondensed );

    mFile << "Energy Balance Table";
    // give a note if we did not include non-calibrated values
    if( !mIncludeNonCalValues ) {
        mFile << "; Warning non-calibrated values are not included in this table.";
    }
    mFile << endl << "Region: " << mRegionName << endl;

    const vector<string> rows = mTable->getRowLabels();
    const vector<string> cols = mTable->getColLabels();

    vector<string>::const_iterator colIt;
    
    // assume that the total column is the last column
    const string totalColName = cols[ cols.size() -1 ];

    // first write heading info
    // we try to line up sectors / subsectors / technology and use --- to fill space
    // note that the corresponding sector to it's subsectors is left justified
    string prevStr = "";
    mFile << "Sectors:,";
    for( colIt = cols.begin(); colIt != cols.end(); ++colIt ) {
        string currSector = objects::searchForValue( mColToSectorMap, *colIt );
        if( *colIt == totalColName ) {
            mFile << totalColName << ',';
        }
        else if( currSector != prevStr ) {
            prevStr = currSector;
            mFile << prevStr << ',';
        }
        else {
            mFile << "---,";
        }
    }
    mFile << "Difference" << endl;

    // note that the corresponding subsector to it's technologies is left justified
    mFile << "Subsectors:,";
    for( colIt = cols.begin(); colIt != cols.end(); ++colIt ) {
        string currSubsector = objects::searchForValue( mColToSubsectorMap, *colIt );
        if( *colIt != totalColName && currSubsector != prevStr ) {
            prevStr = currSubsector;
            mFile << prevStr << ',';
        }
        else {
            mFile << "---,";
        }
    }
    mFile << "---" << endl;

    mFile << "Input,";
    for( colIt = cols.begin(); colIt != cols.end(); ++colIt ) {
        if( *colIt != totalColName ) {
            mFile << objects::searchForValue( mColToTechMap, *colIt );
        }
        else {
            mFile << "---";
        }
        mFile << ',';
    }
    mFile << "---" << endl;

    // write table data
    for( vector<string>::const_iterator rowIt = rows.begin(); rowIt != rows.end(); ++rowIt ) {
        if( *rowIt != "Output" ) {
            mFile << *rowIt << ',';
            for( colIt = cols.begin(); colIt != cols.end(); ++colIt ) {
                mFile << mTable->getValue( *rowIt, *colIt ) << ',';
            }
            
            // calculate the Difference and add it as the last column
            mFile << getTotalSectorOutput( *rowIt ) - mTable->getValue( *rowIt, totalColName ) << endl;
        }
    }

    // output must be the last row
    mFile << "Output,";
    for( colIt = cols.begin(); colIt != cols.end(); ++colIt ) {
        mFile << mTable->getValue( "Output", *colIt ) << ',';
    }
    mFile << endl << endl;
}

/*!
 * \brief Write the condensed version of the table.
 * \details This will print aggregated inputs and outputs at the sector level.
 */
void EnergyBalanceTable::writeCondensedTable() const {
    /*!
     * \pre Assumes that mPrintCondensed is true.
     */
    assert( mPrintCondensed );

    mFile << "Energy Balance Table";
    // give a note if we did not include non-calibrated values
    if( !mIncludeNonCalValues ) {
        mFile << "; Warning non-calibrated values are not included in this table.";
    }
    mFile << endl << "Region: " << mRegionName << endl;

    const vector<string> rows = mTable->getRowLabels();
    const vector<string> cols = mTable->getColLabels();
	
    vector<string>::const_iterator colIt;

    // write heading info
    mFile << "Sectors:,";
    for( colIt = cols.begin(); colIt != cols.end(); ++colIt ) {
        mFile << *colIt << ',';
    }
    mFile << endl;
	
    // write table data
    for( vector<string>::const_iterator rowIt = rows.begin(); rowIt != rows.end(); ++rowIt ) {
        if( *rowIt != "Output" ) {
            mFile << *rowIt << ',';
            for( colIt = cols.begin(); colIt != cols.end(); ++colIt ) {
                mFile << mTable->getValue( *rowIt, *colIt ) << ',';
            }
            mFile << endl;
        }
    }
	
    // output must be the last row
    mFile << "Output,";
    for( colIt = cols.begin(); colIt != cols.end(); ++colIt ) {
        mFile << mTable->getValue( "Output", *colIt ) << ',';
    }
    mFile << endl;
    
    // print a difference row for convienience
    mFile << "Difference,";
    // assume that the total column is the last column
    const string totalColName = cols[ cols.size() -1 ];
    for( colIt = cols.begin(); colIt != cols.end(); ++colIt ) {
        mFile << mTable->getValue( "Output", *colIt ) - mTable->getValue( *colIt, totalColName ) << ',';
    }
    mFile << endl << endl;
}

void EnergyBalanceTable::startVisitSector( const Sector* aSector, const int aPeriod ) {
    mCurrentSector = aSector->getName();
}

void EnergyBalanceTable::endVisitSector( const Sector* aSector, const int aPeriod ) {
    mCurrentSector.clear();
}

void EnergyBalanceTable::startVisitSubsector( const Subsector* aSubsector, const int aPeriod ) {
    mCurrentSubsector = aSubsector->getName();
}

void EnergyBalanceTable::endVisitSubsector( const Subsector* aSubsector, const int aPeriod ) {
    mCurrentSubsector.clear();
}

void EnergyBalanceTable::startVisitTechnology( const Technology* aTechnology, const int aPeriod ) {
    if( aTechnology->mProductionState[ aPeriod ]->isOperating() /*&& aTechnology->getShareWeight() > 0 */) {
        // set the current technology and tell the input to visit
        mIsTechOperating = true;
        mCurrentTech = aTechnology->getName();
        mCalOutput = aTechnology->getCalibrationOutput( false, "", aPeriod );

        // need to add the column labels explicitly
        const string currKey = getKey();
        if( mColToSubsectorMap.find( currKey ) == mColToSubsectorMap.end() ) {
            mTable->addColumn( currKey );
            mColToSectorMap[ currKey ] = mCurrentSector;
            mColToSubsectorMap[ currKey ] = mCurrentSubsector;
            mColToTechMap[ currKey ] = mCurrentTech;
        }
    }
}

void EnergyBalanceTable::endVisitTechnology( const Technology* aTechnology, const int aPeriod ) {
    mIsTechOperating = false;
    mCurrentTech.clear();
    mCalOutput = -1;
}

void EnergyBalanceTable::startVisitMiniCAMInput( const MiniCAMInput* aInput, const int aPeriod ) {
    if( mIsTechOperating && aInput->hasTypeFlag( IInput::ENERGY ) ) {
        const double useForNoCalValue = mIncludeNonCalValues ? aInput->getPhysicalDemand( aPeriod ) : 0;
        mTable->addToType( aInput->getName(), getKey(), aInput->getCalibrationQuantity( aPeriod ) == -1 ?
                          useForNoCalValue : aInput->getCalibrationQuantity( aPeriod ) );
    }
}

void EnergyBalanceTable::startVisitOutput( const IOutput* aOutput, const int aPeriod ) {
    // only visit the output if the technology is operating
    // TODO: is this the way to check primary output and is this the behavior we want?
    //       what to do about secondary outputs?
    // TODO: why is isSameType checking this value, shouldn't it be the xml name atleast?
    if( mIsTechOperating && aOutput->isSameType( "primary-output" ) ) {
        const double useForNoCalValue = mIncludeNonCalValues ? aOutput->getPhysicalOutput( aPeriod ) : 0;
        mTable->addToType( "Output", getKey(), mCalOutput == -1 ? useForNoCalValue : mCalOutput );
    }
    else if( mIsTechOperating && mPrintCondensed && mCalOutput > -1 ) {
        // this is to handle secondary outputs
        // TODO: this currently does not handle mIncludeNonCalValues
        // TODO: this is currently only occuring for the condensed table since there
        //       is nowhere to attribute this secondary output, I could create a bogus
        //       subsector/technology to put these values into
        const string currSectorTemp = mCurrentSector;
        const double primaryOutput = mCalOutput;
        typedef IOutput::OutputList::const_iterator COutputListIterator;
        IOutput::OutputList outputList =
            aOutput->calcPhysicalOutput( primaryOutput,
                                    mRegionName,
                                    0,
                                    aPeriod );
        
        for( COutputListIterator i = outputList.begin(); i != outputList.end(); ++i ) {
            // have to override the current sector so getKey gets the key we really want
            // TODO: what if this column has not been set up yet?
            mCurrentSector = i->first;
            mTable->addToType( "Output", getKey(), i->second );
        }
        mCurrentSector = currSectorTemp;
    }
}

void EnergyBalanceTable::startVisitEnergyFinalDemand( const EnergyFinalDemand* aEnergyFinalDemand, const int aPeriod ) {
    // add the demands for goods by the final demands
    const string key = aEnergyFinalDemand->getName() + "-Final-Demand";
    mTable->addColumn( key );
    mColToSectorMap[ key ] = key;
    mColToSubsectorMap[ key ] = key;
    mColToTechMap[ key ] = key;
    // TODO: putting base service would not be correct after calibration periods
    mTable->addToType( aEnergyFinalDemand->getName(), key, aEnergyFinalDemand->mBaseService[ aPeriod ] );
}

void EnergyBalanceTable::startVisitResource( const AResource* aResource, const int aPeriod )
{
    mCurrentSector = aResource->getName();
}

void EnergyBalanceTable::endVisitResource( const AResource* aResource, const int aPeriod )
{
    mCurrentSector.clear();
}

void EnergyBalanceTable::startVisitSubResource( const SubResource* aSubResource, const int aPeriod )
{
    mCurrentSubsector = mCurrentTech = aSubResource->getName();
    const string key = getKey();
    const double useForNoCalValue = mIncludeNonCalValues ? aSubResource->getAnnualProd( aPeriod ) : 0;
    mTable->addColumn( key );
    mColToSectorMap[ key ] = key;
    mColToSubsectorMap[ key ] = key;
    mColToTechMap[ key ] = key;
    mTable->addToType( "Output", key, aSubResource->mCalProduction[ aPeriod ] == -1 
                                           ? useForNoCalValue : aSubResource->mCalProduction[ aPeriod ] );
}

void EnergyBalanceTable::endVisitSubResource( const SubResource* aSubResource, const int aPeriod )
{
    mCurrentSubsector.clear();
    mCurrentTech.clear();
}

/*!
 * \brief Gets the key, or lowest level column name, which can be used to index into name lookup maps or table.
 * \details If mPrintCondensed is set that it will use the mCurrentSector otherwise it will combine the
 *			tech name + subsector name + sector name to have a complete unique description of the column.
 * \pre The method assumes that the mCurrentSector, mCurrentSubsector, and mCurrentTech have all been set.
 * \return The key that should be used.
 */
string EnergyBalanceTable::getKey() const {
    return  mPrintCondensed ? mCurrentSector : mCurrentTech + mCurrentSubsector + mCurrentSector;
}

/*!
 * \brief Gets the total sector output for the given sector name.
 * \details This is useful when creating an expanded table since multiple
 *          columns will need to be sumed to get the total.  We do this by
 *          iterating over each entry in mColToSectorMap and comparing the
 *          value to see if it matches the given sector name and add the
 *          "Output" row of that column to the sum if it does.
 * \param aSectorName The name of the sector to get the total output for.
 * \return The total output for the sector.
 */
double EnergyBalanceTable::getTotalSectorOutput( const string& aSectorName ) const {
    double outputSum = 0;
    for( map<string, string>::const_iterator it = mColToSectorMap.begin(); it != mColToSectorMap.end(); ++it ) {
        if( (*it).second == aSectorName ) {
            outputSum += mTable->getValue( "Output", (*it).first );
        }
    }
    return outputSum;
}

void EnergyBalanceTable::setRegionName( const string& aRegionName ) {
    mRegionName = aRegionName;
}
