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
* \file batch_csv_outputter.cpp
* \ingroup Objects
* \brief The BatchCSVOutputter class source file for writing batch results to a csv file.
* \details This source file contains the definition for the startVisit and endVisit methods
*          for each class that the visitor visits.  Values will be written directly to the
*          file specified by the configuration paramater batchCSVOutputFile.
* \todo Figure out why I need a string to output doubles/chars to an AutoOutputFile.
* \author Pralit Patel
*/

#include "util/base/include/definitions.h"

#include "util/base/include/configuration.h"
#include "util/base/include/model_time.h"
#include "containers/include/scenario.h"
#include "marketplace/include/market.h"
#include "marketplace/include/imarket_type.h"
#include "climate/include/iclimate_model.h"

#include <string>

#include "reporting/include/batch_csv_outputter.h"


extern Scenario* scenario;

using namespace std;

/*! \brief Constructor
*/
BatchCSVOutputter::BatchCSVOutputter():
mFile( "batchCSVOutputFile", "batch-csv-out.csv" ),
mIsFirstScenario(true)
{
}

/*!
 * \brief Destructor
 */
BatchCSVOutputter::~BatchCSVOutputter(){
}


void BatchCSVOutputter::startVisitScenario( const Scenario* aScenario, const int aPeriod ) {

    // only output the header info one time
    // we can not put this in the constructor because we will not have a model time
    // a that point
    if( mIsFirstScenario ) {
        mFile << "Scenario" << ',';
        const Modeltime* modeltime = aScenario->getModeltime();

        for( int period = 0; period < modeltime->getmaxper(); ++period ) {
            // TODO: hard coding CO2
            const int year = modeltime->getper_to_yr( period );
            mFile << year << ' '<< "CO2 Price" << ',';
        }
        for( int period = 0; period < modeltime->getmaxper(); ++period ) {
            // TODO: hard coding CO2
            const int year = modeltime->getper_to_yr( period );
            mFile << year << ' '<< "CO2 Emissions" << ',';
        }

        int outputInterval = Configuration::getInstance()->getInt( "climateOutputInterval",
                                               scenario->getModeltime()->gettimestep( 0 ) );
        
        // print at least to 2100 if interval is set appropriately
        int endingYear = max( scenario->getModeltime()->getEndYear(), 2100 );
        
        for( int year = scenario->getModeltime()->getStartYear();
             year <= endingYear; year += outputInterval )
        {
            // TODO: hard coding CO2
            mFile << year << ' '<< "CO2 Concentration" << ',';
        }
        for( int year = scenario->getModeltime()->getStartYear();
             year <= endingYear; year += outputInterval )
        {
            // TODO: hard coding CO2
            mFile << year << ' '<< "CO2 Radiative Forcing" << ',';
        }
        for( int year = scenario->getModeltime()->getStartYear();
             year <= endingYear; year += outputInterval )
        {
            // TODO: hard coding CO2
            mFile << year << ' '<< "CO2 Temperature Change" << ',';
        }
        mFile << "Solved" << endl;
    }
    mIsFirstScenario = false;

    mFile << aScenario->getName() << ',';
    // TODO: perhaps write some date/time or something
}

void BatchCSVOutputter::startVisitMarket( const Market* aMarket, const int aPeriod ) {
    if( aMarket->getType() == IMarketType::TAX ) {
        /*!
         * \warninng This is assuming the periods will be visited in appropriate order.
         */
        mFile << "" << aMarket->getPrice() << ',';
        
        // would this be wrong if it didn't solve?
        //mFile << "" << aMarket->getDemand() << ',';
    }
}

void BatchCSVOutputter::startVisitClimateModel( const IClimateModel* aClimateModel, const int aPeriod ) {
    const Modeltime* modeltime = scenario->getModeltime();
    for( int period = 0; period < modeltime->getmaxper(); ++period ) {
        const int year = modeltime->getper_to_yr( period );
        mFile << "" << aClimateModel->getEmissions( "CO2", year ) << ',';
    }
    
    int outputInterval
        = Configuration::getInstance()->getInt( "climateOutputInterval",
                                   modeltime->gettimestep( 0 ) );

    // print at least to 2100 if interval is set appropriately
    int endingYear = max( modeltime->getEndYear(), 2100 );

    // Write the climate variables for all applicable years.
    for( int year = modeltime->getStartYear();
         year <= endingYear; year += outputInterval )
    {
        mFile << "" << aClimateModel->getConcentration( "CO2", year ) << ',';
    }
    for( int year = modeltime->getStartYear();
         year <= endingYear; year += outputInterval )
    {
        mFile << "" << aClimateModel->getForcing( "CO2", year) << ',';
    }
    for( int year = modeltime->getStartYear();
         year <= endingYear; year += outputInterval )
    {
        mFile << "" << aClimateModel->getTemperature( year ) << ',';
    }
}

/*!
 * \brief Writes whether the current scenario had solved or not.
 * \param aDidSolve Whether the current scenario solved.
 */
void BatchCSVOutputter::writeDidScenarioSolve( bool aDidSolve ) {
    mFile << aDidSolve << endl;
}
