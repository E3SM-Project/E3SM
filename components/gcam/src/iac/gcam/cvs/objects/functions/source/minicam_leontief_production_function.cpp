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
 * \file minicam_leontief_production_function.cpp
 * \ingroup Objects
 * \brief The MinicamLeontiefProductionFunction class source file.
 * \author Josh Lurz
 */

#include "util/base/include/definitions.h"
#include <string>
#include <vector>
#include <cmath>
#include <cassert>

#include "functions/include/minicam_leontief_production_function.h"
#include "functions/include/iinput.h"
#include "containers/include/scenario.h"
#include "marketplace/include/marketplace.h"
#include "util/base/include/util.h"
#include "util/base/include/model_time.h"
#include "functions/include/function_utils.h"
#include "util/logger/include/ilogger.h"

using namespace std;

extern Scenario* scenario;

double MinicamLeontiefProductionFunction::calcCosts( const InputSet& aInputs,
                                                     const string& aRegionName,
                                                     const double aAlphaZero,
                                                     int aPeriod ) const
{
    double totalCost = 0;
    for( CInputSetIterator input = aInputs.begin(); input != aInputs.end(); ++input ) {
        totalCost += ( *input )->getCoefficient( aPeriod ) * ( *input )->getPrice( aRegionName, aPeriod );
    }
    return totalCost / aAlphaZero;
}

double MinicamLeontiefProductionFunction::calcProfits( InputSet& aInputs,
                                                       const string& aRegionName,
                                                       const string& aSectorName,
                                                       const double aShutdownCoef,
                                                       int aPeriod,
                                                       double aCapitalStock,
                                                       double aAlphaZero,
                                                       double aSigma ) const
{
    double profits = calcUnscaledProfits( aInputs, aRegionName, aSectorName,
                                          aPeriod, aCapitalStock, aAlphaZero, aSigma );

    // Scale the profits by the scaler which determines how much of the vintage
    // was shutdown.
    return aShutdownCoef * profits;
}

double MinicamLeontiefProductionFunction::calcCoefficient( InputSet& aInputs,
                                                           double aConsumption,
                                                           const string& aRegionName,
                                                           const string& aSectorName,
                                                           int aPeriod,
                                                           double aSigma,
                                                           double aIndBusTax,
                                                           double aCapitalStock ) const
{
    // MiniCAM coeffiicients are read in and do not require calculating.
    return 1;
}

double MinicamLeontiefProductionFunction::changeElasticity( InputSet& input,
                                                            const string& aRegionName,
                                                            double aPriceReceived,
                                                            double aProfits,
                                                            double aCapitalStock,
                                                            const int aPeriod,
                                                            double aAlphaZero,
                                                            double aSigmaNew,
                                                            double aSigmaOld ) const
{
    // MiniCAM coefficients cannot have their elasticity changed.
    return 1;
}

double MinicamLeontiefProductionFunction::calcDemand( InputSet& aInputs,
                                                      double aPersonalIncome,
                                                      const string& aRegionName,
                                                      const string& aSectorName,
                                                      const double aShutdownCoef,
                                                      int aPeriod,
                                                      double aCapitalStock,
                                                      double aAlphaZero,
                                                      double aSigma,
                                                      double aIBT ) const
{
    assert( aAlphaZero >= 1 );

    // PersonalIncome == demand.
    double totalDemand = 0;
    for( CInputSetIterator input = aInputs.begin(); input != aInputs.end(); ++input ) {
        double inputDemand = ( *input )->getCoefficient( aPeriod ) * aPersonalIncome / aAlphaZero;
        ( *input )->setPhysicalDemand( inputDemand, aRegionName, aPeriod );
        totalDemand += inputDemand;
    }
    return totalDemand;
}

double MinicamLeontiefProductionFunction::calcExpProfitRate( const InputSet& aInputs,
                                                             const string& aRegionName,
                                                             const string& aSectorName,
                                                             double aLifetimeYears,
                                                             int aPeriod,
                                                             double aAlphaZero,
                                                             double aSigma ) const
{
    // MiniCAM technologies do not have an expected profit rate.
    return -1;
}

double MinicamLeontiefProductionFunction::calcLevelizedCost( const InputSet& aInputs,
                                                             const string& aRegionName,
                                                             const string& aSectorName,
                                                             int aPeriod,
                                                             double aAlphaZero,
                                                             double aSigma ) const
{
    assert( aAlphaZero >= 1 );

    // Loop through the inputs and add on their costs.
    double levelizedCost = 0;
    for( CInputSetIterator input = aInputs.begin(); input != aInputs.end(); ++input ) {
        levelizedCost += ( *input )->getCoefficient( aPeriod )
                         * ( *input )->getPrice( aRegionName, aPeriod );
    }
    return levelizedCost / aAlphaZero;
}

double MinicamLeontiefProductionFunction::calcOutput( InputSet& aInputs,
                                                      const string& aRegionName,
                                                      const string& aSectorName,
                                                      const double aShutdownCoef,
                                                      int aPeriod,
                                                      double aCapitalStock,
                                                      double aAlphaZero,
                                                      double aSigma ) const
{
    // MiniCAM output is calculated exogenously from the production function.
    return 0;
}

double MinicamLeontiefProductionFunction::getCapitalOutputRatio( const InputSet& aInputs,
                                                                 const string& aRegionName,
                                                                 const string& aSectorName,
                                                                 double aLifetimeYears,
                                                                 int aPeriod,
                                                                 double aAlphaZero,
                                                                 double aSigma ) const
{
    // MiniCAM does not have capital so this is undefined.
    return -1;
}

double MinicamLeontiefProductionFunction::applyTechnicalChange( InputSet& aInputs,
                                                                const TechChange& aTechChange,
                                                                const string& aRegionName,
                                                                const string& aSectorName,
                                                                const int aPeriod,
                                                                double aAlphaZero,
                                                                double aSigma ) const
{
    return FunctionUtils::applyTechnicalChangeInternal( aInputs, aTechChange, aRegionName, aSectorName,
                                                        aPeriod, aAlphaZero, aSigma );
}

double MinicamLeontiefProductionFunction::calcUnscaledProfits( const InputSet& aInputs,
                                                               const string& aRegionName,
                                                               const string& aSectorName,
                                                               const int aPeriod,
                                                               const double aCapitalStock,
                                                               const double aAlphaZero,
                                                               const double aSigma ) const
{
    // Return the price of this good minus the levelized cost.
    return scenario->getMarketplace()->getPrice( aSectorName, aRegionName, aPeriod )
           - calcLevelizedCost( aInputs, aRegionName, aSectorName, aPeriod, aAlphaZero, aSigma );
}
