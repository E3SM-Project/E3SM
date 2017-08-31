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
 * \file minicam_price_elasticity_function.cpp
 * \ingroup Objects
 * \brief The MinicamPriceElasticityFunction class source file.
 * \author Josh Lurz
 * \todo Rename this file to match the header!
 */

#include "util/base/include/definitions.h"
#include <string>
#include <vector>
#include <cmath>
#include <cassert>

#include "functions/include/minicam_price_elasticity_function.h"
#include "functions/include/iinput.h"
#include "containers/include/scenario.h"
#include "marketplace/include/marketplace.h"
#include "util/base/include/util.h"
#include "util/base/include/model_time.h"
#include "functions/include/function_utils.h"
#include "util/logger/include/ilogger.h"
#include "sectors/include/sector_utils.h"

using namespace std;

extern Scenario* scenario;

double MinicamPriceElasticityFunction::calcCosts( const InputSet& aInput,
                                                  const string& aRegionName,
                                                  const double aAlphaZero,
                                                  int aPeriod ) const
{
    int normPeriod = SectorUtils::getDemandNormPeriod( aPeriod );

    double totalCost = 0;
    for( CInputSetIterator input = aInput.begin(); input != aInput.end(); ++input ) {
        double priceRatio = FunctionUtils::calcPriceRatio( aRegionName, *input,
                                                           normPeriod, aPeriod );
        totalCost +=  ( *input )->getCoefficient( aPeriod )
                      * pow( priceRatio, ( *input )->getPriceElasticity() )
                      * ( *input )->getPrice( aRegionName, aPeriod );
    }
    return totalCost / aAlphaZero;
}

double MinicamPriceElasticityFunction::calcProfits( InputSet& input,
                                                    const string& regionName,
                                                    const string& sectorName,
                                                    const double aShutdownCoef,
                                                    int period,
                                                    double capitalStock,
                                                    double alphaZero,
                                                    double sigma ) const
{
    double profits = calcUnscaledProfits( input, regionName, sectorName, period, capitalStock, alphaZero, sigma );
    return aShutdownCoef * profits;
}

double MinicamPriceElasticityFunction::calcCoefficient( InputSet& input,
                                                        double consumption,
                                                        const string& regionName,
                                                        const string& sectorName,
                                                        int period,
                                                        double sigma,
                                                        double indBusTax,
                                                        double capitalStock ) const
{
    // MiniCAM coefficients are read in and do not need to be calculated.
    return 1;
}

double MinicamPriceElasticityFunction::changeElasticity( InputSet& input,
                                                         const string& aRegionName,
                                                         double priceReceived,
                                                         double aProfits,
                                                         double capitalStock,
                                                         const int aPeriod,
                                                         double alphaZero,
                                                         double sigmaNew,
                                                         double sigmaOld ) const
{

    // MiniCAM coefficients cannot have their elasticity changed.
    return 1;
}

double MinicamPriceElasticityFunction::calcDemand( InputSet& aInputs,
                                                   double aConsumption,
                                                   const string& aRegionName,
                                                   const string& aSectorName,
                                                   const double aShutdownCoef,
                                                   int aPeriod,
                                                   double aCapitalStock,
                                                   double aAlphaZero,
                                                   double aSigma,
                                                   double aIBT ) const
{
    int normPeriod = SectorUtils::getDemandNormPeriod( aPeriod );

    // personalIncome == demand.
    double totalDemand = 0;
    for( CInputSetIterator input = aInputs.begin(); input != aInputs.end(); ++input ) {
        double priceRatio = FunctionUtils::calcPriceRatio( aRegionName, *input,
                                                           normPeriod, aPeriod );

        // Calculate the input demand including any changes due to the price ratio.
        double inputDemand = ( *input )->getCoefficient( aPeriod )
                             * aConsumption
                             * pow( priceRatio, ( *input )->getPriceElasticity() ); 

        ( *input )->setPhysicalDemand( inputDemand, aRegionName, aPeriod );
        totalDemand += inputDemand;
    }
    return totalDemand;
}

double MinicamPriceElasticityFunction::calcExpProfitRate( const InputSet& aInputs,
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

double MinicamPriceElasticityFunction::calcLevelizedCost( const InputSet& aInputs,
                                                          const string& aRegionName,
                                                          const string& aSectorName,
                                                          int aPeriod,
                                                          double aAlphaZero,
                                                          double aSigma ) const
{
    int normPeriod = SectorUtils::getDemandNormPeriod( aPeriod );

    // Loop through the inputs and add on their costs.
    double levelizedCost = 0;
    for( CInputSetIterator input = aInputs.begin(); input != aInputs.end(); ++input ) {
        double priceRatio = FunctionUtils::calcPriceRatio( aRegionName, *input,
                                                           normPeriod, aPeriod );
        levelizedCost += ( *input )->getCoefficient( aPeriod )
                         * pow( priceRatio, ( *input )->getPriceElasticity() )
                         * ( *input )->getPrice( aRegionName, aPeriod );
    }
    return levelizedCost / aAlphaZero;
}

double MinicamPriceElasticityFunction::calcOutput( InputSet& aInput,
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

double MinicamPriceElasticityFunction::getCapitalOutputRatio( const InputSet& aInputs,
                                                              const string& aRegionName,
                                                              const string& aSectorName,
                                                              double aLifetimeYears,
                                                              int aPeriod,
                                                              double aAlphaZero,
                                                              double aSigma ) const
{
    // This undefined for MiniCAM.
    return -1;
}

double MinicamPriceElasticityFunction::applyTechnicalChange( InputSet& aInputs,
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

double MinicamPriceElasticityFunction::calcUnscaledProfits( const InputSet& input,
                                                            const string& regionName,
                                                            const string& sectorName,
                                                            const int period,
                                                            const double capitalStock,
                                                            const double alphaZero,
                                                            const double sigma ) const
{
    // Return the price of this good minus the levelized cost.
    return scenario->getMarketplace()->getPrice( sectorName, regionName, period )
           - calcLevelizedCost( input, regionName, sectorName, period, alphaZero, sigma );
}
