#ifndef _LOGIT_PRODUCTION_FUNCTION_H_
#define _LOGIT_PRODUCTION_FUNCTION_H_
#if defined(_MSC_VER)
#pragma once
#endif

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
* \file logit_production_function.h
* \ingroup Objects
* \brief LogitProductionFunction class header file.
* \author Pralit Patel
*/

#include <string>
#include <vector>
#include "functions/include/nested_ces_production_function.h"

class IInput;

/*! 
* \ingroup Objects
* \brief Defines the Logit production function.
* \details A logit production function which could be used with in a nested
*          input structure.  The sigma paramaters represents the logit exponent
*          to use and a coefficient is sysynonymous with a share weight.  This
*          production function would automatically calculate share weights in
*          the base year making them relative to the first input (TODO: it may be 
*          worth it to make them relative to the largest for the sake of consistency
*          with subsector/technology logits) using the calcCoefficient method.
* \todo The changeElasticity method does not do anything at the moment however it
*       probably should so that a user could change the logit exponent and have the
*       share weights get automatically adjusted such that the shares do not change.
* \todo Not sure what to do about applyTechnicalChange since scaling share weights
*       is probably not what the user would be expecting.  The logit must preserve
*       quatities thus there is no way to alter the input/output ratio except by
*       shifting the share between the other inputs at the same level.  Note to
*       get around this we read technical change at the parent node operating as a
*       CES which has the effect of scaling input/output ratios equally for all inputs
*       in this logit nest.
* \note This function inhereted from the NestedCESProductionFunction just because the
*       calcOutput, calcExpProfitRate, getCapitalOutputRatio, applyTechnicalChange, and
*       calcUnscaledProfits would be the same.  The calculations in those methods are not
*       really specific to logit or CES and probably more so just SGM.  I think it might
*       be a good idea to move those methods to the abstract base class and in general
*       clean up the function interface.
* \author Pralit Patel
*/
class LogitProductionFunction : public NestedCESProductionFunction {
public:
    double calcDemand( InputSet& input, double consumption, const std::string& regionName,
                       const std::string& sectorName, const double aShutdownCoef, int period,
                       double capitalStock = 0, double alphaZero = 0, double sigma = 0, double IBT = 0 ) const;
    
    double calcCoefficient( InputSet& input, double consumption, const std::string& regionName,
                            const std::string& sectorName, int period, double sigma = 0, double IBT = 0,
                            double capitalStock = 0 ) const;
    
    double changeElasticity( InputSet& input, const std::string& aRegionName, double priceReceived,
                             double aProfits, double capitalStock, const int aPeriod, double alphaZero = 0,
                             double sigmaNew = 0, double sigmaOld = 0 ) const;
    
    /*
    double calcOutput( InputSet& input, const std::string& regionName,
                       const std::string& sectorName, const double aShutdownCoef,
                       int period, double capitalStock = 0, double alphaZero = 0, double sigma = 0 ) const;
    
    double calcExpProfitRate( const InputSet& input, const std::string& regionName,
                              const std::string& sectorName, double aLifeTimeYears, int period, double alphaZero = 0,
                              double sigma = 0 ) const;
                              */

    double calcLevelizedCost( const InputSet& aInputs, const std::string& aRegionName,
                              const std::string& aSectorName, int aPeriod, double aAlphaZero, double aSigma ) const;

    /*
    double getCapitalOutputRatio( const InputSet& aInputs, const std::string& aRegionName,
                                  const std::string& aSectorName, double aLifeTimeYears, int aPeriod,
                                  double aAlphaZero, double aSigma ) const;
    
    double applyTechnicalChange( InputSet& input, const TechChange& aTechChange,
                                 const std::string& regionName,const std::string& sectorName, const int aPeriod, 
                                 double alphaZero = 0, double sigma = 0 ) const;
    
    double calcUnscaledProfits( const InputSet& aInputs, 
                                const std::string& aRegionName,
                                const std::string& aSectorName,
                                const int aPeriod,
                                const double aCapitalStock,
                                const double aAlphaZero,
                                const double aSigma ) const;
                                */
};

#endif // _LOGIT_PRODUCTION_FUNCTION_H_
