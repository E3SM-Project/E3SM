#ifndef _MINICAM_LEONTIEF_PRODUCTION_FUNCTION_H_
#define _MINICAM_LEONTIEF_PRODUCTION_FUNCTION_H_
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
 * \file minicam_leontief_production_function.h
 * \ingroup Objects
 * \brief MinicamLeontiefProductionFuntion class header file.
 * \author Josh Lurz
 */

#include <string>
#include <vector>
#include "functions/include/ifunction.h"

class IInput;

/*! 
 * \ingroup Objects
 * \brief Defines the Leontief production function for Minicam.
 * \details This function provides a simplified leontief production function for
 *          use in MiniCAM technologies.
 * \author Josh Lurz
 */
class MinicamLeontiefProductionFunction : public IFunction {
public:

    double calcProfits( InputSet& aInputs,
                        const std::string& aRegionName,
                        const std::string& aSectorName,
						const double aShutdownCoef,
                        int aPeriod,
                        double aCapitalStock = 0,
                        double aAlphaZero = 0,
						double aSigma = 0 ) const;
    
    double calcCosts( const InputSet& aInputs,
                      const std::string& aRegionName,
                      const double aAlphaZero,
                      int aPeriod ) const;

	double calcDemand( InputSet& aInputs,
                       double aConsumption,
                       const std::string& aRegionName,
					   const std::string& aSectorName,
                       const double aShutdownCoef,
                       int aPeriod,
					   double aCapitalStock = 0,
                       double aAlphaZero = 0,
                       double aSigma = 0,
                       double aIBT = 0 ) const;
    
    double calcCoefficient( InputSet& aInputs,
                            double aConsumption,
                            const std::string& aRegionName,
                            const std::string& aSectorName,
                            int aPeriod,
                            double aSigma = 0,
                            double aIBT = 0,
							double aCapitalStock = 0 ) const;
	
    double changeElasticity( InputSet& aInputs,
                             const std::string& aRegionName,
                             double aPriceReceived,
							 double aProfits,
                             double aCapitalStock,
                             const int aPeriod,
                             double alphaZero = 0,
							 double aSigmaNew = 0,
                             double aSigmaOld = 0 ) const;
	
    double calcOutput( InputSet& aInputs,
                       const std::string& aRegionName,
                       const std::string& aSectorName,
                       const double aShutdownCoef,
                       int aPeriod,
                       double aCapitalStock = 0,
                       double aAlphaZero = 0,
                       double aSigma = 0 ) const;
    
    double calcExpProfitRate( const InputSet& aInputs,
                              const std::string& aRegionName, 
                              const std::string& aSectorName,
                              double aLifeTimeYears,
                              int aPeriod,
                              double aAlphaZero = 0,
                              double aSigma = 0 ) const;
    
    double calcLevelizedCost( const InputSet& aInputs,
                              const std::string& aRegionName,
                              const std::string& aSectorName,
                              int aPeriod,
                              double aAlphaZero,
                              double aSigma ) const;

    double getCapitalOutputRatio( const InputSet& aInputs,
                                  const std::string& aRegionName,
                                  const std::string& aSectorName,
                                  double aLifeTimeYears, int aPeriod,
								  double aAlphaZero = 0,
                                  double aSigma = 0 ) const;

    double applyTechnicalChange( InputSet& aInputs,
                                 const TechChange& aTechChange,
                                 const std::string& aRegionName,
                                 const std::string& aSectorName,
                                 const int aPeriod, 
                                 double aAlphaZero = 0,
                                 double aSigma = 0 ) const;

    double calcUnscaledProfits( const InputSet& aInputs, 
                                const std::string& aRegionName,
                                const std::string& aSectorName,
                                const int aPeriod,
                                const double aCapitalStock,
                                const double aAlphaZero,
                                const double aSigma ) const;
};

#endif // _MINICAM_LEONTIEF_PRODUCTION_FUNCTION_H_
