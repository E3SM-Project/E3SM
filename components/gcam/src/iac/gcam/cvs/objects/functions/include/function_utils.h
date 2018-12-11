#ifndef _FUNCTION_UTILS_H_
#define _FUNCTION_UTILS_H_
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
* \file function_utils.h
* \ingroup Objects
* \brief The FunctionUtils class header file.
*
* \author Josh Lurz
* \author Sonny Kim
*/

#include <string>
#include <vector>
class IInput;
class IFunction;
struct TechChange;
class INestedInput;

typedef std::vector<IInput*> InputSet;

/*! \brief A structure which contains the information necessary for a
*          production function.
*/
struct ProductionFunctionInfo {
    //! The vector of inputs.
    const InputSet& mInputs;

    //! Pointer to the technology's production function.
    const IFunction* mProductionFunction;

    //! The current sigma the production function is using.
    const double mSigma;

    //! Alpha zero used to scale the output of the production.
    const double mAlphaZeroScaler;

    //! Amount of capital stock the vintage owns.
    const double mCapitalStock;
    
    //! The root of the nested inputs.
    const INestedInput* mNestedInputRoot;
};

/*! 
* \ingroup Objects
* \brief This class is a set of static helper methods which the various
*        production and demand functions use.
* \details TODO
* \author Sonny Kim, Josh Lurz
*/
class FunctionUtils {
public:
    static void scaleCoefficientInputs( InputSet& input,
                                        double scaler, const int aPeriod );
    
    static double getCurrencyDemandSum( const InputSet& aInputs, const int aPeriod );

    static double getPhysicalDemandSum( const InputSet& aInputs, const int aPeriod );
    
    static double getCoefSum( const InputSet& input, const int aPeriod );
    
    static IInput* getInput( const InputSet& aInputs,
                             const std::string& aInputName );
    
    static IInput* getCapitalInput( const InputSet& aInputs );
    
    static IInput* getInputByType( const InputSet& aInputs, const int aTypeFlag );
    
    static IInput* getNumeraireInput( const InputSet& aInputs );
    
    static double getRho( const double aSigma );
    
    static double PMT( double aRate, double aNper, double aPV );

    // **** Begin Unit Conversions
    // TODO: remove when units conversion is working
    static double HOURS_PER_YEAR( void );

    static double HOURS_PER_DAY( void );

    static double DEFLATOR_1975_PER_DEFLATOR_2003( void );

    static double DEFLATOR_1975_PER_DEFLATOR_2005( void );

    static double DEFLATOR_1990_PER_DEFLATOR_1975( void );

    static double GJ_PER_KWH( void );

    static double GJ_PER_MWH( void );

    static double GJ_PER_EJ( void );

    static double EJ_PER_GWH( void );

    static double MWH_PER_GWH( void );

    static double GWH_PER_GJ( void );

    static double MWH_PER_GJ( void );

    static double KG_PER_METRIC_TON( void );

    // **** End Unit Conversions

    static double getNetPresentValueMult( const InputSet& aInputs,
                                          const std::string& aRegionName,
                                          const double aLifetimeYears,
                                          const int aPeriod  );
    
    static double calcNetPresentValueMult( const double aDiscountRate,
                                           const double aLifetime );
    
    static double getExpectedPriceReceived( const InputSet& aInputs,
                                            const std::string& aRegionName,
                                            const std::string& aGoodName,
                                            const double aLifetimeYears,
                                            const int aPeriod );

    static void setPricePaid( const std::string& aRegionName,
                              const std::string& aGoodName,
                              const int aPeriod,
                              const double aPricePaid );

    static double getPricePaid( const std::string& aRegionName,
                                const std::string& aGoodName,
                                const int aPeriod );

    static void setPriceReceived( const std::string& aRegionName,
                                  const std::string& aGoodName,
                                  const int aPeriod,
                                  const double aPriceReceived );

    static double getPriceReceived( const std::string& aRegionName,
                                    const std::string& aGoodName,
                                    const int aPeriod );

    static double applyTechnicalChangeInternal( InputSet& input,
                                                const TechChange& aTechChange, 
                                                const std::string& regionName,
                                                const std::string& sectorName, 
                                                const int aPeriod,
                                                double alphaZero,
                                                double sigma );

    static double getTechChangeForInput( const IInput* aInput,
                                         const TechChange& aTechChange,
                                         const int aPeriod );

    static bool isFixedPrice( const std::string& aRegionName,
                              const std::string& aGoodName,
                              const int aPeriod );

    static double getMarketConversionFactor( const std::string& aRegionName,
                                             const std::string& aGoodName,
                                             const bool aMustExist = true );

    static void copyInputParamsForward( const InputSet& aPrevInputs,
                                        InputSet& aCurrInputs,
                                        const int aPeriod );

    static double getCO2Coef( const std::string& aRegionName,
                              const std::string& aGoodName,
                              const int aPeriod,
                              const bool aMustExist = true );

    static double calcPriceRatio( const std::string& aRegionName,
                                  const IInput* aInput,
                                  const int aBasePeriod,
                                  const int aCurrentPeriod );

    static InputSet getLeafInputs( const INestedInput* aNestedInput );
    
    // TODO: should the following two be in another utility?
    static void setCapitalGoodPrice( const std::string& aRegionName,
                                     const int aPeriod,
                                     const double aCapitalGoodPrice );

    static double getCapitalGoodPrice( const std::string& aRegionName,
                                       const int aPeriod );
};

#endif // _FUNCTION_UTILS_H_
