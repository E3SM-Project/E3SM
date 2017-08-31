#ifndef _IFUNCTION_H_
#define _IFUNCTION_H_
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
* \file ifunction.h
* \ingroup Objects
* \brief IFunction class header file.
* \author Josh Lurz
*/

#include <string>
#include <vector>
#include "functions/include/iinput.h"

/*! \brief A structure which groups the types of technical change a technology
*          may define.
* \details A technology may define none of these types, in which case only input
*          level technical change will apply, it may define some of these types,
*          or it may define all of these types. Input level technical change
*          will override material or energy technical change. Hicks neutral
*          technical change will always be applied.
* \note TechChange is defined here instead of BaseTechnology to avoid IFunction
*       including BaseTechnology.
* \author Josh Lurz
* \ingroup Objects
*/
struct TechChange {
    /*!
     * \brief Constructor.
     * \param aMaterialTechChange Technical change applied to material inputs.
     * \param aEnergyTechChange Technical change applied to energy inputs.
     * \param aHicksTechChange Technical change applied to all inputs.
     */
    TechChange( const double aMaterialTechChange,
                const double aEnergyTechChange,
                const double aHicksTechChange )
    : mMaterialTechChange( aMaterialTechChange ),
      mEnergyTechChange( aEnergyTechChange ),
      mHicksTechChange( aHicksTechChange )
    {}

    /*!
     * \brief Default constructor.
     */
    TechChange()
        : mMaterialTechChange( 0 ),
          mEnergyTechChange( 0 ),
          mHicksTechChange( 0 )
    {}

    //! Technical change rate for energy usage.
    double mEnergyTechChange;
    
    //! Technical change rate for the materials.
    double mMaterialTechChange;

    //! Hicks neutral technical change(applies to all inputs together)
    double mHicksTechChange;
};

typedef std::vector<IInput*> InputSet;
typedef InputSet::const_iterator CInputSetIterator;
typedef InputSet::iterator InputSetIterator;

/*! \brief The interface to a generic production or demand function.
* \details TODO
* \author Sonny Kim, Josh Lurz
* \ingroup Objects
*/
class IFunction {
public:
	virtual double calcDemand( InputSet& aInputs, double aConsumption, 
        const std::string& aRegionName, const std::string& aSectorName, 
        const double aShutdownCoef,
        int aPeriod, 
        double aCapitalStock = 0, double aAlphaZero = 0, double aSigma = 0, double aIBT = 0 ) const = 0;
    
    virtual double calcCoefficient( InputSet& aInput, double aConsumption, 
        const std::string& aRegionName, const std::string& aSectorName, int aPeriod, double aSigma = 0, 
        double aIBT = 0, double aCapitalStock = 0 ) const = 0;

    virtual double changeElasticity( InputSet& aInputs,  const std::string& aRegionName, 
		double aPriceReceived, double aProfits, double aCapitalStock, const int aPeriod, double aAlphaZero = 0,
		double aSigmaNew = 0, double aSigmaOld = 0 ) const = 0;
	
    // TODO: This is not really a feature of the functions and could be removed
    // once demand and production functions are separated.
    virtual double applyTechnicalChange( InputSet& aInputs, const TechChange& aTechChange,
        const std::string& aRegionName, const std::string& aSectorName, 
        const int aPeriod, double aAlphaZero = 0, double aSigma = 0 ) const = 0;

	virtual double calcOutput( InputSet& aInputs, const std::string& aRegionName,
        const std::string& aSectorName, const double aShutdownCoef, int aPeriod,
        double aCapitalStock = 0, double aAlphaZero = 0, 
        double aSigma = 0 ) const = 0;

    virtual double calcProfits( InputSet& aInputs, const std::string& aRegionName,
        const std::string& aSectorName, const double aShutdownCoef,
        int aPeriod, double aCapitalStock = 0, double aAlphaZero = 0, 
        double aSigma = 0 ) const = 0;

    virtual double calcLevelizedCost( const InputSet& aInputs,
        const std::string& aRegionName, const std::string& aSectorName, int aPeriod,
        double aAlphaZero, double aSigma ) const = 0;

    virtual double calcCosts( const InputSet& aInputs,
                              const std::string& aRegionName,
                              const double aAlphaZero,
                              int aPeriod ) const = 0;

	virtual double calcExpProfitRate( const InputSet& aInputs, const std::string& aRegionName,
        const std::string& aSectorName, double aLifeTimeYears, int aPeriod, double aAlphaZero = 0,
        double aSigma = 0 ) const = 0;

    virtual double getCapitalOutputRatio( const InputSet& aInputs,
        const std::string& aRegionName, const std::string& aSectorName, double aLifeTimeYears,
        int aPeriod, double aAlphaZero, double aSigma ) const = 0;
    
    virtual double calcUnscaledProfits( const InputSet& aInputs, 
                                        const std::string& aRegionName,
                                        const std::string& aSectorName,
                                        const int aPeriod,
                                        const double aCapitalStock,
                                        const double aAlphaZero,
                                        const double aSigma ) const = 0;
};

#endif // _IFUNCTION_H_


