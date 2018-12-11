#ifndef _INESTED_INPUT_H_
#define _INESTED_INPUT_H_
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
* \file inested_input.h
* \ingroup Objects
* \brief INestedInput interface header file.
* \author Pralit Patel
* \author Ron Sands
*/

#include <string>
#include <xercesc/dom/DOMNode.hpp>
#include <memory>

#include "functions/include/iinput.h"

class IFunction;
struct TechChange;

/*! 
* \ingroup Objects
* \brief Interface for a nested Input.
*
* \details Interface for being an IInput which can be contained within a nesting
*          stucture.  An input must provide all of the functionality of an IInput
*          as well as methods which will allow a production function to trigger
*          it's calculations through out the nesting structure.
* \author Pralit Patel
* \author Ron Sands
*/
class INestedInput : public IInput
{
public:
    /*!
     * \brief Constructor.
     * \details Inlined constructor to avoid compiler problems with abstract
     *          base classes. 
     */
    INestedInput();

    /*!
     * \brief Destructor.
     * \details Inlined destructor to avoid compiler problems with abstract base
     *          classes. 
     */
    virtual ~INestedInput();

    /*!
     * \brief Remove any empty inputs.
     */
    virtual void removeEmptyInputs() = 0;

    /*!
     * \brief Initialize the value, and price of the nodes.
     */
    virtual void initialize() = 0;

    /*!
     * \brief Calculates the coefficient by calling a production demand
     *   function's calcCoefficient for each subtree in the nesting.
     * \param aRegionName The name of the region this input is in.
     * \param aSectorName The name of the sector this input is in.
     * \param aTechPeriod The period converted from the tech's year.  This would
     *   be something other than zero for techs that come online later.
     */
    virtual void calcCoefficient( const std::string& aRegionName, const std::string& aSectorName,
        const int aTechPeriod ) = 0;

    /*!
     * \brief Changes elasticity for each node in the nesting.
     * \details This method is intended to be used when converting a new
     *          vintage technology to an old vintage technology.
     * \param aRegionName The name of the containing region.
     * \param aPeriod The current model period.
     * \param aAlphaZero  The alpha zero scaler.
     */
    virtual void changeElasticity( const std::string& aRegionName, const int aPeriod,
        const double aAlphaZero ) = 0;

    /*!
     * \brief Changes elasticity for each node in the nesting.
     * \details This method must be called for a new vintage technology
     *          to take care of the case where a nested input gets merged
     *          into another nested input forward in time that had read in
     *          a different sigma to use.
     * \param aRegionName The name of the containing region.
     * \param aPeriod The current model period.
     * \param aAlphaZero  The alpha zero scaler.
     */
    virtual void changeSigma( const std::string& aRegionName, const int aPeriod,
        const double aAlphaZero ) = 0;

    /*!
     * \brief Calculates the prices for each nested input.
     * \param aRegionName The name of the containing region.
     * \param aSectorName The name of the containing sector.
     * \param aPeriod The current model period.
     * \param aAlphaZero  The alpha zero scaler.
     */
    virtual void calcLevelizedCost( const std::string& aRegionName, const std::string& aSectorName,
        const int aPeriod, const double aAlphaZero ) = 0;

    /*!
     * \brief Calculates the input demand by calling a production demand
     *        function's calcDemand for each subtree in the nesting.
     * \param aRegionName The name of the containing region.
     * \param aSectorName The name of the containing sector.
     * \param aPeriod The current model period.
     * \param aPhysicalOutput The physical output at this nest used to drive demands
     *                        of any inputs nested below.
     * \param aUtilityParameterA The A paramater to be used in a UtilityDemandFunction,
     *                           TODO: this is a hack to get this value from the consumer
     *                           to the function where it is needed.  No other function would
     *                           attempt to utilize this paramater.
     * \param aAlphaZero  The alpha zero scaler.
     */
    virtual double calcInputDemand( const std::string& aRegionName, const std::string& aSectorName,
        const int aPeriod, const double aPhysicalOutput, const double aUtilityParameterA,
        const double aAlphaZero ) = 0;

    /*!
     * \brief Calculates the capital output ratio.
     * \details TODO
     * \todo Should I get rid of this completly, it may still be possible to compute
     *       however is not currently used.
     */
    virtual double calcCapitalOutputRatio( const std::string& aRegionName, const std::string& aSectorName,
        const int aPeriod, const double aAlphaZero ) = 0;

    /*!
     * \brief Calculates the levelized cost ignoring capital.
     * \details Calculates the levelized cost the same as calcLevelizedCost however
     *          any nested inputs with the flag IInput::CAPITAL will be omitted since
     *          capital is considered sunk cost and therefore non-variable.
     * \param aRegionName The name of the containing region.
     * \param aSectorName The name of the containing sector.
     * \param aPeriod The current model period.
     * \param aAlphaZero  The alpha zero scaler.
     */
    virtual void calcVariableLevelizedCost( const std::string& aRegionName, const std::string& aSectorName,
        const int aPeriod, const double aAlphaZero ) = 0;

    /*!
     * \breif Gets the function of the nested input.
     * \return A pointer to the funciton used by the nested input.
     */
    virtual const IFunction* getFunction() const = 0;
    
    /*!
     * \brief Gets the levelized cost for the root node.
     */
    virtual double getLevelizedCost( const std::string& aRegionName, const std::string& aSectorName,
        const int aPeriod ) const = 0;

    /*!
     * \brief Applies technical change at all levels of the nesting.
     * \param aRegionName The name of the containing region.
     * \param aSectorName The name of the containing sector.
     * \param aPeriod The current model period.
     * \param aTechChange A structure which can hold various types of technical
     *                    which may be applied.
     */
    virtual void applyTechnicalChange( const std::string& aRegionName, const std::string& aSectorName,
        const int aPeriod, const TechChange& aTechChange ) = 0;

    /*!
     * \brief Resets the flag to signal to recalculate the levelized cost.
     * \details The idea is to gain a performance increase by avoiding excessive
     *          levelized cost calculations.
     */
    virtual void resetCalcLevelizedCostFlag() = 0;
};

// Inline function definitions.
inline INestedInput::INestedInput() {
}

inline INestedInput::~INestedInput() {
}

#endif // _INESTED_INPUT_H_
