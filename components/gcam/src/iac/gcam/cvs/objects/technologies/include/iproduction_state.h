#ifndef _IPRODUCTION_STATE_H_
#define _IPRODUCTION_STATE_H_
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

#include "util/base/include/istandard_component.h"
#include <string>
#include <vector>

class IShutdownDecider;

/*! 
 * \file iproduction_state.h
 * \ingroup Objects
 * \brief The IProductionState interface header file.
 * \author Josh Lurz
 */

class MarginalProfitCalculator;

/*! 
 * \ingroup Objects
 * \brief This interface defines how a MiniCAM technology decides the quantity
 *        of output to produce.
 * \details MiniCAM technologies are always in exactly one production state. The
 *          production state defines how they determine the aggregate level of
 *          output to produce. The production state is valid only for a single
 *          period, a Technology may progress through several production states
 *          over its lifetime.
 * \author Josh Lurz
 */
class IProductionState: public ISimpleComponent
{
public:
	// Clone operator must be declared explicitly even though it is inherited
    // from IStandardComponent so that the return type can be changed. Since
    // this class is a subtype of IStandardComponent, this is legal and referred
    // to as a covariant return type.
	virtual IProductionState* clone() const = 0;

    /*
     * \brief Calculate the amount for the Technology to produce.
     * \details Calculates the amount of production for the Technology based on
     *          the current state of the technology. The current state changes
     *          over the lifetime of the technology.
     * \param aVariableOutput Variable output requested for production by the
     *        subsector.
     * \param aMarginalProfitCalc Calculator of the marginal profit rate for the
     *        Technology.
     * \param aFixedOutputScaleFactor Fraction by which to reduce all fixed
     *        output.
     * \param aShutdownDeciders Set of objects responsible for determining how
     *        much capital to depreciate or shutdown.
     * \param aPeriod Period for which to calculate production.
     * \return Production quantity.
     */
    virtual double calcProduction( const std::string& aRegionName,
                                   const std::string& aSectorName,
                                   const double aVariableOutput,
                                   const MarginalProfitCalculator* aMarginalProfitCalc,
                                   const double aFixedOutputScaleFactor,
                                   const std::vector<IShutdownDecider*>& aShutdownDeciders,
                                   const int aPeriod ) const = 0;

    /*
     * \brief Set the quantity of output in the initial period of the
     *        technology.
     * \param aBaseOutput Base level of output.
     * \param aBaseYear Initial year of the Technology.
     */
    virtual void setBaseOutput( const double aBaseOutput,
                                const int aBaseYear ) = 0;

    /*! 
     * \brief Returns whether the production state represents a Technology in
     *        its initial investment period.
     * \details Returns whether the Technology is considered new investment.
     * \return Whether the production state represents a Technology in
     *        its initial investment period.
     */
    virtual bool isNewInvestment() const = 0;

    /*! 
     * \brief Returns whether the production state represents an operating
     *        technology.
     * \details Returns whether the technology is operating. A technology is
     *          considered operating if it is not retired, it does not have to
     *          produce positive output. Both a fixed output technology with
     *          zero output and a vintaged technology that has performed a
     *          short-term shutdown decision will report that they are
     *          operating.
     * \return Whether the production state represents an operating technology.
     */
    virtual bool isOperating() const = 0;

    /*!
     * \brief Returns a constant representing no fixed output.
     * \return A constant representing no fixed output.
     */
    static double fixedOutputDefault() {
        return -1;
    }
};

#endif // _ISHUTDOWN_DECIDER_H_
