#ifndef _IGROWTH_CALCULATOR_H_
#define _IGROWTH_CALCULATOR_H_
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
 * \file igrowth_calculator.h
 * \ingroup Objects
 * \brief The IGrowthCalculator interface header file.
 * \author Josh Lurz
 */
#include <string>
#include <vector>
#include <xercesc/dom/DOMNode.hpp>

class Tabs;
class IInvestable;
class Demographic;
class NationalAccount;

/*! 
 * \ingroup Objects
 * \brief Interface to an object which calculates the scalar that grows the
 *        overall level of sector investment.
 * \details Defines an interface to an object which decides by what fraction to
 *          increase investment in a ProductionSector from the previous period.
 *          This interface is defined very generally so that this decision can
 *          be calculated in a variety of ways. This class is only used when the
 *          Accelerator type of IInvestor is used, IGrowthCalculators are not
 *          required when investment levels are determined by solving for
 *          marginal profit equal to zero for each ProductionSector.
 * \author Josh Lurz
 */
class IGrowthCalculator
{
public:
    IGrowthCalculator();
    virtual ~IGrowthCalculator();

    // TODO: Fix this class to inherit so that this documentation is inherited.
    virtual void XMLParse( const xercesc::DOMNode* aCurr ) = 0;
    virtual void toDebugXML( const int period, std::ostream& aOut, Tabs* aTabs ) const = 0;
    virtual void toInputXML( std::ostream& aOut, Tabs* aTabs ) const = 0;
    
    /*!
     * \brief Calculate a scalar which will be used to grow investment from the
     *        previous period.
     * \details Determine the rate of growth for the ProductionSector from the
     *          previous period based on a series of parameters and the
     *          subsectors of the ProductionSector.
     * \todo Fix the ordering of the parameters to be more consistent.
     * \param aInvestable Vector of IInvestable children of the
     *        ProductionSector.
     * \param aDemographic Regional demographics container.
     * \param aNationalAccount Regional national accounts container.
     * \param aSectorName Name of the sector in which investment is occurring.
     * \param aRegionName Name of the region in which investment is occurring.
     * \param aInvestmentLogitExp Investment logit exponential.
     * \param aProfitElasExp Profit rate exponential.
     * \param aPeriod Period in which to calculate investment.
     * \return The investment growth scalar.
     */
    virtual double calcInvestmentDependencyScalar( const std::vector<IInvestable*>& aInvestables,
                                                   const Demographic* aDemographic,
                                                   const NationalAccount& aNationalAccount,
                                                   const std::string& aSectorName,
                                                   const std::string& aRegionName,
                                                   const double aInvestmentLogitExp,
                                                   const double aProfitElasExp,
                                                   const int aPeriod ) = 0;
};

// Define empty inline methods.
//! Constructor
inline IGrowthCalculator::IGrowthCalculator(){
}

//! Destructor
inline IGrowthCalculator::~IGrowthCalculator(){
}

#endif // _IGROWTH_CALCULATOR_H_
