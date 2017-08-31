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

#ifndef _RETIRED_PRODUCTION_STATE_H_
#define _RETIRED_PRODUCTION_STATE_H_
#if defined(_MSC_VER)
#pragma once
#endif

/*! 
 * \file retired_production_state.h
 * \ingroup Objects
 * \brief The RetiredProductionState header file.
 * \author Josh Lurz
 */

#include "technologies/include/iproduction_state.h"

/*! 
 * \ingroup Objects
 * \brief The production state for a retired Technology.
 * \details This production state is used by Technology's which have already
 *          fully retired. This state is reached when the Technology's lifetime
 *          was completed before the current year. This does not include
 *          technologies which have made a short-term shutdown decision or which
 *          have deprecated to near zero output. Retired technologies produce
 *          zero output.
 * \author Josh Lurz
 */
class RetiredProductionState: public IProductionState
{
    friend class ProductionStateFactory;
public:
    // ISimpleComponent methods.
    RetiredProductionState* clone() const;
    virtual bool isSameType( const std::string& aType ) const;
    virtual const std::string& getName() const;
	virtual void toDebugXML( const int aPeriod,
                             std::ostream& aOut,
                             Tabs* aTabs ) const;
    
    // IProductionState methods.
    virtual double calcProduction( const std::string& aRegionName,
                                   const std::string& aSectorName,
                                   const double aVariableOutput,
                                   const MarginalProfitCalculator* aMarginalProfitCalc,
                                   const double aFixedOutputScaleFactor,
                                   const std::vector<IShutdownDecider*>& aShutdownDeciders,
                                   const int aPeriod ) const;

    virtual void setBaseOutput( const double aBaseOutput,
                                const int aBaseYear );

    virtual bool isNewInvestment() const;

    virtual bool isOperating() const;
protected:
    /*!
     * \brief Protected constructor which may only be called from the
     *        ProductionStateFactory.
     */
    RetiredProductionState();

    /*
     * \brief Get the static name of this object.
     * \return The static name of this object.
     */
    static const std::string& getXMLNameStatic();
};

#endif // _RETIRED_PRODUCTION_STATE_H_
