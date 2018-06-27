#ifndef _DEMAND_INPUT_H_
#define _DEMAND_INPUT_H_
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
* \file demand_input.h
* \ingroup Objects
* \brief DemandInput class header file.
* \author Pralit Patel
* \author Sonny Kim
*/

#include <string>
#include <cassert>
#include <xercesc/dom/DOMNode.hpp>

#include "functions/include/sgm_input.h"

class Tabs;
class IVisitor;

/*! 
* \ingroup Objects
* \brief Defines a single input to a demand function.
* \note Not all demand functions use DemandInputs, as some are modified
*       production functions.
* \details TODO
* \author Pralit Patel, Sonny Kim
*/
class DemandInput : public SGMInput
{
public:
	DemandInput();
	DemandInput* clone() const;
    void copyParam( const IInput* aInput,
                    const int aPeriod  );

    double getPriceElasticity() const;
    double getIncomeElasticity() const;
    
	static const std::string& getXMLNameStatic();
    void accept( IVisitor* aVisitor ) const;

    void copyParamsInto( DemandInput& aDemandInput,
                         const int aPeriod ) const;

    void copyParamsInto( ProductionInput& aProductionInput,
                         const int aPeriod ) const;

    void accept( IVisitor* aVisitor, const int aPeriod ) const;
protected:
    bool XMLDerivedClassParse( const std::string& nodeName, const xercesc::DOMNode* curr );
    void toInputXMLDerived( std::ostream& out, Tabs* tabs ) const;
    void toDebugXMLDerived( const int period, std::ostream& out, Tabs* tabs ) const;
    const std::string& getXMLName() const;
private:
    //! Income Elasticity
    Value mIncomeElasticity;

    //! Price Elasticity
    Value mPriceElasticity;
};

#endif // _DEMAND_INPUT_H_
