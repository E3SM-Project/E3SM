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

#ifndef _STANDARD_TECHNICAL_CHANGE_CALC_H_
#define _STANDARD_TECHNICAL_CHANGE_CALC_H_
#if defined(_MSC_VER)
#pragma once
#endif



/*! 
 * \file standard_technical_change_calc.h
 * \ingroup Objects
 * \brief StandardTechnicalChangeCalc class header file.
 * \author Josh Lurz
 */

#include <string>
#include <vector>
#include <xercesc/dom/DOMNode.hpp>
#include "technologies/include/itechnical_change_calc.h"
#include "util/base/include/value.h"

/*! 
 * \ingroup Objects
 * \brief This object is added on to Technologies to allow for Hicks-Neutral,
 *        material, and energy technical change to be applied to the inputs.
 * \details This object is responsible for calculating the cumulative
 *          Hicks-Neutral and input specific technical change. Hicks-Neutral is
 *          applied to all inputs through an out front scalar calculated by this
 *          class. Energy and Material technical change is applied to the
 *          coefficients of all energy or material inputs that do not read in an
 *          input specific technical change. An input cannot be both energy and
 *          material.
 *
 *          <b>XML specification for StandardTechnicalChangeCalc</b>
 *          - XML name: \c standard-technical-change-calc
 *          - Contained by: Technology
 *          - Parsing inherited from class: None
 *          - Attributes: None
 *          - Elements:
 *              - \c hicks-neutral StandardTechnicalChangeCalc::mHicksNeutralTechChange
 *              - \c energy-only StandardTechnicalChangeCalc::mEnergyTechChange
 *              - \c material-only StandardCaptureComponent::mMaterialTechChange
 *
 * \author Josh Lurz
 */
class StandardTechnicalChangeCalc: public ITechnicalChangeCalc {
public:
    static const std::string& getXMLNameStatic();

	StandardTechnicalChangeCalc();

    // Documentation inherits.
	virtual StandardTechnicalChangeCalc* clone() const;
	
    virtual bool isSameType( const std::string& aType ) const;
	
    virtual const std::string& getName() const;
    
    virtual bool XMLParse( const xercesc::DOMNode* aNode );
	
    virtual void toInputXML( std::ostream& aOut,
                             Tabs* aTabs ) const;
	
    virtual void toDebugXML( const int aPeriod,
                             std::ostream& aOut,
                             Tabs* aTabs ) const;

    // ITechnicalChangeCalc
    virtual void completeInit();

    virtual double calcAndAdjustForTechChange( std::vector<IInput*>& aInputs,
                                               PreviousPeriodInfo& aPreviousPeriodInfo,
                                               const IFunction* aProductionFunc,
                                               const std::string& aRegionName,
                                               const std::string& aSectorName,
                                               const int aPeriod ) const;

protected:
    //! Annual Hicks-Neutral technical change.
    Value mHicksNeutralTechChange;

    //! Annual energy only technical change.
    Value mEnergyTechChange;

    //! Annual material only technical change.
    Value mMaterialTechChange;
};

#endif // _STANDARD_TECHNICAL_CHANGE_CALC_H_
