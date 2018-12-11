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

#ifndef _RENEWABLE_INPUT_H_
#define _RENEWABLE_INPUT_H_
#if defined(_MSC_VER)
#pragma once
#endif

/*! 
 * \file renewable_input.h
 * \ingroup Objects
 * \brief RenewableInput class header file.
 * \author Josh Lurz
 */

#include <string>
#include <xercesc/dom/DOMNode.hpp>
#include "functions/include/minicam_input.h"
#include <vector>
#include "util/base/include/value.h"

class Tabs;

/*! 
 * \ingroup Objects
 * \brief Defines a single unlimited zero price input to a MiniCAM production
 *        function.
 * \details MiniCAM allows the production of energy from a renewable input that
 *          has no cost or limit to production. This is used in technologies
 *          such as hydro. These technologies theoretically convert an energy
 *          input into output, but the content of the energy input is not known
 *          to the model. While these inputs do require a cost to obtain, this
 *          cost is instead included in the non-energy costs.
 *
 *          <b>XML specification for RenewableInput</b>
 *          - XML name: \c renewable-input
 *          - Contained by: Technology
 *          - Parsing inherited from class: MiniCAMInput
 *          - Attributes:
 *              - \c name MiniCAMInput::mName
 *          - Elements: None.
 *
 * \author Josh Lurz
 */
class RenewableInput: public MiniCAMInput
{
    friend class InputFactory;
    friend class UnmanagedLandTechnology;
public:
    static const std::string& getXMLNameStatic();

    const std::string& getXMLReportingName() const;
    
    virtual RenewableInput* clone() const;

    virtual void XMLParse( const xercesc::DOMNode* aNode );

    virtual bool isSameType( const std::string& aType ) const;
    

    virtual void toInputXML( std::ostream& aOut,
                             Tabs* aTabs ) const;

    virtual void toDebugXML( const int aPeriod,
                             std::ostream& aOut,
                             Tabs* aTabs ) const;

    virtual void copyParam( const IInput* aInput,
                            const int aPeriod );

    virtual void completeInit( const std::string& aRegionName,
                               const std::string& aSectorName,
                               const std::string& aSubsectorName,
                               const std::string& aTechName,
                               DependencyFinder* aDependencyFinder,
                               const IInfo* aTechInfo);

    virtual void initCalc( const std::string& aRegionName,
                           const std::string& aSectorName,
                           const bool aIsNewInvestmentPeriod,
                           const bool aIsTrade,
                           const int aPeriod );

    virtual double getPhysicalDemand( const int aPeriod ) const;
    
    virtual void setPhysicalDemand( const double aPhysicalDemand,
                                    const std::string& aRegionName,
                                    const int aPeriod );

    virtual double getPrice( const std::string& aRegionName,
                             const int aPeriod ) const;
    
    virtual void setPrice( const std::string& aRegionName,
                           const double aPrice,
                           const int aPeriod );

    virtual double getCO2EmissionsCoefficient( const std::string& aGHGName,
                                            const int aPeriod ) const;
    
    virtual double getCoefficient( const int aPeriod ) const;

    virtual void setCoefficient( const double aCoefficient,
                                 const int aPeriod );

    virtual bool hasTypeFlag( const int aTypeFlag ) const;
    
    virtual double getIncomeElasticity() const;

    virtual double getPriceElasticity() const;

    virtual double getTechChange( const int aPeriod ) const;

	virtual double getCalibrationQuantity( const int aPeriod ) const;
protected:
    RenewableInput();

    // Constuctor to allow internal creation of this object
    RenewableInput( const std::string& aName );

    //! Physical Demand.
    std::vector<Value> mPhysicalDemand;
private:
    const static std::string XML_REPORTING_NAME; //!< tag name for reporting xml db 
};

#endif // _RENEWABLE_INPUT_H_
