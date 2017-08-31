/*
 * This software, which is provided in confidence, was prepared by employees of
 * Pacific Northwest National Labratory operated by Battelle Memorial Institute.
 * Battelle has certain unperfected rights in the software which should not be
 * copied or otherwise disseminated outside your organization without the
 * express written authorization from Battelle. All rights to the software are
 * reserved by Battelle. Battelle makes no warranty, express or implied, and
 * assumes no liability or responsibility for the use of this software.
 */

#ifndef _INPUT_OM_VAR_H_
#define _INPUT_OM_VAR_H_
#if defined(_MSC_VER)
#pragma once
#endif

/*! 
 * \file input_om_fixed.h
 * \ingroup Objects
 * \brief InputOMVar class header file.
 * \author Josh Lurz
 */

#include <string>
#include <vector>
#include <xercesc/dom/DOMNode.hpp>
#include "functions/include/minicam_input.h"
#include "util/base/include/value.h"

class Tabs;

/*! 
 * \ingroup Objects
 * \brief Defines a single variable O&M input to a MiniCAM production function.
 * \details A variable O&M input adds to the marginal cost input in the
 *          production function.  The input does not have an associated good in 
 *          the model, and so the quantity is not tracked.
 *
 *          <b>XML specification for InputOMVar</b>
 *          - XML name: \c input-OM-fixed
 *          - Contained by: Technology
 *          - Parsing inherited from class: MiniCAMInput
 *          - Attributes:
 *              - \c name MiniCAMInput::mName
 *          - Elements:
 *              - \c OM-var InputOMVar::mOMVar
 *              - \c tech-change InputOMVar::mTechChange
 *
 * \author Josh Lurz
 */
class InputOMVar: public MiniCAMInput
{
    friend class InputFactory;
public:

    static const std::string& getXMLNameStatic();

    const std::string& getXMLReportingName() const;    

    InputOMVar* clone() const;

    virtual void XMLParse( const xercesc::DOMNode* aNode );

    virtual bool isSameType( const std::string& aType ) const;
    

    virtual void toInputXML( std::ostream& aOut,
                             Tabs* aTabs ) const;

    virtual void toDebugXML( const int aPeriod,
                             std::ostream& aOut,
                             Tabs* aTabs ) const;

    virtual void copyParam( const IInput* aInput,
                            const int aPeriod );

    virtual void copyParamsInto( InputOMVar& aInput,
                                const int aPeriod ) const;

    virtual void completeInit( const std::string& aRegionName,
                               const std::string& aSectorName,
                               const std::string& aSubsectorName,
                               const std::string& aTechName,
                               DependencyFinder* aDependencyFinder,
                               const IInfo* aTechInfo );
    
    virtual void initCalc( const std::string& aRegionName,
                           const std::string& aSectorName,
                           const bool aIsNewInvestmentPeriod,
                           const bool aIsTrade,
                           const int aPeriod );

    double getPhysicalDemand( const int aPeriod ) const;

    void setPhysicalDemand( const double aPhysicalDemand,
                            const std::string& aRegionName,
                            const int aPeriod );

    double getPrice( const std::string& aRegionName,
                     const int aPeriod ) const;
    
    virtual void setPrice( const std::string& aRegionName,
                           const double aPrice,
                           const int aPeriod );

    double getCO2EmissionsCoefficient( const std::string& aGHGName,
                                       const int aPeriod ) const;

    double getCoefficient( const int aPeriod ) const;

    void setCoefficient( const double aCoefficient,
                         const int aPeriod );

    void tabulateFixedQuantity( const std::string& aRegionName,
                                const double aFixedOutput,
                                const bool aIsInvestmentPeriod,
                                const int aPeriod );

    virtual void scaleCalibrationQuantity( const double aScaleFactor );

    virtual double getCalibrationQuantity( const int aPeriod ) const;

    bool hasTypeFlag( const int aTypeFlag ) const;

    virtual double getIncomeElasticity() const;

    virtual double getPriceElasticity() const;

    virtual double getTechChange( const int aPeriod ) const;
    
    virtual void doInterpolations( const int aYear, const int aPreviousYear,
                                   const int aNextYear, const IInput* aPreviousInput,
                                   const IInput* aNextInput );

protected:

    InputOMVar();

    //! Cost of the capital input adjusted for the additional costs of the
    //! capture component.
    std::vector<Value> mAdjustedCosts;

    //! Coefficient for production or demand function. Coefficients are not
    // read in and are initialized to 1, but can increase over time with
    // technical change.
    std::vector<Value> mAdjustedCoefficients;

    //! Input specific technical change.
    Value mTechChange;

private:
    //! Variable O&M cost.
    Value mOMVar;

    // Convert to 1975$/GJ for now.
    double calcOMVarCost( void ) const;

    const static std::string XML_REPORTING_NAME; //!< tag name for reporting xml db 

};

#endif // _INPUT_OM_VAR_H_
