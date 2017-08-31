#ifndef _TRADE_INPUT_H_
#define _TRADE_INPUT_H_
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
* \file trade_input.h
* \ingroup Objects
* \brief TradeInput class header file.
* \author Pralit Patel
*/

#include <string>
#include <cassert>
#include <xercesc/dom/DOMNode.hpp>

#include "functions/include/sgm_input.h"

class Tabs;

/*! 
 * \ingroup Objects
 * \brief Defines an input that is imported from a trading partner.
 * \details Extends an SGMInput to add the functionality that a good would not
 *          be demanded with in the same region however from another region
 *          which is declared through the XML parse.  This input also applies
 *          two additonal taxes one for importing which goes to the region that
 *          contains this input and an export tax which goes to the foriegn trading
 *          partner.  Note transportation costs are not accounted for with in this
 *          input.
 *
 *          <b>XML specification for TradeInput</b>
 *          - XML name: \c trade-input
 *          - Contained by: NodeInput
 *          - Parsing inherited from class: SGMInput
 *          - Elements:
 *              - \c trade-partner TradeInput::mTradingPartner
 *                   The name of the region this good is imported from.
 *              - \c import-tax-rate TradeInput::mImportTaxRate
 *                   The import tax rate price adjustment.
 *              - \c export-tax-rate TradeInput::mExportTaxRate
 *                   The export tax rate price adjustment.
 *
 * \author Pralit Patel
 */
class TradeInput : public SGMInput
{
public:
    TradeInput();
    TradeInput( const TradeInput& aTradeInput );
    virtual ~TradeInput();
    virtual TradeInput* clone() const ;
    virtual void copyParam( const IInput* aInput, const int aPeriod );
    virtual void copyParamsInto( TradeInput& aTradeInput, const int aPeriod ) const;

    virtual void initCalc( const std::string& aRegionName,
                           const std::string& aSectorName,
                           const bool aIsNewInvestmentPeriod,
                           const bool aIsTrade,
                           const int aPeriod );
    virtual void setPhysicalDemand( const double aPhysicalDemand, const std::string& aRegionName, const int aPeriod );
    virtual double getPrice( const std::string& aRegionName, const int aPeriod ) const;
    virtual void calcPricePaid( const std::string& aRegionName,
                                const std::string& aSectorName,
                                const MoreSectorInfo* aMoreSectorInfo,
                                const std::vector<AGHG*>& aGhgs,
                                const ICaptureComponent* aSequestrationDevice,
                                const int aLifetimeYears,
                                const int aPeriod );
    virtual double calcTaxes( const std::string& aRegionName,
                            NationalAccount* aNationalAccount,
                            Expenditure* aExpenditure,
                            const int aPeriod ) const;
    virtual double getPriceReceived( const std::string& aRegionName, const int aPeriod ) const;
    
    static const std::string& getXMLNameStatic();
    virtual const std::string& getXMLReportingName() const;
    virtual void accept( IVisitor* aVisitor, const int period ) const;
    
    // IInput methods that are not implemented
    virtual double getPriceElasticity() const { return 0; }
    virtual double getIncomeElasticity() const { return 0; }
    virtual void copyParamsInto( ProductionInput& aProductionInput, const int aPeriod ) const {}
    virtual void copyParamsInto( DemandInput& aDemandInput, const int aPeriod ) const {}
protected:
    const std::string& getXMLName() const;
    virtual bool XMLDerivedClassParse( const std::string& aNodeName, const xercesc::DOMNode* aCurrNode );
    virtual void toInputXMLDerived( std::ostream& aOut, Tabs* aTabs ) const;
    virtual void toDebugXMLDerived( const int aPeriod, std::ostream& aOut, Tabs* aTabs ) const;
private:
    //! The tax rate at which an importer has to pay taxes to the local government
    Value mImportTaxRate;

    //! The rate at which an exporter taxes the good.  Note that the export tax
    //! amount goes back to the foreign producer of the good, and the foreign government
    //! pays for that tax/subsidy.
    Value mExportTaxRate;

    //! The region from whom this good is imported from.
    std::string mTradingPartner;
};

#endif // _TRADE_INPUT_H_
