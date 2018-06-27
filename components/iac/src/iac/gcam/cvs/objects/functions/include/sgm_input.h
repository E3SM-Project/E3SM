#ifndef _SGM_INPUT_H_
#define _SGM_INPUT_H_
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
* \file SGMInput.h
* \ingroup Objects
* \brief SGMInput class header file.
* \author Pralit Patel, Sonny Kim
*/

#include <string>
#include <vector>
#include <xercesc/dom/DOMNode.hpp>
#include "util/base/include/value.h"
#include "util/base/include/time_vector.h"
#include "functions/include/iinput.h"
#include "functions/include/inested_input.h"

class Tabs;
class DemandInput;
class ProductionInput;
class CachedMarket;

/*! 
 * \ingroup Objects
 * \brief Defines a single SGMInput to a production or demand function.
 * \details This input keeps track of the usual coefficients for a production function,
 *          price paid, and physical demands.  It also has a sales tax price adjustmnet,
 *          and will set type flags to distiguish between energy, material, land, labor,
 *          and capital mostly to help the accounting structures.  An SGMInput is also
 *           considered a nested input so that they can be included in a nested input 
 *          structure however all the INestedInput methods essentially do nothing since
 *          SGMInputs are intended to be leaves in the nesting structure.
 *
 *          <b>XML specification for SGMInput</b>
 *          - XML name: \c N/A (abstract base class)
 *          - Contained by: NodeInput
 *          - Parsing inherited from class: None
 *          - Attributes: \c name SGMInput::mName
 *                           This is a name which should correspond to a market name.
 *                           Note that certain names will imply that a flag type be
 *                           be set for this input see SGMInput::setFlagsByName(string).
 *          - Elements:
 *              - \c coefficient SGMInput::mCoefficient
 *                   The base year coefficient to be used with a production function.
 *              - \c demandCurrency SGMInput::mPhysicalDemand
 *                   If no year attribute is specified this is the base year currency
 *                   demand otherwise it will be placed with in the vector for the 
 *                   corresponding period.  Note that even though this is called currency
 *                   demand it is being placed into the physical demand vector.  This is
 *                   a hack which was necessary because we use currency demands to calibrate
 *                   CES coefficients in init calc of the  base year but for the rest of the
 *                   model these values will really be physical.
 *              - \c priceAdjustFactor SGMInput::mPriceAdjustFactor
 *                   The a constant price adjustment.  Note that this is no longer used.
 *              - \c technicalChange SGMInput::mTechnicalChange
 *                   An input specific technical change paramater which would only get
 *                   applied for the initial year that this input was operated.
 *              - \c flag SGMInput::mTypeFlags
 *                   - Attributes: \c type
 *                                    Attempt to set an input type flag directly rather than
 *                                    relying on hard coded market names.  This will rely on
 *                                    SGMInput::setFlagsByName(string) see that method for 
 *                                    details.
 *              - \c sales-tax-rate SGMInput::mSalesTaxRate
 *                   A price adjustment to account for sales taxes, see SGMInput::calcPricePaid
 *                   and SGMInput::calcTaxes for details.
 *
 * \author Pralit Patel, Sonny Kim, Josh Lurz
 */
class SGMInput: public INestedInput
{
    friend class SocialAccountingMatrix;
    friend class DemandComponentsTable;
    friend class SectorReport;
    friend class SGMGenTable;
    friend class XMLDBOutputter;
public:
    SGMInput();
    SGMInput( const SGMInput& aInput );
    virtual ~SGMInput();
    // TODO: shouldn't this be IInput* clone()?
    virtual SGMInput* clone() const = 0;
    virtual void copyParam( const IInput* aInput, const int aPeriod ) = 0;

    virtual const std::string& getXMLReportingName() const;

    void XMLParse( const xercesc::DOMNode* node );
    bool isSameType( const std::string& aType ) const;

    void completeInit( const std::string& aRegionName,
                       const std::string& aSectorName,
                       const std::string& aSubsectorName,
                       const std::string& aTechName,
                       DependencyFinder* aDependencyFinder ,
                       const IInfo* aTechInfo);

    virtual void initCalc( const std::string& aRegionName,
                           const std::string& aSectorName,
                           const bool aIsNewInvestmentPeriod,
                           const bool aIsTrade,
                           const int aPeriod );

    void toInputXML( std::ostream& out, Tabs* tabs ) const;
    void toDebugXML( const int period, std::ostream& out, Tabs* tabs ) const;
    const std::string& getName() const;
    double getConversionFactor( const int aPeriod ) const;
    double getCO2EmissionsCoefficient( const std::string& aGHGName, const int aPeriod ) const;
    double getPhysicalDemand( const int aPeriod ) const;
    double getCarbonContent( const int aPeriod ) const;
    double getCurrencyDemand( const int aPeriod ) const;
    double getTechChange( const int aPeriod ) const;

    double getPrice( const std::string& aRegionName, const int aPeriod ) const;

    virtual void setPrice( const std::string& aRegionName,
                           const double aPrice,
                           const int aPeriod );

    double getPricePaid( const std::string& aRegionName, const int aPeriod ) const;
    void setPricePaid( const double aPricePaid, const int aPeriod );
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

    double getPriceReceived( const std::string& aRegionName, const int aPeriod ) const;

    void setCurrencyDemand( const double aCurrencyDemand, const std::string& aRegionName, const int aPeriod );

    void setPhysicalDemand( const double aPhysicalDemand, const std::string& aRegionName, const int aPeriod );

    double getCoefficient( const int aPeriod ) const;

    void setCoefficient( const double aCoefficient, const int aPeriod );
    double getPriceAdjustment() const;
    bool hasTypeFlag( const int aTypeFlag ) const;
    virtual double getPriceElasticity() const = 0;
    virtual double getIncomeElasticity() const = 0;

    virtual double getCalibrationQuantity( const int aPeriod ) const;
    
    virtual void doInterpolations( const int aYear, const int aPreviousYear,
                                   const int aNextYear, const IInput* aPreviousInput,
                                   const IInput* aNextInput );

    void csvSGMOutputFile( std::ostream& aFile, const int period ) const;
    virtual void accept( IVisitor* aVistior, const int aPeriod ) const = 0;
    virtual void copyParamsInto( ProductionInput& aInput, const int aPeriod ) const = 0;
    virtual void copyParamsInto( DemandInput& aInput, const int aPeriod ) const = 0;

    // MiniCAM input types which SGM inputs should ignore.
    virtual void copyParamsInto( EnergyInput& aInput, const int aPeriod ) const {}
    virtual void copyParamsInto( NonEnergyInput& aInput, const int aPeriod ) const {}
    virtual void copyParamsInto( BuildingDemandInput& aInput, const int aPeriod ) const {}
    virtual void copyParamsInto( RenewableInput& aInput, const int aPeriod ) const {}
    virtual void copyParamsInto( InputSubsidy& aInput, const int aPeriod ) const {}
    virtual void copyParamsInto( InputTax& aInput, const int aPeriod ) const {}
    virtual void copyParamsInto( NodeInput& aInput, const int aPeriod ) const {}
    virtual void copyParamsInto( TradeInput& aInput, const int aPeriod ) const {}
	virtual void copyParamsInto( InputOMVar& aInput, const int aPeriod ) const {}
    virtual void copyParamsInto( InputOMFixed& aInput, const int aPeriod ) const {}
    virtual void copyParamsInto( InputCapital& aInput, const int aPeriod ) const {}


    // INestedInput methods
    // define them to do nothing since an SGMInput is a leaf in the nesting structure
    // this should be the end point for recursion
    virtual void removeEmptyInputs() {}

    virtual void initialize() {}

    virtual void calcCoefficient( const std::string& aRegionName, const std::string& aSectorName,
        const int aTechPeriod ) {}

    virtual void changeElasticity( const std::string& aRegionName, const int aPeriod,
        const double aAlphaZero ) {}

    virtual void changeSigma( const std::string& aRegionName, const int aPeriod,
        const double aAlphaZero ) {}

    virtual void calcLevelizedCost( const std::string& aRegionName, const std::string& aSectorName,
        const int aPeriod, const double aAlphaZero ) {}

    virtual double calcInputDemand( const std::string& aRegionName, const std::string& aSectorName,
        const int aPeriod, const double aPhysicalOutput, const double aUtilityParameterA,
        const double aAlphaZero ) { return 0; }

    virtual double calcCapitalOutputRatio( const std::string& aRegionName, const std::string& aSectorName,
        const int aPeriod, const double aAlphaZero ) { return 1.0; }

    virtual void calcVariableLevelizedCost( const std::string& aRegionName, const std::string& aSectorName,
        const int aPeriod, const double aAlphaZero ) {}

    virtual const IFunction* getFunction() const { return 0; }
    
    virtual double getLevelizedCost( const std::string& aRegionName, const std::string& aSectorName,
        const int aPeriod ) const { return 0; }

    virtual void applyTechnicalChange( const std::string& aRegionName, const std::string& aSectorName,
        const int aPeriod, const TechChange& aTechChange ) {}

    virtual void resetCalcLevelizedCostFlag() {}
protected:
    //! Name of the SGMInput.
    std::string mName;

    //! Coefficient for production or demand function.
    objects::PeriodVector<Value> mCoefficient;

    //! Currency Demand.
    objects::PeriodVector<Value> mPhysicalDemand;

    //! Price adjustment factor.  It used to be utilized
    //! to implement price wedges for Capital.  It is no
    //! longer used and should likely be removed.
    Value mPriceAdjustFactor;

    //! Price paid for SGMInput, adjusted from market price.
    objects::PeriodVector<Value> mPricePaid;

    //! Technical Change.
    Value mTechnicalChange;

    //! Conversion factor.  Was used to convert from currency
    //! to energy however we now use prices to take care of this.
    //! This value is no longer used and should likely be removed.
    Value mConversionFactor;

    //! CO2 coefficient.
    Value mCO2Coefficient;

    //! Sales tax rate
    Value mSalesTaxRate;

    //! Type flags.
    int mTypeFlags;

    //! A pre-located market which has been cahced from the marketplace to get
    //! the price and add demands to.
    std::auto_ptr<CachedMarket> mCachedMarket;

    void initializeTypeFlags( const std::string& aRegionName );
    void initializeCachedCoefficients( const std::string& aRegionName );
    void setFlagsByName( const std::string& aTypeName );

    bool isEnergyGood( const std::string& aRegionName ) const;

    bool isPrimaryEnergyGood( const std::string& aRegionName ) const;

    bool isSecondaryEnergyGood( const std::string& aRegionName ) const;

    virtual const std::string& getXMLName() const = 0;
    virtual bool XMLDerivedClassParse( const std::string& nodeName, const xercesc::DOMNode* curr ) = 0;
    virtual void toInputXMLDerived( std::ostream& out, Tabs* tabs ) const = 0;
    virtual void toDebugXMLDerived( const int period, std::ostream& out, Tabs* tabs ) const = 0;

private:
    const static std::string XML_REPORTING_NAME; //!< tag name for reporting xml db 

};

#endif // _SGM_INPUT_H_
