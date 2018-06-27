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

#ifndef _MINICAM_INPUT_H_
#define _MINICAM_INPUT_H_
#if defined(_MSC_VER)
#pragma once
#endif

/*! 
 * \file minicam_input.h
 * \ingroup Objects
 * \brief MiniCAMInput class header file.
 * \author Josh Lurz
 */

#include <string>
#include <map>
#include <xercesc/dom/DOMNode.hpp>
#include "functions/include/iinput.h"

class Tabs;
class ICaptureComponent;

/*! 
 * \ingroup Objects
 * \brief Defines a single input to a MiniCAM production function.
 * \details This is the abstract base class for inputs to a MiniCAM production
 *          function. This class has stub implementations for functions that
 *          having meaning in SGM but not yet in MiniCAM. This class contains no
 *          actual functionality itself.
 *
 *          <b>XML specification for MiniCAMInput</b>
 *          - XML name: None.
 *          - Contained by: Technology
 *          - Parsing inherited from class: None.
 *          - Attributes: None.
 *          - Elements: None.
 *
 * \author Josh Lurz
 */
class MiniCAMInput: public IInput
{
friend class XMLDBOutputter;
public:
    virtual ~MiniCAMInput();
    virtual MiniCAMInput* clone() const = 0;

    virtual void copyParam( const IInput* aInput,
                            const int aPeriod ) = 0;

    virtual const std::string& getXMLReportingName() const = 0;

    virtual void XMLParse( const xercesc::DOMNode* aNode ) = 0;

    virtual bool isSameType( const std::string& aType ) const = 0;
    
    virtual void completeInit( const std::string& aRegionName,
                               const std::string& aSectorName,
                               const std::string& aSubsectorName,
                               const std::string& aTechName,
                               DependencyFinder* aDependencyFinder,
                               const IInfo* aTechInfo ) = 0;

    virtual void initCalc( const std::string& aRegionName,
                           const std::string& aSectorName,
                           const bool aIsNewInvestmentPeriod,
                           const bool aIsTrade,
                           const int aPeriod ) = 0;

    virtual void toInputXML( std::ostream& aOut,
                             Tabs* aTabs ) const = 0;

    virtual void toDebugXML( const int aPeriod,
                             std::ostream& aOut,
                             Tabs* aTabs ) const = 0;

    virtual const std::string& getName() const;

    virtual double getConversionFactor( const int aPeriod ) const;

    virtual double getCO2EmissionsCoefficient( const std::string& aGHGName,
                                            const int aPeriod ) const = 0;

    virtual double getCurrencyDemand( const int aPeriod ) const;
    
    virtual double getTechChange( const int aPeriod ) const = 0;

    virtual double getPrice( const std::string& aRegionName,
                             const int aPeriod ) const = 0;

    virtual double getPricePaid( const std::string& aRegionName,
                                 const int aPeriod ) const;

    virtual void setPricePaid( const double aPricePaid,
                               const int aPeriod );

    virtual double getPriceReceived( const std::string& aRegionName,
                                     const int aPeriod ) const;

    virtual void setCurrencyDemand( const double aCurrencyDemand,
                                    const std::string& aRegionName,
                                    const int aPeriod );

    virtual double getPhysicalDemand( const int aPeriod ) const = 0;

    virtual double getCarbonContent( const int aPeriod ) const;

    virtual void setPhysicalDemand( const double aPhysicalDemand,
                                    const std::string& aRegionName,
                                    const int aPeriod ) = 0;

    virtual double getCoefficient( const int aPeriod ) const = 0;

    virtual void setCoefficient( const double aCoefficient,
                                 const int aPeriod ) = 0;

    virtual double getPriceAdjustment() const;
    
    virtual bool hasTypeFlag( const int aTypeFlag ) const = 0;
    
    virtual double getPriceElasticity() const = 0;
    
    virtual double getIncomeElasticity() const = 0;

    virtual double getCalibrationQuantity( const int aPeriod ) const = 0;

    virtual void csvSGMOutputFile( std::ostream& aFile, const int period ) const;
    
    virtual void doInterpolations( const int aYear, const int aPreviousYear,
                                   const int aNextYear, const IInput* aPreviousInput,
                                   const IInput* aNextInput );

    virtual void accept( IVisitor* aVisitor, const int aPeriod ) const;

    virtual void copyParamsInto( ProductionInput& aInput,
                                 const int aPeriod ) const;

    virtual void copyParamsInto( DemandInput& aInput,
                                 const int aPeriod ) const;

    virtual void copyParamsInto( TradeInput& aInput,
                                 const int aPeriod ) const;

    virtual void copyParamsInto( NodeInput& aInput,
                                 const int aPeriod ) const;

    virtual void copyParamsInto( EnergyInput& aInput,
                                 const int aPeriod ) const {}

    virtual void copyParamsInto( NonEnergyInput& aInput,
                                 const int aPeriod ) const {}

    virtual void copyParamsInto( BuildingDemandInput& aInput,
                                 const int aPeriod ) const {}

    virtual void copyParamsInto( RenewableInput& aInput,
                                 const int aPeriod ) const {}

    virtual void copyParamsInto( InputSubsidy& aInput,
                                 const int aPeriod ) const {}

    virtual void copyParamsInto( InputTax& aInput,
                                 const int aPeriod ) const {}
	
	virtual void copyParamsInto( InputOMVar& aInput,
								const int aPeriod ) const {}
	
	virtual void copyParamsInto( InputOMFixed& aInput,
								const int aPeriod ) const {}
	
	virtual void copyParamsInto( InputCapital& aInput,
								const int aPeriod ) const {}

    virtual void calcPricePaid( const std::string& aRegionName,
                                const std::string& aSectorName,
                                const MoreSectorInfo* aMoreSectorInfo,
                                const std::vector<AGHG*>& aGhgs,
                                const ICaptureComponent* aSequestrationDevice,
                                const int aLifetimeYears,
                                const int aPeriod ) {}

    virtual double calcTaxes( const std::string& aRegionName,
                            NationalAccount* aNationalAccount,
                            Expenditure* aExpenditure,
                            const int aPeriod ) const { return 0; }
protected:
    MiniCAMInput();

    //! Name of the Input.
    std::string mName;
    
    //! A map of a keyword to its keyword group
    std::map<std::string, std::string> mKeywordMap;
};

#endif // _MINICAM_INPUT_H_
