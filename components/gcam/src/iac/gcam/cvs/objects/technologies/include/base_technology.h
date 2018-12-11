#ifndef _BASE_TECHNOLOGY_H_
#define _BASE_TECHNOLOGY_H_
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
* \file base_technology.h
* \ingroup Objects
* \brief BaseTechnology class header file.
* \author Pralit Patel
* \author Sonny Kim
*/

#include <string>
#include <vector>
#include <memory>
#include <xercesc/dom/DOMNode.hpp>
#include <cassert>
#include <map>

#include "technologies/include/expenditure.h"
#include "investment/include/iinvestable.h"
#include "util/base/include/iround_trippable.h"
#include "util/base/include/value.h"

class IInput;
class IOutput;
class IFunction;
class MoreSectorInfo;
class Demographic;
class NationalAccount;
class Consumer;
class HouseholdConsumer;
class GovtConsumer;
class TradeConsumer;
class InvestConsumer;
class ProductionTechnology;
class InvestConsumer;
class IVisitor;
class AGHG;
class IExpectedProfitRateCalculator;
class TechnologyType;
class INestedInput;
class ICaptureComponent;

/*! 
 * \ingroup Objects
 * \brief The common base class of ProductionTechnology and Consumer.
 * \details Abstract base class for an SGM technology.  Subclasses will be responsible
 *          for most functionality since consumers and production technologies could be
 *          very different.  For instance consumers would not be invested in and they
 *          would also have not need for vintaging.
 *
 *          <b>XML specification for BaseTechnology</b>
 *          - XML name: \c N/A (abstract base class)
 *          - Contained by: Subsector
 *          - Parsing inherited from class: None
 *          - Attributes: \c name BaseTechnology::name
 *                           The name of the technology.
 *                        \c year BaseTechnology::year
 *                           The vintage year of this technology (i.e. the first year
 *                           this technology will operate).
 *          - Elements:
 *              - \c SGMOutput::getXMLReportingNameStatic() BaseTechnology::mOutputs
 *                   Parse and output object (used to parse output for pre-base year
 *                   vintage technologies).  Note if non are parsed a default SGMOutput
 *                   will be created.
 *              - \c categoryName BaseTechnology::categoryName
 *                   Category name for reporting.
 *              - \c GHGFactory::isGHGNode( nodeName ) BaseTechnology::mGhgs
 *                   Create and/or parse a GHG object.  Note if non are parsed a
 *                   default CO2 emissions is created.
 *              - \c share-weight BaseTechnology::mShareWeight
 *                   A share weight to use to bais a logit function when investing.
 *              - \c initial-year BaseTechnology::mIsInitialYear
 *                   A flag to signal if this technology will be the first technology
 *                   of it's kind to operate.  This flag is useful so that we can have
 *                   it calculate it's coefficients in the base year.
 *              - \c CaptureComponentFactory::isOfType( nodeName ) BaseTechnology::mSequestrationDevice
 *                   Create and/or parse a ICaptureComponent object which will be used
 *                   to capture emissions.
 *              - \c NodeInput::getXMLNameStatic() BaseTechnology::mNestedInputRoot
 *                   Create and/or parse the root of a nested input structure.  Note that
 *                   there may only be one root.
 *
 * \author Pralit Patel, Sonny Kim
 */

class BaseTechnology: public IInvestable, public IRoundTrippable
{
    friend class XMLDBOutputter;
    friend class InvestableCounterVisitor;
    friend class SetShareWeightVisitor;
public:
    BaseTechnology();
    BaseTechnology( const BaseTechnology& baseTechIn );
    BaseTechnology& operator= (const BaseTechnology& baseTechIn );
    virtual BaseTechnology* clone() const = 0;
    virtual ~BaseTechnology();

    virtual void copyParam( const BaseTechnology* baseTechIn,
                            const int aPeriod ) = 0;

    virtual void copyParamsInto( GovtConsumer& govtConsumerIn,
                                 const int aPeriod ) const { assert( false ); }

    virtual void copyParamsInto( TradeConsumer& tradeConsumerIn,
                                 const int aPeriod ) const { assert( false ); }

    virtual void copyParamsInto( InvestConsumer& investConsumerIn,
                                 const int aPeriod ) const { assert( false ); }

    virtual void copyParamsInto( ProductionTechnology& prodTechIn,
                                 const int aPeriod ) const { assert( false ); }

    virtual void copyParamsInto( HouseholdConsumer& householdConsumerIn,
                                 const int aPeriod ) const { assert( false ); }

    virtual void copyParamsInto( Consumer& consumerIn,
                                 const int aPeriod ) const { assert( false ); }

    void XMLParse( const xercesc::DOMNode* node );

    void toInputXML( std::ostream& out, Tabs* tabs ) const;
    void toDebugXML( const int period, std::ostream& out, Tabs* tabs ) const;
    
    virtual void completeInit( const std::string& aRegionName,
                               const std::string& aSectorName,
                               const std::string& aSubsectorName ) = 0;
    
    virtual void initCalc( const MoreSectorInfo* aMoreSectorInfo,
                           const std::string& aRegionName, 
                           const std::string& aSectorName,
                           NationalAccount& nationalAccount, 
                           const Demographic* aDemographics,
                           const double aCapitalStock,
                           const int aPeriod ) = 0;

    virtual bool isNewInvestment( const int period ) const;

    const std::string& getName() const;
    int getYear() const;
    void setYear( int newYear );

    virtual void operate( NationalAccount& aNationalAccount, const Demographic* aDemographics, 
        const MoreSectorInfo* aMoreSectorInfo, const std::string& aRegionName, 
        const std::string& aSectorName, const bool isNewVintageMode, const int aPeriod ) = 0;

    virtual double setInvestment( const std::string& aRegionName, const double aAnnualInvestment,
                                 const double aTotalInvestment, const int aPeriod ) = 0;

    virtual double getOutput( const int aPeriod ) const;
    virtual double getCapitalStock() const = 0;
    void calcPricePaid( const MoreSectorInfo* aMoreSectorInfo, const std::string& aRegionName,
        const std::string& aSectorName, const int aPeriod, const int aLifetimeYears ) const;

    virtual double getAnnualInvestment( const int aPeriod ) const = 0;
    
    virtual double getExpectedProfitRate( const NationalAccount& aNationalAccount,
                                          const std::string& aRegionName,
                                          const std::string& aSectorName,
                                          const IExpectedProfitRateCalculator* aExpProfitRateCalc,
                                          const double aInvestmentLogitExp,
                                          const bool aIsShareCalc,
                                          const bool aIsDistributing,
                                          const int aPeriod ) const = 0;
    
    virtual double getCapitalOutputRatio( const IDistributor* aDistributor,
                                          const IExpectedProfitRateCalculator* aExpProfitRateCalc,
                                          const NationalAccount& aNationalAccount,
                                          const std::string& aRegionName,
                                          const std::string& aSectorName, 
                                          const int aPeriod ) const = 0;

    virtual double getFixedInvestment( const int aPeriod ) const = 0;
    virtual double distributeInvestment( const IDistributor* aDistributor,
                                         NationalAccount& aNationalAccount,
                                         const IExpectedProfitRateCalculator* aExpProfitRateCalc,
                                         const std::string& aRegionName,
                                         const std::string& aSectorName,
                                         const double aNewInvestment,
                                         const int aPeriod ) = 0;
    virtual void updateMarketplace( const std::string& sectorName, const std::string& regionName,
                                    const int period ) = 0;
    virtual void csvSGMOutputFile( std::ostream& aFile, const int period ) const = 0;
    virtual void accept( IVisitor* aVisitor, const int aPeriod ) const = 0;

    const std::string getIdentifier() const;
    static const std::string createIdentifier( const std::string& aName, int aYear );
    void removeEmptyInputs();
    virtual void setTypeHelper( TechnologyType* aTechType ) = 0;
    virtual void postCalc( const std::string& aRegionName, const std::string& aSectorName,
                           const int aPeriod ) = 0;
    bool hasCalibrationMarket() const;
    virtual double getShareWeight( const int aPeriod ) const;
    bool isInitialYear() const;
protected:
    virtual bool isCoefBased() const = 0;
    
    virtual bool isTrade() const {
        return false;
    };

    virtual bool XMLDerivedClassParse( const std::string& nodeName, const xercesc::DOMNode* curr ) = 0;
    virtual const std::string& getXMLName() const = 0;
    virtual void toInputXMLDerived( std::ostream& out, Tabs* tabs ) const = 0;
    virtual void toDebugXMLDerived( const int period, std::ostream& out, Tabs* tabs ) const = 0;

    std::string name; //!< Name
    std::string categoryName; //!< Category name, used for reporting
    int year; //!< Year the technology will initially operate.
    std::map<std::string,int> mGhgNameMap; //!< Mapping of ghg name to number.

    //! A temporary list of the leaf inputs from mNestedInputRoot
    std::vector<IInput*> mLeafInputs;
    std::vector<AGHG*> mGhgs; //!< Green-House gases.
    
    std::vector<IOutput*> mOutputs; //!< Outputs

    std::vector<Expenditure> expenditures; //!< Keep track of expenditures

    //! The share weight used to bias investment
    double mShareWeight;
    
    //! The root of the nested inputs
    // TODO: should I make this an auto_ptr
    INestedInput* mNestedInputRoot;

    //! A flag to determine if this technology is the first of it's kind
    bool mIsInitialYear;

    //! Flag to help reduce calls to calcPricePaid
    mutable bool mPricePaidCached;

    //! A sequestration device to capture emissions
    std::auto_ptr<ICaptureComponent> mSequestrationDevice;
private:
    void clear();
    void copy( const BaseTechnology& baseTechIn );
    bool doCalibration;
};

#endif // _BASE_TECHNOLOGY_H_
