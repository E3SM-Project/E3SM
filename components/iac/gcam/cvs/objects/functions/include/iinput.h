#ifndef _IINPUT_H_
#define _IINPUT_H_
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
 * \file iinput.h
 * \ingroup Objects
 * \brief IInput interface header file.
 * \author Josh Lurz
 */

#include <string>
#include <xercesc/dom/DOMNode.hpp>
#include <iosfwd> // remove when csv output is removed.

class Tabs;
class DependencyFinder;
class ICaptureComponent;
class IInfo;
class MoreSectorInfo;
class AGHG;
class ICaptureComponent;
class NationalAccount;
class Expenditure;

// Until copyParam is fixed.
class DemandInput;
class ProductionInput;
class NodeInput;
class TradeInput;
class BuildingDemandInput;
class EnergyInput;
class NonEnergyInput;
class RenewableInput;
class InputSubsidy;
class InputTax;
class InputOMVar;
class InputOMFixed;
class InputCapital;

#include "util/base/include/ivisitable.h"

/*! 
 * \ingroup Objects
 * \brief Represents a single generic input to a production function.
 * \details
 * \author Josh Lurz
 */
class IInput: public IVisitable { 
public:
    /*!
     * \brief Define different type attributes of inputs. These are not mutually
     *        exclusive.
     * \details The types are represented as bits of an integer to allow testing
     *          of multiple flags at once. For instance, to test if an input has
     *          both the ENERGY and FACTOR flags, the function hasTypeFlag is
     *          called as:
     *
     *          hasTypeFlag( IInput::ENERGY | IInput::FACTOR )
     *
     *          To add additional flags simply increment the bit-shift by one.
     */
    enum Type {
        //! Energy.
        ENERGY = 1 << 0,

        //! Material.
        MATERIAL = 1 << 1,

        //! Factor supply.
        FACTOR = 1 << 2,

        //! Land.
        LAND = 1 << 3,

        //! Labor.
        LABOR = 1 << 4,

        //! Capital.
        CAPITAL = 1 << 5,

        //! Primary energy.
        PRIMARY = 1 << 6,

        //! Secondary energy.
        SECONDARY = 1 << 7,

        //! Numeraire.
        NUMERAIRE = 1 << 8,

        //! Initialized
        INITIALIZED = 1 << 9,
        
        //! Subsidy.
        SUBSIDY = 1 << 10,

        //! Tax.
        TAX = 1 << 11,

        //! Traded Good.
        TRADED = 1 << 12,
		
		//! O&M Input
        OM_VAR = 1 << 13,
		
		//! O&M Input
        OM_FIXED = 1 << 14

    };

    /*!
     * \brief Constructor.
     * \details Inlined constructor to avoid compiler problems with abstract
     *          base classes. 
     */
    IInput();

    /*!
     * \brief Destructor.
     * \details Inlined destructor to avoid compiler problems with abstract base
     *          classes. 
     */
    virtual ~IInput();
    
    /*!
     * \brief Creates an exact copy of the input.
     * \return An exact copy of the capture input. 
     */
    virtual IInput* clone() const = 0;
    
    /*!
     * \brief Copy parameters from another input.
     * \param aInput An input from which to copy.
     * \param aPeriod Period in which the input is being copied.
     */
    virtual void copyParam( const IInput* aInput,
                            const int aPeriod ) = 0;

    /*!
     * \brief Returns whether the type of the object is the same as the passed
     *        in type.
     * \param aType Type to check the object's type against.
     * \return Whether the type of the object is the same as the passed in type.
     */
    virtual bool isSameType( const std::string& aType ) const = 0;
    
    /*!
     * \brief Return the name of the input.
     * \return The name of the input.
     */
    virtual const std::string& getName() const = 0;

    /*!
     * \brief Return the name of the input for reporting.
     * \return The name of the input for reporting.
     */
    virtual const std::string& getXMLReportingName() const = 0;

    /*!
     * \brief Parse the data for this object starting at a given node.
     * \param aNode Root node from which to parse data.
     */
    virtual void XMLParse( const xercesc::DOMNode* aNode ) = 0;
    
    /*!
     * \brief Write data from this object in an XML format so that it can be
     *        read back in later as input.
     * \param aOut Filestream to which to write.
     * \param aTabs Object responsible for writing the correct number of tabs. 
     */
    virtual void toInputXML( std::ostream& aOut,
                             Tabs* aTabs ) const = 0;
    
    /*!
     * \brief Write data from this object in an XML format for debugging.
     * \param aPeriod Period for which to write data.
     * \param aOut Filestream to which to write.
     * \param aTabs Object responsible for writing the correct number of tabs. 
     */
    virtual void toDebugXML( const int aPeriod,
                             std::ostream& aOut,
                             Tabs* aTabs ) const = 0;
    
    /*!
     * \brief Returns whether the input has the specified type flag set.
     * \details
     * \param aTypeFlag A bit mask of flags containing all flags to check for.
     * \return Whether the specified type flag is set.
     */
    virtual bool hasTypeFlag( const int aTypeFlag ) const = 0;

    /*!
     * \brief Complete the initialization of the input.
     * \param aRegionName Name of the region containing the input.
     * \param aSectorName Name of the sector containing the input.
     * \param aSubsectorName Name of the subsector containing the input.
     * \param aTechName Name of the Technology containing the input.
     * \param aDependencyFinder The input dependency finder, which may be null.
     * \param aTechInfo Technology's info object.
     */
    virtual void completeInit( const std::string& aRegionName,
                               const std::string& aSectorName,
                               const std::string& aSubsectorName,
                               const std::string& aTechName,
                               DependencyFinder* aDependencyFinder,
                               const IInfo* aTechInfo ) = 0;

    /*!
     * \brief Initialize an input for a given period.
     * \param aRegionName Name of the containing region.
     * \param aSectorName Name of the containing sector.
     * \param aIsInvestmentPeriod Whether this is the initial investment period
     *        of the Technology.
     * \param aIsTrade Whether this is a trade technology.
     * \param aPeriod Model period.
     */
    virtual void initCalc( const std::string& aRegionName,
                           const std::string& aSectorName,
                           const bool aIsNewInvestmentPeriod,
                           const bool aIsTrade,
                           const int aPeriod ) = 0;

    /*!
     * \brief Get the currency demand for input used.
     * \param aPeriod Model period.
     * \return The currency demand for each input used.
     */
    virtual double getCurrencyDemand( const int aPeriod ) const = 0;

    /*!
     * \brief Set the currency demand for input used.
     * \param aCurrencyDemand Currency demand.
     * \param aRegionName Region name.
     * \param aPeriod Model period.
     */
    virtual void setCurrencyDemand( const double aCurrencyDemand,
                                    const std::string& aRegionName, 
                                    const int aPeriod ) = 0;

    /*!
     * \brief Get the physical demand for input used.
     * \param aPeriod Model period.
     * \return The physical demand for each input used.
     */
    virtual double getPhysicalDemand( const int aPeriod ) const = 0;
    
    /*!
     * \brief Get the carbon content of the input used.
     * \param aPeriod Model period.
     * \return The carbon content of each input used.
     */
    virtual double getCarbonContent( const int aPeriod ) const = 0;
    
    /*!
     * \brief Set the physical demand for input used.
     * \param aPhysicalDemand Currency demand.
     * \param aRegionName Region name.
     * \param aPeriod Model period.
     */
    virtual void setPhysicalDemand( const double aPhysicalDemand,
                                    const std::string& aRegionName, 
                                    const int aPeriod ) = 0;

    /*!
     * \brief Get the price of the input in a given period.
     * \param aRegionName Name of the region containing the input.
     * \param aPeriod Period for which to return price.
     * \return The price in the given period.
     */
    virtual double getPrice( const std::string& aRegionName,
                             const int aPeriod ) const = 0;

    /*!
     * \brief Set the price of an input in a given period.
     * \param aRegionName Name of the region containing the input.
     * \param aPrice The new price of the input.
     * \param aPeriod Model period.
     */
    virtual void setPrice( const std::string& aRegionName,
                           const double aPrice,
                           const int aPeriod ) = 0;

    /*!
     * \brief Get the price adjustment factor.
     * \details
     * \return The price adjustment factor.
     * \todo Remove this if possible.
     */
    virtual double getPriceAdjustment() const = 0;

    /*!
     * \brief Get the price paid of the input in a given period.
     * \param aRegionName Name of the region containing the input.
     * \param aPeriod Period for which to return price paid.
     * \return The price paid in the given period.
     */
    virtual double getPricePaid( const std::string& aRegionName,
                                 const int aPeriod ) const = 0;

    /*!
     * \brief Set the price paid of the input in a given period.
     * \param aPricePaid The price paid to set.
     * \param aPeriod Period for which to set price paid.
     */
    virtual void setPricePaid( const double aPricePaid,
                               const int aPeriod ) = 0;

    /*!
     * \brief Calculate the price paid of the input.
     * \param aRegionName Name of the region containing the input.
     * \param aSectorName Name of the containing sector.
     * \param aMoreSectorInfo The sector info which may contain additional costs.
     * \param aGhgs GHGs which may add to the cost of the input.
     * \param aSequestrationDevice A capture component which may capture some emssions
     *          and thus reduce emissions tax costs.
     * \param aLifetimeYears The number of years the technology will operate for.
     *          Used to calculate depreciation of capital.
     * \param aPeriod Period for which to calculate price paid.
     */
    virtual void calcPricePaid( const std::string& aRegionName,
                                const std::string& aSectorName,
                                const MoreSectorInfo* aMoreSectorInfo,
                                const std::vector<AGHG*>& aGhgs,
                                const ICaptureComponent* aSequestrationDevice,
                                const int aLifetimeYears,
                                const int aPeriod ) = 0;

    /*!
     * \brief Get the coefficient of the input.
     * \param aPeriod Model period.
     * \return Coefficient
     */
    virtual double getCoefficient( const int aPeriod ) const = 0;

    /*!
     * \brief Set the coefficient of the input.
     * \param aPeriod Model period.
     * \param aCoefficient The new coefficient.
     */
    virtual void setCoefficient( const double aCoefficient,
                                 const int aPeriod ) = 0;

    /*! \brief Get the conversion factor.
    * \param aPeriod Model period
    * \details TODO
    * \return The conversion factor.
    */
    virtual double getConversionFactor( const int aPeriod ) const = 0;

    /*! \brief Get the emissions coefficient of the input for a given gas.
    * \param aGHGName The name of the gas.
    * \param aPeriod Model period
    * \return The emissions coefficient for the gas.
    */
    virtual double getCO2EmissionsCoefficient( const std::string& aGHGName,
                                             const int aPeriod ) const = 0;

    /*!
     * \brief Calculate taxes from the input.
     * \details Calculates the taxes and places them into the appropriate
     *          accounting structure.
     * \param aRegionName Name of the region containing the input.
     * \param aNationalAccount The national account to add taxes into if avaiable.
     * \param aExpenditure The current period expenditure to track technology expenses if availabe.
     * \param aPeriod The period in which to calculate taxes.
     * \return The amount of non-emission taxes collected.
     */
    virtual double calcTaxes( const std::string& aRegionName,
                            NationalAccount* aNationalAccount,
                            Expenditure* aExpenditure,
                            const int aPeriod ) const = 0;

    /*!
     * \brief Get the current calibration quantity.
     * \param aPeriod The period for which to get the calibration quantity.
     * \details
     * \return The current calibration quantity.
     */
    virtual double getCalibrationQuantity( const int aPeriod ) const = 0;

    /*!
     * \brief Get the price elasticity of the input.
     * \return The price elasticity.
     */
    virtual double getPriceElasticity() const = 0;

    /*!
     * \brief Get the income elasticity of the input.
     * \return The income elasticity.
     */
    virtual double getIncomeElasticity() const = 0;

    /*!
     * \brief Get the input specific technical change.
     * \param aPeriod Model period.
     * \return The input specific technical change.
     * \author Josh Lurz
     */
    virtual double getTechChange( const int aPeriod ) const = 0;

    /*! \brief Write out the SGM csv output file.
    * \todo Remove this function.
    * \param aFile Output file.
    * \param aPeriod Model period.
    */
    virtual void csvSGMOutputFile( std::ostream& aFile,
                                   const int aPeriod ) const = 0;
    
    /*!
     * \brief Hook for an input to do interpolations to fill in any data that
     *        should be interpolated to a newly created input for the missing
     *        technology.
     * \param aYear the year to be filled in.
     * \param aPreviousYear The year of the last parsed input.
     * \param aNextYear The year of the next closest parsed input.
     * \param aPreviousInput The previous parsed input.
     * \param aNextInput The next parsed input.
     */
    virtual void doInterpolations( const int aYear, const int aPreviousYear,
                                   const int aNextYear, const IInput* aPreviousInput,
                                   const IInput* aNextInput ) = 0;
    
    virtual void copyParamsInto( ProductionInput& aInput,
                                 const int aPeriod ) const = 0;

    virtual void copyParamsInto( DemandInput& aInput,
                                 const int aPeriod ) const = 0;

    virtual void copyParamsInto( NodeInput& aInput,
                                 const int aPeriod ) const = 0;

    virtual void copyParamsInto( TradeInput& aInput,
                                 const int aPeriod ) const = 0;

    virtual void copyParamsInto( EnergyInput& aInput,
                                 const int aPeriod ) const = 0;

    virtual void copyParamsInto( NonEnergyInput& aInput,
                                 const int aPeriod ) const = 0;
	virtual void copyParamsInto( InputCapital& aInput,
								const int aPeriod ) const = 0;
	
    virtual void copyParamsInto( InputOMFixed& aInput,
								const int aPeriod ) const = 0;
	
    virtual void copyParamsInto( InputOMVar& aInput,
								const int aPeriod ) const = 0;

    virtual void copyParamsInto( BuildingDemandInput& aInput,
                                 const int aPeriod ) const = 0;

    virtual void copyParamsInto( RenewableInput& aInput,
                                 const int aPeriod ) const = 0;

    virtual void copyParamsInto( InputSubsidy& aInput,
                                 const int aPeriod ) const = 0;

    virtual void copyParamsInto( InputTax& aInput,
                                 const int aPeriod ) const = 0;
	
    // IVisitable interface.
    virtual void accept( IVisitor* aVisitor,
                        const int aPeriod ) const = 0;
};

// Inline function definitions.
inline IInput::IInput(){
}

inline IInput::~IInput(){
}

#endif // _IINPUT_H_
