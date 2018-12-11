#ifndef _IOUTPUT_H_
#define _IOUTPUT_H_
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
 * \file ioutput.h
 * \ingroup Objects
 * \brief IOutput interface header file.
 * \author Josh Lurz
 */

#include <string>
#include <xercesc/dom/DOMNode.hpp>
#include <list>

class Tabs;
class DependencyFinder;
class ICaptureComponent;
class IInfo;

#include "util/base/include/ivisitable.h"
#include "util/base/include/iparsable.h"
#include "util/base/include/iround_trippable.h"
#include "land_allocator/include/iland_allocator.h"
/*! 
* \ingroup Objects
* \brief Represents a single generic output of a Technology.
* \details The output interface represents a single output of a technology. All
*          MiniCAM technologies must have a primary output that implements this
*          interface. The primary output level is determined by the technology.
*          The output levels of other outputs may be determined by the scale of
*          the primary output, or through other means. Outputs may have positive
*          or negative monetary value, which is incorporated into the costs of
*          operating the technology. Outputs may also have associated carbon
*          content and emissions. These are accounted for in the emissions cost
*          and quantity calculations.
* \author Josh Lurz
*/
class IOutput : public IVisitable, public IParsable, public IRoundTrippable {
public:
    /*! 
     * \brief Constructor.
     * \details Inlined constructor to avoid compiler problems with abstract
     *          base classes. 
     */
    inline IOutput();

    /*!
     * \brief Destructor.
     * \details Inlined destructor to avoid compiler problems with abstract base
     *          classes and allow deletion through the base class pointer. 
     */
    inline virtual ~IOutput();

    /*!
     * \brief Creates an exact copy.
     * \return An exact copy. 
     */
    virtual IOutput* clone() const = 0;

    /*!
     * \brief Returns whether the type of the object is the same as the passed
     *        in type.
     * \param aType Type to check the object's type against.
     * \return Whether the type of the object is the same as the passed in type.
     */
    virtual bool isSameType( const std::string& aType ) const = 0;

    /*!
     * \brief Return the name of the output.
     * \return The name of the output.
     */
    virtual const std::string& getName() const = 0;

    /*!
     * \brief Return the name of the input for reporting.
     * \return The name of the input for reporting.
     */
    virtual const std::string& getXMLReportingName() const = 0;

    // Documentation is inherited.
    virtual bool XMLParse( const xercesc::DOMNode* aNode ) = 0;

    // Documentation is inherited.
    virtual void toInputXML( std::ostream& aOut,
                             Tabs* aTabs ) const = 0;

    // Documentation is inherited.
    virtual void toDebugXML( const int aPeriod,
                             std::ostream& aOut,
                             Tabs* aTabs ) const = 0;

    /*!
     * \brief Complete the initialization of the output.
     * \param aSectorName Sector name.
     * \param aDependencyFinder The input dependency finder, which may be null.
     * \param aTechInfo The technology's information container.
     * \param aIsTechOperating Whether the Technology can operate.
     */
    virtual void completeInit( const std::string& aSectorName,
                               DependencyFinder* aDependencyFinder,
                               const IInfo* aTechInfo,
                               const bool aIsTechOperating ) = 0;

    /*!
     * \brief Initialize an output for a given period.
     * \param aRegionName Name of the containing region.
     * \param aPeriod Model period.
     * \param aSectorName Sector name.
     */
    virtual void initCalc( const std::string& aRegionName,
                           const std::string& aSectorName,
                           const int aPeriod ) = 0;

    /*!
     * \brief Perform any final operations for a output in a given period.
     * \param aRegionName Name of the containing region.
     * \param aPeriod Model period.
     */
    virtual void postCalc( const std::string& aRegionName,
                           const int aPeriod ) = 0;

    /*! 
     * \brief Scale the output coefficient by a specified value.
     * \details Scales the output coefficient by a specified value. Outputs
     *          are not required to support this operation.
     * \note This is currently a workaround for byproducts not having
             independent coefficients.
     * \param aCoefficientScaler Coefficient scaler.
     */
    virtual void scaleCoefficient( const double aScaler ) = 0;

    /*!
     * \brief Calculate and return the physical output determined by the
     *        specified primary output of the Technology.
     * \details Determine the level of output from the known primary output of
     *          the Technology. Does not add the output to the marketplace.
     * \param aPrimaryOutput Output of the primary good.
     * \param aRegionName Region name.
     * \param aPeriod Period.
     * \return A list of pairs of output name and physical output given the
     *         level of primary output.
     * \sa setPhysicalOutput
     */
    typedef std::list<std::pair<std::string, double> > OutputList;
    virtual OutputList calcPhysicalOutput( const double aPrimaryOutput,
                                           const std::string& aRegionName,
                                           const ICaptureComponent* aCaptureComponent,
                                           const int aPeriod ) const = 0;

    /*!
     * \brief Set the physical output determined by the specified primary output
     *        of the Technology.
     * \details Determine the level of output from the known primary output of
     *          the Technology and add the output to the marketplace.
     * \param aPrimaryOutput Output of the primary good.
     * \param aRegionName Region name.
     * \param aPeriod Period.
     * \sa calcPhysicalOutput
     */
    virtual void setPhysicalOutput( const double aPrimaryOutput,
                                    const std::string& aRegionName,
                                    ICaptureComponent* aCaptureComponent,
                                    const int aPeriod ) = 0;

    /*!
     * \brief Get the quantity of physical output.
     * \param aPeriod Model period.
     * \return The physical output.
     */
    virtual double getPhysicalOutput( const int aPeriod ) const = 0;

    /*!
     * \brief Set the output in currency units.
     * \param aRegionName Region name.
     * \param aOutput Currency output.
     * \param aPeriod Period.
     */
    virtual void setCurrencyOutput( const std::string& aRegionName,
                                    const double aOutput,
                                    const int aPeriod ) = 0;

    /*!
     * \brief Get the output in currency units.
     * \param aPeriod Period.
     * \return Output in currency units.
     */
    virtual double getCurrencyOutput( const int aPeriod ) const = 0;

    /*!
     * \brief Get the monetary value of a single unit of output.
     * \param aRegionName Name of the region containing the output.
     * \param aPeriod Period.
     * \return The value in the given period.
     */
    virtual double getValue( const std::string& aRegionName,
                             const ICaptureComponent* aCaptureComponent,
                             const int aPeriod ) const = 0;

    /*!
     * \brief Get the emissions of a given gas per unit of primary output.
     * \details Return the emissions of the gas per a unit of primary output.
     *          For the primary output, this is the same as the emissions
     *          coefficient. For other outputs, this is determined by their
     *          emissions coefficient and their current ratio of output to the
     *          primary output.
     * \param aGHGName The name of the gas.
     * \param aPeriod Model period
     * \return The emissions of the gas per unit of output.
     */
    virtual double getEmissionsPerOutput( const std::string& aGHGName,
                                          const int aPeriod ) const = 0;

    // Documentation is inherited.
    virtual void accept( IVisitor* aVisitor, const int aPeriod ) const = 0;

    /*!
     * \brief Set the land allocator
     * \param aLandAllocator - the new land allocator
     * \param aName - the name of the technology
     */
    virtual void sendLandAllocator( const ILandAllocator* aLandAllocator,
                                   const std::string& aName ) = 0;
    
    /*!
     * \brief Hook for an output to do interpolations to fill in any data that
     *        should be interpolated to a newly created output for the missing
     *        technology.
     * \param aYear the year to be filled in.
     * \param aPreviousYear The year of the last parsed input.
     * \param aNextYear The year of the next closest parsed input.
     * \param aPreviousOutput The previous parsed output.
     * \param aNextOutput The next parsed output.
     */
    virtual void doInterpolations( const int aYear, const int aPreviousYear,
                                   const int aNextYear, const IOutput* aPreviousInput,
                                   const IOutput* aNextInput ) = 0;
};

// Inline function definitions.
IOutput::IOutput()
{
}

IOutput::~IOutput()
{
}

#endif // _IOUTPUT_H_
