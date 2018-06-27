#ifndef _BY_PRODUCT_H_
#define _BY_PRODUCT_H_
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
* \file by_product.h
* \ingroup Objects
* \brief The ByProduct class header file.
* \author Sonny Kim
*/

#include <xercesc/dom/DOMNode.hpp>
#include <string>
#include <vector>
#include "technologies/include/ioutput.h"
#include "util/base/include/value.h"

// Forward declaration
class Tabs;
class ICaptureComponent;

/*! 
* \ingroup Objects
* \brief ByProduct represents a secondary output of a Technology which has a
*        zero or negative value.
* \details ByProduct represents a secondary output of a Technology which can be
*          calculated based on a primary output quantity and a coefficient. The
*          ByProduct may also be sequestered using an ICaptureComponent from the
*          Technology.
* \todo XML documentation
* \author Sonny Kim, Josh Lurz
*/

class ByProduct : public IOutput {
    friend class OutputFactory;
public:
    /*!
     * \brief Get the XML name for the class.
     * \return The XML name for the class.
     */
    static const std::string& getXMLNameStatic();

    virtual const std::string& getXMLReportingName() const;

    virtual ByProduct* clone() const;

    virtual bool isSameType( const std::string& aType ) const;

    virtual const std::string& getName() const;

    virtual bool XMLParse( const xercesc::DOMNode* aNode );

    virtual void toInputXML( std::ostream& aOut,
                             Tabs* aTabs ) const;

    virtual void toDebugXML( const int aPeriod,
                             std::ostream& aOut,
                             Tabs* aTabs ) const;

    virtual void completeInit( const std::string& aSectorName,
                               DependencyFinder* aDependencyFinder,
                               const IInfo* aTechInfo,
                               const bool aIsTechOperating );

    virtual void initCalc( const std::string& aRegionName,
                           const std::string& aSectorName,
                           const int aPeriod );
    
    virtual void postCalc( const std::string& aRegionName,
                           const int aPeriod );

    virtual void scaleCoefficient( const double aScaler );

    virtual OutputList calcPhysicalOutput( const double aPrimaryOutput,
                                           const std::string& aRegionName,
                                           const ICaptureComponent* aCaptureComponent,
                                           const int aPeriod ) const;

    virtual void setPhysicalOutput( const double aPrimaryOutput,
                                    const std::string& aRegionName,
                                    ICaptureComponent* aCaptureComponent,
                                    const int aPeriod );

    virtual double getPhysicalOutput( const int aPeriod ) const;

    virtual void setCurrencyOutput( const std::string& aRegionName,
                                    const double aOutput,
                                    const int aPeriod )
    {
        // TODO: This could work by converting from physical to currency with
        // the market price.
    }

    virtual double getCurrencyOutput( const int aPeriod ) const
    {
        // TODO: This could work by converting from physical to currency with
        // the market price.
        return 0;
    }

    virtual double getValue( const std::string& aRegionName,
                             const ICaptureComponent* aCaptureComponent,
                             const int aPeriod ) const;

    virtual double getEmissionsPerOutput( const std::string& aGHGName,
                                          const int aPeriod ) const;

    virtual void accept( IVisitor* aVisitor,
                         const int aPeriod ) const;

    // Documentation is inherited.
    virtual void sendLandAllocator(
       const ILandAllocator*    aLandAllocator,
       const std::string& aName ) {}
    
    virtual void doInterpolations( const int aYear, const int aPreviousYear,
                                   const int aNextYear, const IOutput* aPreviousInput,
                                   const IOutput* aNextInput );
protected:

    /*!
     * \brief Protected constructor so the class can only be created by the
     *        OutputFactory.
     */
    ByProduct();
    
    double calcPhysicalOutputInternal( const double aPrimaryOutput,
                                       const ICaptureComponent* aCaptureComponent ) const;

    //! Quantity of the byproduct(calculated).
    std::vector<double> mQuantities;    
    
    //! Name of the byproduct.
    std::string mName;

    //! Coefficient of the byproduct from the output quantity.
    Value mCoef;

    //! Coefficient of the byproduct adjusted for the initial mass.
    Value mAdjustedCoef;

private :
    const static std::string XML_REPORTING_NAME; //!< tag name for reporting xml db 
};

#endif // _BY_PRODUCT_H_

