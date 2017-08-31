#ifndef _PRIMARY_OUTPUT_H_
#define _PRIMARY_OUTPUT_H_
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
 * \file primary_output.h
 * \ingroup Objects
 * \brief PrimaryOutput class header file.
 * \author Josh Lurz
 */

#include <string>
#include <xercesc/dom/DOMNode.hpp>

class Tabs;
class DependencyFinder;
class CachedMarket;

#include "technologies/include/ioutput.h"
#include "util/base/include/value.h"

/*! 
 * \ingroup Objects
 * \brief An output representing the primary output of a Technology.
 * \details The primary output of a Technology is the same as the output of the
 *          sector containing a Technology. The primary output should always be
 *          chosen such that the revenue from the primary output always exceeds
 *          any revenue from secondary outputs. Primary outputs have expenses and
 *          revenues calculated directly by the Technology object.
 *          
 *          <b>XML specification for PrimaryOutput</b>
 *          - XML name: None. The primary output is generated automatically
                        and cannot be parsed.
 *
 * \author Josh Lurz
 */
class PrimaryOutput : public IOutput {
public:
    /*
     * \brief Constructor
     * \param aSectorName Name of the sector and primary output.
     */
    PrimaryOutput( const std::string& aSectorName );
    
    PrimaryOutput( const PrimaryOutput& aPrimaryOutput );
    
    virtual ~PrimaryOutput();

    virtual PrimaryOutput* clone() const;

    virtual bool isSameType( const std::string& aType ) const;

    virtual const std::string& getName() const;

    virtual const std::string& getXMLReportingName() const;

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

    virtual void accept( IVisitor* aVisitor, const int aPeriod ) const;

    // Documentation is inherited.
    virtual void sendLandAllocator( const ILandAllocator* aLandAllocator,
                                   const std::string& aName ) {}
    
    virtual void doInterpolations( const int aYear, const int aPreviousYear,
                                   const int aNextYear, const IOutput* aPreviousInput,
                                   const IOutput* aNextInput ) {}

protected:
    //! Physical output by period.
    std::vector<Value> mPhysicalOutputs;

    //! Name of the primary output. This is the same as the sector name.
    const std::string mName;

    //! CO2 emissions coefficient cached from the marketplace.
    Value mCachedCO2Coef;
    
    //! A pre-located market which has been cached from the marketplace to add supply to.
    std::auto_ptr<CachedMarket> mCachedMarket;

private:
    const static std::string XML_REPORTING_NAME; //!< tag name for reporting xml db 
};

#endif // _PRIMARY_OUTPUT_H_
