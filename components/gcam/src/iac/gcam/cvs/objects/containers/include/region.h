#ifndef _REGION_H_
#define _REGION_H_
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
* \file region.h
* \ingroup Objects
* \brief The Region class header file.
* \author Sonny Kim
*/

#include <map>
#include <vector>
#include <memory>
#include <string>
#include <list>
#include <xercesc/dom/DOMNode.hpp>

#include "util/base/include/ivisitable.h"
#include "util/base/include/iround_trippable.h"
#include "util/base/include/summary.h"
#include <boost/noncopyable.hpp>

// Forward declarations.
class Demographic;
class Sector;
class GHGPolicy;
class PolicyPortfolioStandard;
class Curve;
class AResource;
class IInfo;
/*! 
* \ingroup Objects
* \brief This is an abstract base class for Regions.
*
* \author Sonny Kim
*/

class Region: public IVisitable, public IRoundTrippable, protected boost::noncopyable
{
    friend class XMLDBOutputter;
public:
    Region();
    virtual ~Region();
    void XMLParse( const xercesc::DOMNode* node );
    void toInputXML( std::ostream& out, Tabs* tabs ) const;
    void toDebugXML( const int period, std::ostream& out, Tabs* tabs ) const;
    static const std::string& getXMLNameStatic();
    virtual void completeInit();
    const std::string& getName() const;

    virtual void calc( const int period ) = 0;
    
    virtual void initCalc( const int period ) = 0;
    
    virtual void postCalc( const int aPeriod ) = 0;

    virtual void csvOutputFile() const {};
    virtual void dbOutput( const std::list<std::string>& aPrimaryFuelList ) const {};
    virtual void initializeAgMarketPrices( const std::vector<double>& pricesIn ) {};
    virtual void updateSummary( const std::list<std::string>& aPrimaryFuelList, const int period ) {};
    virtual const Summary& getSummary( const int period ) const { static const Summary nullSummary; return nullSummary; };

    void setTax( const GHGPolicy* aTax );
    const Curve* getEmissionsQuantityCurve( const std::string& ghgName ) const;
    const Curve* getEmissionsPriceCurve( const std::string& ghgName ) const;

    virtual bool isAllCalibrated( const int period, double calAccuracy, const bool printWarnings ) const { return true; };
    virtual void setCalSuppliesAndDemands( const int period ) {};
    virtual void initializeCalValues( const int period ) {};
    virtual void updateAllOutputContainers( const int period ) = 0;
    virtual void updateMarketplace( const int period ) {};

    virtual void csvSGMOutputFile( std::ostream& aFile, const int period ) const {};

    virtual void accept( IVisitor* aVisitor, const int aPeriod ) const;
    virtual void csvSGMGenFile( std::ostream& aFile ) const {};
protected:
    std::string name; //!< Region name
    std::auto_ptr<Demographic> demographic; //!< Population object
    std::vector<Sector*> supplySector; //!< vector of pointers to supply sector objects
    //! vector of pointers to ghg market objects, container for constraints and emissions
    std::vector<GHGPolicy*> mGhgPolicies;
    //! vector of pointers to portfolio standard market objects, container for constraints
    std::vector<PolicyPortfolioStandard*> mPolicies;
    std::vector<AResource*> mResources; //!< vector of pointers to resource objects
    //! The region's information store.
    std::auto_ptr<IInfo> mRegionInfo;

    virtual const std::string& getXMLName() const = 0;
    virtual void toInputXMLDerived( std::ostream& out, Tabs* tabs ) const = 0;
    virtual bool XMLDerivedClassParse( const std::string& nodeName, const xercesc::DOMNode* curr ) = 0;
    virtual void toDebugXMLDerived( const int period, std::ostream& out, Tabs* tabs ) const = 0;
private:
    void clear();
};

#endif // _REGION_H_
