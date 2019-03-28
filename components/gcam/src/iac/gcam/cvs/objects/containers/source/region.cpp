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
* \file region.cpp
* \ingroup Objects
* \brief The Region class source file.
* \author Sonny Kim
*/

#include "util/base/include/definitions.h"
#include <fstream>
#include <string>
#include <vector>
#include <cassert>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include <algorithm>
#include <memory>

#include "containers/include/region.h"
#include "containers/include/scenario.h"
#include "util/base/include/model_time.h" 
#include "sectors/include/sector.h"
#include "demographics/include/demographic.h"
#include "util/base/include/ivisitor.h"
#include "marketplace/include/marketplace.h"

#include "policy/include/policy_portfolio_standard.h"
#include "policy/include/policy_ghg.h"
#include "emissions/include/total_sector_emissions.h"
#include "emissions/include/emissions_summer.h"

#include "util/curves/include/curve.h"
#include "util/curves/include/point_set_curve.h"
#include "util/curves/include/xy_data_point.h"
#include "util/curves/include/point_set.h"
#include "util/curves/include/explicit_point_set.h"

#include "resources/include/resource.h"
#include "resources/include/unlimited_resource.h"
#include "resources/include/depleting_fixed_resource.h"

#include "containers/include/iinfo.h"

#include "util/logger/include/ilogger.h"
#include "util/base/include/xml_helper.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;

typedef std::vector<Sector*>::iterator SectorIterator;
typedef std::vector<Sector*>::const_iterator CSectorIterator;
typedef std::vector<GHGPolicy*>::iterator GHGPolicyIterator;
typedef std::vector<GHGPolicy*>::const_iterator CGHGPolicyIterator;
typedef std::vector<PolicyPortfolioStandard*>::iterator PolicyIterator;
typedef std::vector<PolicyPortfolioStandard*>::const_iterator CPolicyIterator;
typedef std::vector<AResource*>::iterator ResourceIterator;
typedef std::vector<AResource*>::const_iterator CResourceIterator;
//! Default constructor
Region::Region() {
}

//! Default destructor destroys sector, demsector, Resource, and
//! population objects.
Region::~Region() {
    clear();
}

//! Clear member variables and initialize elemental members.
void Region::clear(){
    for ( SectorIterator secIter = supplySector.begin(); secIter != supplySector.end(); secIter++ ) {
        delete *secIter;
    }

    for( GHGPolicyIterator policyIter = mGhgPolicies.begin(); policyIter != mGhgPolicies.end(); ++policyIter ){
        delete *policyIter;
    }

    for( PolicyIterator policyIter = mPolicies.begin(); policyIter != mPolicies.end(); ++policyIter ){
        delete *policyIter;
    }
    
    for ( ResourceIterator rescIter = mResources.begin(); rescIter != mResources.end(); ++rescIter ) {
        delete *rescIter;
    }
}

/*! Return the region name.
* \return The string name of the region is returned.
*/
const string& Region::getName() const {
    return name;
}

/*! 
* \brief Sets the data members from the XML input.
* \details This function parses all XML data from Region down to the lowest set
*          of objects. As the XML data is parsed, new objects are continually
*          added to the object container using the push_back routine.
* \param node XML DOM node of the region
* \todo Change the diagnostic "assert( node );" to fail with a more informative
*       error (file, previous node?, location?)
*/
void Region::XMLParse( const DOMNode* node ){
    // make sure we were passed a valid node.
    assert( node );

    // get the name attribute.
    name = XMLHelper<string>::getAttr( node, "name" );

    // get all child nodes.
    DOMNodeList* nodeList = node->getChildNodes();

    // loop through the child nodes.
    for( unsigned int i = 0; i < nodeList->getLength(); i++ ){
        DOMNode* curr = nodeList->item( i );
        string nodeName = XMLHelper<string>::safeTranscode( curr->getNodeName() );

        if( nodeName == "#text" ) {
            continue;
        }

        else if( nodeName == Demographic::getXMLNameStatic() ){
            parseSingleNode( curr, demographic, new Demographic );
        }
        else if( nodeName == GHGPolicy::getXMLNameStatic() ){
            parseContainerNode( curr, mGhgPolicies, new GHGPolicy() );
        }
        else if( nodeName == PolicyPortfolioStandard::getXMLNameStatic() ){
            parseContainerNode( curr, mPolicies, new PolicyPortfolioStandard() );
        }
        // TODO: should we create a factory for resources?
        else if( nodeName == DepletableResource::getXMLNameStatic() ){
            parseContainerNode( curr, mResources, new DepletableResource() );
        }
        else if( nodeName == FixedResource::getXMLNameStatic() ){
            parseContainerNode( curr, mResources, new FixedResource() );
        }
        else if( nodeName == RenewableResource::getXMLNameStatic() ){
            parseContainerNode( curr, mResources, new RenewableResource() );
        }
        else if( nodeName == UnlimitedResource::getXMLNameStatic() ){
            parseContainerNode( curr, mResources, new UnlimitedResource );
        }
        else if( nodeName == DepletingFixedResource::getXMLNameStatic() ) {
            parseContainerNode( curr, mResources, new DepletingFixedResource() );
        }
        else if( XMLDerivedClassParse(nodeName, curr) ){
            // Do nothing but avoid printing the error.
        }
        else {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            mainLog << "Unrecognized text string: " << nodeName << " found while parsing region." << endl;
        }
    }
}

/*! 
* \brief Write datamembers to datastream in XML format. Calls XMLWriteElement
*        function from the XMLHelper class for the actual writing.
* \param out Output file in XML format.
* \param tabs Tabs object used to track the number of tabs to print.
* \ref faqitem1 
*/
void Region::toInputXML( ostream& out, Tabs* tabs ) const {
    XMLWriteOpeningTag ( getXMLName(), out, tabs, name );

    // write the xml for the class members.
    // write out the single population object.
    if( demographic.get() ){
        demographic->toInputXML( out, tabs);
    }

    // write out supply sector objects.
    for( CSectorIterator j = supplySector.begin(); j != supplySector.end(); j++ ){
        ( *j )->toInputXML( out, tabs );
    }

    // write out mGhgPolicies objects.
    for( CGHGPolicyIterator l = mGhgPolicies.begin(); l != mGhgPolicies.end(); l++ ){
        ( *l )->toInputXML( out, tabs );
    }

    // write out mPolicies objects.
    for( CPolicyIterator l = mPolicies.begin(); l != mPolicies.end(); l++ ){
        ( *l )->toInputXML( out, tabs );
    }
    
    // write out the resources objects.
    for( CResourceIterator i = mResources.begin(); i != mResources.end(); i++ ){
        ( *i )->toInputXML( out, tabs );
    }

    toInputXMLDerived( out, tabs );

    // finished writing xml for the class members.
    XMLWriteClosingTag( getXMLName(), out, tabs );
}

/*! \brief Write datamembers to datastream in XML format for debugging purposes.  
* Calls XMLWriteElement function from the XMLHelper class for the actual
* writing. Calls debug functions in other contained objects. 
*
* \param period Model time period
* \param out Output file for debugging purposes in XML format
* \param tabs Tabs object used to track the number of tabs to print.
*/
void Region::toDebugXML( const int period, ostream& out, Tabs* tabs ) const {

    XMLWriteOpeningTag ( getXMLName(), out, tabs, name );

    // write out mGhgPolicies objects.
    for( CGHGPolicyIterator currPolicy = mGhgPolicies.begin(); currPolicy != mGhgPolicies.end(); ++currPolicy ){
        (*currPolicy)->toDebugXML( period, out, tabs );
    }

    if( demographic.get() ){
        demographic->toDebugXML( period, out, tabs );
    }

    // write out supply sector objects.
    for( CSectorIterator j = supplySector.begin(); j != supplySector.end(); j++ ){
        ( *j )->toDebugXML( period, out, tabs );
    }

    // write out mGhgPolicies objects.
    for( CGHGPolicyIterator currPolicy = mGhgPolicies.begin(); currPolicy != mGhgPolicies.end(); ++currPolicy ){
        (*currPolicy)->toDebugXML( period, out, tabs );
    }

    // write out mPolicies objects.
    for( CPolicyIterator currPolicy = mPolicies.begin(); currPolicy != mPolicies.end(); ++currPolicy ){
        (*currPolicy)->toDebugXML( period, out, tabs );
    }
    
    // write out the resources objects.
    for( CResourceIterator currResource = mResources.begin(); currResource != mResources.end(); ++currResource ){
        (*currResource)->toDebugXML( period, out, tabs );
    }

    toDebugXMLDerived( period, out, tabs );
    // Finished writing xml for the class members.

    XMLWriteClosingTag( getXMLName(), out, tabs );
}

/*! \brief Get the XML node name in static form for comparison when parsing XML.
*
* This public function accesses the private constant string, XML_NAME. This way
* the tag is always consistent for both read-in and output and can be easily
* changed. The "==" operator that is used when parsing, required this second
* function to return static.
* \note A function cannot be static and virtual.
* \author Josh Lurz, James Blackwood
* \return The constant XML_NAME as a static.
*/
const std::string& Region::getXMLNameStatic() {
    static const string XML_NAME = "region";
    return XML_NAME;
}

/*! \brief Complete the initialization. 
 *
 * Completes the initialization for this regions data members.
 * \todo Figure out if Demographics and supplysectors calls to completeInit can
 *       be moved down from RegionMiniCAM and RegionCGE.
 * \author Pralit Patel
 */
void Region::completeInit() {    
    for( GHGPolicyIterator ghgPolicy = mGhgPolicies.begin(); ghgPolicy != mGhgPolicies.end(); ++ghgPolicy ){
        (*ghgPolicy)->completeInit( name );
    }
    for( PolicyIterator policy = mPolicies.begin(); policy != mPolicies.end(); ++policy ){
        (*policy)->completeInit( name );
    }
    for( ResourceIterator resourceIter = mResources.begin(); resourceIter != mResources.end(); ++resourceIter ) {
        (*resourceIter)->completeInit( name, mRegionInfo.get() );
    }
}

/*! \brief Function to finalize objects after a period is solved.
* \details This function is used to calculate and store variables which are only
*          needed after the current period is complete. 
* \param aPeriod The period to finalize.
* \todo Finish this function.
* \author Josh Lurz
*/
void Region::postCalc( const int aPeriod ){
    // Finalize sectors.
    for( SectorIterator sector = supplySector.begin(); sector != supplySector.end(); ++sector ){
        (*sector)->postCalc( aPeriod );
    }
    // Post calculation for resource sectors
    for( ResourceIterator currResource = mResources.begin(); currResource != mResources.end(); ++currResource ){
        (*currResource)->postCalc( name, aPeriod );
    }
}

/*! \brief Update a visitor for a Region.
* \param aVisitor Visitor to update.
* \param aPeriod Period to update.
*/
void Region::accept( IVisitor* aVisitor, const int aPeriod ) const {
    aVisitor->startVisitRegion( this, aPeriod );

    // Visit demographics object.
    if( demographic.get() ){
        demographic->accept( aVisitor, aPeriod );
    }

    // loop for supply sectors
    for( CSectorIterator currSec = supplySector.begin(); currSec != supplySector.end(); ++currSec ){
        (*currSec)->accept( aVisitor, aPeriod );
    }
    
    // loop for resources.
    for( CResourceIterator currResource = mResources.begin(); currResource != mResources.end(); ++currResource ){
        (*currResource)->accept( aVisitor, aPeriod );
    }

    aVisitor->endVisitRegion( this, aPeriod );
}

/*! \brief Set a policy for the region.
* \details Searches through the list of the region's taxes for a tax with the
*          same name as aTax. If the tax is found, it is deleted and replaced
*          with aTax. Otherwise aTax is added to the end of the tax vector.
* \param aTax Tax to add.
*/
void Region::setTax( const GHGPolicy* aTax ){
    /*! \pre Tax is not null. */
    assert( aTax );

    // Check if the tax is applicable.
    if( !aTax->isApplicable( name ) ){
        return;
    }

    GHGPolicy* insertedTax = 0;
    // Search for an existing policy to replace.
    for( unsigned int i = 0; i < mGhgPolicies.size(); i++ ){
        if( mGhgPolicies[ i ]->getName() == aTax->getName() ){
            delete mGhgPolicies[ i ];
            // Create a copy of the tax.
            mGhgPolicies[ i ] = insertedTax = aTax->clone();
        }
    }

    // Need to insert the tax.
    if( !insertedTax ){
        insertedTax = aTax->clone();
        mGhgPolicies.push_back( insertedTax );
    }

    // Initialize the tax.
    insertedTax->completeInit( name );
}

/*! \brief A function to generate a ghg emissions quantity curve based on an
*          already performed model run.
* \details This function used the information stored in it to create a curve,
*          with each datapoint containing a time period and an amount of gas
*          emissions. These values are retrieved from the emissions.
* \note The user is responsible for deallocating the memory in the returned
*       Curve.
* \author Josh Lurz
* \param ghgName The name of the ghg to create a curve for.
* \return A Curve object representing ghg emissions quantity by time period.
*/
const Curve* Region::getEmissionsQuantityCurve( const string& ghgName ) const {
    /*! \pre The run has been completed. */
    const Modeltime* modeltime = scenario->getModeltime();

    auto_ptr<ExplicitPointSet> emissionsPoints( new ExplicitPointSet() );

    for( int i = 0; i < scenario->getModeltime()->getmaxper(); i++ ) {
        EmissionsSummer emissionsSummer( ghgName );
        accept( &emissionsSummer, i );
        XYDataPoint* currPoint = new XYDataPoint( modeltime->getper_to_yr( i ),
            emissionsSummer.getEmissions( i ) );
        emissionsPoints->addPoint( currPoint );
    }

    Curve* emissionsCurve = new PointSetCurve( emissionsPoints.release() );
    emissionsCurve->setTitle( ghgName + " emissions curve" );
    emissionsCurve->setXAxisLabel( "year" );
    emissionsCurve->setYAxisLabel( "emissions quantity" );

    return emissionsCurve;
}

/*! \brief A function to generate a ghg emissions price curve based on an
*          already performed model run.
* \details This function used the information stored in it to create a curve,
*          with each datapoint containing a time period and the price gas
*          emissions. These values are retrieved from the marketplace.
* \note The user is responsible for deallocating the memory in the returned
*       Curve.
* \author Josh Lurz
* \param ghgName The name of the ghg to create a curve for.
* \return A Curve object representing the price of ghg emissions by time period.
*/
const Curve* Region::getEmissionsPriceCurve( const string& ghgName ) const {
    /*! \pre The run has been completed. */
    const Modeltime* modeltime = scenario->getModeltime();
    const Marketplace* marketplace = scenario->getMarketplace();

    auto_ptr<ExplicitPointSet> emissionsPoints( new ExplicitPointSet() );

    for( int i = 0; i < modeltime->getmaxper(); i++ ) {
        XYDataPoint* currPoint = new XYDataPoint( modeltime->getper_to_yr( i ), marketplace->getPrice( ghgName, name, i ) );
        emissionsPoints->addPoint( currPoint );
    }

    Curve* emissionsCurve = new PointSetCurve( emissionsPoints.release() );
    emissionsCurve->setTitle( ghgName + " emissions tax curve" );
    emissionsCurve->setXAxisLabel( "year" );
    emissionsCurve->setYAxisLabel( "emissions tax" );

    return emissionsCurve;
}
