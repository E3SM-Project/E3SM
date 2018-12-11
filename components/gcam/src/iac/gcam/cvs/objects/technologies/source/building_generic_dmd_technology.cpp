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
 * \file building_generic_dmd_technology.cpp
 * \ingroup Objects
 * \brief BuildingGenericDmdTechnology technology source file.
 * \author Steve Smith, Josh Lurz
 */

#include "util/base/include/definitions.h"
#include <string>
#include <cassert>
#include <xercesc/dom/DOMNode.hpp>

#include "technologies/include/building_generic_dmd_technology.h"
#include "util/base/include/xml_helper.h"
#include "containers/include/scenario.h"
#include "functions/include/iinput.h"
#include "util/logger/include/ilogger.h"
#include "containers/include/iinfo.h"
#include "sectors/include/sector_utils.h"
#include "containers/include/info_factory.h"
#include "functions/include/function_manager.h"
#include "functions/include/function_utils.h"
#include "technologies/include/iproduction_state.h"
#include "util/base/include/ivisitor.h"
#include "util/base/include/configuration.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;

typedef vector<IInput*>::iterator InputIterator;
typedef vector<IInput*>::const_iterator CInputIterator;

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
const string& BuildingGenericDmdTechnology::getXMLNameStatic() {
	const static string XML_NAME1D = "building-demand-technology";
	return XML_NAME1D;
}

/*!
 * \brief Constructor
 * \param aName Technology name.
 * \param aYear Technology year.
 */
BuildingGenericDmdTechnology::BuildingGenericDmdTechnology( const string& aName, const int aYear )
: Technology( aName, aYear ){
}

/*!
 * \brief Destructor
 */
BuildingGenericDmdTechnology::~BuildingGenericDmdTechnology() {
}

/*!
 * \brief Copy constructor.
 * \param aTechnology Building demand technology to copy.
 */
BuildingGenericDmdTechnology::BuildingGenericDmdTechnology( const BuildingGenericDmdTechnology& aOther ):
Technology( aOther ){
    mAveInsulation = aOther.mAveInsulation;
    mFloorToSurfaceArea = aOther.mFloorToSurfaceArea;

    // TODO: Can't currently copy info objects. This should not be a problem as
    // copying is done before completeInit is called.
}

BuildingGenericDmdTechnology* BuildingGenericDmdTechnology::clone() const {
	return new BuildingGenericDmdTechnology( *this );
}

/*! \brief Get the XML node name for output to XML.
*
* This public function accesses the private constant string, XML_NAME.
* This way the tag is always consistent for both read-in and output and can be easily changed.
* This function may be virtual to be overridden by derived class pointers.
* \author Josh Lurz, James Blackwood
* \return The constant XML_NAME.
*/
const std::string& BuildingGenericDmdTechnology::getXMLName() const {
	return getXMLNameStatic();
}

void BuildingGenericDmdTechnology::completeInit( const string& aRegionName,
                                                 const string& aSectorName,
                                                 const string& aSubsectorName,
                                                 DependencyFinder* aDepFinder,
                                                 const IInfo* aSubsectorInfo,
                                                 ILandAllocator* aLandAllocator )
{
    // Construct the info object before calling Technology::completeInit so that
    // the info object can be used in the Technology.
    mInfo.reset( InfoFactory::constructInfo( aSubsectorInfo, aRegionName + "-" + aSectorName + "-" + 
        aSubsectorName + "-" + mName ) );
    
    if( util::isEqual( mAveInsulation, 0.0 ) ) {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::WARNING );
        mainLog << "Input variable aveInsulationis 0. Reset to 1." << endl;
        mAveInsulation = 1;
    }

    if( util::isEqual( mFloorToSurfaceArea, 0.0 ) ) {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::WARNING );
        mainLog << "Input variable floorToSurfaceArea is 0. Reset to 1." << endl;
        mFloorToSurfaceArea = 1;
    }
    
    // Add the average insulation value and floor to surface area ratio.
    mInfo->setDouble( "average-insulation", mAveInsulation );
    mInfo->setDouble( "floor-to-surface-area", mFloorToSurfaceArea );

    Technology::completeInit( aRegionName, aSectorName, aSubsectorName, 
                              aDepFinder, aSubsectorInfo, aLandAllocator );
}

void BuildingGenericDmdTechnology::initCalc( const string& aRegionName,
                                             const string& aSectorName,
                                             const IInfo* aSubsectorInfo,
                                             const Demographic* aDemographics,
                                             PreviousPeriodInfo& aPrevPeriodInfo,
                                             const int aPeriod )
{
    Technology::initCalc( aRegionName, aSectorName, aSubsectorInfo,
                          aDemographics, aPrevPeriodInfo, aPeriod );

    // Check coefficients in the initial year of the technology.
    // We skip this check if calibration is active since the coefs
    // will be backed out.  If we are missing calibration values upstream
    // that will be flagged elsewhere
    const bool calibrationPeriod = aPeriod > 0 && aPeriod <= scenario->getModeltime()->getFinalCalibrationPeriod();
    static const bool calibrationActive = Configuration::getInstance()->getBool( "CalibrationActive" );
    if( mProductionState[ aPeriod ]->isNewInvestment() && !( calibrationActive && calibrationPeriod ) ){
        checkCoefficients( aSectorName, aPeriod );
    }
}


void BuildingGenericDmdTechnology::postCalc( const string& aRegionName, const int aPeriod ) {
    Technology::postCalc( aRegionName, aPeriod );
}

bool BuildingGenericDmdTechnology::XMLDerivedClassParse( const string& aNodeName, const DOMNode* aCurr ) {
    if( aNodeName == "aveInsulation" ){
        mAveInsulation = XMLHelper<double>::getValue( aCurr );
    } 
    else if( aNodeName == "floorToSurfaceArea" ){
        mFloorToSurfaceArea = XMLHelper<double>::getValue( aCurr );
    }
    else {
        return false;
    }
    return true;
}

void BuildingGenericDmdTechnology::toInputXMLDerived( ostream& out, Tabs* tabs ) const {  
    XMLWriteElementCheckDefault( mAveInsulation, "aveInsulation", out, tabs, 0.0 );
    XMLWriteElementCheckDefault( mFloorToSurfaceArea, "floorToSurfaceArea", out, tabs, 1.0 );
}	

void BuildingGenericDmdTechnology::toDebugXMLDerived( const int period, ostream& out, Tabs* tabs ) const { 
    XMLWriteElement( mAveInsulation, "aveInsulation", out, tabs );
    XMLWriteElement( mFloorToSurfaceArea, "floorToSurfaceArea", out, tabs );
}	

/*!
* \brief calculate technology unnormalized shares
* \details Building technologies are really just calculating demands for
*          specific services so shares are always 1. This ensures that sector
*          price is correctly calculated.
* \author Steve Smith
* \param aRegionName Region name.
* \param aSectorName Sector name.
* \param aGDP Regional GDP container.
* \param aPeriod Model period.
* \return The building demand technology share which is always one.
*/
double BuildingGenericDmdTechnology::calcShare( const string& aRegionName,
                                                const string& aSectorName,
                                                const GDP* aGDP,
                                                const double aLogitExp,
                                                const int aPeriod ) const
{
    return Technology::calcShare( aRegionName, aSectorName, aGDP, aLogitExp, aPeriod );
}

void BuildingGenericDmdTechnology::production( const string& aRegionName,
                                               const string& aSectorName,
											   double aVariableDemand,
                                               double aFixedOutputScaleFactor,
											   const GDP* aGDP, const int aPeriod )
{
    Technology::production( aRegionName, aSectorName, aVariableDemand,
                            aFixedOutputScaleFactor, aGDP, aPeriod );
}

void BuildingGenericDmdTechnology::calcCost( const string& aRegionName,
                                             const string& aSectorName,
											 const int aPeriod )
{
	Technology::calcCost( aRegionName, aSectorName, aPeriod );
}

double BuildingGenericDmdTechnology::getTotalInputCost( const string& aRegionName,
                                                        const string& aSectorName,
														const int aPeriod ) const 
{
	return Technology::getTotalInputCost( aRegionName, aSectorName, aPeriod );
}

const IInfo* BuildingGenericDmdTechnology::getTechInfo() const {
    return mInfo.get();
}

const IFunction* BuildingGenericDmdTechnology::getProductionFunction() const {
    return FunctionManager::getFunction( "minicam-price-elasticity" );
}

/*!
 * \brief Error check the input coefficients to ensure they have been
 *        initialized correctly.
 * \details Prints a warning if the input coefficients have not been
 *          initialized. The input is valid if:
 *              - It read-in a coefficient adjustment
 *              - It copied a coefficient adjustment from the previous period.
 *              - It is a non-energy input.
 * \note The check for coefficient greater than one does not directly check for
 *       a non-unity coefficient adjustment but should be correct most of the
 *       time.
 * \note If calibration is active we would not want to do this check since
 *       consistency checks will be done elsewhere.
 * \param aSectorName Name of the containing sector.
 * \param aPeriod Model period.
 */
void BuildingGenericDmdTechnology::checkCoefficients( const string& aSectorName,
                                                      const int aPeriod ) const
{
    for( CInputIterator i = mInputs.begin(); i != mInputs.end(); ++i ){
        if( ( *i )->hasTypeFlag( IInput::ENERGY ) &&
            ( *i )->getCoefficient( aPeriod ) >= 1 )
        {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            mainLog << "Input " << ( *i )->getName() << " in technology " << mName
                << " with year " << year << " in sector " << aSectorName
                << " has no calibration quantity and a coefficient greater than one."
                << endl;
        }
    }
}

void BuildingGenericDmdTechnology::accept( IVisitor* aVisitor, const int aPeriod ) const {
    aVisitor->startVisitBuildingGenericDmdTechnology( this, aPeriod );
    Technology::accept( aVisitor, aPeriod );
    aVisitor->endVisitBuildingGenericDmdTechnology( this, aPeriod );
}

void BuildingGenericDmdTechnology::doInterpolations( const Technology* aPrevTech, const Technology* aNextTech )
{
    Technology::doInterpolations( aPrevTech, aNextTech );
    
    const BuildingGenericDmdTechnology* prevBldTech = static_cast<const BuildingGenericDmdTechnology*>( aPrevTech );
    const BuildingGenericDmdTechnology* nextBldTech = static_cast<const BuildingGenericDmdTechnology*>( aNextTech );
    
    /*!
     * \pre We are given a valid BuildingGenericDmdTechnology for the previous tech.
     */
    assert( prevBldTech );
    
    /*!
     * \pre We are given a valid BuildingGenericDmdTechnology for the next tech.
     */
    assert( nextBldTech );
    
    mAveInsulation = util::linearInterpolateY( year, prevBldTech->year, nextBldTech->year,
                                               prevBldTech->mAveInsulation, nextBldTech->mAveInsulation );
}
