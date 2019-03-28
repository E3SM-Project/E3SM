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
 * \file empty_technology.cpp
 * \ingroup Objects
 * \brief EmptyTechnology class source file.
 * \author Pralit Patel
 */

// Standard Library headers
#include "util/base/include/definitions.h"
#include <string>
#include <cassert>
#include "technologies/include/empty_technology.h"

using namespace std;
using namespace xercesc;

EmptyTechnology* EmptyTechnology::getInstance() {
    static EmptyTechnology INSTANCE;
    
    return &INSTANCE;
}

bool EmptyTechnology::isSameType( const string& aType ) const {
    return false;
}

ITechnology* EmptyTechnology::clone() const {
    assert( false );
    
    return 0;
}

const string& EmptyTechnology::getXMLName() const {
    const static string XML_NAME = "empty-technology";
    
    return XML_NAME;
}

bool EmptyTechnology::XMLParse( const DOMNode* aNode )
{
    assert( false );
    
    return false;
}

void EmptyTechnology::completeInit( const string& aRegionName,
                              const string& aSectorName,
                              const string& aSubsectorName,
                              DependencyFinder* aDepFinder,
                              const IInfo* aSubsectorInfo,
                              ILandAllocator* aLandAllocator )
{
}

void EmptyTechnology::toInputXML( ostream& aOut, Tabs* aTabs ) const {
    assert( false );
}

void EmptyTechnology::toDebugXML( const int aPeriod, ostream& aOut, Tabs* aTabs ) const {
}

void EmptyTechnology::initCalc( const string& aRegionName,
                          const string& aSectorName,
                          const IInfo* aSubsectorInfo,
                          const Demographic* aDemographics,
                          PreviousPeriodInfo& aPrevPeriodInfo,
                          const int aPeriod )
{
}

void EmptyTechnology::postCalc( const string& aRegionName, const int aPeriod ) {
}

double EmptyTechnology::getTotalGHGCost( const string& aRegionName,
                                   const string& aSectorName,
                                   const int aPeriod ) const
{
    return 0;
}

double EmptyTechnology::calcShare( const string& aRegionName,
                             const string& aSectorName,
                             const GDP* aGDP,
                             const double aLogitExp,
                             const int aPeriod ) const
{
    return 0;
}

double EmptyTechnology::getFixedOutput( const string& aRegionName,
                                  const string& aSectorName,
                                  const bool aHasRequiredInput,
                                  const string& aRequiredInput,
                                  const int aPeriod ) const
{
    return -1;
}

void EmptyTechnology::production( const string& aRegionName,
                            const string& aSectorName,
                            double aVariableDemand,
                            double aFixedOutputScaleFactor,
                            const GDP* aGDP,
                            const int aPeriod )
{
}

const map<string, double> EmptyTechnology::getEmissions( const string& aGoodName,
                                                   const int aPeriod ) const
{
    const map<string, double> emissions;
    
    return emissions;
}

const map<string, double> EmptyTechnology::getEmissionsByFuel( const string& aGoodName,
                                                         const int aPeriod ) const
{
    const map<string, double> emissionsByFuel;
    
    return emissionsByFuel;
}

const string& EmptyTechnology::getName() const
{
    const static string name = "empty";
    return name;
}

Value EmptyTechnology::getShareWeight() const
{
    const Value zeroShareWeight( 0 );
    
    return zeroShareWeight;
}

Value EmptyTechnology::getParsedShareWeight() const {
    const Value zeroShareWeight( 0 );
    
    return zeroShareWeight;
}

void EmptyTechnology::setShareWeight( double shareWeightValue )
{
}

bool EmptyTechnology::isOutputFixed( const bool aHasRequiredInput,
                               const string& aRequiredInput,
                               const int aPeriod ) const
{
    // zero share weight implies output is fixed
    return true;
}

bool EmptyTechnology::isFixedOutputTechnology( const int aPeriod ) const
{
    // empty technologies do not have fixed output values
    return false;
}

bool EmptyTechnology::isAvailable( const int aPeriod ) const
{
    return false;
}

double EmptyTechnology::getOutput( const int aPeriod ) const
{
    return 0;
}

double EmptyTechnology::getTotalInputCost( const string& aRegionName,
                                     const string& aSectorName,
                                     const int aPeriod ) const
{
    return 0;
}

double EmptyTechnology::getEnergyCost( const string& aRegionName,
                                 const string& aSectorName,
                                 const int aPeriod ) const
{
    return 0;
}

double EmptyTechnology::getEnergyInput( const int aPeriod ) const
{
    return 0;
}

double EmptyTechnology::getCalibrationOutput( const bool aHasRequiredInput,
                                        const string& aRequiredInput,
                                        const int aPeriod ) const
{
    return -1;
}

bool EmptyTechnology::hasCalibratedValue( const int aPeriod ) const {
    return false;
}

void EmptyTechnology::calcCost( const string& aRegionName,
                          const string& aSectorName,
                          const int aPeriod )
{
}

double EmptyTechnology::getCost( const int aPeriod ) const
{
    return 0;
}

double EmptyTechnology::calcFuelPrefElasticity( const int aPeriod ) const
{
    return 0;
}

const vector<string> EmptyTechnology::getGHGNames() const
{
    const vector<string> names;

    return names;
}

int EmptyTechnology::getNumbGHGs()  const {
    return 0;
}

void EmptyTechnology::copyGHGParameters( const AGHG* prevGHG ) {
}

const AGHG* EmptyTechnology::getGHGPointer( const string& aGHGName ) const {
    return 0;
}

double EmptyTechnology::getEmissionsByGas( const string& aGasName,
                                     const int aPeriod ) const
{
    return 0;
}

const map<string, double> EmptyTechnology::getFuelMap( const int aPeriod ) const
{
    const map<string, double> inputMap;

    return inputMap;
}

bool EmptyTechnology::isAllCalibrated( const int aPeriod,
                                 double aCalAccuracy,
                                 const string& aRegionName,
                                 const string& aSectorName,
                                 const string& aSubsectorName,
                                 const bool aPrintWarnings ) const
{
    return true;
}

void EmptyTechnology::setYear( const int aYear )
{
}

void EmptyTechnology::doInterpolations( const Technology* aPrevTech,
                                        const Technology* aNextTech )
{
}

void EmptyTechnology::accept( IVisitor* aVisitor, const int aPeriod ) const {
}

void EmptyTechnology::acceptDerived( IVisitor* aVisitor, const int aPeriod ) const {
}
