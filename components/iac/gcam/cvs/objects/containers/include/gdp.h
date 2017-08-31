#ifndef _GDP_H_
#define _GDP_H_
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
* \file gdp.h
* \ingroup Objects
* \brief The GDP class header file.
* \author Josh Lurz
*/
#include <vector>
#include "util/base/include/iround_trippable.h"
#include "util/base/include/ivisitable.h"
#include "util/base/include/value.h"

class Demographic;

/*! 
* \ingroup Objects
* \brief This class defines an object which contains the GDP information and calculations for a single region
* along with function which can be used to access the GDP in various ways.
* \details This class all controls the read-in and initialization of this data along with calibration
* of GDP. The class contains code to check if adjusted GDP values exist when they are requested and 
* prints an error if this is not the case.
*
* \note This class is constructed of code that was formerly in several classes throughout the model. 
* \author Josh Lurz, Sonny Kim, Steve Smith
*/

class GDP: public IRoundTrippable, public IVisitable
{
    friend class XMLDBOutputter;
private:
    std::vector<Value> laborProdGrowthRate; //!< labor productivity growth rate
    std::vector<Value> laborForceParticipationPercent; //!< labor force participation percent
    std::vector<double> laborForce; //!< actual labor force
    std::vector<double> gdpValue; //!< approximate regional gross domestic product in constant dollars, not adjusted for energy price for this period
    std::vector<double> gdpPerCapita; //!< regional gross domestic product per capita in constant dollars ($)
    std::vector<double> gdpValueAdjusted; //!< regional gross domestic product adjusted for energy price feedback
    std::vector<double> gdpValueNotAdjusted; //!< regional gross domestic product without any adjustments for energy price feedback in any period
    std::vector<double> gdpPerCapitaNotAdjusted; //!< regional GDP per cap without any adjustments for energy price feedback in any period
    std::vector<double> gdpValueAdjustedPPP; //!< regional adjusted GDP in PPP terms
    std::vector<double> gdpPerCapitaAdjusted; //!< regional gross domestic product per capita in constant dollars ($)
    std::vector<double> gdpPerCapitaAdjustedPPP; //!< regional gross domestic product per capita in constant dollars ($)
    std::vector<double> gdpPerCapitaApproxPPP; //!< approximate regional GDP per capita PPP terms (before energy price adjustment)
    std::vector<bool> gdpAdjustedFlag; //!< flag to tell if GDPs have been adjusted yet
    std::vector<double> calibrationGDPs; //!< Calibration values for GDP (constant dollars)
    std::string mGDPUnit; //!< Unit for GDP
    double baseGDP; //!< Base-year value (constant dollars) for regional GDP
    double mEnergyGDPElasticity; //!< Energy service price feedback elasticity for GDP
    double PPPConversionFact; //!< 1990 Ratio of PPP to Market GDP
    double PPPDelta; //!< Internal exponent variable for PPP conversion
    bool constRatio; //!< Flag to turn on dynamic ratio of PPP to Market GDP
    static const std::string XML_NAME; //!< node name for toXML methods
    double calculatePPPPerCap( const int period,const double marketGDPperCap ); // function to calculate PPP values
    double getPPPMERRatio( const int period, const double marketGDPperCap ); // function to calculate PPP/MER ratio
    double getTotalLaborProductivity( const int period ) const;
    double getLaborForce( const int per ) const;
    int findNextPeriodWithValue( const int aStartPeriod, const std::vector<Value>& aValueVector ) const;
public:
    GDP();
    void XMLParse( const xercesc::DOMNode* node );
    void toInputXML( std::ostream& out, Tabs* tabs ) const;
    void toDebugXML( const int period, std::ostream& out, Tabs* tabs ) const;
    void initData( const Demographic* regionalPop );
    void initialGDPcalc( const int period, const double population);
    static const std::string& getXMLNameStatic();
    void setupCalibrationMarkets( const std::string& regionName, const std::vector<double> aCalibrationGDPs  );
    void writeBackCalibratedValues( const std::string& regionName, const int period );
    void csvOutputFile( const std::string& regionName ) const;
    void dbOutput( const std::string& regionName ) const;
    void adjustGDP( const int period, const double priceratio );
    double getApproxGDPperCap( const int period ) const;
    double getApproxScaledGDPperCap( const int period ) const;
    double getApproxScaledGDP( const int period ) const;
    double getScaledGDPperCap( const int period ) const;
    double getGDPperCap( const int period ) const;
    double getPPPGDPperCap( const int period ) const;
    double getGDP( const int period ) const;
    double getApproxGDP( const int period ) const;
    double getBestScaledGDPperCap( const int period ) const;
    double getGDPNotAdjusted( const int period ) const;
    double getGDPPerCapitaNotAdjusted( const int period ) const;
    double getApproxPPPperCap( const int period ) const;
    void accept( IVisitor* aVisitor, const int aPeriod ) const;
 };

#endif // _GDP_H_

