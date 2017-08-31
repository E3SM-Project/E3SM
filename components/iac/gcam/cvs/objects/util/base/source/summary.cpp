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
* \file summary.cpp
* \ingroup Objects
* \brief An object which contains variables for reporting purposes only.
*  Several maps are created to store and add values for reporting.
* \author Sonny Kim
*/

#include "util/base/include/definitions.h"
#include "util/base/include/summary.h"
#include "util/base/include/util.h"
#include "containers/include/world.h"
#include "containers/include/scenario.h"
#include <vector>

using namespace std;

extern Scenario* scenario;


//! Default constructor
/*   \todo Redesign using a map of maps or other improved datastructure.
*/
Summary::Summary() {
}

/*! \brief Initialize the fuel consumption map.
* \details Adds fuel consumption to the existing consumption and
*  adds to total.
* \param fname Fuel name.
* \param value Fuel consumption.
*/
void Summary::initfuelcons( const string& fname, const double value ){

    fuelcons[ fname ] += value;
    fuelcons[ "zTotal" ] += value;
}

/*! \brief Initialize the primary fuel production map.
* \details Checks the primary fuel list and adds only primary fuel production to the 
*  existing production and adds to total.
* \param aPrimaryFuelList Primary fuel list.
* \param aFuelName Fuel name.
* \param aValue Primary fuel production.
*/
void Summary::initpeprod( const list<string>& aPrimaryFuelList, 
                         const string& aFuelName, const double aValue  ) 
{
    // Map primary fuel list only.
    // Primary-fuel list is read-in from the output-meta-data.
    // Use find method so that a non-existent key is not created in the map.
    if( !aPrimaryFuelList.empty() )
    {
        for( list<string>::const_iterator iter = aPrimaryFuelList.begin(); 
            iter!= aPrimaryFuelList.end(); ++iter )
        {
            if( *iter == aFuelName ){
                peprod[ aFuelName ] += aValue;
                peprod[ "zTotal" ] += aValue;
            }
        }
    }
    // Add all if primary fuel list is empty.
    else
    {
        peprod[ aFuelName ] += aValue;
        peprod[ "zTotal" ] += aValue;
    }

}

/*! \brief Return the fuel consumption map.
*/
const Summary::SummaryItem& Summary::getfuelcons() const {
    return fuelcons;
}

/*! \brief Return the primary fuel consumption map.
*/
const Summary::SummaryItem& Summary::getpecons() const {
    return pecons;
}

/*! \brief Return the primary fuel production map.
*/
const Summary::SummaryItem& Summary::getpeprod() const {
    return peprod;
}

/*! \brief Return the primary fuel trade map.
*/
const Summary::SummaryItem& Summary::getpetrade() const {
    return petrade;
}

/*! \brief Return the emission map.
*/
const Summary::SummaryItem& Summary::getemission() const {
    return emission;
}

/*! \brief Return the emissions by fuel map.
*/
const Summary::SummaryItem& Summary::getemfuelmap() const {
    return emissfuel;
}

//! return map of sequestered amount of emissions
const Summary::SummaryItem& Summary::getSequesteredAmountMap() const {
    return sequesteredAmount;
}

//! Add the passed fuel map to the summary fuelinfo map 
/* The consumption values in the fuelinfo map that is passed are added 
to the summary object maps fuelcons and pecons.

The iterator fmap is used to traverse the fuelinfo map.
* \param aPrimaryFuelList Primary fuel list.
* \param fuelinfo Map of Fuel consumption.
*/
void Summary::updatefuelcons( const list<string>& aPrimaryFuelList, const SummaryItem& fuelinfo ) {
    // map all primary and secondary fuel consumption
    for (CSummaryIterator fmap=fuelinfo.begin(); fmap!=fuelinfo.end(); ++fmap) {    // iterate to one less than the end
        fuelcons[fmap->first] += fmap->second; // Add values from the passed map to fuelcons
        // Don't need a zTotal b/c the fuels are not comparable. 
    }

    // map primary fuel list only.
    for( list<string>::const_iterator fuelIter = aPrimaryFuelList.begin();
        fuelIter != aPrimaryFuelList.end(); ++fuelIter )
    {
        CSummaryIterator fmap=fuelinfo.find( *fuelIter );
        if( fmap!=fuelinfo.end() ) {
            pecons[fmap->first] += fmap->second;
            pecons["zTotal"] += fmap->second;
        }
    }
}

//! Update primary energy trade by substracting consumption from production.
void Summary::updatepetrade() {
    // map all primary and secondary fuel consumption
    for ( CSummaryIterator fmap = peprod.begin(); fmap != peprod.end(); ++fmap ) {
        petrade[ fmap->first ] = peprod[ fmap->first ] - pecons[ fmap->first ];
    }
}

//! Update and add to GHG emissions from passed in GHG emissions map.
//! param ghginfo Map of GHG emissions.
void Summary::updateemiss( const SummaryItem& ghginfo ) {
    // map all primary and secondary fuel consumption
    for ( CSummaryIterator fmap = ghginfo.begin(); fmap != ghginfo.end(); ++fmap){
        emission[ fmap->first ] += fmap->second;
    }
}

//! Update and add to GHG emissions by fuel from passed in GHG emissions map.
//! param ghginfo Map of GHG emissions.
void Summary::updateemfuelmap( const SummaryItem& ghginfo ) {
    // map all primary and secondary fuel consumption
    for ( CSummaryIterator fmap = ghginfo.begin(); fmap != ghginfo.end(); ++fmap ) {
        emissfuel[ fmap->first ] += fmap->second;
    }
}

//! update the map of sequestered amount of emissions
void Summary::updateSequesteredAmountMap( const SummaryItem& ghginfo ) {
    // map sequestered amount of CO2 for secondary fuels and zTotal
    for ( CSummaryIterator fmap = ghginfo.begin(); fmap != ghginfo.end(); ++fmap ) {
        sequesteredAmount[ fmap->first ] += fmap->second;
    }
}

//! Clear fuel consumption and primary energy consumption maps.
void Summary::clearfuelcons() {
    fuelcons.clear();
    pecons.clear();
}

//! Clear primary energy production and trade maps.
void Summary::clearpeprod() {
    peprod.clear();
    petrade.clear();
}

//! Clear emissions map.
void Summary::clearemiss() {
    emission.clear();
}

//! Clear emissions by fuel map.
void Summary::clearemfuelmap() {
    emissfuel.clear();
}

//! clear out map of sequestered amount
void Summary::clearSequesteredAmountMap() {
    sequesteredAmount.clear();
}
/*! \todo Fix all these names. Should not include 'second'.*/
//! Return consumption of fuel from the fuel consumption map.
double Summary::get_fmap_second( const string& name ) const {
    return util::searchForValue( fuelcons, name );
}

/*! \todo Fix all these names. Should not include 'second'.*/
//! Return consumption of a primary fuel from the primary energy consumption map.
double Summary::get_pemap_second( const string& name ) const {
    return util::searchForValue( pecons, name );
}

/*! \todo Fix all these names. Should not include 'second'.*/
//! Return the trade of a primary fuel from the primary energy trade map.
double Summary::get_petrmap_second( const string& name ) const {
    return util::searchForValue( petrade, name );
}

/*! \todo Fix all these names. Should not include 'second'.*/
//! Return primary energy production of a fuel from the primary energy production map.
double Summary::get_peprodmap_second( const string& name ) const {
    return util::searchForValue( peprod, name );
}

/*! \todo Fix all these names. Should not include 'second'.*/
//! Return an emission amount from the emissions map.
double Summary::get_emissmap_second( const string& name ) const {
    return util::searchForValue( emission, name );
}

//! Return the sequestered amount from the sequestered amount map
double Summary::getSequesteredAmount( const string& name ) const {
    return util::searchForValue( sequesteredAmount, name );
}

//! Return emissions from a fuel from the emissions by fuel map.
double Summary::get_emissfuelmap_second( const string& name ) const {
    return util::searchForValue( emissfuel, name );
}
