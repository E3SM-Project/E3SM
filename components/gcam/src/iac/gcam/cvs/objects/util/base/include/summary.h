#ifndef _SUMMARY_H_
#define _SUMMARY_H_
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
* \file summary.h
* \ingroup Objects
* \brief The Summary class header file.
* \author Sonny Kim
*/

#include <map>
#include <string>
#include <list>

/*! 
* \ingroup Objects
* \brief An object which contains variables for reporting purposes only.
*  Several maps are created to store and add values for reporting.
* \author Sonny Kim
*/

class Summary
{
private:
    typedef std::map<std::string, double> SummaryItem;
    typedef SummaryItem::const_iterator CSummaryIterator;
    typedef SummaryItem::iterator SummaryIterator;
    SummaryItem fuelcons;  //!< map of fuel name and amount consumed
    SummaryItem pecons;  //!< map of primary energy consumption
    SummaryItem peprod;  //!< map of primary energy production
    SummaryItem petrade;  //!< map of primary energy trade
    SummaryItem emission;  //!< map of ghg emissions
    SummaryItem emissfuel;  //!< map of ghg emissions implicit in fuel
    SummaryItem sequesteredAmount;  //!< map of sequestered amount of emissions
public:
    Summary(); // default constructor
    void initfuelcons( const std::string& fname, const double value );
    //void initpeprod( const std::string& fname, const double value );
    void initpeprod( const std::list<std::string>& aPrimaryFuelList, 
        const std::string& aFuelName, const double aValue );
    const SummaryItem& getfuelcons() const;
    const SummaryItem& getpecons() const;
    const SummaryItem& getpeprod() const;
    const SummaryItem& getpetrade() const;
    const SummaryItem& getemission() const;
    const SummaryItem& getemfuelmap() const;
    const SummaryItem& getSequesteredAmountMap() const;
    void updatefuelcons( const std::list<std::string>& aPrimaryFuelList, const SummaryItem& fuelinfo);
    void updatepetrade();
    void updateemiss( const SummaryItem& ghginfo );
    void updateemfuelmap( const SummaryItem& ghginfo );
    void updateSequesteredAmountMap( const SummaryItem& ghginfo );
    void clearfuelcons();
    void clearpeprod();
    void clearemiss();
    void clearemfuelmap();
    void clearSequesteredAmountMap();
    double get_fmap_second( const std::string& name ) const;
    double get_pemap_second( const std::string& name ) const;
    double get_petrmap_second( const std::string& name ) const;
    double get_peprodmap_second( const std::string& name ) const;
    double get_emissmap_second( const std::string& name ) const;
    double getSequesteredAmount( const std::string& name ) const;
    double get_emissfuelmap_second( const std::string& name ) const;
};

#endif // _SUMMARY_H_

