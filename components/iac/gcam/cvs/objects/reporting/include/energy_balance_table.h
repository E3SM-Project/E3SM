#ifndef _ENERGY_BALANCE_TABLE_H_
#define _ENERGY_BALANCE_TABLE_H_
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
* \file energy_balance_table.h
* \ingroup Objects
* \brief EnergyBalanceTable class header file.
* \author Pralit Patel
*/

#include <string>
#include <map>
#include <memory>
#include <iosfwd>
#include "util/base/include/default_visitor.h"

class Sector;
class Subsector;
class Technology;
class MiniCAMInput;
class IOutput;
class EnergyFinalDemand;
class StorageTable;

/*! 
 * \ingroup Objects
 * \brief An object which outputs an energy balance table.
 * \details Sums the calibration values by production and consumption to help
 *          users better understand why the model was or was not balanced for
 *          calibration.  The table can be constructed in the following ways
 *          each with their own benefits.  If mPrintCondensed is set to true
 *          all subsectors and technologies get collapsed which is useful for 
 *          when checking if a good is balanced.  When mPrintCondensed is false
 *          a complete table is printed which is useful to check each of your
 *          calibration values.  You can also set mIncludeNonCalibrated values
 *          which will let you see when something is not balanced when you have
 *          mixed calibrated and non-calibrated values.
 *
 * \todo CalibrationBalanceTable might be a better name for this class.
 * \author Pralit Patel
 */
class EnergyBalanceTable : public DefaultVisitor {
public:
    EnergyBalanceTable( const std::string& aRegionName, std::ostream& aFile, const bool aPrintCondensed,
                        const bool aIncludeNonCalValues );
    virtual ~EnergyBalanceTable();
    void finish() const;
    void startVisitSector( const Sector* aSector, const int aPeriod );
    void endVisitSector( const Sector* aSector, const int aPeriod );
    void startVisitSubsector( const Subsector* aSubsector, const int aPeriod );
    void endVisitSubsector( const Subsector* aSubsector, const int aPeriod );
    void startVisitTechnology( const Technology* aTechnology, const int aPeriod );
    void endVisitTechnology( const Technology* aTechnology, const int aPeriod );
    void startVisitMiniCAMInput( const MiniCAMInput* aInput, const int aPeriod );
    void startVisitOutput( const IOutput* aOutput, const int aPeriod );
    void startVisitEnergyFinalDemand( const EnergyFinalDemand* aEnergyFinalDemand, const int aPeriod );
    void startVisitResource( const AResource* aResource, const int aPeriod );
    void endVisitResource( const AResource* aResource, const int aPeriod );
    void startVisitSubResource( const SubResource* aSubResource, const int aPeriod );
    void endVisitSubResource( const SubResource* aSubResource, const int aPeriod );
    void setRegionName( const std::string& aRegionName );
private:
    //! Table to store our results
    std::auto_ptr<StorageTable> mTable;

    //! Name of the current region
    std::string mRegionName;

    //! Name of the current sector
    std::string mCurrentSector;

    //! Name of the current subsector
    std::string mCurrentSubsector;

    //! Name of the current technology being visited.
    std::string mCurrentTech;

    //! Map of a column name to the sector it belongs to.
    std::map<std::string, std::string> mColToSectorMap;

    //! Map of a column name to the subsector it belongs to.
    std::map<std::string, std::string> mColToSubsectorMap;

    //! Map of a column name to tehcnology name
    std::map<std::string, std::string> mColToTechMap;

    //! Whether to visit the input and output.
    bool mIsTechOperating;

    //! File to which to write.
    std::ostream& mFile;

    //! the current cal output form the technology
    double mCalOutput;

    //! Whether we should print the condensed version of the table
    const bool mPrintCondensed;

    //! Whether we should include non-calibrated supplies and demands
    //! when a calibrated value was not supplied
    const bool mIncludeNonCalValues;

    void writeFullTable() const;

    void writeCondensedTable() const;

    std::string getKey() const;
    
    double getTotalSectorOutput( const std::string& aSectorName ) const;
};

#endif // _ENERGY_BALANCE_TABLE_H_
