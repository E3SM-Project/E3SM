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

#ifndef _CAPACITY_LIMIT_BACKUP_CALCULATOR_H_
#define _CAPACITY_LIMIT_BACKUP_CALCULATOR_H_
#if defined(_MSC_VER)
#pragma once
#endif

/*!
 * \file capacity_limit_backup_calculator.h
 * \ingroup Objects
 * \brief The CapacityLimitBackupCalculator class header file.
 * \author Josh Lurz
 */

#include <string>
#include "sectors/include/ibackup_calculator.h"

// Forward declaration
class IInfo;

/*!
 * \ingroup Objects
 * \brief The backup calculator for a capacity limited resource.
 * \details Calculates the amount of backup required per unit of output. The marginal
 *          backup requirement is computed as an exponential function with parameters that
 *          specify a maximum backup requirement, the point at which half of the maximum 
 *          backup is required, and the steepness of the curve in reaching the maximum
 *          backup requirement<br>
 *
 *          functional form of marginal backup requirement is
 *          
 *          backupCapacity = mFmax / (1.0 + exp( mC * ( xmid - elecShare ) / mTau ))
 *                  where  mFmax is the maximum backup fraction (defaults to 1)
 *                         mC and mTau are shape parameters
 *                                  (ratio of mC to mTau determined steepness to reach fMax)
 *                         xmid is the value of the share at which backupCapacity = 0.5*fMax
 *                                
 *
 *          The integral of this function, which is used to compute total and average bakcup
 *            is given by
 *
 *          totalBackup = mFmax * mTau / mC * ( log ( 1.0 + exp( c * ( xmid - elecShare ) / mTau))
 *                        - mC * ( xmid - elecShare ) / mTau)
 *
 *
 *          <b>XML specification for CapacityLimitBackupCalculator</b>
 *          - XML name: \c capacity-limit-backup-calculator
 *          - Contained by: Technology
 *          - Parsing inherited from class: None
 *          - Attributes: None
 *          - Elements:
 *              - \c capacity-limit CapacityLimitBackupCalculator::mCapacityLimit
 *
 * \author Josh Lurz, Marshall Wise
 */
class CapacityLimitBackupCalculator: public IBackupCalculator {
    friend class BackupCalculatorFactory;
public:
    virtual CapacityLimitBackupCalculator* clone() const;
    virtual bool isSameType( const std::string& aType ) const;
    virtual const std::string& getName() const;
    virtual bool XMLParse( const xercesc::DOMNode* aNode );
    virtual void toInputXML( std::ostream& aOut, Tabs* aTabs ) const;
    virtual void toDebugXML( const int aPeriod, std::ostream& aOut, Tabs* aTabs ) const;
    virtual void initCalc( const IInfo* aTechInfo );
    
    virtual double getMarginalBackupCapacity( const std::string& aSector,
                                              const std::string& aElectricSector,
                                              const std::string& aResource,
                                              const std::string& aRegion,
                                              const double aReserveMargin,
                                              const double aAverageGridCapacityFactor,
                                              const int aPeriod ) const;
    
    virtual double getAverageBackupCapacity( const std::string& aSector,
                                             const std::string& aElectricSector,
                                             const std::string& aResource,
                                             const std::string& aRegion,
                                             const double aReserveMargin,
                                             const double aAverageGridCapacityFactor,
                                             const int aPeriod ) const;
protected:
    static const std::string& getXMLNameStatic();
    CapacityLimitBackupCalculator();

    double getMarginalBackupCapacityFraction( const std::string& aSector,
                                              const std::string& aElectricSector,
                                              const std::string& aResource,
                                              const std::string& aRegion,
                                              const double aReserveMargin,
                                              const double aAverageGridCapacityFactor,
                                              const int aPeriod ) const;

    double calcIntermittentShare( const std::string& aSector,
                                  const std::string& aElectricSector,
                                  const std::string& aResource,
                                  const std::string& aRegion,
                                  const double aReserveMargin,
                                  const double aAverageGridCapacityFactor,
                                  const int aPeriod ) const;

    //! Capacity limit which sets the fraction of total output at which the backup curve
    //! returns a value of 50% of the upper limit for backup (i.e.,50 of mFmax)
    double mCapacityLimit;

    //! Parameter for max rate of back up (e.g. specify 1 for 1-to-1 backup as max)
    double mFmax;

    //! Parameter for steepness of backup curve. Higher number means steeper ascent.
    double mC;
   
    //! Parameter for steepness of backup curve. Lower means steeper ascent. Ascent depends
    //! on the ratio of c/tau
    double mTau;
  
};

#endif // _CAPACITY_LIMIT_BACKUP_CALCULATOR_H_
