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

#ifndef _CSP_BACKUP_CALCULATOR_H_
#define _CSP_BACKUP_CALCULATOR_H_
#if defined(_MSC_VER)
#pragma once
#endif

/*!
 * \file CSP_backup_calculator.h
 * \ingroup Objects
 * \brief The CSPBackupCalculator class header file.
 * \author Marshall Wise
 */

#include <string>
#include "sectors/include/ibackup_calculator.h"
/*!
 * \ingroup Objects
 * \brief The backup calculator for the CSP (concentrated solar trough) technology.
 * \details Calculates the amount of backup energy required per unit of output. This
 *          differs from the capacity_limit_backup_calculator in that the CSP
 *          backup retunrs energy reuqired rather than capacity.<br>
 *
 *
 *          <b>XML specification for CSPBackupCalculator</b>
 *          - XML name: \c CSP-backup-calculator
 *          - Contained by: IntermittentSubsector
 *          - Parsing inherited from class: None
 *          - Attributes: None
 *          - Elements: None
 *
 * \author Marshall Wise
 */
class CSPBackupCalculator: public IBackupCalculator {
    friend class BackupCalculatorFactory;
public:
    virtual CSPBackupCalculator* clone() const;
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
    CSPBackupCalculator();

    double calcIntermittentShare( const std::string& aSector,
                                  const std::string& aElectricSector,
                                  const std::string& aResource,
                                  const std::string& aRegion,
                                  const double aReserveMargin,
                                  const double aAverageGridCapacityFactor,
                                  const int aPeriod ) const;

    
    // Parameters hard-coded in constructor for first step
    //! Maximum Backup Fraction is the maximum backup fraction required.
    double mMaxBackupFraction;
     //! Fraction of electric sector that is intermediate and peak
    double mMaxSectorLoadServed;
    //! Backup function exponent parameter
    double mBackupExponent;
    //! Fraction of year backup mode operates (during no-sun days)
    double mNoSunDayBackup;
    //! Backup fraction -- cached for reporting
    mutable double mBackupFraction;

};

#endif // _CSP_BACKUP_CALCULATOR_H_
