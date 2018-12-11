#ifndef _INDIRECT_EMISS_COEF_H_
#define _INDIRECT_EMISS_COEF_H_
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
* \file indirect_emiss_coef.h
* \ingroup Objects
* \brief The Emcoef_ind class header file.
* \author Sonny Kim
* \date $Date: 2007/01/11 00:12:29 $
* \version $Revision: 1.2.2.2 $
*/

#include <map>
#include <string>

/*! 
* \ingroup Objects
* \brief This class contains the indirect Greenhouse gas emissions coefficients.
*
* Each object contains a map object of emissions
* coefficients, and each emisscoef_ind object is
* used to represent the secondary energy sector.
* \author Sonny Kim
*/

class Emcoef_ind
{
private:
    std::string name; //!< name of secondary good or sector
    std::map<std::string, double> emcoef; //!< contains all coefficients for all gases
public:
    Emcoef_ind( const std::string sectorName );
    void setemcoef( const std::map<std::string,double>& eminfo, const double toutput );
    const std::string& getName() const;
    double getemcoef( const std::string& gasName ) const;
};

#endif // _INDIRECT_EMISS_COEF_H_

