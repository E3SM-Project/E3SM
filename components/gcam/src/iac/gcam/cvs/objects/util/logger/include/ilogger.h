#ifndef _ILOGGER_H_
#define _ILOGGER_H_
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
* \file ilogger.h
* \ingroup objects
* \brief The ILogger interface header file.
* \author Josh Lurz
* \date $Date: 2007/01/11 23:52:34 $
* \version $Revision: 1.2.2.5 $
*/

#include <iostream>
#include <string>
/*! 
* \ingroup objects
* \brief This is an abstract class which defines the interface to a Logger.
*
* \author Josh Lurz
* \date $Date: 2007/01/11 23:52:34 $
* \version $Revision: 1.2.2.5 $
*/

class ILogger: public std::ostream {
public:
    //! Enumeration which describes the possible levels for messages.
    enum WarningLevel 
    {
        DEBUG, //!< Debugging warning level.
        NOTICE, //!< Notice warning level.
        WARNING, //!< Warning warning level.
        ERROR, //!< Error warning level.
        SEVERE //!< Severe warning level.
    };
    ILogger( std::streambuf* aStreamBuf ): std::ostream( aStreamBuf ){}
    virtual ~ILogger(){};
    virtual void open( const char[] = 0 ) = 0;
    virtual int receiveCharFromUnderStream( int ch ) = 0;
    virtual void close() = 0;
    virtual void setLevel( const WarningLevel newLevel ) = 0;
    static ILogger& getLogger( const std::string& aLoggerName );
};

#endif // _ILOGGER_H_

