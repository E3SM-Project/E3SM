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
* \file plain_text_logger.cpp
* \ingroup Objects
* \brief PlainTextLogger class source file.
* \author Josh Lurz
* \date $Date: 2007/01/11 00:13:34 $
* \version $Revision: 1.5.2.4 $
*/

#include "util/base/include/definitions.h"
#include <iostream>
#include <string>
#include "util/logger/include/logger.h"
#include "util/logger/include/plain_text_logger.h"

using namespace std;

//! Constructor
PlainTextLogger::PlainTextLogger( const string& aLoggerName ):Logger( aLoggerName ){
}

//! Tells the logger to begin logging.
void PlainTextLogger::open( const char[] ){
    if( mFileName.empty() ) { // set a default value
        cout << "Using default log file name." << endl;
        mFileName = "log.txt";
    }

    mLogFile.open( mFileName.c_str(), ios::out );

    // Print the header message
    if( !mHeaderMessage.empty() ){
        parseHeader( mHeaderMessage );
        mLogFile << mHeaderMessage << endl << endl;
    }
}

//! Tells the logger to finish logging.
void PlainTextLogger::close(){
    mLogFile.close();
}

//! Logs a single message.
void PlainTextLogger::logCompleteMessage( const string& aMessage ){
    // Decide whether to print the message
    if ( mCurrentWarningLevel >= mMinLogWarningLevel ){
        // Print the warning level
        if ( mPrintLogWarningLevel || mCurrentWarningLevel >= ILogger::ERROR ) {
            mLogFile << convertLevelToString( mCurrentWarningLevel ) << ":";
        }
        mLogFile << aMessage << endl;
    }
}
