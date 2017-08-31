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
* \file logger.cpp
* \ingroup Objects
* \brief Logger class source file.
* \author Josh Lurz
* \date $Date: 2007/01/11 00:13:34 $
* \version $Revision: 1.7.2.7 $
*/

#include "util/base/include/definitions.h"
#include <string>
#include <iostream>
#include <sstream>
#include <cassert>
#include <ctime>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include "util/logger/include/logger.h"
#include "util/base/include/xml_helper.h"

using namespace std;
using namespace xercesc;

//! Default Constructor
PassToParentStreamBuf::PassToParentStreamBuf():
mParent( 0 ){
}

//! Overriding overflow function which passes its argument to its parent.
int PassToParentStreamBuf::overflow( int aChar ){
	/*! \pre Make sure the parent is not null. */
	assert( mParent );
	return mParent->receiveCharFromUnderStream( aChar );
}

//! Overriding underflow function which should not be reached because this is a write-only stream.
int PassToParentStreamBuf::underflow( int aChar ){
	/*! \pre This function should never be called. */
	assert( false );
	return 0;
}

//! Set the parent of the stream to which we will pass all data.
void PassToParentStreamBuf::setParent( Logger* aParent ) {
	/*! \pre Make sure a non-null parent is passed in. */
	assert( aParent );
	mParent = aParent;
}

//! Constructor which sets default values.
Logger::Logger( const string& aFileName ):
ILogger( &mUnderStream ),
// Initialize all variables which are not set by Configuration values.
mCurrentWarningLevel( ILogger::DEBUG ),
mFileName( aFileName ),
mMinLogWarningLevel( ILogger::DEBUG ),
mMinToScreenWarningLevel( ILogger::SEVERE ),
mPrintLogWarningLevel( false ){
    // Set the understream's parent to this Logger.
	mUnderStream.setParent( this );
}

//! Virtual destructor
Logger::~Logger() {
}

//! Set the current warning level.
void Logger::setLevel( const ILogger::WarningLevel aLevel ){
	mCurrentWarningLevel = aLevel;
}

//! Receive a single character from the underlying stream and buffer it, printing the buffer it is a newline.
int Logger::receiveCharFromUnderStream( int ch ) {
    // Only receive the character or print to the screen if it needed.
    if( mCurrentWarningLevel >= mMinLogWarningLevel || mCurrentWarningLevel >= mMinToScreenWarningLevel ){
        if( ch == '\n' ){
            char tempBuf[ MAX_LINE_SIZE ];
            mBuf.get( tempBuf, MAX_LINE_SIZE );
            string tempString( tempBuf, mBuf.gcount() );
            mBuf.clear(); // The get hits an EOF, so we need to clear the failure flag.
            logCompleteMessage( tempString );
            printToScreenIfConfigured( tempString );
       }
        else {
            mBuf << (char)ch;
        }
    }
    return ch;
}

//! Print the message to the screen if the Logger is configured to.
void Logger::printToScreenIfConfigured( const string& aMessage ){
	// Decide whether to print the message
	if ( mCurrentWarningLevel >= mMinToScreenWarningLevel ) {
		// Print the warning level
		if ( mPrintLogWarningLevel || mCurrentWarningLevel >= ILogger::ERROR ) {
            cout << convertLevelToString( mCurrentWarningLevel ) << ":";
		}
		cout << aMessage << endl;
	}
}

void Logger::XMLParse( const DOMNode* aNode ){
	// assume we were passed a valid node.
	assert( aNode );
	
	mName = XMLHelper<string>::getAttr( aNode, "name" );
	mType = XMLHelper<string>::getAttr( aNode, "type" );
	// get the children of the node.
	DOMNodeList* nodeList = aNode->getChildNodes();

	// loop through the children
	for ( unsigned int i = 0; i < nodeList->getLength(); ++i ){
		DOMNode* curr = nodeList->item( i );
		string nodeName = XMLHelper<string>::safeTranscode( curr->getNodeName() );
		
		if ( nodeName == "FileName" ){
			mFileName = XMLHelper<string>::getValue( curr );
		}
		else if ( nodeName == "printLogWarningLevel" ) {
			mPrintLogWarningLevel = XMLHelper<bool>::getValue( curr );
		}
		// These casts are unsafe.
		else if ( nodeName == "minLogWarningLevel" ) {
			mMinLogWarningLevel = static_cast<WarningLevel>( XMLHelper<int>::getValue( curr ) );
		}
		else if ( nodeName == "minToScreenWarningLevel" ) {
			mMinToScreenWarningLevel = static_cast<WarningLevel>( XMLHelper<int>::getValue( curr ) );
		}
		else if ( nodeName == "headerMessage" ) {
			mHeaderMessage = XMLHelper<string>::getValue( curr );
		}
	}
}

void Logger::toDebugXML( ostream& out, Tabs* tabs ) const {
	out << "<Logger name=\"" << mName << "\" type=\"" << mType << "\" >" << endl;
	tabs->increaseIndent();
	XMLWriteElement( mFileName, "fileName", out, tabs );
	XMLWriteElement( mMinLogWarningLevel, "minLogWarningLevel", out, tabs );
	XMLWriteElement( mMinToScreenWarningLevel, "minToScreenWarningLevel", out, tabs );
	XMLWriteElement( mPrintLogWarningLevel, "printLogWarningLevel", out, tabs );
	XMLWriteClosingTag( "Logger", out, tabs );
}

//! Parses the header of a log file replacing special strings.
void Logger::parseHeader( string& aHeader ) {
	
	static const basic_string <char>::size_type npos = static_cast<char>( -1 );
	int offset = 0;

	// Loop through the string.
	while( offset < aHeader.size() && offset != npos ){
		
		// Find the first left bracket.
		int leftBracket = static_cast<unsigned int>( aHeader.find_first_of( "{", offset ) );
		offset = leftBracket;

		// Exit if we do not find it.
		if( leftBracket == npos ){
			break;
		}
		
		int rightBracket = rightBracket = static_cast<unsigned int>( aHeader.find_first_of( "}", offset ) );
		offset = rightBracket;

		// Exit if we do not find it.
		if( rightBracket == npos ){
			break;
		}

		const string command = aHeader.substr( leftBracket + 1, rightBracket - leftBracket - 1 );
		
		string replaceWithString;
		if( command == "date" ) {
			replaceWithString = getDateString();
		}
		else if( command == "time" ) {
			replaceWithString = getTimeString();
		}
		else {
			continue;
		}

		aHeader.replace( leftBracket, rightBracket - leftBracket + 1, replaceWithString );
	}
}

//! Get the current date as a string
const string Logger::getDateString(){
	time_t currTime;
	time( &currTime );
    tm* timeInfo = util::getLocalTime( currTime );

	// Create the string
	stringstream buffer;
	buffer << ( timeInfo->tm_year + 1900 ); // Set the year
	buffer << "-";
	
	if( timeInfo->tm_mday < 10 ){
		buffer << "0";
	}
	
	buffer << timeInfo->tm_mday; // Set the day
	buffer << "-";

	if( timeInfo->tm_mon + 1 < 10 ) {
		buffer << "0";
	}
	buffer << ( timeInfo->tm_mon + 1 ); // Month's in ctime range from 0-11
	
	string retString;
	buffer >> retString;

	return retString;
}

//! Get the current time as a string.
const string Logger::getTimeString() {
	time_t currTime;	
	time( &currTime );
	tm* timeInfo = util::getLocalTime( currTime );

	stringstream buffer;	
	if( timeInfo->tm_hour < 10 ){
		buffer << "0";
	}
	buffer << timeInfo->tm_hour;
	buffer << ":";
	
	if( timeInfo->tm_min < 10 ){
		buffer << "0";
	}
	buffer << timeInfo->tm_min;
	buffer << ":";
	
	if( timeInfo->tm_sec < 10 ){
		buffer << "0";
	}
	buffer << timeInfo->tm_sec;
	
	string retString;	
	buffer >> retString;
	return retString;
}

/*! \brief Convert a warning level enum into the corresponding string.
* \param aLevel Warning level to convert.
* \return A string describing the warning level.
*/
const string& Logger::convertLevelToString( WarningLevel aLevel ){
	// Create a static array of warning level strings in the order
	// of the enum.
	const static string levelDescriptions[] = 
	{ "Debug", "Notice", "Warning", "ERROR", "SEVERE ERROR" };

	// Use the enum to index into the array.
	return levelDescriptions[ aLevel ];
}
