/*
 *  MAGICC_IO_helpers.cpp
 *  magicc++
 *
 *  Created by d3x290-local on 10/9/09.
 *  Copyright 2009 DOE Pacific Northwest Lab. All rights reserved.
 *
 */

#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>

using namespace std;

// Small helper functions related to file I/O

void openfile_read( ifstream* infile, const string& f, bool echo )
{
    (*infile).open( f.c_str(), ios::in );
    if ( !infile ) {
        cerr << "Unable to open file " << f << " for read\n";
        exit( 1 ); 
    }
    if ( echo ) cout << "Opened file " << f << " for read OK\n";
}

void skipline( ifstream* infile, bool echo )
{
    string line;
    getline( *infile, line );
    
    if( echo ) cout << "Skipping line: " << line << endl;
}

float read_csv_value( ifstream* infile, bool echo )
{
    float f;
    char comma;
    (*infile) >> f >> comma;
    
    if ( echo ) cout << "Read " << f << ",";
    return f;
}

float read_and_discard( ifstream* infile, bool echo )
{
    float f;
    string line;
    
    getline( *infile, line );
    if (EOF == sscanf( line.c_str(), "%f", &f ) ) {
        cerr << "Error in reading data from line: " << line << "\n";
        exit( 1 );
    }
    if ( echo ) cout << "Read value " << f << " from line: " << line << endl;
    return f;
}

void openfile_write( ofstream* outfile, const string& f, bool echo )
{
    (*outfile).open( f.c_str(), ios::out );
    if ( !outfile ) {
        cerr << "Unable to open file " << f << " for write\n";
        exit( 1 ); 
    }
    if ( echo ) cout << "Opened file " << f << " for write OK\n";
}


