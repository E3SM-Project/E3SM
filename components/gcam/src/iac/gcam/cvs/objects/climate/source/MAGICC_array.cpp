/*
 *  MAGICC_array.cpp
 *  magicc++
 *
 *  Created by d3x290-local on 5/6/10.
 *  Copyright 2010 DOE Pacific Northwest Lab. All rights reserved.
 *
 */
#include <iostream>
#include <stdlib.h>
#include <string.h>

#include "climate/include/MAGICC_array.h"

using namespace std;

magicc_array::magicc_array():data( 0 )
{
    initialized = 0;    // Users must call init, below, to set up array(s)
}

magicc_array::magicc_array( const magicc_array& array )
{
    copy( array );
}

magicc_array::~magicc_array() {
    delete[] data;
}

magicc_array& magicc_array::operator=( const magicc_array& array ) {
    if( this != &array ) {
        delete[] data;
        copy( array );
    }
    return *this;
}

void magicc_array::copy( const magicc_array& array ) {
    initialized = array.initialized;
    low1 = array.low1;
    low2 = array.low2;
    high1 = array.high1;
    high2 = array.high2;
    strncpy( name, array.name, MA_NAMELEN );

    const int size = ( high1-low1+1 ) * ( high2-low2+1 );
    data = new float[ size ];
    for( int i = 0; i < size; ++i ) {
        data[ i ] = array.data[ i ];
    }
}

void magicc_array::init( const char* s, int l1, int h1, int l2, int h2 )
{
    if( l1 > h1 || l2 > h2 || initialized )    // Could get fancier and allow this
    {
        cout << "init error! " << s << " " << initialized << " " << l1 << " " << h1 << " " << l2 << " " << h2 << endl;
        exit( -1 );        
    }
    else
    {
        strncpy( name, s, MA_NAMELEN );
        low1 = l1; high1 = h1; low2 = l2; high2 = h2;
        data = new float[ ( high1-low1+1 ) * ( high2-low2+1 ) ];
        initialized = 1;
        for( int i=low1; i <=high1; i++ )
            for( int j=low2; j<=high2; j++)
                setval( 0.0, i, j );
    }
}

int magicc_array::computepos( int i1, int i2 )
{
    return i1-low1 + ( i2-low2 )*( high1-low1+1 );
}

void magicc_array::setval( float v, int i1, int i2 )
{
    if( !initialized || i1 < low1 || i1 > high1 || i2 < low2 || i2 > high2 )
    {
        cout << "setval error! " << name << " " << initialized << " " << i1 << " " << i2 << " " << v << endl;
        exit( -1 );
    }
    else
        data[ computepos( i1, i2 ) ] = v;
//    cout << name << ": writing " << v << " to " << i1 << " " << i2 << endl;
}

float magicc_array::getval( int i1, int i2 )
{
    if( !initialized || i1 < low1 || i1 > high1 || i2 < low2 || i2 > high2 )
    {
        cout << "getval error! " << i1 << " " << i2 << endl;
        exit( -1 );
    }
    else
        return data[ computepos( i1, i2 ) ];
}

float* magicc_array::getptr( int i1, int i2 )
{
    if( !initialized || i1 < low1 || i1 > high1 || i2 < low2 || i2 > high2 )
    {
        cout << "getptr error! " << i1 << " " << i2 << endl;
        exit( -1 );
    }
    else
        return &data[ computepos( i1, i2 ) ];
}

void magicc_array::print()
{
    cout << "-- DATA FOR " << name << " --" << endl;
    cout << "initialized=" << initialized << endl;
    cout << "low1=" << low1 << endl;
    cout << "high1=" << high1 << endl;
    cout << "low2=" << low2 << endl;
    cout << "high2=" << high2 << endl;
    
    for( int i=low1; i <=high1; i++ )
    {
        for( int j=low2; j<=high2; j++)
        {
            cout << getval( i, j ) << "  ";
        }
        cout << endl;
    }
}
