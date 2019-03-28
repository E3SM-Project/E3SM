#ifndef _MAGICC_array_H_
#define _MAGICC_array_H_

/*
 *  MAGICC_array.h
 *  magicc++
 *
 *  Created by d3x290-local on 5/6/10.
 *  Copyright 2010 DOE Pacific Northwest Lab. All rights reserved.
 *
 */
 

/*  MAGICC uses a whole series of arrays to track data with variable indices
    (typically from 226 to iTp). The magicc_array class mimics this functionality
    for us in C++, allowing a one- or two-dimensional array with arbitrary
    index ranges.
*/

#define MA_NAMELEN 20

class magicc_array {
    int initialized, low1, high1, low2, high2;
    float* data;
    char name[ MA_NAMELEN ];
private:
    int computepos( int, int );
    void copy( const magicc_array& array );
public:
    magicc_array();
    magicc_array( const magicc_array& array );
    ~magicc_array();
    magicc_array& operator=( const magicc_array& array );

    void init( const char*, int, int, int=0, int=0 );
    void setval( float, int, int=0 );
    float getval( int, int=0 );    
    float* getptr( int, int=0 );    
    void print();
};


#endif // _MAGICC_array_H_
