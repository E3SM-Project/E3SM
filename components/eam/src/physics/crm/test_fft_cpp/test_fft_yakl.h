#pragma once
#include <iostream>
#include <stdio.h>
#include "YAKL.h"
#include "Array.h"
#include "YAKL_fft.h"

using yakl::c::Bounds;
using yakl::c::SimpleBounds;
using yakl::min;
using yakl::max;
using yakl::abs;
using yakl::c::parallel_for;
using yakl::SArray;
using yakl::ScalarLiveOut;

typedef double real;

typedef yakl::Array<real,1,yakl::memDevice,yakl::styleC> real1d;
typedef yakl::Array<real,2,yakl::memDevice,yakl::styleC> real2d;
typedef yakl::Array<real,3,yakl::memDevice,yakl::styleC> real3d;
typedef yakl::Array<real,4,yakl::memDevice,yakl::styleC> real4d;
typedef yakl::Array<real,5,yakl::memDevice,yakl::styleC> real5d;

typedef yakl::Array<int,1,yakl::memDevice,yakl::styleC> int1d;

int  constexpr nx    = 32;
int  constexpr ny    = 1;
int  constexpr nzm   = 1;
int  constexpr ncrms = 1;
int  constexpr YES3D = 0; // Domain dimensionality: 1 - 3D, 0 - 2D
int  constexpr RUN3D = ny > 1;

int filter_wn_max = 4; // MMF_VT_KMAX

void VT_filter(real4d &f_in, real4d &f_out);

