#ifndef REARR_OPTIONS_H
#define REARR_OPTIONS_H

! communication algorithm options
#define COLLECTIVE 0
#define POINT_TO_POINT 1
#define FLOW_CONTROL 2

! Default values for POINT_TO_POINT and FLOW_CONTROL
#define DEF_P2P_HANDSHAKE .true.
#define DEF_P2P_ISEND .false.
#define DEF_P2P_MAXREQ 64

#ifndef _MPISERIAL
#ifndef _NO_FLOW_CONTROL
#define _USE_COMP2IO_FC 1
#define _USE_IO2COMP_FC 1
#endif  
#endif
!
! The completely unreadable nature of the following lines is required by some compilers
!
#ifdef _USE_ALLTOALLW
#define DEF_COMP2IO_OPTION 0
#define DEF_IO2COMP_OPTION 0
#else
#ifdef _USE_COMP2IO_FC
#define DEF_COMP2IO_OPTION 2
#else
#define DEF_COMP2IO_OPTION 1
#endif
#ifdef _USE_IO2COMP_FC
#define DEF_IO2COMP_OPTION 2
#else
#define DEF_IO2COMP_OPTION 1
#endif
#endif

#endif
