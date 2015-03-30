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

#endif
