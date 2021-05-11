/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_CONNECTIVITY_HELPERS_HPP
#define HOMMEXX_CONNECTIVITY_HELPERS_HPP

#include "Dimensions.hpp"
#include "Types.hpp"
#include "HommexxEnums.hpp"

#include <Kokkos_Array.hpp>

namespace Homme
{

// +--------------------------------------------------------------------------------------------------------------+
// |                                                                                                              |
// |                                  REMARKS ON BOUNDARY EXCHANGE                                                |
// |                                                                                                              |
// | Each element has 8 neighbors: west(W), east(E), south(S), north(N),                                          |
// |                               south-west(SW), south-east(SE), north-west(NW), north-east(NE)                 |
// | The first 4 correspond to the neighbors sharing a full edge with this element,                               |
// | while the latters correspond to the neighbors sharing only a corner with this element.                       |
// | NOTE: if the # of elements on each face of the cube sphere (ne*ne in homme) is greater than one, then        |
// | there are 24 elements that miss one corner neighbor. These are the element that touch one                    |
// | of the cube vertices (if ne=1, then all elements miss all the corner neighbors).                             |
// | The numeration of dofs on each element is the following                                                      |
// |                                                                                                              |
// |  (NW)     -(N)->   (NE)                                                                                      |
// |      12--13--14--15                                                                                          |
// |       |   |   |   |                                                                                          |
// |       |   |   |   |                                                                                          |
// |    ^  8---9--10--11  ^                                                                                       |
// |    |  |   |   |   |  |                                                                                       |
// |   (W) |   |   |   | (E)                                                                                      |
// |    |  4---5---6---7  |                                                                                       |
// |       |   |   |   |                                                                                          |
// |       |   |   |   |                                                                                          |
// |       0---1---2---3                                                                                          |
// |  (SW)     -(S)->   (SE)                                                                                      |
// |                                                                                                              |
// | The arrows around an edge neighbor refer to the 'natural' ordering of the dofs within an edge.               |
// | From the picture we can immediately state the following:                                                     |
// |                                                                                                              |
// |  1) edge neighbors contain 4 points/dofs, while corners contain only 1                                       |
// |  2) the S/N edges store the points contiguously, while the W/E points are strided (stride is NP)             |
// |     (Note: this is relevant only if we decide to switch to RMA, avoiding buffers, and copying                |
// |            data directly to/from host views)                                                                 |
// |                                                                                                              |
// | In addition, we need to understand that the local ordering of edge points may differ on two                  |
// | neighboring elements. For instance, consider two elements sharing their S edge. This is depicted             |
// | in the following:                                                                                            |
// |                                                                                                              |
// |          elem 1                                                                                              |
// |       0---1---2---3                                                                                          |
// |                                                                                                              |
// |       3---2---1---0                                                                                          |
// |          elem 2                                                                                              |
// |                                                                                                              |
// | So the 1st dof on the S edge of elem 1 does not correspond to the 1st dof on the S edge of elem 2.           |
// | This always happen for same-edge neighbors (W/W, E/E, S/S/, N/N), and also for W/N, E/S. For these           |
// | neighbors we need to store a flag marking the neighbor as 'backward' ordered. The other neighbors            |
// | (S/W, S/N, W/E, E/N) are marked as 'forward', meaning that the 1st dof on elem1's edge matches the           |
// | 1st dof on elem2's edge.                                                                                     |
// | NOTE: in F90, in this context, 'edge' means an edge of the dual mesh, i.e., a connection between elements.   |
// |       Here, we reserve the word 'edge' for the mesh edges, and we use 'connection' to refer to the 8         |
// |       possible connections with neighboring elements.                                                        |
// |                                                                                                              |
// | The W/E/S/N/SW/SE/NW/NE values from a neighboring element are stored on a buffer, to allow all MPI           |
// | operations to be over before the local view is updated. The values are then summed into the local field      |
// | view in a pre-established order, to guarantee reproducibility of the accumulation.                           |
// |                                                                                                              |
// +--------------------------------------------------------------------------------------------------------------+

// ============ Constexpr counters =========== //

constexpr int NUM_CONNECTION_KINDS     = 3;
constexpr int NUM_CONNECTION_SHARINGS  = 3;
constexpr int NUM_CONNECTIONS_PER_KIND = 4;

constexpr int NUM_CORNERS     = NUM_CONNECTIONS_PER_KIND;
constexpr int NUM_EDGES       = NUM_CONNECTIONS_PER_KIND;
constexpr int NUM_CONNECTIONS = NUM_CORNERS + NUM_EDGES;

// =========== A simple type for a Gauss Point =========== //

// A simple struct to store i,j indices of a gauss point. This is much like an std::pair,
// but with shorter and more meaningful member names than 'first' and 'second'.
// Note: we want to allow aggregate initialization, so no explitit constructors (and no non-static methods)!
struct GaussPoint
{
  int ip;   // i
  int jp;   // j
};
using ArrayGP = Kokkos::Array<GaussPoint,NP>;

// =========== A container struct for the information about connections =========== //

// Here we define a bunch of conxtexpr int's and arrays (of arrays (of arrays)) of ints, which we can
// use to easily retrieve information about a connection, such as the kind (corner or edge), the ordering
// on the remote (only relevant for edges), the (i,j) coordinates of the Gauss point(s) in the connection,
// and more.
struct ConnectionHelpers {

  ConnectionHelpers () {}

  ConnectionHelpers& operator= (const ConnectionHelpers& src) { return *this; }

  // Unpacking edges in the following order: S, N, W, E. For corners, order doesn't really matter
  const int UNPACK_EDGES_ORDER  [NUM_EDGES]   = { etoi(ConnectionName::SOUTH), etoi(ConnectionName::NORTH), etoi(ConnectionName::WEST),  etoi(ConnectionName::EAST) };
  const int UNPACK_CORNERS_ORDER[NUM_CORNERS] = { etoi(ConnectionName::SWEST), etoi(ConnectionName::SEAST), etoi(ConnectionName::NWEST), etoi(ConnectionName::NEAST)};

  const int CONNECTION_SIZE[NUM_CONNECTION_KINDS] = {
    NP,   // EDGE
    1,    // CORNER
    0     // MISSING (for completeness, but probably never used)
  };

  const ConnectionKind CONNECTION_KIND[NUM_CONNECTIONS] = {
      ConnectionKind::EDGE,     // S
      ConnectionKind::EDGE,     // N
      ConnectionKind::EDGE,     // W
      ConnectionKind::EDGE,     // E
      ConnectionKind::CORNER,   // SW
      ConnectionKind::CORNER,   // SE
      ConnectionKind::CORNER,   // NW
      ConnectionKind::CORNER    // NE
  };

  const Direction CONNECTION_DIRECTION[NUM_CONNECTIONS][NUM_CONNECTIONS] = {
    {Direction::BACKWARD, Direction::FORWARD , Direction::FORWARD,  Direction::BACKWARD, Direction::INVALID, Direction::INVALID, Direction::INVALID, Direction::INVALID}, // S/(S-N-W-E)
    {Direction::FORWARD,  Direction::BACKWARD, Direction::BACKWARD, Direction::FORWARD,  Direction::INVALID, Direction::INVALID, Direction::INVALID, Direction::INVALID}, // N/(S-N-W-E)
    {Direction::FORWARD,  Direction::BACKWARD, Direction::BACKWARD, Direction::FORWARD,  Direction::INVALID, Direction::INVALID, Direction::INVALID, Direction::INVALID}, // W/(S-N-W-E)
    {Direction::BACKWARD, Direction::FORWARD , Direction::FORWARD,  Direction::BACKWARD, Direction::INVALID, Direction::INVALID, Direction::INVALID, Direction::INVALID}, // E/(S-N-W-E)
    {Direction::INVALID,  Direction::INVALID,  Direction::INVALID,  Direction::INVALID,  Direction::FORWARD, Direction::FORWARD, Direction::FORWARD, Direction::FORWARD},
    {Direction::INVALID,  Direction::INVALID,  Direction::INVALID,  Direction::INVALID,  Direction::FORWARD, Direction::FORWARD, Direction::FORWARD, Direction::FORWARD},
    {Direction::INVALID,  Direction::INVALID,  Direction::INVALID,  Direction::INVALID,  Direction::FORWARD, Direction::FORWARD, Direction::FORWARD, Direction::FORWARD},
    {Direction::INVALID,  Direction::INVALID,  Direction::INVALID,  Direction::INVALID,  Direction::FORWARD, Direction::FORWARD, Direction::FORWARD, Direction::FORWARD}
  };

  // We only need 12 out of these 16, but for clarity, we define them all, plus an invalid one
  const GaussPoint GP_0       {  0,  0 };
  const GaussPoint GP_1       {  0,  1 };
  const GaussPoint GP_2       {  0,  2 };
  const GaussPoint GP_3       {  0,  3 };
  const GaussPoint GP_4       {  1,  0 };
  const GaussPoint GP_5       {  1,  1 };
  const GaussPoint GP_6       {  1,  2 };
  const GaussPoint GP_7       {  1,  3 };
  const GaussPoint GP_8       {  2,  0 };
  const GaussPoint GP_9       {  2,  1 };
  const GaussPoint GP_10      {  2,  2 };
  const GaussPoint GP_11      {  2,  3 };
  const GaussPoint GP_12      {  3,  0 };
  const GaussPoint GP_13      {  3,  1 };
  const GaussPoint GP_14      {  3,  2 };
  const GaussPoint GP_15      {  3,  3 };
  const GaussPoint GP_INVALID { -1, -1 };

  const ArrayGP SOUTH_PTS_FWD = {{ GP_0 , GP_1 , GP_2 , GP_3  }};
  const ArrayGP NORTH_PTS_FWD = {{ GP_12, GP_13, GP_14, GP_15 }};
  const ArrayGP WEST_PTS_FWD  = {{ GP_0 , GP_4 , GP_8 , GP_12 }};
  const ArrayGP EAST_PTS_FWD  = {{ GP_3 , GP_7 , GP_11, GP_15 }};

  const ArrayGP SOUTH_PTS_BWD = {{ GP_3 , GP_2 , GP_1 , GP_0  }};
  const ArrayGP NORTH_PTS_BWD = {{ GP_15, GP_14, GP_13, GP_12 }};
  const ArrayGP WEST_PTS_BWD  = {{ GP_12, GP_8 , GP_4 , GP_0  }};
  const ArrayGP EAST_PTS_BWD  = {{ GP_15, GP_11, GP_7 , GP_3  }};

  const ArrayGP SWEST_PTS = {{ GP_0 , GP_INVALID, GP_INVALID, GP_INVALID }};
  const ArrayGP SEAST_PTS = {{ GP_3 , GP_INVALID, GP_INVALID, GP_INVALID }};
  const ArrayGP NWEST_PTS = {{ GP_12, GP_INVALID, GP_INVALID, GP_INVALID }};
  const ArrayGP NEAST_PTS = {{ GP_15, GP_INVALID, GP_INVALID, GP_INVALID }};

  const ArrayGP NO_PTS = {{ }}; // Used as a placeholder later on

  // Now we pack all the connection points

  // Connections fwd
  const ArrayGP CONNECTION_PTS_FWD [NUM_CONNECTIONS] =
    { SOUTH_PTS_FWD, NORTH_PTS_FWD, WEST_PTS_FWD, EAST_PTS_FWD, SWEST_PTS, SEAST_PTS, NWEST_PTS, NEAST_PTS };

  // Connections bwd
  const ArrayGP CONNECTION_PTS_BWD [NUM_CONNECTIONS] =
    { SOUTH_PTS_BWD, NORTH_PTS_BWD, WEST_PTS_BWD, EAST_PTS_BWD, SWEST_PTS, SEAST_PTS, NWEST_PTS, NEAST_PTS };

  // All connections
  // You should never access CONNECTIONS_PTS with Direction=INVALID
  const ArrayGP CONNECTION_PTS[NUM_DIRECTIONS][NUM_CONNECTIONS] =
    {
      { SOUTH_PTS_FWD, NORTH_PTS_FWD, WEST_PTS_FWD, EAST_PTS_FWD, SWEST_PTS, SEAST_PTS, NWEST_PTS, NEAST_PTS },
      { SOUTH_PTS_BWD, NORTH_PTS_BWD, WEST_PTS_BWD, EAST_PTS_BWD, SWEST_PTS, SEAST_PTS, NWEST_PTS, NEAST_PTS },
      { NO_PTS }
    };

  // Edges and corners (fwd), used in the unpacking
  const ArrayGP EDGE_PTS_FWD [NUM_CONNECTIONS_PER_KIND] =
    { SOUTH_PTS_FWD, NORTH_PTS_FWD, WEST_PTS_FWD, EAST_PTS_FWD };

  const ArrayGP CORNER_PTS_FWD [NUM_CONNECTIONS_PER_KIND] =
    { SWEST_PTS, SEAST_PTS, NWEST_PTS, NEAST_PTS};
};

} // namespace Homme

#endif // HOMMEXX_CONNECTIVITY_HELPERS_HPP
