(omega-user-horz-operators)=

# Horizontal Operators

The horizontal discretization in Omega is a staggered numerical scheme known as
TRiSK. It defines discrete versions of basic differential operators (divergence,
gradient, and curl) as well as other operators that are needed by the scheme,
such as reconstruction of tangential velocity from its normal components. Omega
provides reference implementations of these operators, each in a separate C++
class. The class name describes the operator and the mesh element associated
with its result.

The following operators are currently implemented:
- `DivergenceOnCell`
- `GradientOnEdge`
- `CurlOnVertex`
- `TangentialReconOnEdge`
