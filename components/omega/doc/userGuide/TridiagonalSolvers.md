(omega-user-tridiagonal)=

## Tridiagonal Solvers

OMEGA provides two types of solvers for solving many tridiagonal systems of equations
in one go:
- `TriDiagSolver` for general tridiagonal systems
- `TriDiagDiffSolver` for diffusion-type tridiagonal systems

Different solver algorithms are used on CPU and GPU platforms.
The optimal solver algorithm for a given platform is chosen automatically.
There are no user-configurable parameters.
The [dev guide](#omega-dev-tridiagonal) provides instructions on how to use the solvers.
