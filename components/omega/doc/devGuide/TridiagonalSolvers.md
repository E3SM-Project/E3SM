(omega-dev-tridiagonal)=

# Tridiagonal Solvers

Omega provides batched solvers for tridiagonal systems of equations.
Here, batched means that many systems of the same size need to be solved at the same time.
Typically, this is because one system needs to be solved per ocean column.
In addition to providing solvers for general tridiagonal systems,
specialized solvers for diffusion-type systems with better stability
properties are available.
To obtain portable performance, different solution algorithms are utilized
on CPUs and GPUs.
The solvers implement the same user interface for CPUs and GPUs, so users
should be able to ignore which algorithm is used.
Two type aliases are provided that resolve to the optimal solver for a given
architecture:
- `TriDiagSolver` for general systems
- `TriDiagDiffSolver` for diffusion-type systems

## Representations of tridiagonal systems

### General tridiagonal system
A general tridiagonal system is represented by three vectors
`DL`, `D`, `DU`, representing the lower, main, and upper diagonal, respectively, and a right-hand-side vector `X`
$$
(DL)_i Y_{i - 1} + D_i Y_i + (DU)_i Y_{i + 1} = X_i.
$$

### Diffusion-type tridiagonal system

A diffusion-type system is represented by two vectors `G` and `H`, and a right-hand-side vector `X`. It has the form
$$
-G_{i - 1} Y_{i - 1} + (G_{i - 1} + G_i + H_i) Y_i - G_{i} Y_{i + 1} = X_i.
$$

## Using tridiagonal solvers

There are two main ways to use Omega batched tridiagonal solvers:
- Create global arrays storing the coefficients of all systems and call the top-level `solve` method
- Use a team-level solver inside a parallel loop using Kokkos team policy

While the first method is simpler, for computational performance it is better to use
the second method. This is because the second method allows one to define a system, solve it,
and do some post-processing of the solution in one parallel loop.
With the first method this would require multiple parallel loops, and allocation of storage
for the coefficients of all systems, which are typically no longer needed after the solution is obtained.

### Top-level method

To solve `NBatch` systems of size `NRow` using the general solver
four arrays `DL`, `D`, `DU`, and `X` need to be defined.
All arrays are of size (`NBatch`, `NRow`) and contain the
coefficients for each of the systems.
```c++
Array2DReal DL("DL", NBatch, NRow); // DL(n, :) = lower diagonal of system n
Array2DReal D("D", NBatch, NRow);   // D(n, :)  = main diagonal of system n
Array2DReal DU("DU", NBatch, NRow); // DU(n, :) = upper diagonal of system n
Array2DReal X("X", NBatch, NRow);   // X(n, :)  = rhs of system n
```
Once these arrays are filled with the coefficients, to solve the systems the
static `TriDiagSolver:solve` method is called
```
TriDiagSolver::solve(DL, D, DU, X);
```
After this call, the rhs array `X` contains the solution.

Solving diffusion-type systems is similar, except the input arrays contain coefficients in the
form presented [above](#diffusion-type-tridiagonal-system)
and the static `TriDiagDiffSolver::solve` method is called
```c++
Array2DReal G("G", NBatch, NRow); // G(n, :) = coefficients G_i of system n
Array2DReal H("H", NBatch, NRow); // H(n, :) = coefficients H_i of system n
Array2DReal X("X", NBatch, NRow); // X(n, :) = rhs of system n

// fill G, H, X

TriDiagDiffSolver::solve(G, H, X);
```

### Team-level method

Note: since the steps to use the team-level general and specialized solvers are similar,
this subsection first shows how to use the general solver.
At the end, an example of using the specialized solver is presented.

The team-level solvers need to be used inside a parallel loop using Kokkos team
policy. To create a team policy for solving `NBatch` systems of size `NRow`
the static member function `makeTeamPolicy` is used
```c++
TeamPolicy Policy = TriDiagSolver::makeTeamPolicy(NBatch, NRow);
```
Every team of threads except the last is responsible for solving `VecLength` systems.
Since the total number of systems likely doesn't evenly divide `VecLength`, the last
team solves fewer than `VecLength` systems. Hence, it is necessary to handle that with
if statements.

Typically, the parallel loop consists of three main parts:
- defining the system
- solving the system
- using the solution

At the beginning of the loop, temporary scratch storage for the coefficients of `VecLength` systems
needs to be allocated.
This is done by creating a `TriDiagScratch` struct, which contains
four member arrays `DL`, `D`, `DU`, and `X` of size (`NRow`, `VecLength`).
Note that, contrary to the global arrays approach, the system dimension comes first.

These arrays need to be filled with the coefficients of `VecLength` systems solved by the team.
This is done in a parallel loop over system size (using `TeamThreadRange`)
and a serial loop over `VecLength`.
Overall, the code that defines the systems coefficients looks typically like this:
```c++
parallel_for(Policy, KOKKOS_LAMBDA (TeamMember &member) {
   // create scratch data
   TriDiagScratch Scratch(Member, NRow);

   int IStart = Member.league_rank() * VecLength;

   // define the systems coefficients
   parallel_for(TeamThreadRange(Member, NRow), [=] (int K) {
      for (int IVec = 0; IVec < VecLength; ++IVec) {
         const int I = IStart + IVec;
         // handle last team having fewer systems
         if (I < NBatch) {
            Scratch.DL(K, IVec) = ...;
            Scratch.D(K, IVec)  = ...;
            Scratch.DU(K, IVec) = ...;
            Scratch.X(K, IVec)  = ...;
         }
      }
      ...
   });
```
To perform a team-level solve a different overload of the static member function `TriDiagSolver::solve`
needs to be called.
This overload takes the team member `Member` and the filled scratch struct `Scratch`.
This call needs to be surrounded by barriers to ensure that the
parallel loop setting the coefficients has finished and to make the computed solution available to all threads.
```c++
parallel_for(Policy, KOKKOS_LAMBDA (TeamMember &member) {
   // previous code

   // solve the system
   Member.team_barrier();
   TriDiagSolver::solve(Member, Scratch);
   Member.team_barrier();

   ...
});
```
After the call to `solve` the solution is available in `Scratch.X`. It can be manipulated further or simply stored.
```cpp
parallel_for(Policy, KOKKOS_LAMBDA (TeamMember &member) {
   // previous code

   // multiply the solution by 2 (just an example)
   parallel_for(TeamThreadRange(Member, NRow), [=] (int K) {
      for (int IVec = 0; IVec < VecLength; ++IVec) {
        Scratch.X(K, IVec) *= 2;
      }
   });

   Member.team_barrier();

   // store the solution
   parallel_for(TeamThreadRange(Member, NRow), [=] (int K) {
      for (int IVec = 0; IVec < VecLength; ++IVec) {
         const int I = IStart + IVec;
         // handle last team having fewer systems
         if (I < NBatch) {
            X(I, K) = Scratch.X(K, IVec);
         }
      }
   });
});
```
The team-level specialized diffusion solver can be used similarly, the only substantial difference is in the form
of the coefficients. A minimal complete example of using it is shown below.
```cpp
// create team policy
TeamPolicy Policy = TriDiagDiffSolver::makeTeamPolicy(NBatch, NRow);

parallel_for(Policy, KOKKOS_LAMBDA (TeamMember &member) {
   // create scratch data
   TriDiagDiffScratch Scratch(Member, NRow);

   int IStart = Member.league_rank() * VecLength;

   // define the systems coefficients
   parallel_for(TeamThreadRange(Member, NRow), [=] (int K) {
      for (int IVec = 0; IVec < VecLength; ++IVec) {
         const int I = IStart + IVec;
         // handle last team having fewer systems
         if (I < NBatch) {
            Scratch.G(K, IVec) = ...;
            Scratch.H(K, IVec)  = ...;
            Scratch.X(K, IVec)  = ...;
         }
      }

   // solve the system
   Member.team_barrier();
   TriDiagDiffSolver::solve(Member, Scratch);
   Member.team_barrier();

   // store the solution
   parallel_for(TeamThreadRange(Member, NRow), [=] (int K) {
      for (int IVec = 0; IVec < VecLength; ++IVec) {
         const int I = IStart + IVec;
         // handle last team having fewer systems
         if (I < NBatch) {
            X(I, K) = Scratch.X(K, IVec);
         }
      }
   });
});
