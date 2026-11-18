# Tridiagonal Solver

<!--- use table of contents if desired for longer documents  -->
**Table of Contents**
1. [Overview](#1-overview)
2. [Requirements](#2-requirements)
3. [Algorithmic Formulation](#3-algorithmic-formulation)
4. [Design](#4-design)
5. [Verification and Testing](#5-verification-and-testing)

## 1 Overview

In geophysical models, implicit time integration of vertical terms often requires solution to tridiagonal systems of equations.
Typically, a separate system needs to be solved in each vertical column, which requires an efficient batched tridiagonal solver.
One common situation where a tridiagonal system arises is implicit treatment of vertical diffusion/mixing.
In principle, this problem results in a symmetric diagonally-dominant tridiagonal system.
However, in ocean modeling, the coefficients of vertical mixing can vary by orders of magnitude and
can become very large in some layers.
This requires specialized algorithms that can handle this situation stably.

## 2 Requirements

### 2.1 Requirement: Modularity

There should be one module that provides tridiagonal solvers that are sufficiently general to handle every need in Omega.

### 2.2 Requirement: Stability for vertical mixing problems

A specialized solver that is stable for vertical mixing problems with large variations in mixing coefficients must be
provided.

### 2.3 Requirement: Top and bottom boundary conditions

When applied to the implicit vertical mixing of momentum, the specialized diffusion solver should be
sufficiently general to be able to incorporate various bottom drag and wind stress formulations.

### 2.4 Requirement: Performance

All solvers must be performant on CPU and GPU architectures. On CPUs this means supporting vectorization.

### 2.5 Desired: Ability to fuse kernels

Implicit vertical mixing will require some pre-processing (e.g. setup of the system) and post-processing work.
It is desirable to handle all of that in one computational kernel. This requires the ability to call the
solvers inside a `parallelFor`.

## 3 Algorithmic Formulation
A general tridiagonal system has the form:
$$
a_i x_{i - 1} + b_i x_i + c_i x_{i + 1} = y_i.
$$
The standard second-order discretization of a diffusion problem leads to a system of the form
$$
-g_{i - 1} x_{i - 1} + (g_{i - 1} + g_i + h_i) x_i - g_{i} x_{i + 1} = y_i,
$$
where $g_i$, $h_i$ are positive and $g_i$ is proportional to the mixing coefficient
and can be much greater than $h_i$, which is the layer thickness.


### 3.1 Thomas algorithm
The Thomas algorithm is a simplified form of Gaussian elimination for tridiagonal systems of equations.
In its typical implementation, the forward elimination phase proceeds as follows
$$
\begin{aligned}
c'_i &= \frac{c_i}{b_i - a_i c'_{i - 1}}, \\
y'_i &= \frac{y_i - a_i y'_{i - 1}}{b_i - a_i c'_{i - 1}}.
\end{aligned}
$$
When applied to the diffusion system this leads to
$$
\begin{aligned}
c'_i &= \frac{-g_i}{h_i + g_{i - 1} (1 + c'_{i - 1}) + g_{i}}, \\
y'_i &= \frac{y_i + g_{i - 1} y'_{i - 1}}{h_i + g_{i - 1} (1 + c'_{i - 1}) + g_i}.
\end{aligned}
$$
Let's consider what happens when the mixing coefficient, and hence $g_i$, abruptly changes.
Suppose that $h_i$ is small, $g_{i - 1}$ is small, but $g_i$ is very large.
Then $c'_i \approx -1$ and in the next iteration the term $g_i (1 + c'_i)$ in
the denominator of both expression will multiply a very small number by a very large number.

Following [Appendix E in Schopf and Loughe](https://journals.ametsoc.org/view/journals/mwre/123/9/1520-0493_1995_123_2839_argiom_2_0_co_2.xml), to remedy that we can introduce
$$
\alpha'_i = g_i (1 + c'_i),
$$
which satisfies the following recursion relation
$$
\alpha'_i = \frac{g_i (h_i + \alpha'_{i - 1})}{h_i + \alpha'_{i - 1} + g_i}.
$$
The above equation together with
$$
\begin{aligned}
c'_i &= \frac{-g_i}{h_i + \alpha'_{i - 1} + g_{i}}, \\
y'_i &= \frac{y_i + g_{i - 1} y'_{i - 1}}{h_i + \alpha'_{i - 1} + g_i},
\end{aligned}
$$
forms the modifed stable algorithm.

### 3.2 (Parallel) cyclic reduction

The Thomas algorithm is work-efficient, but inherently serial.
While systems in different columns can
be solved in parallel, this might not expose enough parallelism on modern GPUs.
There are parallel tridiagonal algorithms that perform better on modern GPUs,
see [Zhang, Cohen, and Owens](https://doi.org/10.1145/1837853.1693472).
The two algorithms best suited for small systems are cyclic reduction and parallel cyclic reduction.

The basic idea of both cyclic reduction algorithms is as follows. Let's consider three consecutive equations corresponding to $y_{i- 1}$, $y_i$, and $y_{i+1}$
$$
\begin{aligned}
a_{i - 1} x_{i - 2} + b_{i - 1} x_{i - 1} + c_{i - 1} x_i &= y_{i - 1}, \\
a_i x_{i - 1} + b_i x_i + c_i x_{i + 1} &= y_i, \\
a_{i + 1} x_{i} + b_{i + 1} x_{i + 1} + c_{i + 1} x_{i + 2} &= y_{i + 1}.
\end{aligned}
$$
We can eliminate $x_{i - 1}$ and $x_{i + 1}$ to obtain a system of the form
$$
\hat{a}_i x_{i - 2} + \hat{b}_i x_i + \hat{c}_i x_{i + 2} = \hat{y}_i,
$$
where the modified coefficients $\hat{a}_i$, $\hat{b}_i$, $\hat{c}_i$ and the modified rhs $\hat{y}_i$ are
$$
\begin{aligned}
\hat{a}_i &= -\frac{a_{i - 1} a_i}{b_{i - 1}}, \\
\hat{b}_i &= b_i - \frac{c_{i - 1} a_i}{b_{i - 1}} - \frac{c_{i} a_{i + 1}}{b_{i + 1}}, \\
\hat{c}_i &= -\frac{c_i c_{i + 1}}{b_{i + 1}}, \\
\hat{y}_i &= y_i - \frac{a_i y_{i - 1}}{b_{i - 1}} - \frac{c_i y_{i + 1}}{b_{i + 1}}.
\end{aligned}
$$
The resulting system of equations for $x_{2j}$ is still tridiagonal and has roughly half the size of the original.

The cyclic reduction algorithm has two phases.
In the first phase, the above elimination step is iterated until the system is reduced to either one or two equations, which can then be directly solved.
In each iteration the computation of modified coefficients can be done in parallel.
The second phase involves finding the rest of the solution by using the final coefficients.
The second phase is also iterative, where at each iteration the number of know solution values increases by a factor of two. A drawback of this algorithm is that the amount of parallel computations available at each iteration is not constant.

The parallel cyclic reduction is based on the same idea, but has only one phase. In the first iteration it reduces the
original system to two systems of half the size. The second iteration reduces it to four systems of quarter the size
and so on. In the final iteration systems of size one or two are solved to obtain the whole solution at once.
In contrast to the cyclic reduction, this algorithm has constant amount of parallelism available at the cost of performing
more redundant work.

A naive application of the cyclic reduction to the diffusion system would result
in the following equation for the modified main diagonal
$$
\hat{b}_i = h_i
             + g_{i - 1} - \frac{g_{i-1}^2}{h_{i - 1} + g_{i - 2} + g_{i - 1}}
             + g_i - \frac{g_i^2}{h_{i + 1} + g_{i} + g_{i + 1}}.
$$
Using this expression can potentially result in catastrophic cancellation errors and overflows if
$g_{i - 1}$ or $g_{i + 1}$ are very large.
To improve its stability, this expression can be rewritten as
$$
\hat{b}_i = h_i
             + h_{i - 1} \frac{g_{i - 1}}{h_{i - 1} + g_{i - 2} + g_{i - 1}}
             + h_{i + 1} \frac{g_i}{h_{i + 1} + g_i + g_{i + 1}}
             + g_{i - 1} \frac{g_{i - 2}}{h_{i - 1} + g_{i - 2} + g_{i - 1}}
             + g_{i + 1} \frac{g_i}{h_{i + 1} + g_i + g_{i + 1}}.
$$
This equation can be shown to be in the form
$$
\hat{b}_i = \hat{h}_i + \hat{g}_{i - 2} + \hat{g}_i,
$$
where
$$
\begin{aligned}
\hat{h}_i &= h_i
             + h_{i - 1} \frac{g_{i - 1}}{h_{i - 1} + g_{i - 2} + g_{i - 1}}
             + h_{i + 1} \frac{g_i}{h_{i + 1} + g_i + g_{i + 1}}, \\
\hat{g}_i &= g_{i + 1} \frac{g_i}{h_{i + 1} + g_i + g_{i + 1}}.
\end{aligned}
$$
These two equations form the basis of stable (parallel) cyclic reduction for diffusion problems.


## 4 Design

Four different algorithms will be implemented:
- Thomas algorithm for general tridiagonal systems on CPUs
- PCR algorithm for general tridiagonal systems on GPUs
- Thomas algorithm for diffusion systems with improved stability properties on CPUs
- PCR algorithm for diffusion systems with improved stability properties on GPUs

The algorithms will be designed to work within Kokkos team policies, with each team of threads
solving one column system on GPUs and `VecLength` column systems on CPUs.
The user interface for CPU and GPU solvers will be the same. There will be
two type aliases `TriDiagSolver` and `TriDiagDiffSolver` that will resolve to the
optimal solver class based on the chosen architecture.

### 4.1 Data types and parameters

#### 4.1.1 Parameters

No parameters are required.

#### 4.1.2 Class/structs/data types

There will be a solver struct for each algorithm. The solvers inputs and outputs will be encapsulated in
two scratch data structs, one for the general solver and one for the diffusion solver.

#### 4.1.2.1 Scratch data structs

To facilitate constructing the systems on the fly and for performance reasons the solver inputs and outputs will be using the Kokkos scratch memory space. The scratch data for the general tridiagonal solver will be encapsulated in a struct called `TriDiagScratch`
```c++
struct TriDiagScratch {
   ScratchArray2DReal DL;
   ScratchArray2DReal D;
   ScratchArray2DReal DU;
   ScratchArray2DReal X;
};
```
where `DL`, `D`, `DU`, `X` are views of size (`NRow`, `VecLength`) in the scratch memory space. The views
`DL`, `D`, `DU` are inputs denoting the lower, main, and upper diagonal, respectively. The view `X` should contain
the rhs at input and will be overwritten with the solution after `solve` is called.

The scratch data for the specialized diffusion tridiagonal solver will be encapsulated in a struct called `TriDiagDiffScratch`
```c++
struct TriDiagDiffScratch {
   ScratchArray2DReal G;
   ScratchArray2DReal H;
   ScratchArray2DReal X;
   ScratchArray2DReal Alpha;
};
```
where `G`, `H`, `X`, `Alpha` are views of size (`NRow`, `VecLength`) in the scratch memory space. The views
`G` and `H` are inputs corresponding to the variables `g` and `h` introduced in [Section 3](#3-algorithmic-formulation).
The view `X` has the same meaning as in the general case.
The view `Alpha` is an internal workspace used by the algorithm.

#### 4.1.2.2 Solver structs

The four solver algorithms will be implemented as four structs `ThomasSolver`, `PCRSolver`, `ThomasDiffusionSolver`, and `PCRDiffusionSolver`. Currently, there is no plan for those structs to have any data members and they will only provide static methods, acting essentially as namespaces.

### 4.2 Methods

#### 4.2.1 Scratch Constructors
The constructors of scratch spaces take a team member and system size (`NRow`)
```c++
    TriDiagScratch(const TeamMember &Member, int NRow);
    TriDiagDiffScratch(const TeamMember &Member, int NRow);
```

#### 4.2.2 Policy creation
Every solver will provide a static method
```c++
   static TeamPolicy makeTeamPolicy(int NBatch, int NRow);
```
that creates an appropriate team policy for solving `NBatch` systems of size `NRow`.

#### 4.2.3 Solve Methods
The general solvers `ThomasSolver` and `PCRSolver` will have a static solve method
```c++
   static void solve(const TeamMember &Member, const TriDiagScratch &Scratch);
```
that takes a team member and an initialized general scratch space. After calling this method `Scratch.X`
will contain the solution. There will also be a convenience method
```c++
static void solve(const TeamMember &Member,
                  const Array2DReal &DL, const Array2DReal &D, const Array2DReal &DU, const Array2DReal &X);
```
that loads the inputs from global arrays.

The diffusion solvers `ThomasDiffusionSolver` and `PCRDiffusionSolver` will provide a similar method
```c++
   static void solve(const TeamMember &Member, const TriDiagDiffScratch &Scratch);
```
differing only in the type of scratch space. Similarly, there will be a convenience method
```c++
   static void solve(const TeamMember &Member,
                     const Array2DReal &G, const Array2DReal &H, const Array2DReal &X);
```
which loads the inputs from global arrays.

## 5 Verification and Testing

### 5.1 Test solvers correctness using prescribed matrix

Given analytically prescribed matrix `A` and vector `y` the solution to the problem `A x = z` with `z = A y` will be
checked to see if the resulting `x` is equal to the prescribed vector `y`. This will be done for all solvers for a variety of (batch size, system size) combinations.

### 5.2 Test diffusion solvers convergence using manufactured solution

The convergence of diffusion solvers will be tested using a manufactured solution.

### 5.3 Test stability

The diffusion solvers stability will be tested on an idealized vertical mixing problem with abrupt changes in
the diffusion coefficient.
