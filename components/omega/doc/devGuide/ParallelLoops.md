(omega-dev-parallel-loops)=

# Parallel loops

Omega adopts the Kokkos programming model to express on-node parallelism. To provide
simplified syntax for the most frequently used computational patterns, Omega provides
wrappers funtions that internally handle creating and setting-up Kokkos policies.

## Flat multi-dimensional parallelism

### parallelFor

To perform parallel iteration over a multi-dimensional index range Omega provides the
`parallelFor` wrapper. For example, the following code shows how to set every element of
a 3D array in parallel.
```c++
   Array3DReal A("A", N1, N2, N3);
   parallelFor(
       {N1, N2, N3},
       KOKKOS_LAMBDA(int J1, int J2, int J3) {
          A(J1, J2, J3) = J1 + J2 + J3;
       });
```
Ranges with up to five dimensions are supported.
Optionally, a label can be provided as the first argument of `parallelFor`.
```c++
   parallelFor("Set A",
       {N1, N2, N3},
       KOKKOS_LAMBDA(int J1, int J2, int J3) {
          A(J1, J2, J3) = J1 + J2 + J3;
       });
```
Adding labels can result in more informative messages when
Kokkos debug variables are defined.

### parallelReduce

To perform parallel reductions over a multi-dimensional index range the
`parallelReduce` wrapper is available. The following code sums
every element of `A`.
```c++
   Real SumA;
   parallelReduce(
       {N1, N2, N3},
       KOKKOS_LAMBDA(int J1, int J2, int J3, Real &Accum) {
          Accum += A(J1, J2, J3);
       },
       SumA);
```
Note the presence of an accumulator variable `Accum` in the `KOKKOS_LAMBDA` arguments.
You can use `parallelReduce` to perform other types of reductions.
As an example, the following snippet finds the maximum of `A`.
```c++
   Real MaxA;
   parallelReduce(
       {N1, N2, N3},
       KOKKOS_LAMBDA(int J1, int J2, int J3, Real &Accum) {
          Accum = Kokkos::max(Accum, A(J1, J2, J3));
       },
       Kokkos::Max<Real>(MaxA));
```
To perform reductions that are not sums, in addition to modifying the lambda body,
the final reduction variable needs to be cast to the appropriate type. In the above example,
`MaxA` is cast to `Kokkos::Max<Real>` to perform a max reduction.
The `parallelReduce` wrapper supports performing multiple reductions at the same time.
You can compute `SumA` and `MaxA` in one pass over the data:
```c++
   parallelReduce(
       {N1, N2, N3},
       KOKKOS_LAMBDA(int J1, int J2, int J3, Real &AccumSum, Real &AccumMax) {
          AccumSum += A(J1, J2, J3);
          AccumMax = Kokkos::max(AccumMax, A(J1, J2, J3));
       },
       SumA, Kokkos::Max<Real>(MaxA));
```
There is no limit to how many reductions can be done at the same time.
It is usually beneficial for performance to group a small number (2-8) of simple
reductions together.
Similarly to `parallelFor`, `parallelReduce` supports labels and up to five dimensions.

## Hierarchical parallelism

A more flexible alternative to flat parallelism is hierarchical parallelism.
In general, Kokkos supports a three-level parallelism hierarchy of
- teams of a league
- threads of a team
- vector lanes of a thread

In Omega, we simplify this design to a two-level hierarchy (outer and inner) that naturally
corresponds to the  model parallelism in the horizontal and vertical dimensions.
The outer level of parallelism in Omega corresponds to parallelism over Kokkos teams,
whereas the inner level may be a combination of parallelism over team threads and thread vector lanes.
The simplest example of a hierarchical parallel loop using Omega wrappers looks as follows
```c++
   parallelForOuter(
       {NCells},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
           parallelForInner(Team, NVertLayers, INNER_LAMBDA (int KLayer) {
               // Inner body
           });
   });
```
Note the appearance of a variable named `Team` of type `TeamMember` in the outer lambda
argument list. This variable handles launching of inner loops and team synchronization.
The flexibility of hierarchical parallelism comes from:
- the possibility of the inner loop range to depend on the outer index
- the option to use different patterns at the outer and inner levels
- the option to have multiple parallel loops at the inner level

The following outer iteration patterns are supported in Omega:
- `parallelForOuter`
- `parallelReduceOuter`

The following inner iteration patterns are supported in Omega:
- `parallelForInner`
- `parallelReduceInner`
- `parallelScanInner`
- `parallelSearchInner`

To provide even more flexibility, the outer loops support iterating over a multi-dimensional range.
Currently, the inner loops are limited to one dimension.

### INNER_LAMBDA

When using hierarchical parallelism the outer loop lambda has to be annotated with `KOKKOS_LAMBDA` to
be device accessible. In the inner lambda it must be possible to obtain correct results using capture by value (`[=]`).
However, in some cases, better performance might be possible by using
capture by reference (`[&]`). For that reason, Omega provides the
`INNER_LAMBDA` macro that expends to the best performing capture clause depending on the chosen
architecture.

### teamBarrier

When multiple parallel inner loops are present within one outer loop, synchronization within
a thread team might be necessary. Omega provides the `teamBarrier` function that
synchronizes the threads in a team, which can be called as follows.
```c++
teamBarrier(Team);
```
This function may only be called at the outer level and needs to be called by every thread in a team.
Failure to follow these rules might result in deadlocks.

### Kokkos::single

When using hierarchical parallelism, it might be necessary to restrict execution to just one thread of a team.
To do that Kokkos provides the `single` function. To execute a statement once per team with a single thread do
```c++
  Kokkos::single(PerTeam(Team), [=](){
      // executed once per team (by one of the team threads)
  });
```

### parallelForOuter
To start outer iterations over a multidimensional index range the `parallelForOuter` wrapper is available.
A call to `parallelForOuter` might look as follows.
```c++
   parallelForOuter(
       {N1},
       KOKKOS_LAMBDA(int J1, const TeamMember &Team) {
           // launch inner loops here
       });
```
Labels and multi-dimensional iteration ranges of up to five dimensions are supported by `parallelForOuter`.


### parallelReduceOuter
To start outer reductions over a multi-dimensional index range the `parallelReduceOuter` wrapper is available.
A call to `parallelReduceOuter` might look as follows.
```c++
   Real OuterSum;
   parallelReduceOuter(
       {N1},
       KOKKOS_LAMBDA(int J1, const TeamMember &Team, Real &Accum) {
           Real InnerContribution;
           // compute InnerContribution
           Kokkos::single(PerTeam(Team), [&](){
               Accum += InnerContribution;
           });
       }, OuterSum);
```
Note how in this example the addition of `InnerContribution` is done inside `Kokkos::single`
with just one thread from a team. This shows how to add one contribution per team to the total sum.
The reduction accumulator comes after the team member variable in the lambda argument list.
Labels and multi-dimensional iteration ranges of up to five dimensions are supported in `parallelReduceOuter`.
Different types of reductions (min, max) and multiple reducers are also supported.

### parallelForInner
To launch inner parallel iterations Omega provides the
`parallelForInner` wrapper. For example, the following code shows how to set every element of
a 3D array in parallel using hierarchical parallelism.
```c++
   Array3DReal A("A", N1, N2, N3);
   parallelForOuter(
       {N1, N2}, KOKKOS_LAMBDA(int J1, int J2, const TeamMember &Team) {
        parallelForInner(Team, N3, INNER_LAMBDA(Int J3) {
          A(J1, J2, J3) = J1 + J2 + J3;
        });
    });
```
Labels are not supported by `parallelForInner` and only one-dimensional index range can be used.
The inner loop range can depend on the outer loop index. For example, to set the elements below the main
diagonal of a square matrix one can do:
```c++
   Array2DReal M("M", N, N);
   parallelForOuter(
       {N}, KOKKOS_LAMBDA(int J1, const TeamMember &Team) {
        parallelForInner(Team, J1, INNER_LAMBDA(Int J2) {
          M(J1, J2) = J1 + J2;
        });
    });
```

### parallelReduceInner
To launch inner parallel reductions Omega provides the `parallelReduceInner` wrapper.
For example, computing sums along the third dimension of a 3D array in parallel and storing them
in a 2D array might be done as follows.
```c++
   Array3DReal A("A", N1, N2, N3);
   Array2DReal B("B", N1, N2);
   parallelForOuter(
       {N1, N2}, KOKKOS_LAMBDA(int J1, int J2, const TeamMember &Team) {
        Real SumD3;
        parallelReduceInner(Team, N3, INNER_LAMBDA(Int J3, Real &Accum) {
            Accum += A(J1, J2, J3);
        }, SumD3);
        B(J1, J2) = SumD3;
    });
```
Labels are not supported by `parallelReduceInner` and only one-dimensional index range can be used.
Different types of reductions (min, max) and multiple reducers are supported.
For example, to additionally compute and store maxima along the third dimension of A do
```c++
   Array2DReal C("C", N1, N2);
   parallelForOuter(
       {N1, N2}, KOKKOS_LAMBDA(int J1, int J2, const TeamMember &Team) {
        Real SumD3, MaxD3;
        parallelReduceInner(Team, N3, INNER_LAMBDA(Int J3, Real &AccumSum, Real &AccumMax) {
            AccumSum += A(J1, J2, J3);
            AccumMax = Kokkos::Max(AccumMax, A(J1, J2, J3));
        }, SumN3, MaxN3);
        B(J1, J2) = SumD3;
        C(J1, J2) = MaxD3;
    });
```
Inner reductions are collective operations within a team. This means that the reduction results
are available to every thread of a team, without any need for synchronization

### parallelScanInner
To launch inner parallel scans Omega provides the `parallelScanInner` wrapper.
For example, computing prefix sums (partial sums) along the third dimension of a 3D array might
be done as follows.
```c++
   Array3DReal A("A", N1, N2, N3);
   Array3DReal D("D", N1, N2, N3);
   parallelForOuter(
       {N1, N2}, KOKKOS_LAMBDA(int J1, int J2, const TeamMember &Team) {
       parallelScanInner(Team, N1, INNER_LAMBDA(Int J3, Real &Accum, bool IsFinal) {
            Accum += A(J1, J2, J3);
            if (IsFinal) {
              D(J1, J2, J3) = Accum;
            }
        });
       });
```
This example computes partial sums up to and including `A(J1, J2, J3)` because the accumulator is updated
before the `if` statement. That is, it performs an inclusive scan. To compute an exclusive scan
simply move the addition after the `if` statement.
```c++
  Real FinalScanValue;
  parallelScanInner(Team, N1, INNER_LAMBDA(Int J3, Real &Accum, bool IsFinal) {
       if (IsFinal) {
         D(J1, J2, J3) = Accum;
       }
       Accum += A(J1, J2, J3);
  }, FinalScanValue);
```
Moreover, this example illustrates that the final scan value can be obtained by providing
an additional argument `FinalScanValue`. Labels are not supported by `parallelScanInner`
and only one-dimensional index range can be used. In contrast to `parallelReduceInner`,
`parallelScanInner` supports only sum-based scans and only one scan variable.

### parallelSearchInner
To search an index range in parallel for the first index where a given condition occurs Omega
provides the `parallelSearchInner` function.
For example, the following code finds, for each row of a matrix, the first column index where
the matrix element is above a certain threshold. If no element matches the condition then
`parallelSearchInner` returns `-1`.
```c++
   Array2DReal M("M", N1, N2);
   Array1DI3 ThresholdIdx("ThresholdIdx", N1);
   parallelForOuter(
       {N1}, KOKKOS_LAMBDA(int J1, const TeamMember &Team) {

       int Idx;
       parallelSearchInner(Team, N2, INNER_LAMBDA(Int J2) {
            return M(J1, J2) > Threshold;
       }, Idx);

       ThresholdIdx(J1) = Idx;
   });
```
Labels are not supported by `parallelSearchInner` and only one-dimensional index range can be used.
