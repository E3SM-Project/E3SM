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

### Inner Iteration Ranges

There are two ways of specifying the iteration range of an inner loop.
The first takes the total number of iterations `N` as the second argument
```c++
   parallelForInner(Team, N, INNER_LAMBDA (int K) {
   });
```
and the loop index `K` takes values from `0` up to and including `N - 1`.
The second way uses a helper struct `Range` to provide a range of valid indices
```c++
   parallelForInner(Team, Range{N1, N2}, INNER_LAMBDA (int K) {
   });
```
Note that this range is inclusive, i.e. the loop index `K` takes values from `N1` up to and including `N2`.
This means that `Range{0, N}` specifies a different range than the first example.
For simplicity, most examples in this document use the first way of specifying the range,
but a `Range` argument can be passed to all inner iteration patterns.


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
        parallelForInner(Team, N3, INNER_LAMBDA(int J3) {
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
        parallelForInner(Team, J1, INNER_LAMBDA(int J2) {
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
        parallelReduceInner(Team, N3, INNER_LAMBDA(int J3, Real &Accum) {
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
        parallelReduceInner(Team, N3, INNER_LAMBDA(int J3, Real &AccumSum, Real &AccumMax) {
            AccumSum += A(J1, J2, J3);
            AccumMax = Kokkos::max(AccumMax, A(J1, J2, J3));
        }, SumD3, Kokkos::Max<Real>(MaxD3));
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
       parallelScanInner(Team, N3, INNER_LAMBDA(int J3, Real &Accum, bool IsFinal) {
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
  parallelScanInner(Team, N3, INNER_LAMBDA(int J3, Real &Accum, bool IsFinal) {
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
To search an index range in parallel for the first index at which a given condition occurs,
Omega provides the `parallelSearchInner` function.
For example, the following code finds, for each row of a matrix, the first column index where
the matrix element is above a certain threshold. If no element matches the condition then
`parallelSearchInner` returns `-1`.
```c++
   Array2DReal M("M", N1, N2);
   Array1DI4 ThresholdIdx("ThresholdIdx", N1);
   const Real Threshold = 0.5;
   parallelForOuter(
       {N1}, KOKKOS_LAMBDA(int J1, const TeamMember &Team) {

       int Idx;
       parallelSearchInner(Team, N2, INNER_LAMBDA(int J2) {
            return M(J1, J2) > Threshold;
       }, Idx);

       ThresholdIdx(J1) = Idx;
   });
```
Labels are not supported by `parallelSearchInner` and only one-dimensional index range can be used.

### Launch Config

While specifying loop bounds is enough to start an outer parallel loop, sometimes more control over the underlying
Kokkos `TeamPolicy` is desired. The most common use case is utilizing scratch memory, a concept discussed more
thoroughly in the next sub-section. To enable more control, outer loops can be launched by providing
a `LaunchConfig` struct as the first argument, which is composed of three parts:
- loop bounds,
- team size,
- amount of scratch memory.

For example, the following snippet launches a loop iterating over a two-dimensional index range
with team size of 32 and enough scratch memory for 8 `Real` values and 4 `I4` values per team.
```c++
   auto LConfig = LaunchConfig({N1, N2}, 32, TeamScratch<Real, I4>(8, 4));
   parallelForOuter(LConfig,
       KOKKOS_LAMBDA(int J1, int J2, const TeamMember &Team) {
   });
```
It is not necessary to provide all three arguments to `LaunchConfig`. If you want the default team size,
or you don't need any scratch memory, you can use the following constructors.
```c++
   auto LConfig1 = LaunchConfig({N1, N2}, TeamScratch<Real, I4>(8, 4));
   auto LConfig2 = LaunchConfig({N1, N2}, 32);
```
For simplicity, most examples in this document use the simple form of launching outer loops with just the bounds,
but `LaunchConfig` can be used for all types of outer parallel loops.
Inner parallel loops cannot use `LaunchConfig`.

### Team Scratch Memory

In hierarchical code, it is often useful to have some amount of scratch memory private to each team.
Scratch memory enables reuse of expensive to compute data in inner loops.
To enable scratch memory, the outer loops needs to be launched with the `LaunchConfig` parameter described above,
configured with the requested number of scratch values.
Inside the outer loop, unmanaged scratch arrays can be created from a pool of memory accessible
by calling the `teamScratch(Team)` function.
Scratch arrays have a different type than normal Omega arrays, for example `ArrayScratch1DReal` is the
type of a 1D scratch array of Reals. They also cannot have labels.

As an example, the following code uses scratch memory to compute an expensive function on elements of a 2D array `A`.
It then computes finite differences along the second dimension of the scratch array, and stores them in `A`.
By using scratch memory, the expensive function is only computed once for every element, and there is no need for global memory allocation.
```c++
   Array2DReal A("A", N1, N2);
   parallelForOuter(
       LaunchConfig({N1}, TeamScratch<Real>(N2)),
       KOKKOS_LAMBDA(int J1, const TeamMember &Team) {

        ArrayScratch1DReal SA(teamScratch(Team), N2);

        parallelForInner(Team, N2, INNER_LAMBDA (int J2) {
            SA(J2) = expensiveFunc(A(J1, J2));
        });

        teamBarrier(Team);

        parallelForInner(Team, N2, INNER_LAMBDA (int J2) {

            const int J2M1 = Kokkos::max(J2 - 1, 0);
            const int J2P1 = Kokkos::min(J2 + 1, N2 - 1);

            A(J1, J2) = SA(J2P1) - SA(J2M1);
        });
   });
```
You can create multiple scratch arrays of different types, as in the following code.
```c++
   parallelForOuter(
       LaunchConfig({N1}, TeamScratch<Real, I4>(4, 8)),
       KOKKOS_LAMBDA(int J1, const TeamMember &Team) {
        ArrayScratch1DI4 ScratchI4(teamScratch(Team), 8);
        ArrayScratch1DReal ScratchReal(teamScratch(Team), 4);
   });
```
As the above example illustrates, the order in which the arrays are created inside the outer region
doesn't need to match the order of arguments to `TeamScratch`.
