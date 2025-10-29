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
The `parallelReduce` wrapper supports performing multiple reduction at the same time.
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
Similarly to `parallelFor`, `parallelReduce` supports labels and up to five dimensions.
