# AMRVW

[![Build Status](https://travis-ci.com/jverzani/AMRVW.jl.svg?branch=master)](https://travis-ci.com/jverzani/AMRVW.jl)


Implementation of core-chasing algorithms for finding eigenvalues of factored matrices.  Fortan code for such methods is provided in the [eiscor](https://github.com/eiscor/eiscor) repository.

This repository provides a `Julia` package implementing the methods,
as applied to the problem of finding the roots of polynomials through
the computation of the eigenvalues of a sparse factorization of the
companion matrix in:

* Fast and Backward Stable Computation of Roots of Polynomials.
Jared L. Aurentz, Thomas Mach, Raf Vandebril, and David S. Watkins
SIAM J. Matrix Anal. Appl., 36(3), 942–973. (2015)
https://doi.org/10.1137/140983434


* Fast and backward stable comptation of roots of polynomials, Part II: backward error analysis; companion matrix and companion pencil. By
Jared L. Aurentz, Thomas Mach, utation of roots of polynomials, Part II: backward error analysis; companion matrix and companion pencil. By
Jared L. Aurentz, Thomas Mach, Leonardo Robol, Raf Vandebril, David S. Watkins; arXiv:1611.02435

The methods are summarized in monograph format:

Core-Chasing Algorithms for the Eigenvalue Problem; by Jared L. Aurentz, Thomas Mach, Leonardo Robol, Raf Vandebril, and David S. Watkins; https://doi.org/10.1137/1.9781611975345

As well, the twisted algorithm from "A generalization of the multishift QR algorithm" by Raf Vandebril and David S. Watkins; https://doi.org/10.1137/11085219X is implemented here.

The core-chasing algorithms utilize Francis's QR algorithm on sparse factorizations of the respected companion matrix. For polynomials with real coefficients, the storage requirements are O(n) and the algorithm requires O(n) flops per iteration, or O(n^2) flops overall. The basic QR algorithm applied to the full companion matrix would require O(n^2) storage and O(n^3) flops overall.


## Examples

```julia
julia> import AMRVW; const A = AMRVW
julia> p4 = [24.0, -50.0, 35.0, -10.0, 1.0]  # (x-1) * (x-2) * (x-3) * (x-4)

5-element Array{Float64,1}:
  24.0
 -50.0
  35.0
 -10.0
   1.0

julia> A.roots(p4)
4-element Array{Complex{Float64},1}:
 0.9999999999999996 + 0.0im
 2.0000000000000027 + 0.0im
 2.9999999999999876 + 0.0im
 4.000000000000012 + 0.0im
```

By means of comparison, using the `Polynomials` package:

```julia
julia> using Polynomials
julia> roots(Polynomial(p4))

4-element Array{Float64,1}:
 1.000000000000002
 1.9999999999999805
 3.0000000000000386
 3.9999999999999822
```


The advantage of the methods comes from the fact that this algorithm
can be used for much larger polynomials.

* Compared to the `roots` function of the
  [Polynomials](https://github.com/JuliaMath/Polynomials.jl) package,
  the methods are faster once the degree is around 75, and much faster
  as this grows. These methods are O(n) in storage and O(n^2) in time,
  whereas `roots` is O(n^2) in storage (the full companion matrix is
  created) and O(n^3) in time (the running time of a generic
  eigenvalue solver). As well, the `roots` function only computes over
  `Float64` values, not generic floating point values.

* Compared to the `roots` function of the
  [PolynomialsRoots](https://github.com/giordano/PolynomialRoots.jl)
  package, these methods are a bit slower, and perhaps a bit less
  accurate. This `roots` function is O(n) in storage and also appears
  to be O(n^2) in time. This `roots` function works over generic
  floating point values. However, this `roots` method will run into
  numeric issues for polynomials of degree n larger than 350 or so.


The backward stability of the methods is shown in the second paper to
grow linearly in the norm of the coefficients, so the following should
be quite accurate and is computable in a reasonable time:


```julia
## by DOI:	10.1142/S0219199715500522, this should have expected value ~ 2/pi*log(n) + .625738072 + 2/(pi*n) ~ 6.48
julia> rs = rand(Float64, 10_000) .- 1/2
julia> @time rts  = A.roots(rs)
julia> sum(isreal, rts)
 15.955615 seconds (35 allocations: 1017.297 KiB)
5
```

As this is relatively speedy, statistics can be generated, albeit the following will take some time to  finish:

```julia
julia> xs = [sum(isreal, A.roots(randn(3000))) for _ in 1:3000]
julia> using StatsBase
julia> xbar, s = mean_and_std(xs)
julia> n = 3000
julia> xbar .+ 1.96*s/sqrt(n) * [-1,1], 2/pi*log(n) + .625738072 + 2/(pi*n)
 ([5.67865426156726, 5.820012405099407], 5.7229621769994745)
```


----

There are no exported functions, as of now. But the internal functions may be of interest. For example, the paper [Fast and stable unitary QR algorithm](http://etna.mcs.kent.edu/volumes/2011-2020/vol44/abstract.php?vol=44&pages=327-341) discusses a situation where a matrix `A` is unitary Hessenberg, and so is factored in terms of a descending chain of rotatorrs. To fit this matrix into the current framework, we have, for example:

```
using LinearAlgebra
T = Float64
const A =  AMRVW
Qs = A.random_rotator.(T, 1:10)
Q = A.DescendingChain(Qs)
QF = A.q_factorization(Q)
F = A.QRFactorization(QF)
eigvals(F)
```

Which can be compared with:

```
MI = diagm(0 => ones(T, 11))
Qs * MI |> eigvals
```
