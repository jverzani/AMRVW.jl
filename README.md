# AMRVW


Implementation of core-chasing algorithms for finding eigenvalues of factored matrices.



These are used to find the roots of polynomials through a sparse factorization of the companion matrix in:

* Fast and backward stable computation of the eigenvalues of matrix polynomials. By Aurentz, Jared & Mach, Thomas & Robol, Leonardo & Vandebril, Raf & Watkins, David. (2016). Mathematics of Computation. 88. DOI: 10.1090/mcom/3338.


* Fast and backward stable computation of roots of polynomials, Part II: backward error analysis; companion matrix and companion pencil
Jared L. Aurentz, Thomas Mach, Leonardo Robol, Raf Vandebril, David S. Watkins; arXiv:1611.02435

These papers are summarized in monograph format:

Core-Chasing Algorithms for the Eigenvalue Problem; by Jared L. Aurentz, Thomas Mach, Leonardo Robol, Raf Vandebril, and David S. Watkins; https://doi.org/10.1137/1.9781611975345

## Examples

```julia
julia> import AMRVW: roots

julia> p4 = [24.0, -50.0, 35.0, -10.0, 1.0]  # (x-1) * (x-2) * (x-3) * (x-4)
5-element Array{Float64,1}:
  24.0
 -50.0
  35.0
 -10.0
   1.0

julia> roots(p4)
4-element Array{Complex{Float64},1}:
 2.0000000000000226 + 0.0im
 0.9999999999999978 + 0.0im
  4.000000000000043 + 0.0im
 2.9999999999999423 + 0.0im
```

By means of comparison

```julia
julia> roots(Poly(p4))
4-element Array{Float64,1}:
 1.000000000000002
 1.9999999999999805
 3.0000000000000386
 3.9999999999999822
```


The advantage comes in the fact that this algorithm can be used for much larger polynomials. For degree 4 polynomials, the `roots` method is competitive, but for larger degree polynomials (~ 75) the sparse method is faster. The second paper shows that the backward error grows linearly in the norm of the coefficients, so the following should be quite accurate and computable

```julia
## by DOI:	10.1142/S0219199715500522, this should have expected value ~ 2/pi* log(10_000) = 5.86...
julia> rs = rand(Float64, 10_000) .- 1/2
julia> @time rts  = A.roots(rs)
julia> rts  .|> isreal |> sum
 18.205041 seconds (31 allocations: 1.146 MiB)
5
```




----

The core-chasing algorithms utilize Francis's QR algorithm on sparse factorizations of the respected companion matrix. For polynomials with real coefficents, the storage requirements are O(n) and the algorithm requires O(n) flops per iteration, or O(n^2) flops overall. The basic QR algorithm applied to the full companion matrix would require O(n^2) storage and O(n^3) flops overall.
