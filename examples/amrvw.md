# An overview of `AMRVW`

The `AMRVW` package implements some numerical linear algebra algorithms of Jared L. Aurentz, Thomas Mach, Leonardo Robol, Raf Vandebril, and David S. Watkins for finding eigenvalues of matrices through Francis's method.

An inspiration is to find  the roots of a polynomial through the eigenvalues of a companion matrix. This is implemented in `AMRVW` through the `roots` function.

To illustrate,

```
using AMRVW
const A = AMRVW  # currently there are no exports
ps = [24.0, -50.0, 35.0, -10.0, 1.0]  #  (x-1)(x-2)(x-3)(x-4) = 24 -50x + 25x^2  -10x^3 +  x^4
A.roots(ps)
```

The example shows

* the interface expects polynomials specified through a vector of coefficients, `[a0, a1, ..., an]`
* the  4 roots, always as complex numbers
* the fact that the roots are numeric approximations due to accumulated round off errors.

Similarly, roots of polynomials over the  complex numbers can be found

```
ps  =  [0.0 + 1.0im, -1.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im, 0.0 - 1.0im, 1.0 + 0.0im] # (x-1)(x+1)(x-i)^2)(x+i)
A.roots(ps)
```

There are other ways  to numerically find roots  of polynomials  in `Julia`, notably the `roots` function  of  the `Polynomials` package and the `roots` function of the `PolynomialRoots` package:

* unlike `Polynomials.roots` (but similar to `PolynomialRoots.roots`) this  `roots` function can work with  big floats and other floating point types.

* For moderate sized polynomials ($n \approx 50$), `PolynomialRoots.roots` is faster than `Polynomials.roots` which is faster than `roots`, though all are fast. When $n \approx 75$, `roots` is faster (much so for large $n$ than `Polynomials.roots`).

* Unlike both `Polynomials.roots` and `PolynomialRoots.roots` this `roots` function can accurately identify roots of polynomials of high degree.  (For a polynomial with $n$ random coefficients, e.g., `ps = rand(n)`) `Polynomials.roots` will have troubles for n around 50; `PolynomialRoots.roots` will have issues for `n` around 300; `roots` can quickly handle degree 3000, and still be accurate for higher degrees.)

* The `roots` function is shown by  the authors  to be  backward stable. The same isn't the case  for the other  two.

* In the first example, the residual errors are  similar in size to `Polynomials.roots`, but  `PolynomialRoots.roots` the residual errors seem to be generally a bit smaller.

## The companion matrix

Both `roots` and `Polynomials.roots` use a companion matrix representation using the  eigenvalues of this matrix to identify the roots of the polynomial. (The  `PolynomialRoots.roots` function relies on a different method following a paper by [Skowron and Gould](https://arxiv.org/abs/1203.1034).

Using some functions within `AMRVW` we can see the companion matrix:

```
ps = [24.0, -50.0, 35.0, -10.0, 1.0]  #  (x-1)(x-2)(x-3)(x-4) = 24 -50x + 25x^2  -10x^3 +  x^4
F = A.amrvw(ps)
M = Matrix(F) |> round2   # round2 is just M -> round.(M, digits=2)
```


```
using LinearAlgebra
eigvals(F)
```

The eigvenvalues of this matrix indeed are the roots of the polynomials. Internally, `amrvw` actually uses an enlarged matrix, with an extra dimension that is not shown here.

## Francis's Algorithm

Francis's Algorithm begins with a QR decomposition `M`. For example,

```
LinearAlgebra.qr(M)
```

the decomposition in `AMRVW` is slightly different, though similar

```
Matrix(F.QF)
```

```
Matrix(F.RF) |> round2
```

(Here the `RF` matrix shows the extra size used internally in the algorithm.)

The idea of Francis's shifted algorithm is to identify shifts $\rho_1$, $\rho_2$, $\dots$, $\rho_m$ and generate a *unitary* matrix $V_0 = \alpha (A-\rho_1 I)(A-\rho_2 I)\cdots(A-\rho_m)I \cdot e_1$, $e_1$ being a unit vector with $1$ in the $1$ entry and $0$ elsewhere. As $V_0$ is unitary, the product $V_0^{-1}F V_0$ will have the same eigenvalues.  When $F$ is upper Hessenberg (upper triangular starting with the subdiagonal), as is the case with the companion matrix, then this product will be almost upper Hessenberg, save for a bulge.

In the real-coefficient case,  $m=2$ is used to allow the calculations to be done over the real numbers. For the complex-coefficient case, $m=1$ is possible.

```
ps  =  [0.0 + 1.0im, -1.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im, 0.0 - 1.0im, 1.0 + 0.0im]
F = A.amrvw(ps)
M = Matrix(F); n = size(M, 1)
MI = diagm(0 => ones(Complex{Float64}, 5)) # identity matrix
storage, ctr, m = A.make_storage(F), A.make_counter(F), 1
A.create_bulge(F.QF, F.RF, storage, ctr) # finds shifts and creates V_0
(V0 = storage.VU[1] * MI) |> round2
```

Up to rounding, $V_0$ is unitary:

```
isapprox(V0 * V0', MI, atol=1e-8)
```


The matrix $V_0' M V_0$ has a bulge below the subdiagonal (the `[3,1]` position, illustrated with a `1` below):

```
V0' * M * V0 |> round2 .|> !iszero
```

The algorithm finds $V_1$ to chase the bulge downward:

```
A.absorb_Ut(F.QF, F.RF, storage, ctr)
A.passthrough_triu(F.QF, F.RF, storage, ctr, Val(:right))
A.passthrough_Q(F.QF, F.RF, storage, ctr, Val(:right))
(V1 = storage.VU[1] * MI) |> round2
```

And this produce will have a bulge in `[4,2]` position:

```
V1' * (V0' * M * V0) * V1 |> round2 .|> !iszero
```

And again:

```
A.passthrough_triu(F.QF, F.RF, storage, ctr, Val(:right))
A.passthrough_Q(F.QF, F.RF, storage, ctr, Val(:right))
(V2 = storage.VU[1] * MI) |> round2
```

```
V2' * (V1' * (V0' * M * V0) * V1) * V2 |> round2 .|> !iszero
```

Once pushed to the bottom, the bulge is absorbed into the matrix, leaving an upper Hessenberg form.

If the shifts are appropriately chosen, after a few iterations this resulting matrix can have an eigenvalue immediately read off and after deflation subsequent eigenvalues can be found.

### Shifts

Above the bulge is created with a single rotator. As mentioned, for
the real variable case, two rotators are used, so that the
computations can be kept using real numbers. In general, the AMRVW
algorithm can be defined for $m$ rotators. These rotators are produced
by `create_bulge`, as illustrated above, and stored in `storage.UV`.

The choice of shifts is *essentially* the eigenvalues of the lower $2
\times 2$ submatrix (In the $m$-shift case, the lower $m\times m$ submatrix). In the Vandrebril and Watkins paper it is
mentioned that this choice *normally* yields quaratic
convergence. That is, one of the rotators will become a diagonal
rotator with `s` part $0$. When that happens, deflation can occur. The
algorithm is applied on this smaller, deflated matrix, until the
deflated matrix is comprised of no more than $m$ rotators, at least one of which is a
diagonal rotator. At this point one or more eigenvalues can be found. This
quadratic convergence implies that generally there are
$\mathcal{O}(n)$ steps taken.

## The AMRVW decomposition of the companion matrix

The main result of the two papers on "Fast and Backward Stable Computation of Roots of Polynomials" is a proof of backward stability of a method that utilizes a sparse factorization of both the `Q` and `R` parts of the QR decomposition of a companion matrix.

Returning to the real case, and digging into some structures, we can illustrate:

```
ps = [24.0, -50.0, 35.0, -10.0, 1.0]
F = A.amrvw(ps)
F.QF.Q
```

To explain, this is a "chain" of real rotators, more clearly seen with:

```
Vector(F.QF.Q)
```

A rotator is a matrix which is identical to the identity matrix except in the `[i,i+1] × [i, i+1]` block, in which case it takes the form of a rotator: `[c s; -s c]`. (Our rotators are in the different direction than those in the papers.) Here `c` and `s` are the cosine and sine of some angle. These rotators are indexed by `i` and we use the notation $U_i$ to indicate a rotator of this form for a given $i$. In the above, we  can see  with inspection that there are 3 rotators with $i$ being 1, 2, and 3. This set of rotators is "descending" due to their order (1 then 2 then 3); ascending would be 3 then 2 then 1. The product of descending rotators will be upper Hessenberg:

```
Matrix(F.QF.Q)
```

A rotator at level $i$ will commute with a rotator at level $j$ unless $|i-j| \leq 1$. In the case where $i-j = \pm 1$, a key computation is the "turnover", which represents $U_i V_j W_i$ as $VV_j WW_i UU_j$. With the turnover, we can easily pass a rotator through an ascending or descending chain without disturbing those patterns.

In the  above illustration of Francis's  algorithm, the matrices  $V_0$, $V_1$,  etc. can be seen to be  rotators of this type. More generally, a unitary matrix with $m$ shifts can be viewed as a product of $m$  such rotators.

The $R$ decomposition  is trickier.  In the initial QR decomposition, $R$ has a simple structure plus a rank one part (coming from the coefficients). The  decomposition has two chains, an ascending one, `Ct`, and a descending one, `B`:


```
Ct = F.RF.Ct
B = F.RF.B
```

These almost begin as inverses:

```
MI = diagm(0 => ones(5))
(Z = Ct * (B * MI)) |> round2
```

However, `Ct` is cleverly chosen to encode the rank 1 part. This can be uncovered through the following:

```
e1 = vcat(1, zeros(4))
en1 = vcat(zeros(4), 1)
rho = (en1' * (Ct * MI) * e1)
```

and

```
yt = -(1/rho * en1' * (Ct*(B*MI)))
Ct * (e1 * yt) |> round2
```

Leading to `R`:

```
Z + Ct * (e1 * yt) |> round2
```

The algorithm passes a rotator through this decomposition, which in turn relies on passing a rotator through the two chains `B` and `Ct`, which, with the turnover computation, is easily computed.


With these decompositions in mind, the computation above $V_0' M V_0$, can be seen as $U_1' Q R U_1$.  The product $U_1' Q = U_1' Q_1 Q_2 \cdots Q_k$ is just a product of two rotators at level 1, so $U_1' Q_1$ can be fused to give a new $\tilde{Q}_1$ in the descending chain factorization of $Q$.  The passthrough just mentioned allows $\tilde{Q} R U_1$ to have this form $\tilde{Q} \tilde{U}_1 \tilde{R}$ and by passing through the descending chain, we have this form $U_2 \hat{Q} \tilde{R}$. The matrix $V_1$ (of Francis's algorithm above) is seen to be $U_2$, as the similarity transform using $Q_1=U_2$ leaves the product $\hat{Q} \tilde{R} U_2$ having the same eigen values, but with the bulge shifted down one level. This basic idea forms the algorithm to chase the bulge. In the $m > 1$ case, some other details are included.

The decomposition of the companion matrix is sparse. Rather than require $O(n^2)$ storage, it only needs $O(n)$. The iterative algorithm is $O(n)$ per iteration  and  $O(n^2)$ overall, as compared to the $O(n^3)$ required in general. This reduction allows the method to be practical for large $n$, unlike `Polynomial.roots` which uses an $O(n^3)$ algorithm.


### Pencil decompositions

In "*Fast and backward stable computation of roots of polynomials, Part II*" the method is extended to the pencil decomposition of a polynomial.  A pencil decomposition of a polynomial, is a specification where if $p = a_0 + a_1x^1 + \cdots + a_n x^n$ then $v_1 = a_0$, $v_{i+1} + w_i = a_i$, and $w_n = a_n$. This has some advantages in cases where the polynomial has a particularly small leading coefficient, since division by a tiny $a_n$ will result in very large entries. The algorithm uses two upper triangular matrices.

The `roots` function allows a pencil decomposition to be passed in as two vectors:

```
ps = [24.0, -50.0, 35.0, -10.0, 1.0]
vs, ws = A.basic_pencil(ps)
A.roots(vs, ws)
```


#### The Wilkinson polynomial

The Wilkinson polynomial, $(x-1)(x-2)\cdots(x-20)$, poses a challenge for `roots`:

```
import Polynomials
ps = Polynomials.coeffs(Polynomials.poly(1.0:20.0))
A.roots(ps)
```

The answer involves complex-valued roots, even though the roots are
clearly integer valued. (The `Polynomials.roots` function and
`PolynomialRoots.roots` do get  real answers for  this polynomial). As an aside, the
exact implementation of the fundamental `turnover` operation will
effect the number of such roots.  The issue of spurious complex-valued
roots might be addressed by separating out the smaller coefficients
from the large ones. Here is a function to split a polynomial into two
pieces.

```
function pencil_split(ps, n)
    vs = zeros(eltype(ps), length(ps)-1)
    ws = zeros(eltype(ps), length(ps)-1)

    vs[1:n] = ps[1:n]
    ws[n:end] = ps[n+1:end]

    vs, ws
end
```

It turns out that with `n=13` only real roots are identified:

```
A.roots(pencil_split(ps, 13)...)
```

The pencil here does a better job with this polynomial, but the choice of `13` was made with hindsight, not foresight.


## Other uses

The `AMRVW.jl` package allows some other applications, though the exact interface needs to be settled on.

If we take rotators $Q_1, Q_2, \dots, Q_k$ their product will be upper Hessenberg:

```
Qs = A.random_rotator.(Float64, [1,2,3,4])
```

```
MI = diagm(0 => ones(5))
(M = Qs * MI) |> round2
```

Their eigenvalues can be found:

```
eigvals(M)
```

But the sparse representation can be used to also find such eigenvalues:

```
QF = A.q_factorization(A.DescendingChain(Qs))
F = A.QRFactorization(QF)  # defaults to identify R factorization
Matrix(F) |> round2 # same as M
```



```
eigvals(F)
```


The `qr_factorization` method can take a Hessenberg matrix and complete the factorization. For this case, we have:

```
F = A.qr_factorization(M, unitary=true)
eigvals(F)
```

When `unitary=true` is the case, this will outperform `eigvals` once the matrix size is around 50 by 50 and is significantly more performant for the 150 by 150 case.



----

Not all upper Hessenberg matrices can he expressed as a descending chain of rotators, as the latter is unitary. However, any upper Hessenberg matrix can easily be seen to be represented as a descending chain of rotators times an upper triangular matrix.

The Givens rotation is a rotator, $U$, chosen so that if $x= [a,b]$, then $Ux = [r,0]$. This allows, for example, the following:

```
M = triu(rand(5,5), -1)  # upper Hessenberg
```

```
c,s,r = A.givensrot(M[1,1], M[2,1])
U1 = A.Rotator(c,s,1)
U1 * M |> round2
```

That is the subdiagonal in column 1 is 0
Similarly, a `U2` could then be found so that subdiagonal in column 2 is 0, etc. That is a choice of rotators is available for $U_k U_{k-1} U_{k-2} \cdots U_2 U_1 M = R$. Setting $V_i = U_i'$, we have then $M = V_1 V_2 \cdots V_k R = QR$, where $Q$ is unitary and  in decomposed form, and $R$ is upper triangular.


For example:

```
Us = Any[]
G = copy(M)
for i in 1:4
  c,s,r = A.givensrot(G[i,i], G[i+1,i])
  Ui =  A.Rotator(c,s,i)
  pushfirst!(Us, Ui)
  lmul!(Ui, G)
end
R = G
Qs = reverse(adjoint.(Us))
```

With this, we can do the following:

```
QF = A.q_factorization(A.DescendingChain(Qs))
RF = A.RFactorizationUpperTriangular(R)

F = A.QRFactorization(QF, RF)

Matrix(F)  - M |> round2# same as F
```

And

```
[eigvals(F) eigvals(M)]
```



This patttern is encoded in the `qr_factorization` function, mentioned above. For any Hessenberg matrix it can be employed:


```
H = hessenberg(rand(5,5))
F = A.qr_factorization(H.H)
eigvals(F)
```

This approach is in the ballpark with `eigvals(Matrix(H.H))`, though
slower by a factor performance-wise, as the `R` part is not factored
into rotators, so uses $\mathcal{O}(n)$ operations -- not
$\mathcal{O}(1)$ -- in each step of the algorithm. (When
`unitary=true` is the case, the `R` part is much faster, so is
competitive.


----

The product of a descending chain of rotators is upper Hessenberg and the product of an  ascending chain or rotators is lower Hessenberg. With modifications, described in "A generalization of the multishift QR algorithm" a related algorithm can be employed. A "twisted chain" is one where a descending chain of rotators is permuted, an ascending chain being one possible twisting among many others.

```
N = 4
T = Float64; S = T # real case here
MI = diagm(0 => ones(N+1))
Qs = A.random_rotator.(T, [4,3,2,1])
M = Qs * MI

QF = A.q_factorization(A.TwistedChain(Qs))
F = A.QRFactorization(QF)   # use default identify R factorization

[eigvals(M) eigvals(F)]
```

Here is another example with a CMV pattern to the twisted

```
N = 5
MI = diagm(0 => ones(N+1))
Qs = A.random_rotator.(T, [1,3,5,2,4])
M = Qs * MI

QF = A.q_factorization(A.TwistedChain(Qs))

F = A.QRFactorization(QF)

[eigvals(M) eigvals(F)]
```

The implementation for twisted chains is not nearly as efficient as that for descending chains.
