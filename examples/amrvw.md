# An overview of `AMRVW`

The `AMRVW` package implements some numerical linear algebra algorithms of Jared L. Aurentz, Thomas Mach, Leonardo Robol, Raf Vandebril, David S. Watkins for finding eigenvalues of matrices through Francis's method.

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

Both `roots` and `Polynomials.roots` use a companion matrix representation and use the  eigenvalues of this matrix to identify the roots of the polynomial. (The  `PolynomialRoots.roots` function relies on a different method following a paper by [Skowron and Gould](https://arxiv.org/abs/1203.1034).

Using some functions within `AMRVW` we can see the companion matrix:

```
ps = [24.0, -50.0, 35.0, -10.0, 1.0]  #  (x-1)(x-2)(x-3)(x-4) = 24 -50x + 25x^2  -10x^3 +  x^4
state = A.amrvw(ps)
F = Matrix(state) |> round2   # round2 is just M -> round.(M, digits=2)
```

This isn't quite the classic decomposition, where the coefficients are in the last column, but we see this has the proper eigenvalues:

```
using LinearAlgebra
eigvals(F)
```

Well, *almost*. This companion matrix has an extra row and column added, introducing an eigenvalue of $0$.

## Francis's Algorithm

Francis's Algorithm begins with a QR decomposition `F`. For example,

```
LinearAlgebra.qr(F)
```

the decomposition in `AMRVW` is slightly different, though similar

```
Matrix(state.QF)
```

```
Matrix(state.RF) |> round2
```

The idea of Francis's shifted algorithm is to identify shifts $\rho_1$, $\rho_2$, $\dots$, $\rho_m$ and generate a *unitary* matrix $V_0 = \alpha (A-\rho_1 I)(A-\rho_2 I)\cdots(A-\rho_m)I \cdot e_1$, $e_1$ being a unit vector with $1$ in the $1$ entry and $0$ elsewhere. As $V_0$ is unitary, the product $V_0^{-1}F V_0$ will have the same eigenvalues.  When $F$ is upper Hessenberg (upper triangular starting with the subdiagonal), as is the case with the companion matrix, then this product will be almost upper Hessenberg, save for a bulge.

In the real-coefficient case,  $m=2$ is used to allow the calculations to be done over the real numbers. For the complex-coefficient case, $m=1$ is possible.

```
ps  =  [0.0 + 1.0im, -1.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im, 0.0 - 1.0im, 1.0 + 0.0im]
state = A.amrvw(ps)
F = Matrix(state)
M = diagm(0 => ones(Complex{Float64}, 6)) # identity matrix
A.create_bulge(state)   # finds shifts and creates V_0
(V0 = state.UV[1] * M) |> round2
```

Up to rounding, $V_0$ is unitary:

```
isapprox(V0 * V0', M, atol=1e-8)
```


The matrix $V_0' F V_0$ has a bulge below the subdiagonal (the `[3,1]` position):

```
V0' * F * V0 |> round2
```

The algorithm finds $V_1$ to chase the bulge downward:

```
A.absorb_Ut(state)
A.passthrough_triu(state, Val(:right))
A.passthrough_Q(state, Val(:right))
(V1 = state.UV[1] * M) |> round2
```

And this produce will have a bulge in `[4,2]` position:

```
V1' * (V0' * F * V0) * V1 |> round2
```

And again:

```
A.passthrough_triu(state, Val(:right))
A.passthrough_Q(state, Val(:right))
(V2 = state.UV[1] * M) |> round2
```

```
V2' * (V1' * (V0' * F * V0) * V1) * V2 |> round2
```

Once pushed to the bottom, the bulge is absorbed into the matrix, leaving an upper Hessenberg form.

If the shifts are appropriately chosen, after a few iterations this resulting matrix can have an eigenvalue immediately read off and after deflation subsequent eigenvalues can be found.

### Shifts

XXX 1,2,m

## The AMRVW decomposition of the companion matrix

The main result of the two papers on "Fast and Backward Stable Computation of Roots of Polynomials" is a proof of backward stability of a method that utilizes a sparse factorization of both the `Q` and `R` parts of the QR decomposition of a companion matrix.

Returning to the real case, and digging into some structures, we can illustrate:

```
ps = [24.0, -50.0, 35.0, -10.0, 1.0]
state = A.amrvw(ps)
state.QF.Q
```

To explain, this is a "chain" of real rotators. A rotator is a matrix which is identical to the identity matrix except in the `[i,i+1] × [i, i+1]` block, in which case it takes the form of a rotator: `[c s; -s c]`. (Our rotators are in the different direction than those in the papers.) Here `c` and `s` are the cosine and sine of some angle. These rotators are indexed by `i` and we use the notation $U_i$ to indicate a rotator of this form for a given $i$. In the above, we  can see  with inspection that there are 3 rotators with $i$ being 1, 2, and 3. This set of rotators is "descending" due to their order (1 then 2 then 3); ascending would be 3 then 2 then 1. The product of descending rotators will be upper Hessenberg:

```
M = diagm(0 => ones(Float64, 4))
state.QF.Q * M
```

A rotator at level $i$ will commute with a rotator at level $j$ unless $|i-j| \leq 1$. In the case where $i-j = \pm 1$, a key computation is the "turnover", which represents $U_i V_j W_i$ as $VV_j WW_i UU_j$. With the turnover, we can easily pass a rotator through an ascending or descending chain without disturbing those patterns.

In the  above illustration of Francis's  algorithm, the matrices  $V_0$, $V_1$,  etc. can be seen to be  rotators of this type. More generally, a unitary matrix with $m$ shifts can be viewed as a product of $m$  such rotators.

The $R$ decomposition  is trickier.  In the initial QR decomposition, $R$ has a simple structure plus a rank one part (coming from the coefficients). The  decomposition has two chains, an ascending one and a descending one,


```
Ct = state.RF.Ct
B = state.RF.B
```

These almost begin as inverses:

```
M = diagm(0 => ones(5))
(Z = Ct * (B * M)) |> round2
```

However, `Ct` is cleverly chosen to encode the rank 1 part. This can be uncovered through the following:

```
e1 = vcat(1, zeros(4))
en1 = vcat(zeros(4), 1)
rho = (en1' * (Ct * M) * e1)
```

and

```
yt = -(1/rho * en1' * (Ct*(B*M)))
Ct * (e1 * yt) |> round2
```

Leading to:

```
Z + Ct * (e1 * yt) |> round2
```

The algorithm passes a rotator through this decomposition, which in turn relies on passing a rotator through the two chains `B` and `Ct`, which, with the turnover computation, is easily computed.


With these decompositions in mind, the computation above $V_0' F V_0$, can be seen as $U_1' Q R U_1$.  The product $U_1' Q = U_1' Q_1 Q_2 \cdots Q_k$ is just a product of two rotators at level 1, so $U_1' Q_1$ can be fused to give a new $\tilde{Q}_1$ in the descending chain factorization of $Q$.  The passthrough just mentioned allows $\tilde{Q} R U_1$ to have this form $\tilde{Q} \tilde{U}_1 \tilde{R}$ and by passing through the descending chain, we have this form $U_2 \hat{Q} \tilde{R}$. The matrix $V_1$ (of Francis's algorithm above) is seen to be $U_2$, as the similarity transform using $Q_1=U_2$ leaves the product $\hat{Q} \tilde{R} U_2$ having the same eigen values, but with the bulge shifted down one level. This basic idea forms the algorithm to chase the bulge. In the $m > 1$ case, some other details are included.

The decomposition of the companion matrix is sparse. Rather than require $O(n^2)$ storage, it only needs $O(n)$. The iterative algorithm is $O(n)$ per iteration  and  $O(n^2)$ overall, as compared to the $O(n^3)$ required in general. This reduction allows the method to be practical for large $n$, unlike `Polynomial.roots` which uses an $O(n^3)$ algorithm.


### Pencil decompositions

In "Fast and backward stable computation of roots of polynomials, Part II" the method is extended to the pencil decomposition of a polynomial.  A pencil decomposition of a polynomial, is a specification where if $p = a_0 + a_1x^1 + \cdots + a_n x^n$ then $v_1 = a_0$, $v_{i+1} + w_i = a_i$, and $w_n = a_n$. This has some advantages in cases where the polynomial has a particularly small leading coefficient, since division by a tiny $a_n$ will result in very large entries. The algorithm uses two upper triangular matrices.

The `roots` function allows a pencil decomposition to be passed in as two vectors:

```
ps = [24.0, -50.0, 35.0, -10.0, 1.0]
vs, ws = A.basic_pencil(ps)
A.roots(vs, ws)
```


## Other uses

The `AMRVW.jl` package allows some other applications, though the exact interface needs to be settled on.

If we take rotators $Q_1, Q_2, \dots, Q_k$ their product will be upper Hessenberg:

```
Qs = A.random_rotator.(Float64, [1,2,3,4])
```

```
M = diagm(0 => ones(5))
(F = Qs * M) |> round2
```

Their eigenvalues can be found:

```
eigvals(F)
```

But the sparse representation can be used to also find such eigenvalues:

```
T = Float64
N = 5
D = A.SparseDiagonal(T, N)
QF = A.QFactorization(A.DescendingChain(Qs), D)

RF = A.RFactorizationIdentity{T,T}()
state = A.QRFactorization(QF, RF)
Matrix(state) |> round2 # same as F
```



```
eigvals(state)
```


----

Not all upper Hessenberg matrices can he expressed as a descending chain of rotators, as the latter is unitary. However, any upper Hessenberg matrix can easily be seen to be represented as a descending chain of rotators times an upper triangular matrix.

The Givens rotation is a rotator, $U$, chosen so that if $x= [a,b]$, then $Ux = [r,0]$. This allows, for example, the following:

```
F = triu(rand(5,5), -1)  # upper Hessenberg
```

```
c,s,r = A.givensrot(F[1,1], F[2,1])
U1 = A.Rotator(c,s,1)
U1 * F |> round2
```

That is the subdiagonal in column 1 is 0
Similarly, a `U2` could then be found so that subdiagonal in column 2 is 0, etc. That is a choice of rotators is available for $U_k U_{k-1} U_{k-2} \cdots U_2 U_1 F = R$. Setting $V_i = U_i'$, we have then $F = V_1 V_2 \cdots V_k R = QR$, where $Q$ is unitary and  in decomposed form, and $R$ is upper triangular.


For example:

```
Us = Any[]
G = copy(F)
for i in 1:4
  c,s,r = A.givensrot(G[i,i], G[i+1,i])
  Ui =  A.Rotator(c,s,i)
  pushfirst!(Us, Ui)
  G .= Ui * G
end
R = G
Qs = reverse(adjoint.(Us))
```

With this, we can do the following:

```
N = 5
D = A.SparseDiagonal(T, N)
QF = A.QFactorization(A.DescendingChain(Qs), D)

RF = A.RFactorizationUpperTriangular(R)
state = A.QRFactorization(QF, RF)
Matrix(state)  - F |> round2# same as F
```

And

```
[eigvals(state) eigvals(F)]
```



With this patttern, we might provide an `eigvals` for `Hessenberg` matrices:

```
function LinearAlgebra.eigvals(H::Hessenberg)
R = Matrix(H.H) # need this for R
S = eltype(R)
T = real(S)
N = size(R)[1]
Qs = Vector{A.Rotator{T,S}}(undef, N-1)
for i in 1:N-1
  c,s,r = A.givensrot(R[i,i], R[i+1,i])
  Ui =  A.Rotator(c,s,i)
  Qs[i] = Ui'
  for k in i:N
    rik, rjk =  R[i,k],  R[i+1,k] # non-allocating multiplication Ui*R
    R[i,k]    = c * rik + s * rjk
    R[i+1,k]  = -conj(s) * rik + conj(c) * rjk
	end
end
D = A.SparseDiagonal(S, N)
QF = A.QFactorization(A.DescendingChain(Qs), D)
RF = A.RFactorizationUpperTriangular(R)
state = A.QRFactorization(QF, RF)
eigvals(state)
end
```

```
H = hessenberg(rand(5,5))
[eigvals(H) eigvals(Matrix(H.H))]
```

Sadly this is not competitive performance-wise with `eigvals(Matrix(H.H))`, as the `R` part is not factored. Though, this approach can be made to work with other floating point types, such as `BigFloat`, that `LinearAlgebra.eigvals` can not currently handle.


----

The product of a descending chain of rotators is upper Hessenberg and the product of an  ascending chain or rotators is lower Hessenberg. With modifications, described in "A generalization of the multishift QR algorithm" a related algorithm can be employed. A "twisted chain" is one where a descending chain of rotators is permuted, an ascending chain being one possible twisting among many others.

```
N = 4
T = Float64; S = T # real case here
M = diagm(0 => ones(N+1))
Qs = A.random_rotator.(T, [4,3,2,1])
F = Qs * M

D = A.SparseDiagonal(T,N+1)
QF = A.QFactorizationTwisted(A.TwistedChain(Qs), D)
RF = A.RFactorizationIdentity{T,S}()
state = A.QRFactorizationTwisted(QF, RF)

[eigvals(state) eigvals(F)]
```

Here is another example with a CMV pattern to the twisted

```
N = 5
M = diagm(0 => ones(N+1))
Qs = A.random_rotator.(T, [1,3,5,2,4])
F = Qs * M

D = A.SparseDiagonal(T,N)
QF = A.QFactorizationTwisted(A.TwistedChain(Qs), D)
RF = A.RFactorizationIdentity{T, S}()
state = A.QRFactorization(QF, RF)

[eigvals(state) eigvals(F)]
```

The implementation for twisted chains is not nealy as efficient as that  for descending chains.
