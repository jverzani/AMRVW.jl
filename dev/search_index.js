var documenterSearchIndex = {"docs":
[{"location":"#AMRVW","page":"AMRVW","title":"AMRVW","text":"","category":"section"},{"location":"","page":"AMRVW","title":"AMRVW","text":"Implementation of core-chasing algorithms for finding eigenvalues of factored matrices.","category":"page"},{"location":"","page":"AMRVW","title":"AMRVW","text":"This repository provides a Julia package implementing the methods, as applied to the problem of finding the roots of polynomials through the computation of the eigenvalues of a sparse factorization of the companion matrix in:","category":"page"},{"location":"","page":"AMRVW","title":"AMRVW","text":"Fast and Backward Stable Computation of Roots of Polynomials. Jared L. Aurentz, Thomas Mach, Raf Vandebril, and David S. Watkins SIAM J. Matrix Anal. Appl., 36(3), 942–973. (2015) https://doi.org/10.1137/140983434","category":"page"},{"location":"","page":"AMRVW","title":"AMRVW","text":"Fast and backward stable comptation of roots of polynomials, Part II: backward error analysis; companion matrix and companion pencil. By Jared L. Aurentz, Thomas Mach, utation of roots of polynomials, Part II: backward error analysis; companion matrix and companion pencil. By Jared L. Aurentz, Thomas Mach, Leonardo Robol, Raf Vandebril, David S. Watkins; arXiv:1611.02435","category":"page"},{"location":"","page":"AMRVW","title":"AMRVW","text":"The methods are summarized in monograph format:","category":"page"},{"location":"","page":"AMRVW","title":"AMRVW","text":"Core-Chasing Algorithms for the Eigenvalue Problem; by Jared L. Aurentz, Thomas Mach, Leonardo Robol, Raf Vandebril, and David S. Watkins; https://doi.org/10.1137/1.9781611975345","category":"page"},{"location":"","page":"AMRVW","title":"AMRVW","text":"As well, the twisted algorithm from \"A generalization of the multishift QR algorithm\" by Raf Vandebril and David S. Watkins; https://doi.org/10.1137/11085219X is implemented here.","category":"page"},{"location":"","page":"AMRVW","title":"AMRVW","text":"The core-chasing algorithms utilize Francis's QR algorithm on sparse factorizations of the respected companion matrix. For polynomials with real coefficients, the storage requirements are O(n) and the algorithm requires O(n) flops per iteration, or O(n^2) flops overall. The basic QR algorithm applied to the full companion matrix would require O(n^2) storage and O(n^3) flops overall.","category":"page"},{"location":"#Reference","page":"AMRVW","title":"Reference","text":"","category":"section"},{"location":"","page":"AMRVW","title":"AMRVW","text":"Modules = [AMRVW]","category":"page"},{"location":"#AMRVW.AbstractRFactorization","page":"AMRVW","title":"AMRVW.AbstractRFactorization","text":"An R factorization encapsupulate the R in a QR factorization. For the companion matrix case, this is a sparse factorization in terms of rotators. Other scenarios are different sub-types.\n\n\n\n\n\n","category":"type"},{"location":"#AMRVW.RFactorizationRankOne","page":"AMRVW","title":"AMRVW.RFactorizationRankOne","text":"A companion matrix will have a QR decomposition wherein R is essentially an identy plus a rank one matrix\n\n\n\n\n\n","category":"type"},{"location":"#AMRVW.RFactorizationUnitaryDiagonal","page":"AMRVW","title":"AMRVW.RFactorizationUnitaryDiagonal","text":"For the case where the QR decomposion has R as a diagonal matrix that is unitary\n\n\n\n\n\n","category":"type"},{"location":"#AMRVW.RFactorizationUpperTriangular","page":"AMRVW","title":"AMRVW.RFactorizationUpperTriangular","text":"For the case where the QR decomposition has R as a full, upper-triangular matrix\n\n\n\n\n\n","category":"type"},{"location":"#AMRVW.amrvw-Union{Tuple{S}, Tuple{Vector{S}, Vector{S}}} where S","page":"AMRVW","title":"AMRVW.amrvw","text":"amrvw(vs, ws)\n\nComputes a sparse factorization of of the companion matrix of a polynomial specified througha  pencil decomposition.\n\nA pencil decomposition of a polynomial, is a specification where if p = a0 + a1x^1 + ... + xn x^n then vs[1] = a0, vs[i+1] + ws[i] = ai, and ws[n] = an.\n\n\n\n\n\n","category":"method"},{"location":"#AMRVW.amrvw-Union{Tuple{Vector{S}}, Tuple{S}} where S","page":"AMRVW","title":"AMRVW.amrvw","text":"amrvw(ps)\n\nFor a polynomial specified by ps computes a sparse factorization of its companion matrix.\n\n\n\n\n\n","category":"method"},{"location":"#AMRVW.create_bulge-Union{Tuple{VV}, Tuple{S}, Tuple{T}, Tuple{AMRVW.QFactorization{T, S, VV}, Any, Any, Any}} where {T, S<:Real, VV}","page":"AMRVW","title":"AMRVW.create_bulge","text":"create_bulge(state)\n\nFinds m=1 or 2 shifts (m=1 in the CSS case, m=2 in the RDS case) based on the eigenvalues of the lower 2x2 block (using stop_index). The vector x = alpha (A - rho_1 I) e_1 or x = alpha (A-rho_1 I) (A-rho_2 I) e1 is found. From this, one or two core transforms are found so that U_1' x = gamma e_1 or U_1' U_2' x = gamma e_1. The values U_1 or U_1, U_2 are store in state.\n\n\n\n\n\n","category":"method"},{"location":"#AMRVW.givensrot-Union{Tuple{T}, Tuple{T, T}} where T<:Real","page":"AMRVW","title":"AMRVW.givensrot","text":"c,s,r  = givensrot(a,b)\n\nwhere `[c s; -s conj(c)] * [a,b] = [r, 0]\n\n\n\n\n\n","category":"method"},{"location":"#AMRVW.qr_factorization-Tuple{AbstractMatrix}","page":"AMRVW","title":"AMRVW.qr_factorization","text":"qr_factorization(H::Matrix; unitary=false)\n\nFor a Hessenberg matrix H return a factorization,Q⋅R, where Q is a QFactorization object of descending rotators and R is either a AMRVW.RFactorizationUnitaryDiagonal object (if unitary=true) or a RFactorizationUpperTriangular object.\n\nH may be a full matrix or the .H compontent  of a hessenberg factorization.\n\n\n\n\n\n","category":"method"},{"location":"#AMRVW.roots-Union{Tuple{Vector{S}}, Tuple{S}} where S","page":"AMRVW","title":"AMRVW.roots","text":"roots(ps)\nroots(vs, ws)\n\nUse an algorithm of AMRVW to find roots of a polynomial over the reals of complex numbers.\n\nps: The coefficients, [p0, p1, p2, ...., pn], of the polynomial p0 + p1*x + p2*x^2 + ... + pn*x^n\n\n[vs, ws]: If a second set of coefficients is passed, a pencil decomposition is used. A pencil decomposition satisfies vs[1]=p0; vs[i+1]+ws[i] = pi and ws[n] = pn. A pencil can be used when there is a wide range in the coefficients of the polynomial, as though slower, it can be more stable. The non-exported function basic_pencil implements the pencil which pulls off the leading coefficient a_n.\n\nReturns a complex vector of roots. When the algorithm fails, a warning is issued about the number of non-identified roots.\n\nExamples:\n\nusing AMRVW; const A = AMRVW\nrs = rand(10)\nA.roots(rs)      # uses Real double shift algorithm\n\nrs = rand(Complex{Float64}, 10)\nA.roots(rs)      # uses Complex Single Shift algorithm\n\nrs = rand(10)\nvs, ws = A.basic_pencil(rs)\nA.roots(vs, ws)  # uses QZ pencil factorization\n\nReferences:\n\nFast and backward stable computation of the eigenvalues of matrix polynomials. By Aurentz, Jared & Mach, Thomas & Robol, Leonardo & Vandebril, Raf & Watkins, David. (2016). Mathematics of Computation. 88. DOI: 10.1090/mcom/3338.\n\nFast and backward stable computation of roots of polynomials, Part II: backward error analysis; companion matrix and companion pencil\n\nJared L. Aurentz, Thomas Mach, Leonardo Robol, Raf Vandebril, David S. Watkins; arXiv:1611.02435\n\n\n\n\n\n","category":"method"}]
}
