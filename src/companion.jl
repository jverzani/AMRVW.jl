##  Factor a  companion matrix
##  so  that  eigenvalue are the roots

## We have ps = [a0, a1, ..., an]
## we factor -(QR) so we use p(-x) coefficients here
function basic_decompose(ps::Vector{T}) where {T}

    #    ps, k = deflate_leading_zeros(ps)  ## ASSUMED
    n = length(ps) - 1
    p0, pNi = ps[1], 1/ps[end]
    qs = Vector{T}(undef, n+1)
    par = iseven(n) ? one(T) : -one(T)

    @inbounds for i in 1:n-1
        qs[i] =  par * (-1)^(i+1) * ps[i+1] * pNi
    end


    qs[n] =  par * p0 * pNi
    qs[n+1] = -one(T)

    qs
end


## A pencil function takes
## ps = [a0, a1, a2, ..., an] and splits into vs, ws with:
## v1 = a0,
## wn = an, and
## v_{i+1} + w_i = a_i for i-1,2,...,n-1.
##
## This gives matrices:
## V = [0 ....     -v1   W = [1  .... w1
##      1 0 .....  -v2        0 1 ... w2
##      0 1 .....  -v3        0 0 1 . w3
##          .....               .....
##      0 0 .... 1 -vn]       0 ....0 wn]
## V is Hessenberg, W upper triangular.
##
## The default, `basic_pencil`, is
## vs = [a0, a1, ..., a_{n-1}], ws = [0,0,..., an]
##
function basic_pencil(ps::Vector{T}) where {T}

    N = length(ps)-1

    qs = zeros(T, N)
    qs[N]  = ps[end]

    ps = ps[1:end-1]
    ps, qs

end



##################################################
##
## Factor R: (Z +  yt)
## xs are decompose(ps)

## Constructor for Real coefficients; no pencil
function r_factorization(xs::Vector{S}) where {S}
    N = length(xs) - 1
    T = S <: Real ? S : real(S)
    Ct = AscendingChain(Vector{Rotator{T,S}}(undef, N))
    B =  DescendingChain(Vector{Rotator{T,S}}(undef, N))
    D = SparseDiagonal(S, N+1)

    _r_factorization(xs, Ct, B, D)

    RFactorizationRankOne(Ct, B, D)
end



# constructor populating Ct, B, D for Real Case
function _r_factorization(xs::Vector{T}, Ct, B, D) where {T <: Real}
    N = length(xs) - 1
    c, s, tmp = givensrot(xs[N], -one(T))
    C = Rotator(c, s, N)
    Ct[1] = C' #Rotator(c, -s, N)
    B[N] = Rotator(s, -c, N)  # not the adjoint, so that Ct_1 * B_n = [0 -1; 1 0]

    @inbounds for i in (N-1):-1:1
        c,s,tmp = givensrot(xs[i], tmp)
        C = Rotator(c, s, i)
        Ct[N + 1 - i], B[i] = C', C
    end

    nothing
end


# constructor populating Ct, B, D for Complex Case (phase is impt)
function _r_factorization(xs::Vector{S}, Ct, B, D) where {S <: Complex}
    N = length(xs) - 1

    c, s, tmp = givensrot(xs[N], -one(S))

    nrm = norm(c)
    alpha = c/nrm

    C = Rotator(c, s, N)
    Ct[1] = C'
    B[N] = Rotator(s * alpha, -nrm, N)
    D[N] = conj(alpha)  # Ct*B*D gives [0 -1; 1 0] block at tail
    D[N+1] = alpha


    @inbounds for i in (N-1):-1:1
        c,s,tmp = givensrot(xs[i], tmp)
        C = Rotator(c, s, i)
        Ct[N + 1 - i], B[i] = C', C
    end


    nothing
end

## Factor Q

function q_factorization(xs::Vector{S}) where {S}
    N = length(xs) - 1
    T = S <: Real ? S : real(S)

    Q =  DescendingChain(Vector{Rotator{T,S}}(undef, N-1))
    zt,ot,zs,os = _zero_one(xs)

    @inbounds for ii = 1:(N-1)
        Q[ii] = Rotator(zero(S), one(T), ii)
    end

    D = SparseDiagonal(S, N)

    return QFactorization(Q, D)

end

##################################################
##
##

"""

    amrvw(ps)

For a polynomial specified by `ps` computes a sparse factorization of its companion matrix.

"""
function amrvw(ps::Vector{S}) where {S}

    xs = basic_decompose(ps)

    QF = q_factorization(xs)
    RF = r_factorization(xs)

    QRFactorization(QF, RF)

end


## decomposition of ps into vs, ws needs to be expanded to work with framework
function adjust_pencil(vs::Vector{T}, ws::Vector{T}) where {T}
    ps = basic_decompose(vcat(vs, one(T)))
    qs = vcat(ws, -one(T))  # or +1?

    ps, qs
end

"""

    amrvw(vs, ws)

Computes a sparse factorization of of the companion matrix of a polynomial specified througha  pencil decomposition.

A pencil decomposition of a polynomial, is a specification where
if `p = a0 + a1x^1 + ... + xn x^n`
then
`vs[1] = a0`,
`vs[i+1] + ws[i] = ai`, and
`ws[n] = an`.

"""
function amrvw(vs::Vector{S}, ws::Vector{S}) where {S}

    N = length(vs)

    ps, qs = adjust_pencil(vs, ws)
    QF = q_factorization(Vector{S}(undef, N+1))
    ZF = pencil_factorization(ps, qs)

    QRFactorization(QF, ZF)

end



"""

    roots(ps)
    roots(vs, ws)

Use an algorithm of AMRVW to find roots of a polynomial over the reals of complex numbers.

ps: The coefficients, `[p0, p1, p2, ...., pn]`, of the polynomial `p0 + p1*x + p2*x^2 + ... + pn*x^n`

[vs, ws]: If a second set of coefficients is passed, a pencil decomposition is used.
A pencil decomposition satisfies `vs[1]=p0; vs[i+1]+ws[i] = pi and ws[n] = pn`.
A pencil can be used when there is a wide range in the coefficients of the polynomial, as though slower, it can be more stable.
The non-exported function `basic_pencil` implements the pencil which pulls off the leading coefficient `a_n`.

Returns a complex vector of roots. When the algorithm fails, a warning is issued about the number of non-identified roots.

Examples:

```
using AMRVW
rs = rand(10)
roots(rs)      # uses Real double shift algorithm

rs = rand(Complex{Float64}, 10)
roots(rs)      # uses Complex Single Shift algorithm

rs = rand(10)
vs, ws = A.basic_pencil(rs)
roots(vs, ws)  # uses QZ pencil factorization

```

References:

* Fast and backward stable computation of the eigenvalues of matrix polynomials. By Aurentz, Jared & Mach, Thomas & Robol, Leonardo & Vandebril, Raf & Watkins, David. (2016). Mathematics of Computation. 88. DOI: 10.1090/mcom/3338.


* Fast and backward stable computation of roots of polynomials, Part II: backward error analysis; companion matrix and companion pencil
Jared L. Aurentz, Thomas Mach, Leonardo Robol, Raf Vandebril, David S. Watkins; arXiv:1611.02435


"""
function roots(ps::Vector{S}) where {S}
    # deflate 0s
    ps, K = deflate_leading_zeros(ps)

    if length(ps) <= 3
        rts = solve_simple_cases(ps)
        append!(rts, zeros(Complex{real(S)}, K))
        return rts
    end

    state = amrvw(ps)

    AMRVW_algorithm(state)

    if K > 0 # put back in leading zeros
        ZS = zeros(real(S), K)
        append!(state.REIGS, ZS)
        append!(state.IEIGS, ZS)
    end
#@show "done"
    rts = complex.(state.REIGS, state.IEIGS)
    try
        _sorteig!(rts)
    catch
    end
    rts
end

## A Pencil algorithm
## if p = a0 + a1x^1 + ... + xn x^n
## then
## vs[1] = a0
## vs[i+1] + ws[i] = ai
## ws[n] = an
## is a pencil decompostion.
## `basic_pencil` will break [a0, a1, ..., an] into vs=[a0,a1, ..., a_{n-1}], ws=[0,0,..., 0, an]
function roots(vs::Vector{S}, ws::Vector{S}) where {S}

    ps, qs, K = deflate_leading_zeros(vs, ws)

    ## solve simple cases
    N = length(ps)
    if N <= 2
        rs = vcat(ps, zero(S)) + vcat(zero(S), qs)
        rts = roots(rs)
        append!(rts, zeros(Complex{real(S)}, K))
        return rts
    end

    state = amrvw(ps, qs)
    AMRVW_algorithm(state)

    if K > 0 # put back in leading zeros
        ZS = zeros(real(S), K)
        append!(state.REIGS, ZS)
        append!(state.IEIGS, ZS)
    end

    complex.(state.REIGS, state.IEIGS)
end
