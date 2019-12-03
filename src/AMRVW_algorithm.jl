
## Main algorithm of AMRVW
## This follows that given in the paper very closely
## should trace algorithm better with `verbose`
function AMRVW_algorithm(state::QRFactorization)

    it_max = 20 * length(state)
    kk = 0

    while kk <= it_max

        kk += 1
        state.ctrs.stop_index <= 0 && return     ## finished up!
        state.ctrs.it_count += 1

        ## show_status(state)


        check_deflation(state)

        k = state.ctrs.stop_index
        delta =  state.ctrs.stop_index - state.ctrs.zero_index

        if delta >= 2

            bulge_step(state)

        elseif delta == 1

            diagonal_block(state,  k + 1)

            e1r,e1i, e2r,e2i = eigen_values(state)
            state.REIGS[k], state.IEIGS[k] = e1r, e1i
            state.REIGS[k+1], state.IEIGS[k+1] = e2r, e2i

            # can finish up if near end
            if state.ctrs.stop_index == 2
                diagonal_block(state, 2)
                e1 = state.A[1,1]
                state.REIGS[1] = real(e1)
                state.IEIGS[1] = imag(e1)
            end

            state.ctrs.zero_index = 0
            state.ctrs.start_index = 1
            state.ctrs.stop_index -= 2

        elseif delta == 0

            diagonal_block(state, state.ctrs.stop_index + 1)
            e1, e2 = state.A[1,1], state.A[2,2]
            if state.ctrs.stop_index == 1
                state.REIGS[1], state.IEIGS[1] = real(e1), imag(e1)
                state.REIGS[2], state.IEIGS[2] = real(e2), imag(e2)
                state.ctrs.stop_index = 0

            else

                state.REIGS[k+1], state.IEIGS[k+1] = real(e2), imag(e2)
                state.ctrs.zero_index = 0
                state.ctrs.start_index = 1
                state.ctrs.stop_index = k - 1

            end

        end

    end


    @warn "Not all eigenvalues were found. The first $(state.ctrs.stop_index-1) are missing."
end

## Find eigenvalues from factorization
function LinearAlgebra.eigvals(state::QRFactorization, args...)
    AMRVW_algorithm(state)
    complex.(state.REIGS, state.IEIGS)
end


##################################################
## Deflation
## when a Q[k] matrix becomes a diagonal matrix, we deflate.
## This is checked by the sine term being basically 0.
function check_deflation(state::AbstractFactorizationState{T,S,Rt, Pt, Twt}, tol=eps(T)) where {T, S,Rt, Pt, Twt}
    QF = state.QF
    for k in state.ctrs.stop_index:-1:state.ctrs.start_index
        c, s = vals(QF.Q[k])
        if abs(s) <= tol
            deflate(QF, k)
            state.ctrs.zero_index = k      # points to a matrix Q[k] either
            state.ctrs.start_index = k + 1
            state.ctrs.it_count = 1        # reset counter
            return
        end
    end
end

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



"""

    amrvw(ps)

For a polynomial specified by `ps` computes a sparse factorization of its companion matrix.

"""
function amrvw(ps::Vector{S}) where {S}

    N = length(ps) - 1
    xs = basic_decompose(ps)

    QF = q_factorization(xs)
    RF = r_factorization(xs)

    qrfactorization(N, QF, RF)

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
A pencil decomposition of a polynomial, is a specificaiton where
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
    ZF = z_factorization(ps, qs)

    qrfactorization(N, QF, ZF)

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
    complex.(state.REIGS, state.IEIGS)
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
