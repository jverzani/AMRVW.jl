##
## Implementation of the algorithm from
## A GENERALIZATION OF THE MULTISHIFT QR ALGORITHM by RAF VANDEBRIL AND DAVID S. WATKINS.
## https://doi.org/10.1137/11085219X
## pdf file found at http://www.sci.wsu.edu/math/faculty/watkins/pdfiles/vw12_SIMAX.pdf
##
## Q is neither descending or ascending, but rather twisted
## the unitary transformation is not limited to 1 or 2 rotators, but rather m
##
## TODO:
##
## * use a D matrix for real case with twisted
##
## * move storage, ctrs into algorithm, not factorization
##
## * Factorization should just include terms to reconstruct matrix
##
## *
##
## * using linearAlgebra.Diagonal
##
## * Clean up allocations
##
## DONE * store Decoupled and Ms in same space to avoid allocations; work with start_index, stop_index
##
## * get shifts by paper (m -> eigenvalues of mxm matrix, continue...)
##
## This algorithm encompasses those in CSS and RDS, but those are more efficient space wise
## This is more general and perhaps of interest as building blocks for experimentation

##################################################


##
## Twisted refers to the QFactorization
##
struct QFactorizationTwisted{T, S, Vt, PVt} <: AbstractQFactorization{T, S}
  Q::TwistedChain{T,S, Vt, PVt} ## Twisted
  D::SparseDiagonal{S}
end

## unlike faster case, here instead of checking for parity in the diagonal rotators, we move -1 terms into  D
function deflate(QF::QFactorizationTwisted{T, S, Vt, Pvt}, k) where {T,S <: Real, Vt, Pvt}

    c,s = vals(QF.Q[k])
    i = idx(QF.Q[k])
    QF.Q[k] = Rotator(one(T), zero(T), i) # ± 1, not just 1
    if sign(c) < 0
        ## move into D term so no dflip concerns
        Asc = ascending_part(QF.Q, i)
        passthrough_phase(DiagonalRotator(c, i), Asc, QF.D)
        Des = descending_part(QF.Q, i)
        passthrough_phase(DiagonalRotator(c, i), Des, QF.D)
    end

end

## # XXX  This is  in CSS
## function deflateXXX(QF::QFactorizationTwisted{T, S, Vt, Pvt}, k) where {T, S <: Complex, Vt, Pvt}

##     alpha, s = vals(QF.Q[k])
##     i = idx(QF.Q[k])

##     ## Make Q[k] an identity rotator
##     QF.Q[k] = Rotator(one(Complex{T}), zero(T), i)

## #    @show k, QF.Q[k]

##     # absorb Di into D
##     Di =  DiagonalRotator(alpha, i)
##     ## passthrough Ascending and Descdending parts of QF.Q merge with D
##     passthrough_phase!(Di, QF)

## end

function passthrough_phase!(Di::DiagonalRotator, QF::QFactorizationTwisted)
    i = idx(Di)

    #ainds = iget(QF.Q, i, Val(:Asc))
    #dinds = iget(QF.Q, i, Val(:Des))
    asc = ascending_part(QF.Q, i)
    des = descending_part(QF.Q, i)

    passthrough_phase!(Di, (des, asc), QF.D)
#    passthrough_phase!(Di,(DescendingChain(dMMs), AscendingChain(aMMs)), QF.D)
#    for (i,j) in enumerate(ainds)
#        QF.Q[j] = aMMs[i]
#    end
#    for (i,j) in enumerate(dinds)
#        QF.Q[j] = dMMs[i]
#    end
end

Base.eltype(QF::QFactorizationTwisted) = eltype(QF.Q.x)
## return q[i:k, j:k]
function Base.getindex(QF::QFactorizationTwisted, j, k)
    ## This is wrong as QF.Q[j,k] is not correct
    QF.Q[j,k] * (k == 0 ? 0 : QF.D[k])
end

function Base.Matrix(QF::QFactorizationTwisted{T, S}) where {T, S}

    n = length(QF) + 1
    M = diagm(0 => ones(S, n))
    D = Matrix(QF.D) * M
    return Vector(QF.Q) * D

end


##################################################
# Twisting *would* require a different bulge chasing algorithm, so
# we hold it in the type for dispatch
## XXX This is likely not what we want...
struct QRFactorizationTwisted{T, S, Vt, PVt, Rt<:AbstractRFactorization{T,S}} <: AbstractQRFactorizationState{T, S, Val{:twisted}}
N::Int
m::Int
  QF::QFactorizationTwisted{T,S, Vt, PVt}
  RF::Rt
  UV::Vector{Rotator{T,S}}    # Ascending chain for cfreating bulge
  W::Vector{Rotator{T,S}}     # the limb when m > 1, a twisted chain
  A::Matrix{S}
  REIGS::Vector{T}
  IEIGS::Vector{T}
  ctrs::AMRVW_Counter
#  QRFactorizationTwisted{T,S,Vt, PVt, Rt}(N, QF, RF, UV,W, A, REIGS, IEIGS, ctrs) where {T, S, Vt, PVt, Rt} = new(N, QF, RF, UV,W, A, REIGS, IEIGS, ctrs)
#  QRFactorizationTwisted(N::Int, QF::QFactorizationTwisted, RF, UV,W, A, REIGS, IEIGS, ctrs) = new(N, QF, RF, UV, A, REIGS, IEIGS, ctrs)
end

# XXX should pass in m, so that A=m x m matrix
function QRFactorizationTwisted(
                         QF::QFactorizationTwisted{T, S, Vt, PVt},
                         RF::AbstractRFactorization{T, S},
                         m::Int = 0  # specify or use 1 for Complex, 2 for Real, larger if desired
                         ) where {T, S, Vt, PVt}

    if iszero(m)
        m = S <: Complex ? 1 : 2
    end
    N = length(QF) + 1
    M = max(2,m)
    A = zeros(S, M, M)
    reigs = zeros(T, N)
    ieigs = zeros(T, N)
    ctr = AMRVW_Counter(0, 1, N-1, 0, N-2)
    UV = Vector{Rotator{T, S}}(undef, m)   # Ascending Chain
    W = Vector{Rotator{T, S}}(undef, m-1)  # the lim
    QRFactorizationTwisted(length(QF), m,  QF, RF, UV, W, A, reigs, ieigs, ctr)
end


## struct QRFactorizationTwisted{T, S, Rt, QFt, RFt} <: AbstractQRFactorizationState{T, S, Rt, QFt, RFt, Val{:twisted}}
##   N::Int
##   QF::QFt
##   RF::RFt
##   UV::Vector{Rt}    # temp storage for U[V] Vector{Rt}, ,
##   A::Matrix{S}
##   REIGS::Vector{T}
##   IEIGS::Vector{T}
##   ctrs::AMRVW_Counter
## end

## This is not correct!
## This should be matrix multiplication
## as Twisted Q is no long Hessenberg
function diagonal_block(state::QRFactorizationTwisted{T, S}, k) where {T,  S}

    A  = state.A
    QF,  RF =  state.QF, state.RF

    i, j = k-2, k-1

    qji, qjj, qjk, qkj, qkk = QF[j,i], QF[j,j], QF[j,k], QF[k,j], QF[k,k]

    rjj, rjk = RF[j,j], RF[j,k]
    rkk = RF[k,k]
    rij, rik = RF[i,j], RF[i,k]

    ## This is matrix multiplication of a Hessenberg matrix times a upper triangular
    ## (triu(Q,-1) * triu(R))[k-1:k, k-1:k]
    ## could be just sum(QF[i,l] * RF[l,j] for l in i-1:j)

    A[1,1] = qji * rij + qjj * rjj
    A[1,2] = qji * rik + qjj * rjk + qjk * rkk
    A[2,1] = qkj * rjj
    A[2,2] = qkj * rjk + qkk * rkk

    return nothing
end

# Identify m shifts by taking the eigenvalues of the lower  m x m matrix, where k specifies the bottom corner of this matrix
# approximate=true only constructs an approximate  part of the m x m matrix using just the lower rotators
function find_shifts(state::QRFactorizationTwisted{T, S}, m, k; approximate=true) where {T,  S}

    A = state.A # storage space for lower mxm block


end

## return state,m rotators to create a bulge
function create_bulge(state::QRFactorizationTwisted{T, S, Vt, PVt, Rt}; approximate=false) where {T, S, Vt, PVt, Rt}

    m = min(state.m, state.ctrs.stop_index - state.ctrs.zero_index)
    if mod(state.ctrs.it_count, 15) == 0
        @show :randomXXX
        for i in 1:m
            state.UV[i] = random_rotator.(S, state.ctrs.start_index + m - i)
        end
    else
        create_bulge(Val(approximate), Val(S <: Real), state)
    end

    return nothing

end

## approximate, Real
function create_bulge(::Val{true}, ::Val{true}, state::QRFactorizationTwisted{T, S, Vt, PVt, Rt}) where   {T, S, Vt, PVt, Rt}

end

## approximate, Complex
## TODO
## * remove eigvals call for m <= 2
## * make approximate
##
function create_bulge(a::Val{true}, r::Val{false}, state::QRFactorizationTwisted{T, S, Vt, PVt, Rt}) where   {T, S, Vt, PVt, Rt}
#    @show :hi_approx_complex
    ## testing
    ## idea: use Am for eigenvalues
    ## idea: use (A-Q)x iterations to solve inverses
    n = state.ctrs.stop_index - state.ctrs.zero_index
    m = min(state.m, n)

    alpha, beta = state.ctrs.start_index, state.ctrs.stop_index
    offset = alpha
    ## XXX get Am more efficiently
    ## We have state.RF[i,j] and can generate matrix for state.QF
    ## from last m rotators....
    A = Matrix(state.QF) * Matrix(state.RF) # full matrix for exact case
    n = size(A)[1]
    if m > 1
        Am = A[beta+1-m:beta+1, beta+1-m:beta+1]
        rhos = eigvals(Am)
        #@show rhos
    else
        Am = A[beta:beta+1, beta:beta+1]
        #@show Am
        e1, e2 = eigvals(Am)
        if norm(Am[2,2] - e1) < norm(Am[2,2] - e2)
            rhos = [e1]
        else
            rhos = [e2]
        end
    end
    e_alpha = zeros(S, m+1); e_alpha[1] = 1

    Aalpha = A[offset:offset+m, offset:offset+m]
    x = e_alpha
    for i in 1:m
        x = (Aalpha - rhos[i] * I) * x
        pi = (i==1) ? :left : state.QF.Q.pv[i-1]
        if  pi == :right
            x = Aalpha \ x
        end
    end

    # now roll up
    _rollup!(state, x, offset-1)

    return nothing
end

## exact, Real
function create_bulge(approximate::Val{false}, real::Val{true}, state::QRFactorizationTwisted{T, S, Vt, PVt, Rt}) where   {T, S, Vt, PVt, Rt}
#    @show :exact_real
    m = state.m
    @assert m >= 2  # for even we need atleast a quadratic poly
    alpha, beta = state.ctrs.start_index, state.ctrs.stop_index
    offset = alpha - 1
    ## XXX get Am in more efficient manner (use last m rotators...)
    A = Matrix(state.QF) * Matrix(state.RF) # full matrix for exact case
    n = size(A)[1]
    Am = A[beta-m+2:beta+1, beta-m+2:beta+1]
    rhos = _sorteig!(eigvals(Am))  ## sorted by LinearAlgebra.sorteig!

    e_alpha = zeros(S, n); e_alpha[alpha] = 1
    x = e_alpha
    i = 1
    while i <= m
        rhoi = rhos[i]
        pi = (i==1) ? :left : state.QF.Q.pv[i-1]
        if i < m && isapprox(rhoi, conj(rhos[i+1]))
            ## XXX deal with left right...
            #@show :complex_roots
            x = (A^2 - 2real(rhos[i])*A + norm(rhos[i])^2*I) * x
            i += 1
        else
            x = (A - rhos[i]*I) * x
        end

        ## XXX fix me
        if  pi == :right
            x = A \ x
        end
        i += 1
    end

    # now roll up
    _rollup!(state, x[offset+1:offset+m+1])

    return nothing

end

## exact, complex
function create_bulge(approximate::Val{false}, real::Val{false}, state::QRFactorizationTwisted{T, S, Vt, PVt, Rt}) where   {T, S, Vt, PVt, Rt}

    m = min(state.m, state.ctrs.stop_index - state.ctrs.zero_index)
    alpha, beta = state.ctrs.start_index, state.ctrs.stop_index
    offset = alpha - 1
    A = Matrix(state.QF) * Matrix(state.RF) # full matrix for exact case
    n = size(A)[1]
    if m > 1
        Am = A[beta-m+1:beta+1, beta-m+1:beta+1]
        rhos = eigvals(Am)
    else
        Am = A[beta:beta+1, beta:beta+1]
        e1, e2 = eigvals(Am)
        if norm(Am[2,2] - e1) < norm(Am[2,2] - e2)
            rhos = [e1]
        else
            rhos = [e2]
        end
    end

    e_alpha = zeros(S, n); e_alpha[alpha] = 1

    x = e_alpha
    for i in 1:m
        x = (A - rhos[i] * I) * x
        pi = (i==1) ? :left : state.QF.Q.pv[i-1]
        if  pi == :right
            x = A \ x
        end
    end

    # now roll up
    _rollup!(state, view(x, offset+1:offset+m+1), offset)

    return nothing
end

function _rollup!(state, x, offset=0)

    r = x[end]
    m = length(x)
    for i in (m-1):-1:1
        c, s, r = givensrot(x[i], r)
        U = Rotator(c, s, i+offset)
        state.UV[m-i] = U'
    end
end


##################################################

function eigvals(state::QRFactorizationTwisted)

    AMRVW_algorithm(state)
    es = complex.(state.REIGS, state.IEIGS)
    LinearAlgebra.sorteig!(es)
    es

end

## We don't go for the fastest method
## We have the first n-m  steps through, a straightening out of the twisted factorization in QF,
## This should be  O(n^2) steps
## Then, once untwisted, we pass down to the more efficient O(n^2) algorithm (in time) to solve
##  XXX  work in deflation, etc...
function AMRVW_algorithm(state::QRFactorizationTwisted{T, S}) where  {T, S}

    n = length(state.QF.Q)

    N0 = findlast(x->x==:right, state.QF.Q.pv)
    N = N0 == nothing ? 1 : N0


    it_max = 20 * length(state)
    kk = 0

    while kk <= it_max

        kk += 1
        state.ctrs.stop_index <= 0 && return     ## finished up!
        state.ctrs.it_count += 1


        check_deflation(state)
        show_status(state)



        ## Delta is number of rotators in piece being considered
        ## m must be less than delta
        k = state.ctrs.stop_index
        delta =  state.ctrs.stop_index - state.ctrs.zero_index  ## number of rotators in Q factorization

        if delta == 1

            diagonal_block(state,  k + 1)
            a11,a12,a21,a22 = state.A[1,1], state.A[1,2], state.A[2,1], state.A[2,2]
            e1r,e1i, e2r,e2i = eigen_values(a11, a12, a21, a22)

            state.REIGS[k], state.IEIGS[k] = e1r, e1i
            state.REIGS[k+1], state.IEIGS[k+1] = e2r, e2i
            @show :add_2, complex(e1r, e1i), complex(e2r, e2i)

            # can finish up if near end
            if state.ctrs.stop_index == 2
                diagonal_block(state, 2)
                e1 = state.A[1,1]

                @show :add_1_all_done, e1

                state.REIGS[1] = real(e1)
                state.IEIGS[1] = imag(e1)
            end

            state.ctrs.zero_index = 0
            state.ctrs.start_index = 1
            state.ctrs.stop_index -= 2

        elseif delta <= 0

            diagonal_block(state, k + 1)
            a11,a12,a21,a22 = state.A[1,1], state.A[1,2], state.A[2,1], state.A[2,2]
            #e1, e2 = a11, a22
            e1r,e1i, e2r,e2i = eigen_values(a11, a12, a21, a22)
            e1 = complex(e1r, e1i)
            e2 = complex(e2r, e2i)

            if state.ctrs.stop_index == 1
                @show :add_2, e1, e2
                state.REIGS[1], state.IEIGS[1] = real(e1), imag(e1)
                state.REIGS[2], state.IEIGS[2] = real(e2), imag(e2)
                state.ctrs.stop_index = 0
            else
                @show :add_1, e2, e1
                state.REIGS[k+1], state.IEIGS[k+1] = real(e2), imag(e2)
                state.ctrs.zero_index = 0
                state.ctrs.start_index = 1
                state.ctrs.stop_index -= 1
            end

        ## elseif m > delta

        ##     δ, Δ = state.ctrs.start_index, state.ctrs.stop_index
        ##     @show  δ, Δ, m, delta
        ##     @show "Need to take $delta $m eigemvalues here, ..."
        ##     es = eigvals(Matrix(state)[δ:Δ+1,δ:Δ+1])
        ##     for (ind, ev) in enumerate(es)
        ##         state.REIGS[δ+ind-1] = real(ev)
        ##         state.IEIGS[δ+ind-1] = imag(ev)
        ##     end
        ##     @show k
        ##     @show state.QF.Q.x
        ##     state.ctrs.zero_index = 0
        ##     state.ctrs.start_index = 1
        ##     state.ctrs.stop_index -= (m-1)

        else

            bulge_step(state)

        end

    end


    ## while N >= 0

    ##     bulge_step(state)

    ##     check for deflation of last one
    ##     c,s = vals(Ms[state.ctrs.stop_index])
    ##     if abs(s) <= eps(T)
    ##         state.ctrs.stop_index -= 1
    ##     end

    ##     N -= m
    ## end

    ## now untwisted, so we can change over
    ## QF = state.QF
    ## Ms = QF.Q.x
    ## QF_new = QFactorization(DescendingChain(Ms), QF.D)
    ## new_state = QRFactorization(QF_new, state.RF)
    ## AMRVW_algorithm(new_state)
    ## new_state

end




## We use `bulge_step!` for more general usage; this is
## tied to AbstractQRFactorizationState{T, S, Rt, QFt, RFt, Val{:twisted}}
function bulge_step(state::QRFactorizationTwisted{T, S, V, Rt})  where {T,S, V, Rt}

    # create bulge is only right when  all rotators in QFt are straightened out,
    # but  this should step in the correct direction.
    create_bulge(state)
    n = state.ctrs.stop_index - state.ctrs.zero_index + 1
    m = min(state.m, n-1)
    m < state.m && @show :m_less

    Asc = copy(state.UV[1:m])  # stored as U, V ## XXX copy!
    Ms = state.QF.Q
    D = state.QF.D
    RF = state.RF
    δ, Δ = state.ctrs.start_index, state.ctrs.stop_index

    # modify Ms, D, RF
    bulge_step!(n, view(Ms, δ:Δ), D, RF, Asc)

    return nothing

end

## Can enter here for more general usage
## Note: Ms is a Twisted Chain, but need not stay that way
## We could reuse the storage space from Ms to hold the values stored in Decoupled
## at the cost of keeping track of the respective boundaries (k in the case of Decoupled, and k+m)
## for Ms. This could save some allocations, as this step allocates a new Decoupled vector
## each pass through.
function bulge_step!(n, Ms::TwistedChain, D, RF::AbstractRFactorization, Asc, choice = :left)


    Rt = eltype(Ms)
    #Decoupled = Rt[]

    Des = reverse(adjoint.(Asc))


    m = length(Asc)

    # psd[k] = ps[m+k] with padding soecufued by choice
    ps = copy(Ms.pv)
    psd = ps[m+1:end]
    if isa(choice, Symbol)
        append!(psd, repeat([choice], min(m, length(Ms.pv))))
    else
        append!(psd, choice)
    end

    ## M0a = Des * (Ms * (Matrix(D) * (Matrix(RF)* Asc)))
    ## @show eigvals(M0a)[1]

    ## Must put AscendingChain into place
    ## py passing through RF and D <--
    passthrough!(RF, AscendingChain(Asc))
    passthrough!(D, AscendingChain(Asc))

#    M0 = Des * (Ms * (Asc * (Matrix(D) * Matrix(RF))))
#    @show eigvals(M0)[1]

    limb, limb_side = step_0!(m, ps, psd, Ms, Des, Asc, D, RF)


    ## _Ms = view(Ms, (m+1):length(Ms))
    ## if limb_side == :left
    ##     M1 = limb * (Des * (_Ms * (Asc * (Matrix(D) * Matrix(RF)))))
    ##     @show eigvals(M1)[1]
    ## else
    ##     M1 = Des * (_Ms * (Asc * (limb *(Matrix(D) * Matrix(RF)))))
    ##     @show eigvals(M1)[1]
    ## end

    for k in 1:(n-m-2)
        limb_side = step_k!(k, n, m, psd, limb_side, limb, Des, Asc, Ms, D, RF)
    end

    ## @show :step_k
    ## _Decoupled = TwistedChain(Decoupled, psd[1:length(Decoupled)-1])
    ## if limb_side == :left
    ##     if psd[n-m-2] == :left
    ##         Mkll = _Decoupled * (limb * (Des * (Ms * (Asc * (Matrix(D) * Matrix(RF))))))
    ##         @show eigvals(Mkll)[1]
    ##     else
    ##         Mklr = limb * (Des * (Ms * (Asc * (_Decoupled * (Matrix(D) * Matrix(RF))))))
    ##         @show eigvals(Mklr)[1]
    ##     end
    ## else
    ##     if psd[n-m-2] == :left
    ##         Mkrl = _Decoupled * (Des * (Ms * (Asc * (limb *(Matrix(D) * Matrix(RF))))))
    ##         @show eigvals(Mkrl)[1]
    ##     else
    ##         Mkrr = Des * (Ms * (Asc * (limb * (_Decoupled *(Matrix(D) * Matrix(RF))))))
    ##         @show eigvals(Mkrr)[1]
    ##     end
    ## end

    # use new choice
    if n-m-1 > 0 # case n=m skips this
        limb_side = step_k!(n-m-1, n, m, psd, limb_side, limb, Des, Asc, Ms, D, RF)
    end

    ## @show :step_kk
    ## _Decoupled = TwistedChain(Decoupled, psd[1:length(Decoupled)-1])
    ## if limb_side == :left
    ##     if psd[n-m-1] == :left
    ##         Mkll = _Decoupled * (limb * (Des * (Ms * (Asc * (Matrix(D) * Matrix(RF))))))
    ##         @show eigvals(Mkll)[1]
    ##     else
    ##         Mklr = limb * (Des * (Ms * (Asc * (_Decoupled * (Matrix(D) * Matrix(RF))))))
    ##         @show eigvals(Mklr)[1]
    ##     end
    ## else
    ##     if psd[n-m-1] == :left
    ##         Mkrl = _Decoupled * (Des * (Ms * (Asc * (limb *(Matrix(D) * Matrix(RF))))))
    ##         @show eigvals(Mkrl)[1]
    ##     else
    ##         Mkrr = Des * (Ms * (Asc * (limb * (_Decoupled *(Matrix(D) * Matrix(RF))))))
    ##         @show eigvals(Mkrr)[1]
    ##     end
    ## end



    # now knit in  limb, Des, Asc
    step_knit!(n, m, psd, limb_side, limb, Des,  Asc, Ms, D, RF)

    ## Decoupled is now Ms, assign and update position vector
    #    append!(Ms.x, Decoupled)
    Ms.pv[:] = psd
#    empty!(Ms.pv)
#    append!(Ms.pv, psd)

    return nothing
end


############################################################
##
## Implement thee different steps for one pass of Francis' algorithm


function step_0!(m,  ps, psd, Ms::TwistedChain, Des, Asc, D, RF) where {T}

    has_limb = ifelse(m > 1, true, false)

    # grab limb, but keep order
    # limb is a copy, not a view
    limb = TwistedChain(copy(Ms.x[1:(m-1)]), copy(Ms.pv[1:(m-2)]))
    limb_side = :nothing

    #@show idx.(Des), idx.(Asc), idx.(limb), idx.(Ms), m

    if has_limb
        if ps[m-1] == :left
            passthrough!(DescendingChain(Des), limb)  # pass limb through left
            limb_side = :left
        else
            passthrough!(limb, AscendingChain(Asc))   # pass limb to right
            limb_side = :right
        end
    end



    U = Ms.x[m]
    psm = m <= length(ps) ? ps[m] : :right # doesn't matter in this case

    if psm == :left
        Des[end], Di = fuse(Des[end], U)  # fuse with descending; aka Ms[m]

        ## need to passthrough Ms too now
        i = idx(Di)
        dMs = descending_part(Ms, i) # stop_index

        if limb_side == :right
            passthrough_phase!(Di, dMs,  (AscendingChain(Asc), limb), D)
        else
            passthrough_phase!(Di, dMs,  (AscendingChain(Asc), ), D)
        end

    else

        AA = popfirst!(Asc)
        AA, Di = fuse(U, AA)

        M2 =  Des * (view(Ms, (m+1):length(Ms)) * (AA * (Di * (Asc * (limb *(Matrix(D) * Matrix(RF)))))))
        #@show eigvals(M2)[1]
        #@show idx.(Des), idx.(Asc)
        if limb_side == :right
            passthrough_phase!(Di, (AscendingChain(Asc), limb), D)
        else
            passthrough_phase!(Di, (AscendingChain(Asc), ), D)
        end
        pushfirst!(Asc, AA)

        ## Asc[1], Di = fuse(U, Asc[1])



        ## if limb_side == :right
        ##     @show :right
        ##     passthrough_phase!(Di, (AscendingChain(view(Asc,2:length(Asc))), limb), D)
        ## else
        ##     passthrough_phase!(Di, (AscendingChain(view(Asc,2:length(Asc))), ), D)
        ## end

    end

    limb, limb_side

end


function step_k!(k, n, m, psd, limb_side, limb, Des, Asc, Ms, D, RF) where {T}

    phatk = psd[k]

    ## 4 cases based on:
    ## phatk/limb_side
    ## A bit redundandant, but easier to verify

    U = Ms.x[m+k]

    if phatk == :left && limb_side == :left

        push!(Des, U)  ## Augment Descending
        passthrough!(DescendingChain(Des), AscendingChain(Asc)) ## translate ascending

        passthrough!(limb, AscendingChain(Asc))

        ## similarity transform Asc to right side
        passthrough!(RF, AscendingChain(Asc)) # <-- pass Asc through RF
        passthrough!(D, AscendingChain(Asc))

        U = popfirst!(Des)

        if m > 1
            limb_side = :left
        end


    elseif phatk == :left && (limb_side == :right || limb_side == :nothing)

        push!(Des, U)  ## Augment Descending
        passthrough!(DescendingChain(Des), AscendingChain(Asc)) ## translate ascending
        passthrough!(DescendingChain(Des), limb)

        ## similarity transform Asc to right side
        passthrough!(RF, AscendingChain(Asc)) # <-- pass Asc through RF
        passthrough!(D, AscendingChain(Asc))

        U = popfirst!(Des)

        if m > 1
            limb_side = :left
        end

    elseif phatk == :right && limb_side == :left

        pushfirst!(Asc, U)
        passthrough!(DescendingChain(Des), AscendingChain(Asc))

        passthrough!(limb, AscendingChain(Asc))

        ## similarity transform descending to left side
        passthrough!(DescendingChain(Des), D)
        passthrough!(DescendingChain(Des), RF) #--> pass Des through Rf;  similarity transform

        U = pop!(Asc)

        if m > 1
            limb_side = :right
        end

    elseif phatk == :right && (limb_side == :right || limb_side == :nothing)

        pushfirst!(Asc, U)
        passthrough!(DescendingChain(Des), AscendingChain(Asc))

        passthrough!(DescendingChain(Des), limb)

        ## similarity transform descending to left side
        passthrough!(DescendingChain(Des), D)
        passthrough!(DescendingChain(Des), RF) #--> pass Des through Rf;  similarity transform

        U = pop!(Asc) # k = idx(U)

        if m > 1
            limb_side = :right
        end

    end

    if k == 1 || psd[k-1] == :left
        Ms.x[k] = U # Decoupled is in the k the position
    else
        Ms.x[k] = U
    end

    return limb_side

end

## steps k=n-m to  n-2
function step_knit!(n, m, psd, limb_side, limb, Des,  Asc, Ms, D, RF) where {T}
    # make bottom
    U = pop!(Des)
    V = popfirst!(Asc)
    bottom, Di = fuse(U,V)

    ## Decoupled too! Decoupled on right? as we augmented Ascending
    ## have a complicted check, as m can be big
    psdm = length(psd)-m+1 > 0 ? psd[end-m+1] : :left

    if psdm == :right
        cDecoupled = view(Ms, 1:n-m)
        if limb_side == :right
            passthrough_phase!(Di, (AscendingChain(Asc), Decoupled, limb), D)
        else
            passthrough_phase!(Di, (AscendingChain(Asc),  Decoupled), D)
        end
    else
        if limb_side == :right
            passthrough_phase!(Di, (AscendingChain(Asc), limb), D)
        else
            passthrough_phase!(Di, (AscendingChain(Asc),), D)
        end
    end


    ## Final Steps, knit in
    ## setps k=n-m to  n-2

    for k in (n-m):(n-2)
        phatk = k > 1 ? psd[k-1] : :left
        if k > (n-m) && k > 2 && psd[k-2] != phatk
            ## need to reposition V structure to other side
            if phatk == :left
                passthrough!(RF, DescendingChain(Des))
                passthrough!(D, DescendingChain(Des))
                bottom = passthrough!(RF, bottom)
                bottom = passthrough!(D, bottom)
                passthrough!(RF, AscendingChain(Asc))
                passthrough!(D, AscendingChain(Asc))
            else
                passthrough!(AscendingChain(Asc),D)
                passthrough!(AscendingChain(Asc),RF)
                bottom = passthrough!(bottom, D)
                bottom = passthrough!(bottom, RF)
                passthrough!(DescendingChain(Des),D)
                passthrough!(DescendingChain(Des),RF)
            end
            # flip limb side
            if limb_side == :right
                limb_side = :left
            elseif limb_side == :left
                limb_side = :right
            end
        end

        ## lengthen one side of V
        if phatk == :left
            Des = push!(Des, bottom)
        else
            Asc = pushfirst!(Asc, bottom)
        end
        passthrough!(DescendingChain(Des), AscendingChain(Asc))

        # remove bottom rotator from limb; keep track of if it is right or left
        lpk = length(limb.pv) > 0 ? limb.pv[end] : :nothing
        L, lpk =  pop!(limb)
        ## we have
        ## 4 cases in terms of limb_side and lpk:
        ## :l, :l -> fuse (Asc, Des), translate
        ## :l, :r -> translate, fuse (Asc, limb, Des)
        ## :r, :l -> translate, fuse (limb)
        ## :r, :r -> fuse (.), translate
        if limb_side == :left || limb_side == :nothing
            if lpk == :left || lpk == :nothing

                # fuse then translate
                Asc[1], Di = fuse(L, Asc[1])
                if length(Asc) > 1
                    AA = Asc[2:end]
                    passthrough_phase!(Di, (AscendingChain(AA), DescendingChain(Des)), D)
                    Asc[2:end] = AA
                else
                    passthrough_phase!(Di, (DescendingChain(Des), ), D)
                end

                ## translate
                length(limb) > 0 && passthrough!(limb, AscendingChain(Asc))

            elseif lpk == :right

                # translate then fuse
                length(limb) > 0 && passthrough!(limb, AscendingChain(Asc))

                # then fuse
                Asc[1], Di = fuse(L, Asc[1])
                if length(Asc) > 1
                    AA = Asc[2:end]
                    passthrough_phase!(Di, (AscendingChain(AA), limb, DescendingChain(Des)), D)
                    Asc[2:end] = AA
                else
                    passthrough_phase!(Di, (limb, DescendingChain(Des)), D)
                end

            end

        else # limb on right

            if  lpk == :left

                # translate then fuse
                length(limb) > 0 && passthrough!(DescendingChain(Des), limb)

                Des[end], Di = fuse(Des[end], L)
                passthrough_phase!(Di, (), D)
            else
                # fuse  then translate
                Des[end], Di = fuse(Des[end], L)
                passthrough_phase!(Di, limb, (), D)

                length(limb) > 0 && passthrough!(DescendingChain(Des), limb)
            end
        end

        if phatk == :left
            # similarity to move Asc to right side
            passthrough!(RF, AscendingChain(Asc))  # <--
            passthrough!(D, AscendingChain(Asc))

            limb_side = :left
        else
            # similarity to move Des to left side
            passthrough!(DescendingChain(Des), D)
            passthrough!(DescendingChain(Des), RF)

            limb_side = :right

        end

        ## pop off top, and add to Decopuled
        ## we have idx(U) ==   k; so psd[k-1] = phatk indicates direction to  add
        if phatk == :left
            U  = popfirst!(Des)
            Ms.x[k] = U
            #push!(Decoupled, U)
        else
            U = pop!(Asc)
            #push!(Decoupled, U)
            Ms.x[k] = U
            #pushfirst!(Decoupled, U)
        end

        U = pop!(Des)
        V = popfirst!(Asc)
        bottom, Di = fuse(U,V)

        if phatk == :right
            Decoupled = view(Ms, 1:k)
            if limb_side == :right
                passthrough_phase!(Di, (AscendingChain(Asc), Decoupled, limb), D)
            else
                passthrough_phase!(Di, (AscendingChain(Asc), Decoupled), D)
            end
        else
            if limb_side == :right
                passthrough_phase!(Di, (AscendingChain(Asc), limb), D)
            else
                passthrough_phase!(Di, (AscendingChain(Asc),), D)
            end
        end

    end

    # just the bottom left, but may be on wrong side if m > 1

    if psd[end] == :left
        if m > 1 && (length(psd) == 1 || psd[end-1] != psd[end])
            bottom = passthrough!(RF, bottom)  # move to other side
            bottom = passthrough!(D, bottom)
        end
        Ms.x[end] = bottom

    else
        if m > 1 && (length(psd) == 1 || psd[end-1] != psd[end])
            bottom = passthrough!(bottom, D)   # move bottom to other side
            bottom = passthrough!(bottom, RF)
        end
        Mx.x[end] = bottom
        #push!(Decoupled, bottom)
    end

    return nothing

end
