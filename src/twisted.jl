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
## ## XXX This is likely not what we want...
## struct QRFactorizationTwisted{T, S, Vt, PVt, Rt<:AbstractRFactorization{T,S}} <: AbstractQRFactorizationState{T, S, Val{:twisted}}
## N::Int
## m::Int
##   QF::QFactorizationTwisted{T,S, Vt, PVt}
##   RF::Rt
##   UV::Vector{Rotator{T,S}}    # Ascending chain for cfreating bulge
##   W::Vector{Rotator{T,S}}     # the limb when m > 1, a twisted chain
##   A::Matrix{S}
##   REIGS::Vector{T}
##   IEIGS::Vector{T}
##   ctrs::AMRVW_Counter
## #  QRFactorizationTwisted{T,S,Vt, PVt, Rt}(N, QF, RF, UV,W, A, REIGS, IEIGS, ctrs) where {T, S, Vt, PVt, Rt} = new(N, QF, RF, UV,W, A, REIGS, IEIGS, ctrs)
## #  QRFactorizationTwisted(N::Int, QF::QFactorizationTwisted, RF, UV,W, A, REIGS, IEIGS, ctrs) = new(N, QF, RF, UV, A, REIGS, IEIGS, ctrs)
## end

## # XXX should pass in m, so that A=m x m matrix
## function QRFactorizationTwisted(
##                          QF::QFactorizationTwisted{T, S, Vt, PVt},
##                          RF::AbstractRFactorization{T, S},
##                          m::Int = 0  # specify or use 1 for Complex, 2 for Real, larger if desired
##                          ) where {T, S, Vt, PVt}

##     if iszero(m)
##         m = S <: Complex ? 1 : 2
##     end
##     N = length(QF) + 1
##     M = max(2,m)
##     A = zeros(S, M, M)
##     reigs = zeros(T, N)
##     ieigs = zeros(T, N)
##     ctr = AMRVW_Counter(0, 1, N-1, 0, N-2)
##     UV = Vector{Rotator{T, S}}(undef, m)   # Ascending Chain
##     W = Vector{Rotator{T, S}}(undef, m-1)  # the lim
##     QRFactorizationTwisted(length(QF), m,  QF, RF, UV, W, A, reigs, ieigs, ctr)
## end

struct QRFactorizationTwisted{T, S, Vt, PVt, Rt<:AbstractRFactorization{T,S}} <: AbstractQRFactorizationState{T, S, Val{:twisted}}
  QF::QFactorizationTwisted{T,S, Vt, PVt}
  RF::Rt
end

# XXX should pass in m, so that A=m x m matrix
## function QRFactorizationTwisted(
##                          QF::QFactorizationTwisted{T, S, Vt, PVt},
##                          RF::AbstractRFactorization{T, S}
##                          ) where {T, S, Vt, PVt}

##     QRFactorizationTwisted(QF, RF)
## end


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
## function diagonal_block(state::QRFactorizationTwisted{T, S}, k) where {T,  S}
##     A  = state.A
##     QF,  RF =  state.QF, state.RF

##     i, j = k-2, k-1

##     qji, qjj, qjk, qkj, qkk = QF[j,i], QF[j,j], QF[j,k], QF[k,j], QF[k,k]

##     rjj, rjk = RF[j,j], RF[j,k]
##     rkk = RF[k,k]
##     rij, rik = RF[i,j], RF[i,k]

##     ## This is matrix multiplication of a Hessenberg matrix times a upper triangular
##     ## (triu(Q,-1) * triu(R))[k-1:k, k-1:k]
##     ## could be just sum(QF[i,l] * RF[l,j] for l in i-1:j)

##     A[1,1] = qji * rij + qjj * rjj
##     A[1,2] = qji * rik + qjj * rjk + qjk * rkk
##     A[2,1] = qkj * rjj
##     A[2,2] = qkj * rjk + qkk * rkk

##     return nothing
## end

## fill in A[j:k+1, j:k+1]
## using rotators Q.x[j], ..., Q.x[k]
## Approximate in general, as not all rotators are used, but
##
function diagonal_block!(A, QF::QFactorizationTwisted, RF, j, k)
    approx_diagonal_block!(A, QF, RF, j, k)
end

function approx_diagonal_block!(A, QF, RF, j, k) where {T, S}
    R = RF
    D = QF.D
    Q = QF.Q
    k = min(k, length(Q))
    # make 1:(k-j+1) identity
    A .= zero(eltype(A))
    if j == 1 #|| iszero(Q.x[j].s)
        _approx_diagonal_block!(A, Q, D, R, j, k)
    else
        _approx_diagonal_block!(A, Q, D, R, j-1, k)
        # shift up and left by 1
        m,n = size(A)
        for j in 1:n-1  # across columns, then down rows
            for i in 1:m-1
                A[i,j] = A[i+1, j+1]
            end
        end
    end

    nothing
end

function _approx_diagonal_block!(A, Q, D, R, j, k)

    S = eltype(A)
    for i in j:k+1
        for ii in j:(i-1)
            A[i-j+1,ii-j+1] = zero(S)
        end
        for ii in i:k+1
            A[i-j+1,ii-j+1] = D[i] *  R[i,ii]
        end
    end

    # we get first one
    Tw = view(Q, j:k)
    for i in reverse(position_vector_indices(Tw.pv))  # allocates here
        U = Tw.x[i]
        V = Rotator(vals(U)..., idx(U) - j + 1)
        mul!(V, A)
    end

    return  nothing
end




## return state,m rotators to create a bulge
function create_bulge(state::QRFactorizationTwisted{T, S, Vt, PVt, Rt}, storage, ctr, m) where {T, S, Vt, PVt, Rt}


    if mod(ctr.it_count, 15) == 0
        #@show :randomXXX
        for i in 1:m
            storage.VU[i] = random_rotator.(S, ctr.start_index + m - i)
        end
    else
        create_bulge(Val(S <: Real), state, storage, ctr, m)
    end

    return nothing

end

## Real
function create_bulge(::Val{true}, state::QRFactorizationTwisted{T, S, Vt, PVt, Rt}, storage, ctr, m) where   {T, S, Vt, PVt, Rt}
    #    @show :exact_real

    @assert m >= 2  # for even we need atleast a quadratic poly

    alpha, beta = ctr.start_index, ctr.stop_index
    n = beta - alpha
    A = storage.A


    offset = alpha - 1

    diagonal_block!(A, state.QF, state.RF, beta-m+1, beta)

    Am = view(A, 2:m+1, 2:m+1)

    if m == 2

        rhos = eigen_values(Am)

    else

        rhos = _eigvals(Am)  ## sorted by LinearAlgebra.sorteig!

    end

    e_alpha = zeros(S, m+1)
    e_alpha[1] = one(S)

    diagonal_block!(A, state.QF, state.RF, alpha, alpha + m - 1)
    Am = view(A, 1:m+1, 1:m+1)

    x = e_alpha
    i = 1
    while i <= m
        rhoi = rhos[i]
        pi = (i==1) ? :left : state.QF.Q.pv[i-1]
        pj = state.QF.Q.pv[i]

        if i < m && isapprox(rhoi, conj(rhos[i+1]))
            ## XXX deal with left right...
            #@show :complex_roots
            x = (Am^2 - 2real(rhos[i])*Am + norm(rhos[i])^2*I) * x
            if pi == :right
                x = Am \ x
            end
            if pj == :right
                x = Am \ x
            end
            i += 1
        else
            x = (Am - real(rhos[i])*I) * x
        end

        i += 1
    end

    # now roll up
    _rollup!(storage.VU, x, alpha-1)
    return nothing

end

## complex
function create_bulge(real::Val{false}, state::QRFactorizationTwisted{T, S, Vt, PVt, Rt}, storage, ctr, m) where   {T, S, Vt, PVt, Rt}

    alpha, beta = ctr.start_index, ctr.stop_index
    n = beta -  alpha
    A = storage.A

    diagonal_block!(A, state.QF, state.RF, beta-m+1, beta)

    if m == 1

        Am = view(A, 1:2, 1:2)
        e1, e2 = eigen_values(Am)

        if norm(Am[2,2] - e1) < norm(Am[2,2] - e2)
            rhos = [e1]
        else
            rhos = [e2]
        end

    elseif m == 2

        Am = view(A, 2:3, 2:3) #A[beta:beta+1, beta:beta+1]
        rhos = eigen_values(Am)

    else

        Am = view(A, 2:m+1, 2:m+1)
        rhos = _eigvals(Am)

    end

    e_alpha = zeros(S, m+1)
    e_alpha[1] = one(S)

    diagonal_block!(A, state.QF, state.RF, alpha, alpha + m - 1)
    Am = view(A, 1:m+1, 1:m+1)

    x = e_alpha
    for i in 1:m
        x = (Am - rhos[i] * I) * x
        pi = (i==1) ? :left : state.QF.Q.pv[i-1]
        if  pi == :right
            x = Am \ x
        end
    end

    # now roll up
    _rollup!(storage.VU, x, alpha - 1)

    return nothing
end

function _rollup!(VU, x, offset=0)

    r = x[end]
    m = length(x)
    for i in (m-1):-1:1
        c, s, r = givensrot(x[i], r)
        U = Rotator(c, s, i+offset)
        VU[m-i] = U'
    end
end


##################################################

function eigvals(state::QRFactorizationTwisted)

    es = AMRVW_algorithm!(state)
    LinearAlgebra.sorteig!(es)
    es

end

function create_storage(QF::QFactorization{T,S}, m) where {T,  S}
    A = zeros(S, 4, 4)
    VU = Vector{Rotator{T, S}}(undef, m)   # Ascending Chain
    limb = DescendingChain(copy(QF.Q.x[1:(m-1)]))
    (A=A, VU=VU, limb=limb)
end

function create_storage(QF::QFactorizationTwisted{T,S}, m) where {T,  S}
    A = zeros(S, max(4, m+2), max(4, m+2))
    VU = Vector{Rotator{T, S}}(undef, m)   # Ascending Chain
    limb = TwistedChain(copy(QF.Q.x[1:(m-1)]), copy(QF.Q.pv[1:(m-2)]))
    (A=A, VU=VU, limb=limb)
end

function AMRVW_algorithm!(state::AbstractQRFactorizationState{T, S,  Twt}) where  {T, S, Twt}
    m = S <: Real ? 2 : 1
    AMRVW_algorithm!(state, m)
end


function AMRVW_algorithm!(state::AbstractQRFactorizationState{T, S, Twt}, m) where  {T, S, Twt}

    QF, RF = state.QF, state.RF
    # set counters
    n = length(QF.Q)
    ctr = AMRVW_Counter(0, 1, n, 0, n-1)

    # set up storage
    storage = create_storage(QF, m)
    A = storage.A  # cleanup
    EIGS = zeros(Complex{T}, n+1)



    it_max = 20*n
    kk = 0

    while kk <= it_max

        if ctr.stop_index <= 0
            break
        end


        #show_status(state, ctr)

        kk += 1
        ctr.it_count += 1


        check_deflation(state, ctr)


        ## Delta is number of rotators in piece being considered
        ## m must be less than delta tto run a bulge step
        k = ctr.stop_index
        delta = ctr.stop_index - ctr.zero_index  ## number of rotators in Q factorization

        if delta == 1


            diagonal_block!(A, QF, RF, k, k+1)
            e1::Complex{T}, e2::Complex{T} = eigen_values(A)

            EIGS[k] = e1
            EIGS[k+1] = e2

            # can finish up if near end
            if ctr.stop_index == 2

                diagonal_block!(A, QF, RF, 1, 2)
                e1 = A[1,1]

                EIGS[1] = e1

                break
            end

            ctr.zero_index = 0
            ctr.start_index = 1
            ctr.stop_index -= 2

        elseif delta <= 0

            diagonal_block!(A, QF, RF, k, k+1)
            e1, e2 = A[1,1], A[2,2]  # 0s off diagonal

            if ctr.stop_index == 1
                # finish up

                EIGS[1], EIGS[2] = e1, e2
                ctr.stop_index = 0

                break
            else

                EIGS[k+1] = e2

                ctr.zero_index = 0
                ctr.start_index = 1
                ctr.stop_index -= 1

            end

        elseif m > 2 && m >= delta - 1
            ## When finding shifts we find m eigenvalues from some other means
            ## so here we find the eigvals

            δ, Δ = ctr.start_index, ctr.stop_index
            diagonal_block!(A, QF, RF, δ, Δ+1)

            es = _eigvals(view(A,1:Δ-δ+2,1:Δ-δ+2))
            for (ind, ev) in enumerate(es)
                EIGS[δ+ind-1] = ev
            end

            if  δ == 2
                # finish up
                diagonal_block!(A, QF, RF, 1, 2)

                e1 = A[1,1] # δ == 2 means first rotator a diagonal, so product will be diagonal in 1:2x1:2 block
                EIGS[1] = e1

                break
            end

            ctr.zero_index = 0
            ctr.start_index = 1
            ctr.stop_index = δ - 2

        else

            bulge_step(state, storage, ctr, m)

        end

    end

    LinearAlgebra.sorteig!(EIGS)

    EIGS

end




## We use `bulge_step!` for more general usage; this is
## tied to AbstractQRFactorizationState{T, S, Rt, QFt, RFt, Val{:twisted}}
function bulge_step(state::QRFactorizationTwisted{T, S, V, Rt}, storage, ctr, m)  where {T,S, V, Rt}

    create_bulge(state, storage, ctr, m)

    n = ctr.stop_index - ctr.zero_index + 1

    Asc = copy(storage.VU[1:m])  # stored as V,U ## XXX copy! Could eliminate this
    Ms = state.QF.Q
    D = state.QF.D
    RF = state.RF
    δ, Δ = ctr.start_index, ctr.stop_index

    bulge_step!(n, m, view(Ms, δ:Δ), D, RF, Asc)

    return nothing

end

## Can enter here for more general usage
## n is the size of matrix, not the number of rotators
function bulge_step!(n, m, Ms::TwistedChain, D, RF::AbstractRFactorization, Asc, choice = :left)


    ## XXX
    ## We *can* save memory if we
    ## * didn't copy VU to Asc
    ## * put Asc, Des in storage and didn't adjust sizes in step_knit
    ## * used the limb defined in create storage
    ## This would require some reworking of step_knit!

    Des = reverse(adjoint.(Asc))

    # psd[k] = ps[m+k] with padding soecufued by choice
    ps = copy(Ms.pv)
    psd = ps[m+1:end]
    if isa(choice, Symbol)
        append!(psd, repeat([choice], min(m, length(Ms.pv))))
    else
        append!(psd, choice)
    end


    ## Must put AscendingChain into place
    ## py passing through RF and D <--
    passthrough!(RF, AscendingChain(Asc))
    passthrough!(D, AscendingChain(Asc))

    limb, limb_side = step_0!(m, ps, psd, Ms, Des, Asc, D, RF)


    for k in 1:(n-m-2)

        limb_side = step_k!(k, n, m, psd, limb_side, limb, Des, Asc, Ms, D, RF)
    end

    # use new choice
    if n-m-1 > 0

        limb_side = step_k!(n-m-1, n, m, psd, limb_side, limb, Des, Asc, Ms, D, RF)
    end



    # now knit in  limb, Des, Asc
    step_knit!(n, m, psd, limb_side, limb, Des,  Asc, Ms, D, RF)

    Ms.pv[:] = psd

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

        if limb_side == :right
            passthrough_phase!(Di, (AscendingChain(Asc), limb), D)
        else
            passthrough_phase!(Di, (AscendingChain(Asc), ), D)
        end
        pushfirst!(Asc, AA)

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

    Ms.x[k] = U

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
        Decoupled = view(Ms, 1:n-m)
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
#            @show :reposition, phatk
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
#        @show limb_side, lpk
        if limb_side == :left || limb_side == :nothing
            if lpk == :left || lpk == :nothing

                # fuse then translate
                Asc[1], Di = fuse(L, Asc[1])
                if length(Asc) > 1
                    AA = Asc[2:end]
                    passthrough_phase!(Di,  (AscendingChain(AA), DescendingChain(Des)), D)
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
        else
            U = pop!(Asc)
            Ms.x[k] = U
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
        if m > 1 && (length(psd) > 1 && psd[end-1] != psd[end])
            bottom = passthrough!(RF, bottom)  # move to other side
            bottom = passthrough!(D, bottom)
        end
        Ms.x[end] = bottom

    else
        if m > 1 && (length(psd) > 1 && psd[end-1] != psd[end])
            bottom = passthrough!(bottom, D)   # move bottom to other side
            bottom = passthrough!(bottom, RF)
        end
        Mx.x[end] = bottom
    end

    return nothing

end
