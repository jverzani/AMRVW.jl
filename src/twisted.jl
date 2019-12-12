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
## * step_knit! for phatk = :right
## * twisted limbs (need passthrough functions)
## * write some tests
##
## This algorithm encompasses those in CSS and RDS, but those are more efficient space wise
## This is more general and perhaps of interest as building blocks for experimentation

## Structure to hold a twisted represenation
## pv is thee position vector of length  n-1
## m is lowest index of rotators
struct TwistedChain{T} <: AbstractRotatorChain{T}
  x::Vector{T}
  pv::Vector{Symbol}
  m::Int
end

## Constructor
function TwistedChain(xs::Vector{T}) where {T}
    sigma = idx.(xs)
    m = minimum(sigma)
    ps = position_vector(sigma)
    TwistedChain(xs, ps, m)
end

## Constructor of a chain
function Chain(xs::Vector{T}) where  {T}
    if length(xs)  <=  1
        return DescendingChain(xs)
    else
        inds =  idx.(xs)
        pv = position_vector(inds)
        if all(pv .==  :left)
            return DescendingChain(xs)
        elseif all(pv .== :right)
            return AscendingChain(xs)
        else
            return TwistedChain(xs,  pv, minimum(inds))
        end
    end
end


## Find position vector from a permuation of 1...n
## If i is to the left of i+1, set  ps[i] = :left
## if i is to the right if i+1, set ps[i] = :right
## e.g. 4,3,2,1 -> r,r,r (descending)
##      1,2,3,4 -> l,l,l (ascending)
##      1,3,2,4 -> rlr  (CMV)
function position_vector(sigma)
    ## sigma a perm of 1....n
    sigma = sigma .- minimum(sigma) .+ 1
    n = length(sigma)
    n <= 1 && return Symbol[]
    ps =  repeat([:nothing],n-1)
    for (i,v) in enumerate(sigma)
       if v == 1

            if ps[1] == :nothing
              ps[1] = :left
            end
       elseif  v == n
            if ps[end] == :nothing
                ps[end] = :right
            end
       else
          if ps[v-1] == :nothing
                ps[v-1] =  :right
            end
            if ps[v]  == :nothing
                ps[v] = :left
            end
        end
    end
    ps
end

Base.adjoint(A::TwistedChain) = TwistedChain(reverse(adjoint.(A.x)))

## Get i,j entry of twisted chain
## XXX This is not correct XXX
## The general formula is complicated, so we use
## the formula for a descending chain, as our algorithm is set up to turn twisted
## chains into descending chains
function Base.getindex(A::TwistedChain{T}, i, j) where {T}
    DescendingChain(A.x)[i,j]
end



# Fish out of M the rotator with idx i
function iget(Ms, i)
    for (j,M) in enumerate(Ms)
        if idx(M) == i
            return (j, M)
        end
    end
end

# fish out and remove
function iget!(Ms,i)
    j, M = iget(Ms, i)
    deleteat!(Ms, j)
    M
end

## return inds and rotators for rotators with indexes i through j
## function iget(Ms, i, j)
##     if i <= j
##         MMs = iget.(Ref(Ms, ), i:j)
##         [u[1] for u in MMs], [u[2] for u in MMs]
##     else
##         Int[], eltype(MMs)[]
##     end
## end

## get ascending part from Un ... U_{i-1}
function iget(Ms::TwistedChain, i, ::Val{:Asc})
    n = Ms.m
    pv = Ms.pv
    N = n + length(pv)
    inds = Int[]
    Asc = eltype(Ms)[]
    ii = i - 1
    (i <  n || i > N) && return (inds, Asc)
    while true
        if ii >= n && pv[ii - n +  1] == :right
           j,U = iget(Ms, ii)
            push!(inds, j)
            push!(Asc, U)
        else
            break
        end
        ii -= 1
    end

    inds, Asc
end


## get descending part from U_{i+1} ... U_N
function iget(Ms::TwistedChain, i, ::Val{:Des})
    n = Ms.m
    pv = Ms.pv
    N = n + length(pv)
    inds = Int[]
    Asc = eltype(Ms)[]
    (i <  n || i > N) && return (inds, Asc)
    while true
        if i  <= N - 1 && pv[i] == :left
            j,U = iget(Ms, i+1)
            push!(inds, j)
            push!(Asc, U)
        else
            break
        end
        i += 1
    end

    inds, Asc
end



## Need passthroughs
## of a vector through an ascending chain  ->
## of a  vector through a descending chain <-
function passthrough!(L::Vector, A::AscendingChain)
    n = idx(A.x[end])
    N = idx(A.x[1]) # extrema without search
    ##  passing left to right, so we need to go  from  end
    for i in length(L):-1:1
        j = idx(L[i])
        l = N + 1 - j
        if n == N
        else
            A[l-1], A[l], L[i] = turnover(L[i], A[l-1], A[l])
        end
    end
end

function passthrough!(A::DescendingChain, L::Vector)
    length(L) == 0 && return nothing
    n = idx(A.x[1])
    N = idx(A.x[end])

    for i in 1:length(L)
        j  = idx(L[i])
        l  = j - n  +  1
        L[i], A[l],  A[l+1] = turnover(A[l], A[l+1],  L[i])
    end
end


function passthrough!(L::Vector, A::DescendingChain)
    length(L) == 0 && return nothing
    n = idx(A.x[1])
    N = idx(A.x[end])

    for i in length(L):-1:1
        j  = idx(L[i])
        l  = j - n  +  1
        A[l-1],  A[l], L[i] = turnover(L[i], A[l-1], A[l])
    end
end

# need to specilize for diagonal rotator
function turnover(Q1::DiagonalRotator{T},
                  Q2::ComplexRealRotator,
                  Q3::ComplexRealRotator) where {T}

    c1, s1 = vals(Q1); c2, s2 = vals(Q2); c3,s3 = vals(Q3)
    i,j,k = idx(Q1), idx(Q2), idx(Q3)
    # @assert i == k && (abs(j-i) == 1)

    c4,s4,c5,s5,c6,s6 = _turnover(c1,s1, c2,s2, c3,s3)
    R1 = ComplexRealRotator(c4, s4, j)
    R2 = ComplexRealRotator(c5, s5, i)
    R3 = DiagonalRotator(c6, j)

    # we have Q1*Q2*Q3 = R1*R2*R3
    R1, R2, R3

end

# This is used to pass vector of diagonal rotators through the limb
function passthrough!(D::Vector, A::TwistedChain)
    error("Implement me")
end


## always return a diagonal rotator, perhaps an identity one
## this makes the bulge step algorithm generic
function _fuse(U::RealRotator{T}, V::RealRotator{T}) where {T}
    fuse(U,V), IdentityDiagonalRotator{T}()
end
function _fuse(U::ComplexRealRotator{T}, V::ComplexRealRotator{T}) where {T}
    fuse(U,V)
end

# pass a diagonal rotator from a fuse operation through one or more rotator chains, merge with D
# a diagonal matrix
# fused to M, so Ui is really Ui*Di
# needs to move right through Vs
# Set M == nothing if no more interaction with M, as in fusing to descending chain at bottom
passthrough_phase(M, Di::IdentityDiagonalRotator, args...) = nothing

# This case only hits with initial fuse *and* when Ms is descending at top
# so first case is a :left
function passthrough_phase(M::TwistedChain, Di::DiagonalRotator{T}, Vs::Tuple, D) where {T}

    ## We have two rules
    ##
    ##   U1     V1 E1  where if U1=(c,s), D2=(alpha,0), then Vi=(c/alpha,s) Ei=(alpha,0)
    ## D2   -->    E2
    ##   U3     V3 E3
    ## and just
    ##   U1     V1 E1 where V1=(c/alpha,s), Ei = (alpha,0)
    ## D2   -->    E2
    ## Similarly
    ## D2    -->    E2
    ##    U3     V3 E3

    ## Deal with initial Descending part of M, then passthrough

    j = idx(Di) + 1
    _, Uj = iget(M, j)
    DesM = [Uj]
    for pj  in M.pv
        if pj == :right
            break
        else
            j = j + 1
            push!(DesM, iget(M.x,j)[2])
        end
    end


    passthrough_phase(nothing, Di, (DescendingChain(DesM), Vs...), D)
    for U in DesM
        j,V = iget(M.x, idx(U))
        M.x[j] = U
    end

    return nothing
end

passthrough_phase(M::Nothing, Di::DiagonalRotator{T}, Vs::Tuple, D) where {T} = passthrough_phase(Di, Vs, D)

passthrough_phase(Di, Vs::Tuple, D::IdentityDiagonal) = nothing
function passthrough_phase(Di::DiagonalRotator{T}, Vs::Tuple, D) where {T}
    if length(Vs) > 0
        Vhead, Vtail = Vs[1], Vs[2:end]
        passthrough_phase(Di, Vhead, Vtail, D)
    else
        fuse!(Di, D)
    end
    return nothing
end

# peeled one off
function passthrough_phase(Di::DiagonalRotator{T}, V::DescendingChain, Vs::Tuple, D) where {T}

    i = idx(Di)
    alpha, _ = vals(Di)

    if length(V) == 0
        return passthrough_phase(Di, Vs, D)
    end

    n, N = extrema(V)

    if i == n-1

        # Di * Uj = Ujalpha * Di * Dj
        j  =  i  + 1
        Uj = V.x[1]
        @assert idx(Uj) == j
        c, s = vals(Uj)
        Ujalpha = Rotator(c*conj(alpha), s, j)
        Dj = DiagonalRotator(alpha, j)
        V.x[1] = Ujalpha

        ## Di can be passed on to Vs
        passthrough_phase(Di, Vs, D)

        ##  Dj
        ## pass Dj through
        ## we don't have views so we must modify
        VV = V.x[2:end]
        passthrough_phase(Dj,  DescendingChain(VV), Vs, D)
        V.x[2:end] = VV

    elseif i == n

        # D_i(alpha, i) * R(c,s,i) = R(c alpha/conj(alpha), s, i) * D(conj(alpha),i)
        # fuse with U leaves a new Di, called Dii  below

        U = V.x[1]
        c,s = vals(U)
        Ualpha = Rotator(c*alpha/conj(alpha), s, i)
        V.x[1] = Ualpha
        Dii = DiagonalRotator(conj(alpha), i)

        VV  = V.x[2:end]
        passthrough_phase(Dii, DescendingChain(VV), Vs, D)
        V.x[2:end] = VV

    elseif n < i <= N

        # turnover
        V.x[i-n], V.x[i-n+1], Di = turnover(Di, V.x[i-n], V.x[i-n+1])
        passthrough_phase(Di, Vs, D)

    elseif i == N+1

        U = V.x[end]
        @assert idx(U) == N
        c,s = vals(U)
        Ualpha = Rotator(c*conj(alpha), s, N)
        Dh = DiagonalRotator(alpha, N)
        V.x[end] = Ualpha

        passthrough_phase(Dh, Vs, D)
        passthrough_phase(Di, Vs, D)

    else

        passthrough_phase(Di, Vs, D)

    end


end

function passthrough_phase(Di::DiagonalRotator{T}, V::AscendingChain, Vs::Tuple, D) where {T}

    i = idx(Di)
    alpha, _ = vals(Di)

    if length(V) == 0
        return passthrough_phase(Di, Vs, D)
    end

    n, N = extrema(V)

    if i == N +  1

        # Di * Uj = Ujalpha * Di * Dj
        j  =  i  - 1
        Uj = V.x[1]
        c, s = vals(Uj)
        Ujalpha = Rotator(c*conj(alpha), s, j)
        Dj = DiagonalRotator(alpha, j)
        V.x[1] = Ujalpha

        ## Di
        passthrough_phase(Di, Vs, D)

        ## pass Dj through rest of V
        VV = V.x[2:end]
        passthrough_phase(Dj,  AscendingChain(VV), Vs, D)
        V.x[2:end] = VV

    elseif i == N

        # D_i(alpha, i) * R(c,s,i) = R(c alpha/conj(alpha), s, i) * D(conj(alpha),i)
        # fuse with U
        U = V.x[1]
        c,s = vals(U)
        Ualpha = Rotator(c*alpha/conj(alpha), s, i)
        V.x[1] = Ualpha
        Dii = DiagonalRotator(conj(alpha), i)

        VV  = V.x[2:end]
        passthrough_phase(Dii, AscendingChain(VV), Vs, D)
        V.x[2:end] = VV

    elseif n <= i < N

        # turnover
        V.x[N+1-i-1], V.x[N+1-i], Di = turnover(Di, V.x[N+1-i-1], V.x[N+1-i])
        passthrough_phase(Di, Vs, D)

    elseif i == n-1

        U = V.x[end]
        c,s = vals(U)
        Ualpha = Rotator(c*conj(alpha), s, n)
        Dh = DiagonalRotator(alpha, i+1)
        V.x[end] = Ualpha

        passthrough_phase(Dh, Vs, D)
        passthrough_phase(Di, Vs, D)

    else
        passthrough_phase(Di, Vs, D)
    end


end


function passthrough_phase(Di::DiagonalRotator{T}, V::TwistedChain, Vs::Tuple, D) where {T}
    i = idx(Di)
    alpha, _ = vals(Di)
    n = V.m
    N = n + length(V.pv)
    ps = V

    if i < n - 1
        passthrough_phase(Di, Vs, D)
    elseif i == n - 1
        idx, Ui = iget(V, i)
        Uialpha = Rotator(c*conj(alpha), s, i)
        V[idx] = Uialpha
        Dj = DiagonalRotator(alpha, i+1)
        passthrough_phase(Di, Vs, D)

        ## bottom
        inds, MMs = iget(Ms, i, Val{:Des})
        passthrough_phase(Dj, TwistedChain(MMs), Vs, D)
        for (i,j) in enumerate(inds)
            V[j] = MMs[i]
        end

    elseif i == n
        if ps[1] == :right
            #turnover, move on with bottom
        else
            #fuse, move on with bottom
        end
    elseif i == N
        if ps[N-1] == :left
            #turnover, move on
        else
            #fuse, move on with top
        end
    elseif i == N + 1
        idx, Ui = iget(V, i)
        Uialpha = Rotator(c*conj(alpha), s, i)
        V[idx] = Uialpha
        Dh = DiagonalRotator(alpha, i-1)
        passthrough_phase(Di, Vs, D)

        ## top
        inds, MMs = iget((Ms,),i, Val(:Asc))
        passthrough_phase(Dj, TwistedChain(MMs), Vs, D)
        for (i,j) in enumerate(inds)
            V[j] = MMs[i]
        end

    elseif i > N + 1
        passthrough_phase!(Di, Vs, D)
    else
        ## in between. Depends on pattern...
        pattern = (ps[i-n], ps[i-n+1])
        if pattern == (:left, :left)
            ##         ?
            ##    ⌈    ?
            ## D  ⌊  ⌈
            ## D     ⌊ ⌈
            ##         ⌊
            ##
            #turnover, move on
            j1,U1 = iget(Vs, i-1)
            j2,U2 = iget(Vs, i)
            U1, U2,  Di = turnover(Di, U1,  U2)
            Vs[j1] = U1
            Vs[j2] = U2
            inds, MMs = iget(Vs, i-1)
            passthrough_phase!(Di, TwistedChain(MMs), Vs, D)
            for  (i,j) in enumerate(inds)
                Vs[j] = MMs[i]
            end
        elseif pattern == (:right, :right)
            ##         ⌈
            ## D    ⌈  ⌊
            ## D ⌈  ⌊
            ##   ⌊     ?
            ##         ?
            ##  turnover, move on
            j1,U1 = iget(Vs, i)
            j2,U2 = iget(Vs, i+1)
            U2, U1,  Di = turnover(Di, U2,  U1)
            Vs[j1] = U1
            Vs[j2] = U2
            inds, MMs = iget(Vs, i+1, Val(:Des))
            passthrough_phase!(Di, TwistedChain(MMs), Vs, D)
            for  (i,j) in enumerate(inds)
                Vs[j] = MMs[i]
            end
        elseif pattern == (:right, :left)
            ##     ⌈
            ## D ⌈ ⌊
            ## D ⌊ ⌈
            ##     ⌊
            ##
            ## passing through gives Di, Dh, and Dj

            j, Ui = iget(Vs,  i)
            c, s = vals(Ui)
            Uialpha = Rotator(c * alpha / conj(alpha), s, i)
            Vs[j] = Uialpha
            Dh, Dj = DiagonalRotator.(conj(alpha), (i-1, i+1))
            Dii = DiagonalRotator(conj(2alpha), i)
            # passthrough_top and bottttom
            inds, MMs = iget(Vs, i-1, Val(:Asc))
            passthrough_phase!(Di, AscendingChain(MMs), Vs, D)
            for  (i,j) in enumerate(inds)
                Vs[j] = MMs[i]
            end

            passthrough_phase!(Dii, Vs, D)

            inds, MMs = iget(Vs, i+1, Val(:Des))
            passthrough_phase!(Di, DescendingChain(MMs), Vs, D)
            for  (i,j) in enumerate(inds)
                Vs[j] = MMs[i]
            end

        elseif pattern == (:left, :right)
            ##        ?
            ##    ⌈   ?
            ## D  ⌊ ⌈
            ## D  ⌈ ⌊
            ##    ⌊   ?
            ##        ?
            Uh = iget(Vs, i-1)
            c, s = vals(Uh)
            Uhalpha = Rotator(c*conj(alpha), s, i-1)

            Ui = iget(Vs, i)
            c, s = vals(Ui)
            Uialpha = Rotator(c*conj(2*alpha), s, i)  # note this  gets hit twice

            Uj = iget(Vs, i+1)
            c, s = vals(Uj)
            Ujalpha = Rotator(c*conj(alpha), s, i+1)

            Dg, Dh, Di, Dj, Dk  = DiagonalRotator.(alpha, i-2:i+2)

            #... realease the hounds...
            inds, MMs = iget(Vs, i-2, Val(:Asc))
            passthrough_phase!(Dg, AscendingChain(MMs), Vs, D)
            for  (i,j) in enumerate(inds)
                Vs[j] = MMs[i]
            end
            passthrough_phase!(Dh, AscendingChain(MMs), Vs, D)
            for  (i,j) in enumerate(inds)
                Vs[j] = MMs[i]
            end
            passthrough_phase!(Di, Vs, D)
            inds, MMs = iget(Vs, i+2, Val(:Des))
            passthrough_phase!(Dj, DescendingChain(MMs), Vs, D)
            for  (i,j) in enumerate(inds)
                Vs[j] = MMs[i]
            end
            passthrough_phase!(Dk, DescendingChain(MMs), Vs, D)
            for  (i,j) in enumerate(inds)
                Vs[j] = MMs[i]
            end
        end
    end






end

##################################################


##
## Twisted refers to the QFactorization
##
struct QFactorizationTwisted{T, Rt} <: AbstractQFactorization{T, Rt, Val{:twisted}}
  Q::TwistedChain{Rt}
  D::AbstractSparseDiagonalMatrix{T, Rt}
end


## return q[i:k, j:k]
function Base.getindex(QF::QFactorizationTwisted, j, k)
    ## This is wrong as QF.Q[j,k] is not correct
    QF.Q[j,k] * (k == 0 ? 0 : QF.D[k])
end

##################################################
# Twisting *would* require a different bulge chasing algorithm, so
# we hold it in the type for dispatch
struct QRFactorizationTwisted{T, S, Rt, QFt, RFt} <: AbstractFactorizationState{T, S, Rt, QFt, RFt, Val{:twisted}}
  N::Int
  QF::QFt
  RF::RFt
  UV::Vector{Rt}    # temp storage for U[V] Vector{Rt}, ,
  A::Matrix{S}
  REIGS::Vector{T}
  IEIGS::Vector{T}
  ctrs::AMRVW_Counter
end


function eigvals(state::QRFactorizationTwisted)

    new_state = AMRVW_algorithm(state)
    complex.(new_state.REIGS, new_state.IEIGS)

end

## We don't go for the fastest method
## We have the first n-m  steps through, a straightening out of the twisted factorization in QF,
## This should be  O(n^2) steps
## Then, once untwisted, we pass down to the more efficient O(n^2) algorithm (in time) to solve
function AMRVW_algorithm(state::QRFactorizationTwisted{T, Rt}) where  {T, Rt}

    Ms = state.QF.Q.x # the vector
    n = length(Ms)
    stop_ctr = n
    m = isa(Rt, RealRotator{T})  ? 2 : 1
    sigma = [idx(R) for R in Ms]
    ps = position_vector(sigma)
    N0 = findlast(x->x==:right, ps)
    N = N0 == nothing ? 1 : N0
    for _ in 1:(N+2-m)
        bulge_step(state)

        # check for deflation of last one
        c,s = vals(Ms[state.ctrs.stop_index])
        if abs(s) <= eps(T)
            state.ctrs.stop_index -= 1
        end
    end

    # now untwisted, so we can change over
    QF = state.QF
    Ms = QF.Q.x
    QF_new = QFactorization(DescendingChain(Ms), QF.D, Vector{eltype(Ms)}(undef, 1))

    new_state = qrfactorization(state.N, QF_new, state.RF)
    AMRVW_algorithm(new_state)
    new_state

end





##  We  break bulge step into two cases
## * for real rotators, we don't worry about the fuse operation, so leave in the generality of the paper (m >=1)
## * for complex/real rotator, we use the m=1 case

## Real Rotator Case

function bulge_step(state::AbstractFactorizationState{T, S, Rt, QFt, RFt, Val{:twisted}})  where {T,S, Rt, QFt, RFt}

    # create bulge is only right when  all rotators in QFt are straightened out,
    # but  this should step in the correct direction.
    create_bulge(state)
    Asc = reverse(state.UV)  # stored as U, V
    Ms = state.QF.Q.x # the vector part without the TwistedChain wrapper
    RF = state.RF
    D = state.QF.D
    # modify Ms, D, RF
    bulge_step!(Ms, D, RF, Asc)

end

## Can enter here for more general usage
function bulge_step!(Ms, D, RF, Asc)

    Rt = eltype(Ms)
    Decoupled = Rt[]

    Des = reverse(adjoint.(Asc))
    n = length(Ms) + 1
    m = length(Asc)


    sigma = [idx(R) for R in Ms]
    ps = position_vector(sigma)


    ## Must put AscendingChain into place
    ## py passing through RF and D <--
    passthrough!(RF, Asc)
    passthrough!(D, AscendingChain(Asc))

    limb, limb_side = step_0!(ps, m, Ms, Des, Asc, D)

    for k in 1:(n-m-2)

        limb_side = step_k!(k, n, m, ps, limb_side, limb, Des, Asc, Decoupled, Ms, D, RF)

    end


    # now we choose :left
    limb_side = step_k!(n-m-1, n, m, :left, limb_side, limb, Des, Asc, Decoupled, Ms, D, RF)

    # now knit in  limb, Des, Asc
    step_knit!(n, m, :left, limb_side, limb, Des,  Asc, Decoupled, Ms, D, RF)

    ## Decoupled is now Ms, assign
    append!(Ms, Decoupled)

    return nothing
end



function step_0!(ps, m,  Ms, Des, Asc, D) where {T}


    has_limb = ifelse(m > 1, true, false)
    limb = has_limb ? iget!.(Ref(Ms,), (1:m-1)) : eltype(Ms)[]

    limb_side = :nothing

    if has_limb
        if ps[m-1] == :left
            passthrough!(DescendingChain(Des), limb)  # pass limb through left
            limb_side = :left
        else
            passthrough!(limb, AscendingChain(Asc))    # pass limb to right
            limb_side = :right
        end
    end

    U = iget!(Ms, m)

    if ps[m] == :left

        Des[end], Di = _fuse(Des[end], U)  # fuse with descending; aka Ms[m]

        ## need to passthrough Ms too now
        if limb_side == :right
            passthrough_phase(TwistedChain(Ms), Di, (AscendingChain(Asc), Chain(limb)), D) ### need to pass through limb as well
        else
            passthrough_phase(TwistedChain(Ms), Di, (AscendingChain(Asc), ), D)
        end


    elseif ps[m] == :right

        AA = popfirst!(Asc)
        AA, Di = _fuse(U, AA)

        if limb_side == :right
            passthrough_phase(Di, (AscendingChain(Asc), Chain(limb)), D)
        else
            passthrough_phase(Di, (AscendingChain(Asc), ), D)
        end
        pushfirst!(Asc, AA)


    end



    limb, limb_side



end

function step_k!( k, n, m, ps, limb_side, limb, Des, Asc, Decoupled, Ms, D, RF) where {T}

    phatk = isa(ps, Symbol) ? ps :  ps[k]

    if phatk == :left

        U = iget!(Ms, m + k )
        push!(Des, U)  ## Augment Descending
        passthrough!(DescendingChain(Des), Asc) ## translate ascending

    else

        U = iget!(Ms, m + k)
        pushfirst!(Asc, U)
        passthrough!(Des, AscendingChain(Asc))

    end


    if limb_side == :left
        passthrough!(limb, AscendingChain(Asc))
    else
        passthrough!(DescendingChain(Des), limb)
    end


    ## similarity tranform
    if phatk == :left

        ## similarity transform Asc to right side
        passthrough!(RF, Asc) # <-- pass Asc through RF
        passthrough!(D, AscendingChain(Asc))
        limb_side = :left

    else

        ## similarity transform descending to left sider
        passthrough!(DescendingChain(Des), D)
        passthrough!(Des, RF) #--> pass Des through Rf;  similarity transform
        limb_side = :right

    end

    ## Decouple
    if phatk == :left

        U = popfirst!(Des)
        push!(Decoupled, U)  # position consistent with phatk

    else

        U = pop!(Asc)
        pushfirst!(Decoupled, U) # posoition consistent with phatk

    end



    return limb_side
end



function step_knit!( n, m, phatk::Symbol, limb_side, limb, Des,  Asc, Decoupled, Ms, D, RF) where {T}

    # make bottom
    U = pop!(Des)
    V = popfirst!(Asc)
    bottom, Di = _fuse(U,V)

    if limb_side == :right
        passthrough_phase(nothing, Di, (AscendingChain(Asc), Chain(limb)), D)
    else
        passthrough_phase(nothing, Di, (AscendingChain(Asc),), D)
    end


    ## Final Steps, knit in
    ## setps k=n-m to  n-2


    for k in (n-m):(n-2)

        if phatk == :left
            Des = push!(Des, bottom)
            passthrough!(DescendingChain(Des), Asc)
        else
            Asc = pushfirst!(Asc, bottom)
            passthrough(Des, AscendingChain(Asc))
        end

        if length(limb) > 0
            lps = position_vector(sortperm([idx(l) for l in limb]))

            L =  pop!(limb)
            if limb_side == :left
                if length(lps) > 0
                    lpk = pop!(lps)
                else
                    lpk = :nothing
                end
                if lpk == :right
                    # translate then fuse
                    passthrough!(limb, AscendingChain(Asc))
                    Asc[1], Di = _fuse(L, Asc[1])
                    if length(Asc) > 1
                        AA = Asc[2:end]
                        passthrough_phase(nothing, Di, (AscendingChain(AA),), D)  # XXX And something prior?
                        Asc[2:end] = AA
                    end
                else
                    # fuse then translate
                    Asc[1], Di = _fuse(L, Asc[1])
                    if length(Asc) > 1
                        AA = Asc[2:end]
                        passthrough_phase(nothing, Di, (AscendingChain(AA),), D)  # XXX and something prior?
                        Asc[2:end] = AA
                        passthrough!(limb, AscendingChain(Asc))
                    end
                end
            else
                if length(lps) >  0
                    lpk = pop!(lps)
                else
                    lpk = :nothing
                end
                if  lpk == :left
                    # translate then fuse
                    passthrough(DescendingChain(Des), limb)
                    Des[end], Di = _fuse(Des[end], L)
                else
                    Des[end], Di = _fuse(Des[end], L)
                    passthrough(DescendingChain(Des), limb)
                end
            end
        end


        if phatk == :left
            # similarity to move Asc to right side
            passthrough!(RF, Asc)  # <--
            passthrough!(D, AscendingChain(Asc))
            push!(Decoupled,  popfirst!(Des))  # take Des and  move top to decoupled
            limb_side = :left

        else
            passthrough!(DescendingChain(Des), D)
            passthrough!(Des, RF)
            push!(Decoupled, pop!(Asc))
            limb_side = :right
        end

        U = pop!(Des)
        V = popfirst!(Asc)
        bottom, Di = _fuse(U,V)

        if limb_side == :right
            passthrough_phase(nothing, Di, (AscendingChain(Asc), Chain(limb)), D)
        else
            passthrough_phase(nothing, Di, (AscendingChain(Asc),), D)
        end

    end

    push!(Decoupled, bottom)

    return nothing

end


## For fun....
## This should be matrix multiplication
## as Twiste Q is no long Hessenberg
function diagonal_block(state::QRFactorizationTwisted{T, S, Rt, QFt, RFt}, k) where {T,  S, Rt, QFt,  RFt}

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
