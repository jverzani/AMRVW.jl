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
    n = Ms.m[]
    pv = Ms.pv
    N = n + length(pv)


    inds = Int[]
    Asc = eltype(Ms)[]
    ii = i - 1
    (i <  n || i > N + 1) && return (inds, Asc)
    while true
        if ii >= n && (ii == N || pv[ii - n +  1] == :right)
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
    n = Ms.m[]
    pv = Ms.pv
    N = n + length(pv)


    inds = Int[]
    Asc = eltype(Ms)[]
    (i <  n - 1 || i >= N) && return (inds, Asc)
    while true
        if i  < N && (i < n || pv[i-n+1] == :left)
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


## Move twisted through down R <- L
function  passthrough!(A::Union{DescendingChain, AscendingChain}, B::TwistedChain)
    length(B) == 0 && return nothing
    ##  startt at bottom of B
    ## must calibratet A and B
    n = B.m[]
    N = idx(A[1])

    pv = B.pv
    inds = [iget(B, l)[1] for  l  in n:n+length(pv)]  # Can use j's from iget, as we  modify B  along the way

    for (i, pi) in enumerate(reverse(pv)) # work from bottom

        #i = 2 -> length(pv) - 2 + 1 idx in pv
        #n + length(pv) - i + 1 idx in A of lower one
        #j = inds[length(pv) - 1 + 1] index of lower one

        if pi  == :right
            j = inds[length(pv) - i + 1 + 1]
            B[j] = passthrough!(A, B[j])
        end
    end

    # move first over
    j =  inds[1]
    B[j] = passthrough!(A, B[j])

    for (i, pi) in enumerate(pv) # work down
        if pi == :left
            j = inds[i + 1]
            B[j] = passthrough!(A, B[j])
        end
    end

    if isa(A, DescendingChain)
        B.m[] += 1
    else
        B.m[] -= 1
    end

end

## Move twisted through  up R --> L
function  passthrough!(B::TwistedChain, A::Union{DescendingChain,AscendingChain})
    length(B) == 0 && return

    n = B.m[]
    N = idx(A.x[end])
    pv = B.pv

    inds = [iget(B, l)[1] for  l  in n:n+length(pv)]  # Can't use j's from iget, as we  modify B  along the way

    for (i, pi) in enumerate(reverse(pv)) # work from bottom

        if pi  == :left
            j = inds[length(pv) - i + 1 + 1]
            U = B[j]
            B[j] = passthrough!(U, A)
        end
    end

    # move first over
    j =  inds[1]
    B[j] = passthrough!(B[j], A)

    for (i, pi) in enumerate(pv) # work down
        if pi == :right
            j = inds[i+1]
            B[j] = passthrough!(B[j], A)
        end
    end


    if isa(A, DescendingChain)
        B.m[] -= 1
    else
        B.m[] += 1
    end
end


# This is used to pass vector of diagonal rotators through the limb
function passthrough!(D::Vector, A::TwistedChain)
    error("Implement me")
end


## always return a diagonal rotator, perhaps an identity one
## this makes the bulge step algorithm generic
function _fuse(U::RealRotator{T}, V::RealRotator{T}) where {T}
    D = IdentityDiagonalRotator{T}(idx(V))
    fuse(U,V), D
end

function _fuse(U::ComplexRealRotator{T}, V::ComplexRealRotator{T}) where {T}
    fuse(U,V)
end


# pass a diagonal rotator from a fuse operation through one or more rotator chains
# merge with D a diagonal matrix
passthrough_phase!(Di, Vs::Tuple, D::IdentityDiagonal) = nothing
passthrough_phase!(Di, C::AbstractRotatorChain, Vs::Tuple, D::IdentityDiagonal) = nothing

function passthrough_phase!(Di::DiagonalRotator{T}, Vs::Tuple, D) where {T}
    if length(Vs) > 0
        Vhead, Vtail = Vs[1], Vs[2:end]
        passthrough_phase!(Di, Vhead, Vtail, D)
    else
        fuse!(Di, D)
    end
    return nothing
end

# peeled one off
function passthrough_phase!(Di::DiagonalRotator{T}, V::DescendingChain, Vs::Tuple, D) where {T}

    i = idx(Di)
    alpha, _ = vals(Di)

    if length(V) == 0
        return passthrough_phase!(Di, Vs, D)
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
        passthrough_phase!(Di, Vs, D)

        ##  Dj
        ## pass Dj through
        ## we don't have views so we must modify
        VV = V.x[2:end]
        passthrough_phase!(Dj,  DescendingChain(VV), Vs, D)
        V.x[2:end] = VV

    elseif i == n

        V.x[1], Dii = passthrough!(Di, V.x[1])
        VV  = V.x[2:end]
        passthrough_phase!(Dii, DescendingChain(VV), Vs, D)
        V.x[2:end] = VV

    elseif n < i <= N

        # turnover
        V.x[i-n], V.x[i-n+1], Di = turnover(Di, V.x[i-n], V.x[i-n+1])
        passthrough_phase!(Di, Vs, D)

    elseif i == N+1

        U = V.x[end]
        @assert idx(U) == N
        c,s = vals(U)
        Ualpha = Rotator(c*conj(alpha), s, N)
        Dh = DiagonalRotator(alpha, N)
        V.x[end] = Ualpha

        passthrough_phase!(Dh, Vs, D)
        passthrough_phase!(Di, Vs, D)

    else

        passthrough_phase!(Di, Vs, D)

    end


end



function passthrough_phase!(Di::DiagonalRotator{T}, V::AscendingChain, Vs::Tuple, D) where {T}

    i = idx(Di)
    alpha, _ = vals(Di)

    if length(V) == 0
        return passthrough_phase!(Di, Vs, D)
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
        passthrough_phase!(Di, Vs, D)

        ## pass Dj through rest of V
        VV = V.x[2:end]
        passthrough_phase!(Dj,  AscendingChain(VV), Vs, D)
        V.x[2:end] = VV

    elseif i == N

        V.x[1], Dii = passthrough!(Di, V.x[1])

        VV  = V.x[2:end]
        passthrough_phase!(Dii, AscendingChain(VV), Vs, D)
        V.x[2:end] = VV

    elseif n <= i < N

        # turnover
        V.x[N+1-i-1], V.x[N+1-i], Di = turnover(Di, V.x[N+1-i-1], V.x[N+1-i])
        passthrough_phase!(Di, Vs, D)

    elseif i == n-1

        U = V.x[end]
        c,s = vals(U)
        Ualpha = Rotator(c*conj(alpha), s, n)
        Dh = DiagonalRotator(alpha, i+1)
        V.x[end] = Ualpha

        passthrough_phase!(Dh, Vs, D)
        passthrough_phase!(Di, Vs, D)

    else
        passthrough_phase!(Di, Vs, D)
    end


end

# two special functions useful to simplify the below:
function passthrough_phase!(Di::DiagonalRotator{T}, V::TwistedChain, Vs::Tuple, D, j, ::Val{:Des}) where {T}
    inds, Des = iget(V, j, Val(:Des))
    passthrough_phase!(Di, DescendingChain(Des), Vs, D)
    for (i,j) in enumerate(inds)
        V[j] = Des[i]
    end
end

function passthrough_phase!(Di::DiagonalRotator{T}, V::TwistedChain, Vs::Tuple, D, j, ::Val{:Asc}) where {T}
    inds, Asc = iget(V, j, Val(:Asc))
    passthrough_phase!(Di, AscendingChain(Asc), Vs, D)
    for (i,j) in enumerate(inds)
        V[j] = Asc[i]
    end
end

function passthrough_phase!(Di::DiagonalRotator{T}, V::TwistedChain, Vs::Tuple, D) where {T}

    length(V.x) <= 2 && return passthrough_phase!(Di, Chain(V.x), Vs, D)

    n = V.m[]
    pv = V.pv
    N = n + length(pv)
    i = idx(Di)  # pv[i-n+1] informs if Ui+1 is left or right of Ui
    if i < n - 1
        passthrough_phase!(Di, Vs, D) # misses chain
    elseif i == n - 1
        passthrough_phase!(Di, V, Vs, D,  n-1, Val(:Des))
    elseif i == n
        if pv[1] == :left
            ind, Ui = iget(V, i)
            Ui, Di = passthrough!(Di, Ui)
            V[ind] = Ui
            passthrough_phase!(Di, V, Vs, D, i, Val(:Des))
        else ## right, so turnover
            if pv[2] == :right
                ii, Ui = iget(V, i)
                ij, Uj = iget(V, i+1)
                Uj, Ui, Dj = turnover(Di, Uj, Ui)
                V[ii] = Ui
                V[ij] = Uj
                passthrough_phase!(Dj, Vs, D)
            else
                # D1, 2,1,3 pattern
                ii,  Ui  =  iget(V,  i)
                ij,  Uj =  iget(V, i+1)
                Uj, Ui, Dj = turnover(Di, Uj, Ui)
                V[ii] = Ui
                V[ij] = Uj
                passthrough_phase!(Dj, V, Vs,  D,  i+1, Val(:Des))
            end
        end
    elseif i == N
        if pv[end] == :right
            # flip, get Ascending
            ind, Ui = iget(V, i)
            Ui, Di = passthrough!(Di, Ui)
            V[ind] = Ui
            passthrough_phase!(Di, V, Vs, D,   i, Val(:Asc))
        elseif pv[end] == :left && pv[end-1] == :left
            ii, Ui = iget(V, i)
            ih, Uh = iget(V, i-1)
            Uh, Ui, Dh = turnover(Di, Uh, Ui)
            V[ii] = Ui
            V[ih] = Uh
            passthrough_phase!(Dh, Vs, D)
        else
            # D3, 2,1,3 pattern
            ih, Uh =  iget(V, i-1)
            ii, Ui = iget(V, i)
            Uh, Ui,  Dh =  turnover(Di, Uh, Ui)
            V[ii] =  Ui
            V[ih] = Uh
            passthrough_phase!(Dh, V, Vs,  D,  i-1, Val(:Asc))
        end
    elseif i == N + 1
        passthrough_phase!(Di, V, Vs, D, N+1, Val(:Asc))
    elseif i > N + 1
        passthrough_phase!(Di, Vs, D)
    else # n < i < N
        ## 4  cases
        if  (pv[i-n+1] == :left) &&  (pv[i-n] == :left)
            # [
            #   [
            #     [.
            ih, Uh = iget(V, i-1)
            ii, Ui = iget(V, i)
            Uh, Ui, Dh = turnover(Di, Uh, Ui)
            V[ih] = Uh
            V[ii] = Ui
            # We might have  to pass  through i-2, but *only* if Dh is on the right
            if i-n-1 > 0 && pv[i-n-1] == :right
                passthrough_phase!(Dh, V, Vs, D, i-1, Val(:Asc))
            else
                passthrough_phase!(Dh, Vs, D)
            end

        elseif (pv[i-n] == :right) && (pv[i-n+1] == :right)
            #     [.
            #   [
            # [
            # turnover
            ii, Ui = iget(V, i)
            ij, Uj = iget(V, i + 1)
            Uj, Ui, Dj = turnover(Di, Uj, Ui)
            V[ii] = Ui
            V[ij] = Uj
            if i-n+2 <= length(pv)  && pv[i-n+2] == :left
                passthrough_phase!(Dj, V, Vs, D, i+1, Val(:Des))
            else
                passthrough_phase!(Dj, Vs, D)
            end
        elseif pv[i-n] == :right && pv[i-n+1] == :left
            #    [         [          [ Dh        [ Dj
            # Di[ --   [ Di  --> [ Di   Di  --> [   Di
            #    [         [               [      [ Dj
            #
            ih, Uh = iget(V, i-1)
            ii, Ui = iget(V, i)
            ij, Uj = iget(V, i+1)
            Ui, Di = passthrough!(Di, Ui)
            Uh, Dh, = tip(Di, Uh)
            Uj, Dj = tip(Di, Uj)

            V[ih] = Uh
            V[ii] = Ui
            V[ij] = Uj
            ## pass Di through
            passthrough_phase!(Di, Vs, D)
            ## Move Dh pass  Asc, if present
            if i-n-1 > 0 && pv[i-n-1] == :right
                passthrough_phase!(Dh, V, Vs, D, i-1, Val(:Asc))
            else
                passthrough_phase!(Dh, Vs, D)
            end

            ## Move Dj passed Des, if present
            if i-n + 2 <= length(pv)  && pv[i-n+2] == :left
                passthrough_phase!(Dj, V, Vs, D, i+1, Val(:Des))
            else
                passthrough_phase!(Dj, Vs, D)
            end


        elseif  pv[i-n] == :left && pv[i-n+1] == :right
            #  [     [ Dh       [ Dh       [   Dh  if rig
            # D [ ->   Di  [ ->      [  ->   [ Di
            #  [         [          [ Dj   [   Dj
            #
            ih, Uh = iget(V, i-1)
            ii, Ui = iget(V, i)
            ij, Uj = iget(V, i+1)
            ##  pass through Dh:
            Uh, Dh = tip(Di, Uh)
            ## turnover
            Uj, Ui, Dj = turnover(Di, Uj, Ui)
            ## tip with Dh
            Ui, Di = tip(Dh, Ui)
            V[ih] = Uh
            V[ii] = Ui
            V[ij] = Uj
            ## pass Di through
            passthrough_phase!(Di, Vs, D)
            ## Move Dh pass  Asc, if present
            if i-n-1 > 0 && pv[i-n-1] == :right
                passthrough_phase!(Dh, V, Vs, D, i-1, Val(:Asc))
            else
                passthrough_phase!(Dh, Vs, D)
            end

            ## Move Dj passed Des, if present

            if i-n+2 <= length(pv)  && pv[i-n+2] == :left
                passthrough_phase!(Dj, V, Vs, D, i+1, Val(:Des))
            else
                passthrough_phase!(Dj, Vs, D)
            end
        end

    end

    return nothing
end


##################################################


##
## Twisted refers to the QFactorization
##
struct QFactorizationTwisted{T, Rt} <: AbstractQFactorization{T, Rt, Val{:twisted}}
  Q::TwistedChain{Rt}
  D::AbstractSparseDiagonalMatrix{T, Rt}
end

Base.eltype(QF::QFactorizationTwisted) = eltype(QF.Q.x)
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


    sigma = idx.(Ms)
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

# limb is first m-1 rotators keeping their positions
function _get_limb!(Ms::TwistedChain, m::Int)

    if m > 1
        pv = Ms.pv #position_vector(idx.(Ms))
        Us =  iget!.(Ref(Ms,), (1:m-1))
        limb = [popfirst!(Us)]
        while length(Us) > 0
            U = popfirst!(Us)
            dir = popfirst!(pv)
            dir == :left ? push!(limb, U) : pushfirst!(limb, U)
        end
        Ms.m[] = Ms.m[] + (m-1)
        popfirst!(Ms.pv)
    else
        limb = eltype(Ms)[]
    end

    TwistedChain(limb)
end

function step_0!(ps, m,  Ms, Des, Asc, D) where {T}

    has_limb = ifelse(m > 1, true, false)

    # grab limb, but keep order
    MMs = TwistedChain(Ms)
    limb = _get_limb!(MMs, m)

    # TODO: wrap limb=TwistedChain(...), as o/w position vecgtor is computed each time

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

    U = iget!(MMs, m)
    if ps[m] == :left

        Des[end], Di = _fuse(Des[end], U)  # fuse with descending; aka Ms[m]

        ## need to passthrough Ms too now
        i = idx(Di)
        #        inds, MMs = iget(TwistedChain(Ms), i+1, Val(:Des))
        inds, _MMs = iget(MMs, i+1, Val(:Des))
        if limb_side == :right
            passthrough_phase!(Di, DescendingChain(_MMs),  (AscendingChain(Asc), limb), D)
        else
            passthrough_phase!(Di, DescendingChain(_MMs),  (AscendingChain(Asc), ), D)
        end
        for (i,j) in enumerate(inds)
            MMs[j] = _MMs[i]
        end


    else

        AA = popfirst!(Asc)
        AA, Di = _fuse(U, AA)

        if limb_side == :right
            passthrough_phase!(Di, (AscendingChain(Asc), limb), D)
        else
            passthrough_phase!(Di, (AscendingChain(Asc), ), D)
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
        passthrough!(DescendingChain(Des), AscendingChain(Asc)) ## translate ascending

    else

        U = iget!(Ms, m + k)
        pushfirst!(Asc, U)
        passthrough!(DescendingChain(Des), AscendingChain(Asc))

    end

    if limb_side == :right
        passthrough!(DescendingChain(Des), limb)
    elseif limb_side == :left
        passthrough!(limb, AscendingChain(Asc))
    else
        nothing
    end


    ## similarity tranform
    if phatk == :left

        ## similarity transform Asc to right side
        passthrough!(RF, Asc) # <-- pass Asc through RF
        passthrough!(D, AscendingChain(Asc))
        if m > 1
            limb_side = :left
        end

    else

        ## similarity transform descending to left sider
        passthrough!(DescendingChain(Des), D)
        passthrough!(Des, RF) #--> pass Des through Rf;  similarity transform
        if m > 1
            limb_side = :right
        end

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



function step_knit!(n, m, phatk::Symbol, limb_side, limb, Des,  Asc, Decoupled, Ms, D, RF) where {T}
    # make bottom
    U = pop!(Des)
    V = popfirst!(Asc)
    bottom, Di = _fuse(U,V)

    if limb_side == :right
        passthrough_phase!(Di, (AscendingChain(Asc), TwistedChain(limb)), D)
    else
        passthrough_phase!(Di, (AscendingChain(Asc),), D)
    end

    ## Final Steps, knit in
    ## setps k=n-m to  n-2

    for k in (n-m):(n-2)

        if phatk == :left
            Des = push!(Des, bottom)
            passthrough!(DescendingChain(Des), AscendingChain(Asc))
        else
            Asc = pushfirst!(Asc, bottom)
            passthrough!(DescendingChain(Des), AscendingChain(Asc))
        end

        if length(limb) > 0

            lpk = length(limb.pv) > 0 ? limb.pv[end] : :nothing
            L =  pop!(limb)

            ## we have
            ## 4 cases in terms of limb_side and lpk:
            ## :l, :l -> fuse (Asc, Des), translate
            ## :l, :r -> translate, fuse (Asc, limb, Des)
            ## :r, :l -> translate, fuse (limb)
            ## :r, :r -> fuse (.), translate
            if limb_side == :left
                if lpk == :left || lpk == :nothing
                    # fuse then translate

                    Asc[1], Di = _fuse(L, Asc[1])
                    if length(Asc) > 1
                        AA = Asc[2:end]
                        passthrough_phase!(Di, (AscendingChain(AA), DescendingChain(Des)), D)  # XXX and something prior?
                        Asc[2:end] = AA
                    else
                        passthrough_phase!(Di, (DescendingChain(Des),), D)
                    end
                    ## translate
                    passthrough!(limb, AscendingChain(Asc))

                elseif lpk == :right
                    # translate then fuse
                    passthrough!(limb, AscendingChain(Asc))
                    Asc[1], Di = _fuse(L, Asc[1])

                    if length(Asc) > 1
                        AA = Asc[2:end]
                        passthrough_phase!(Di, (AscendingChain(AA), limb, DescendingChain(Des)), D)
                        Asc[2:end] = AA
                    else
                        passthrough_phase!(Di, (limb, DescendingChain(Des)), D)
                    end

                end
            else # limb on right

                if length(lps) >  0
                    lpk = pop!(lps)
                else
                    lpk = :nothing
                end

                if  lpk == :left
                    # translate then fuse
                    passthrough!(DescendingChain(Des), limb)
                    Des[end], Di = _fuse(Des[end], L)
                    passthrough_phase!(Di, (), D)
                else
                    Des[end], Di = _fuse(Des[end], L)
                    passthrough_phase!(Di, limb, (), D)

                    passthrough!(DescendingChain(Des), limb)
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
            passthrough_phase!(Di, (AscendingChain(Asc), limb), D)
        else
            passthrough_phase!(Di, (AscendingChain(Asc),), D)
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
