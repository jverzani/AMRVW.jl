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
## * Only works for unitary matrices (R an indentity!!!)
## * get shifts by paper (m -> eigenvalues of mxm matrix, continue...)
##
## This algorithm encompasses those in CSS and RDS, but those are more efficient space wise
## This is more general and perhaps of interest as building blocks for experimentation



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

            V = iget(Ms, i+1)
            if V != nothing
                j,U = V
                push!(inds, j)
                push!(Asc, U)
            else
                break
            end

        else
            break
        end
        i += 1
    end

    inds, Asc
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

## This is not correct!
## This should be matrix multiplication
## as Twisted Q is no long Hessenberg
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
    pv = state.QF.Q.pv
    n = length(Ms)
    stop_ctr = n

    m = isa(Rt, RealRotator{T})  ? 2 : 1

    N0 = findlast(x->x==:right, pv)
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



## We use `bulge_step!` for more general usage; this is
## tied to AbstractFactorizationState{T, S, Rt, QFt, RFt, Val{:twisted}}
function bulge_step(state::AbstractFactorizationState{T, S, Rt, QFt, RFt, Val{:twisted}})  where {T,S, Rt, QFt, RFt}

    # create bulge is only right when  all rotators in QFt are straightened out,
    # but  this should step in the correct direction.
    create_bulge(state)
    Asc = reverse(state.UV)  # stored as U, V
    Ms = state.QF.Q
    D = state.QF.D
    RF = state.RF
    # modify Ms, D, RF
    bulge_step!(Ms, D, RF, Asc)

end

## Can enter here for more general usage
function bulge_step!(Ms::TwistedChain, D, RF::AbstractRFactorization, Asc, choice = :left)


    Rt = eltype(Ms)
    Decoupled = Rt[]

    Des = reverse(adjoint.(Asc))


    n = length(Ms) + 1
    m = length(Asc)

    ps = copy(Ms.pv)
    psd = ps[m+1:end]
    if isa(choice, Symbol)
        append!(psd, repeat([choice],m))
    else
        append!(psd, choice)
    end

    M0 = Des * (Ms * ((Matrix(D) * Matrix(RF)) * Asc))
    @show eigvals(M0)[1]

    ## Must put AscendingChain into place
    ## py passing through RF and D <--
    passthrough!(RF, AscendingChain(Asc))
    passthrough!(D, AscendingChain(Asc))


    limb, limb_side = step_0!(m, ps, Ms, Des, Asc, D)

    if limb_side == :left
        M1 = limb * (Des * (Ms * ( Asc * (Matrix(D) * Matrix(RF)))))
        @show eigvals(M1)[1]
    else
        M1 = Des * (Ms * ( Asc * (limb * (Matrix(D) * Matrix(RF)))))
        @show eigvals(M1)[1]
    end

    for k in 1:(n-m-2)
        @show :step_k
        limb_side = step_k!(k, n, m, ps, psd, limb_side, limb, Des, Asc, Decoupled, Ms, D, RF)




    end


    # now we choose :left
    @show :last_step_k
    limb_side = step_k!(n-m-1, n, m, ps, psd, limb_side, limb, Des, Asc, Decoupled, Ms, D, RF)

    @show :step_knit
    # now knit in  limb, Des, Asc
    step_knit!(n, m, psd, limb_side, limb, Des,  Asc, Decoupled, Ms, D, RF)

    ## Decoupled is now Ms, assign and update position vector
    append!(Ms.x, Decoupled)
    empty!(Ms.pv)
    append!(Ms.pv, psd) #position_vector(idx.(Ms.x)))
    Ms.m[] = 1

    M7 = Ms * (Matrix(D) * Matrix(RF))
    @show eigvals(M7)[1]

    return nothing
end


############################################################
##
## Implement thee different steps for one pass of Francis' algorithm


function step_0!(m,  ps, Ms::TwistedChain, Des, Asc, D) where {T}

    has_limb = ifelse(m > 1, true, false)


    # grab limb, but keep order
    limb = _get_limb!(Ms, m)
    # TODO: wrap limb=TwistedChain(...), as o/w position vecgtor is computed each time

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


    U = iget!(Ms, m)
    popfirst!(Ms.pv)
    Ms.m[] += 1

    if ps[m] == :left
        Des[end], Di = _fuse(Des[end], U)  # fuse with descending; aka Ms[m]

        ## need to passthrough Ms too now
        i = idx(Di)

        inds, MMs = iget(Ms, i, Val(:Des))
        if limb_side == :right
            passthrough_phase!(Di, DescendingChain(MMs),  (AscendingChain(Asc), limb), D)
        else
            passthrough_phase!(Di, DescendingChain(MMs),  (AscendingChain(Asc), ), D)
        end
        for (i,j) in enumerate(inds)
            Ms[j] = MMs[i]
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

function step_k!(k, n, m, ps, psd, limb_side, limb, Des, Asc, Decoupled, Ms, D, RF) where {T}

    if m+k <= length(ps)
        phatk = ps[m+k]
    else
        @show k, n-m-1, n,m,psd
        phatk = psd[k+1]
    end
@show n,k,m, ps, phatk, limb_side
    ## 4 cases based on:
    ## phatk/limb_side
    ## A bit redundandant, but easier to verify

    U = iget!(Ms, m + k )

    if phatk == :left && limb_side == :left

        push!(Des, U)  ## Augment Descending
        passthrough!(DescendingChain(Des), AscendingChain(Asc)) ## translate ascending

        passthrough!(limb, AscendingChain(Asc))

        ## similarity transform Asc to right side
        passthrough!(RF, AscendingChain(Asc)) # <-- pass Asc through RF
        passthrough!(D, AscendingChain(Asc))

        U = popfirst!(Des)
        if k == 1 || isa(ps, Symbol) || ps[m+k-1] == :left
            push!(Decoupled, U)
        else
            pushfirst!(Decoupled, U)
        end

        M4ll = limb * (Decoupled * (Des * (Ms * (Asc * (Matrix(D) * Matrix(RF))))))
        @show eigvals(M4ll)[1]

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
        if k == 1 || isa(ps, Symbol) || ps[m+k-1] == :left
            push!(Decoupled, U)
        else
            pushfirst!(Decoupled, U)
        end


        M4lr = limb * (Decoupled * (Des * (Ms * (Asc * (Matrix(D) * Matrix(RF))))))
        @show eigvals(M4lr)[1]

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
        if k == 1 || isa(ps, Symbol) || ps[m+k-1] == :left
            push!(Decoupled, U)
        else
            pushfirst!(Decoupled, U)
        end
        #push!(Decoupled, U)

        M4rl = Des * (Ms * (Asc * (Decoupled * (limb * (Matrix(D) * Matrix(RF))))))
        @show eigvals(M4rl)[1]

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

        U = pop!(Asc)
        if k == 1 || isa(ps, Symbol) || ps[m+k-1] == :left
            push!(Decoupled, U)
        else
            pushfirst!(Decoupled, U)
        end
        #push!(Decoupled, U)

        M4rr = Des * (Ms * (Asc * (Decoupled * (limb * (Matrix(D) * Matrix(RF))))))
        @show eigvals(M4rr)[1]

        if m > 1
            limb_side = :right
        end

    end

    return limb_side

end

## steps k=n-m to  n-2
function step_knit!(n, m, psd, limb_side, limb, Des,  Asc, Decoupled, Ms, D, RF) where {T}
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
        phatk = psd[k]
        if phatk == :left
            Des = push!(Des, bottom)
        else
            Asc = pushfirst!(Asc, bottom)
        end
        passthrough!(DescendingChain(Des), AscendingChain(Asc))

        if length(limb) > 0

            # remove bottom rotator from limb
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

                    # then fuse
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
                    # fuse  then translate
                    Des[end], Di = _fuse(Des[end], L)
                    passthrough_phase!(Di, limb, (), D)

                    passthrough!(DescendingChain(Des), limb)
                end
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
        @show :push_decoupled
        if phatk == :left
            push!(Decoupled, popfirst!(Des))
        else
            pushfirst!(Decoupled, pop!(Asc))
        end


        M5 = limb * (Decoupled * (Des * (Asc * (Matrix(D) * Matrix(RF)))))
        @show eigvals(M5)[1]


        U = pop!(Des)
        V = popfirst!(Asc)
        bottom, Di = _fuse(U,V)

        if limb_side == :right
            passthrough_phase!(Di, (AscendingChain(Asc), limb), D)
        else
            passthrough_phase!(Di, (AscendingChain(Asc),), D)
        end

    end

    @show :push_decoupled_last
    if psd[end] == :left
        push!(Decoupled, bottom)
    else
        pushfirst!(Decoupled, bottom)
    end

    M6 = Decoupled * (Matrix(D) * Matrix(RF))

    return nothing

end
