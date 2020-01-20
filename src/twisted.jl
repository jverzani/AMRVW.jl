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
## * Clean up allocations
##
## * store Decoupled and Ms in same space to avoid allocations
##
## * get shifts by paper (m -> eigenvalues of mxm matrix, continue...)
##
## This algorithm encompasses those in CSS and RDS, but those are more efficient space wise
## This is more general and perhaps of interest as building blocks for experimentation

##################################################


##
## Twisted refers to the QFactorization
##
struct QFactorizationTwisted{T, S} <: AbstractQFactorization{T, S}
  Q::TwistedChain{T,S} ## Twisted
  D::SparseDiagonal{S}
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
struct QRFactorizationTwisted{T, S,  Qt<:QFactorizationTwisted{T,S}, Rt<:AbstractRFactorization{T,S}} <: AbstractFactorizationState{T, S, Val{:twisted}}
  N::Int
  QF::Qt
  RF::Rt
  UV::Vector{Rotator{T,S}}    # Ascending chain for cfreating bulge
  W::Vector{Rotator{T,S}}     # the limb when m > 1, a twisted chain
  A::Matrix{S}
  REIGS::Vector{T}
  IEIGS::Vector{T}
  ctrs::AMRVW_Counter
  QRFactorizationTwisted{T,S,Qt, Rt}(N, QF, RF, UV,W, A, REIGS, IEIGS, ctrs) where {T, S, Qt, Rt} = new(N, QF, RF, UV,W, A, REIGS, IEIGS, ctrs)
  QRFactorizationTwisted(N::Int, QF::QFactorizationTwisted{T,S}, RF::Rt, UV,W, A, REIGS, IEIGS, ctrs) where {T, S, Rt} = QRFactorizationTwisted{T, S, QFactorizationTwisted{T,S}, Rt}(N,QF,RF, UV,W, A, REIGS, IEIGS, ctrs)
end

# XXX should pass in m, so that A=m x m matrix
function QRFactorization(
                         QF::QFactorizationTwisted{T, S},
                         RF::AbstractRFactorization{T, S},
                         m = 2
                         ) where {T, S}

    N = length(QF) + 1
    A = zeros(S, m, m)
    reigs = zeros(T, N)
    ieigs = zeros(T, N)
    ctr = AMRVW_Counter(0, 1, N-1, 0, N-2)
    m = S == T ? 2 : 1
    UV = Vector{Rotator{T, S}}(undef, m)
    W = Vector{Rotator{T, S}}(undef, m-1)
    QRFactorizationTwisted(length(QF), QF, RF, UV, W, A, reigs, ieigs, ctr)
end


## struct QRFactorizationTwisted{T, S, Rt, QFt, RFt} <: AbstractFactorizationState{T, S, Rt, QFt, RFt, Val{:twisted}}
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
end

function _create_bulge(state, shifts)

end

##################################################

function eigvals(state::QRFactorizationTwisted)

    new_state = AMRVW_algorithm(state)
    es = complex.(new_state.REIGS, new_state.IEIGS)
    LinearAlgebra.sorteig!(es)
    es

end

## We don't go for the fastest method
## We have the first n-m  steps through, a straightening out of the twisted factorization in QF,
## This should be  O(n^2) steps
## Then, once untwisted, we pass down to the more efficient O(n^2) algorithm (in time) to solve
function AMRVW_algorithm(state::QRFactorizationTwisted{T, S}) where  {T, S}

    Ms = state.QF.Q.x # the vector
    pv = state.QF.Q.pv
    n = length(Ms)
    stop_ctr = n
    m = (T == S) ? 2 : 1

    N0 = findlast(x->x==:right, pv)
    N = N0 == nothing ? 1 : N0


    while N >= 0

        bulge_step(state)

        # check for deflation of last one
        c,s = vals(Ms[state.ctrs.stop_index])
        if abs(s) <= eps(T)
            state.ctrs.stop_index -= 1
        end

        N -= m
    end

    # now untwisted, so we can change over
    QF = state.QF
    Ms = QF.Q.x
    QF_new = QFactorization(DescendingChain(Ms), QF.D)
    new_state = QRFactorization(QF_new, state.RF)
    AMRVW_algorithm(new_state)
    new_state

end




## We use `bulge_step!` for more general usage; this is
## tied to AbstractFactorizationState{T, S, Rt, QFt, RFt, Val{:twisted}}
function bulge_step(state::QRFactorizationTwisted{T, S})  where {T,S}

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
## Note: Ms is a Twisted Chain, but need not stay that way
## We could reuse the storage space from Ms to hold the values stored in Decoupled
## at the cost of keeping track of the respective boundaries (k in the case of Decoupled, and k+m)
## for Ms. This could save some allocations, as this step allocates a new Decoupled vector
## each pass through.
function bulge_step!(Ms::TwistedChain, D, RF::AbstractRFactorization, Asc, choice = :left)


    Rt = eltype(Ms)
    Decoupled = Rt[]

    Des = reverse(adjoint.(Asc))


    n = length(Ms) + 1
    m = length(Asc)

    # psd[k] = ps[m+k] with padding soecufued by choice
    ps = copy(Ms.pv)
    psd = ps[m+1:end]
    if isa(choice, Symbol)
        append!(psd, repeat([choice],m))
    else
        append!(psd, choice)
    end


    ## M0a = Des * (Ms * (Matrix(D) * (Matrix(RF)* Asc)))
    ## @show eigvals(M0a)[1]

    ## Must put AscendingChain into place
    ## py passing through RF and D <--
    passthrough!(RF, AscendingChain(Asc))
    passthrough!(D, AscendingChain(Asc))

    ## M0 = Des * (Ms * (Asc * (Matrix(D) * Matrix(RF))))
    ## @show eigvals(M0)[1]


    limb, limb_side = step_0!(m, ps, Ms, Des, Asc, D)


## if limb_side == :left
##     M1 = limb * (Des * (Ms * (Asc * (Matrix(D) * Matrix(RF)))))
##     @show eigvals(M1)[1]
## else
##     M1 = Des * (Ms * (Asc * (limb *(Matrix(D) * Matrix(RF)))))
##     @show eigvals(M1)[1]
## end

    for k in 1:(n-m-2)
        limb_side = step_k!(k, n, m, psd, limb_side, limb, Des, Asc, Decoupled, Ms, D, RF)
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
    limb_side = step_k!(n-m-1, n, m, psd, limb_side, limb, Des, Asc, Decoupled, Ms, D, RF)

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
    step_knit!(n, m, psd, limb_side, limb, Des,  Asc, Decoupled, D, RF)

    ## Decoupled is now Ms, assign and update position vector
    append!(Ms.x, Decoupled)
    empty!(Ms.pv)
    append!(Ms.pv, psd)

    return nothing
end


############################################################
##
## Implement thee different steps for one pass of Francis' algorithm


function step_0!(m,  ps, Ms::TwistedChain, Des, Asc, D) where {T}
    has_limb = ifelse(m > 1, true, false)

    # grab limb, but keep order
    limb = _get_limb!(Ms, m)

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

    if ps[m] == :left
        Des[end], Di = fuse(Des[end], U)  # fuse with descending; aka Ms[m]
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

# limb is first m-1 rotators keeping their positions
function _get_limb!(Ms::TwistedChain, m::Int)
    ## XXX for now, we shorten twisted chain
    lx =  Ms.x[1:m-1]
    for i in 1:m-1
        popfirst!(Ms.x)
    end
    ps = m > 2 ? Ms.pv[1:m-2] : Symbol[]
    for i in 1:(m-1)
        !isempty(Ms.pv) && popfirst!(Ms.pv)
    end
    return TwistedChain(lx, ps)
end

function step_k!(k, n, m, psd, limb_side, limb, Des, Asc, Decoupled, Ms, D, RF) where {T}
    phatk = psd[k]

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
        push!(Decoupled, U)
    else
        push!(Decoupled, U)
        #pushfirst!(Decoupled, U)
    end

    return limb_side

end

## steps k=n-m to  n-2
function step_knit!(n, m, psd, limb_side, limb, Des,  Asc, Decoupled, D, RF) where {T}


    # make bottom
    U = pop!(Des)
    V = popfirst!(Asc)
    bottom, Di = fuse(U,V)

    ## Decoupled too! Decoupled on right? as we augmented Ascending
    if psd[end-m+1] == :right
        _Decoupled = TwistedChain(Decoupled, psd[1:length(Decoupled)-1])
        if limb_side == :right
            passthrough_phase!(Di, (AscendingChain(Asc), _Decoupled, limb), D)
        else
            passthrough_phase!(Di, (AscendingChain(Asc),  _Decoupled), D)
        end
    else
        if limb_side == :right
            passthrough_phase!(Di, (AscendingChain(Asc), limb), D)
        else
            passthrough_phase!(Di, (AscendingChain(Asc),), D)
        end
    end



    ## if limb_side == :right
    ##     if psd[end-m+1] == :right
    ##         passthrough_phase!(Di, (AscendingChain(Asc), TwistedChain(Decoupled, psd[1:length(Decoupled)-1]), limb), D)
    ##     else
    ##         passthrough_phase!(Di, (AscendingChain(Asc), limb), D)
    ##     end
    ## else
    ##     if psd[end-m+1] == :right
    ##         passthrough_phase!(Di, (AscendingChain(Asc),  TwistedChain(Decoupled, psd[1:length(Decoupled)-1])), D)
    ##     else
    ##         passthrough_phase!(Di, (AscendingChain(Asc),), D)
    ##     end
    ## end

    ## Final Steps, knit in
    ## setps k=n-m to  n-2

    for k in (n-m):(n-2)
        phatk = psd[k-1]
        if k > (n-m) && psd[k-2] != phatk
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
            push!(Decoupled, U)
        else
            U = pop!(Asc)
            push!(Decoupled, U)
            #pushfirst!(Decoupled, U)
        end

        ## M = Matrix(D) * Matrix(RF)
        ## @show limb_side, phatk
        ## if limb_side == :left && phatk == :left
        ##     Mnll = limb * (Decoupled * (Des * (Asc * M)))
        ##     @show eigvals(Mnll)[1]
        ## elseif limb_side == :left && phatk == :right
        ##     Mnlr = limb * (Des * (Asc * (Decoupled * M)))
        ##     @show eigvals(Mnlr)[1]
        ## elseif limb_side == :right && phatk == :left
        ##     Mnrl = Decoupled * (Des * (Asc * (limb * M)))
        ##     @show eigvals(Mnrl)[1]
        ## elseif limb_side == :right && phatk == :right
        ##     Mnrr = Des * (Asc * (limb * (Decoupled * M)))
        ##     @show eigvals(Mnrr)[1]
        ## end


        U = pop!(Des)
        V = popfirst!(Asc)
        bottom, Di = fuse(U,V)

        if phatk == :right
            _Decoupled =  TwistedChain(Decoupled, psd[1:length(Decoupled)-1])
            if limb_side == :right
                passthrough_phase!(Di, (AscendingChain(Asc), _Decoupled, limb), D)
            else
                passthrough_phase!(Di, (AscendingChain(Asc), _Decoupled), D)
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
        if m > 1 && psd[end-1] != psd[end]
            bottom = passthrough!(RF, bottom)  # move to other side
            bottom = passthrough!(D, bottom)
        end
        push!(Decoupled, bottom)
    else
        if m > 1 && psd[end-1] != psd[end]
            bottom = passthrough!(bottom, D)   # move bottom to other side
            bottom = passthrough!(bottom, RF)
        end
        push!(Decoupled, bottom)
    end

    return nothing

end
