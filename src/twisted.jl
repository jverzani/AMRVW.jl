struct TwistedChain{T} <: AbstractRotatorChain{T}
x::Vector{T}
pv::Vector{Symbol}
m::Int
end
function TwistedChain(xs::Vector{T}) where {T}
    sigma = idx.(xs)
    m = minimum(sigma)
    ps = position_vector(sigma .- (m-1))
    TwistedChain(xs, ps, m)
end
Base.adjoint(A::TwistedChain) = TwistedChain(reverse(adjoint.(A.x)))
## XXX this is serious, now
function Base.getindex(A::TwistedChain{T}, i, j) where {T}
    DescendingChain(A.x)[i,j]
end


## Twisted refers to the QFactorization
##
struct QFactorizationTwisted{T, Rt} <: AbstractQFactorization{T, Rt, Val{:twisted}}
  Q::TwistedChain{Rt}
  D::AbstractSparseDiagonalMatrix{T, Rt}
end


## return q[i:k, j:k]
function Base.getindex(QF::QFactorizationTwisted, j, k)
    ## This is wrong...
    DescendingChain(QF.Q.x)[j,k]
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


# Fish out M with idx i
function iget(Ms, i)
    for (j,M) in enumerate(Ms)
        if idx(M) == i
            return (j, M)
        end
    end
end

function iget!(Ms,i)
    for (j,M) in enumerate(Ms)
        if idx(M) == i
            deleteat!(Ms, j)
            return M
        end
    end
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
    n = idx(A.x[1])
    N = idx(A.x[end])

    for i in 1:length(L)
        j  = idx(L[i])
        l  = j - n  +  1
        L[i], A[l],  A[l+1] = turnover(A[l], A[l+1],  L[i])
    end
end


function passthrough!(L::Vector, A::DescendingChain)
    n = idx(A.x[1])
    N = idx(A.x[end])

    for i in length(L):-1:1
        j  = idx(L[i])
        l  = j - n  +  1
        A[l-1],  A[l], L[i] = turnover(L[i], A[l-1], A[l])
    end
end

# need to specilize for diagonal
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
    @show :twisted_chain
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
    DesM = [M.x[1]]
    ps = M.pv
    for i in 1:length(ps)
        if ps[i] == :right
            break
        else
            push!(DesM, M.x[i+1])
        end
    end

    passthrough_phase(Di, DescendingChain(DesM), Vs, D)
    M.x[1:length(DesM)] = DesM

    passthrough_phase(nothing, Di, Vs, D)
    return nothing
end

function passthrough_phase(M::Nothing, Di::DiagonalRotator{T}, Vs::Tuple, D) where {T}

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
        if length(Vs) > 0
            Vhead, Vtail = Vs[1], Vs[2:end]
            passthrough_phase(Di, Vhead, Vtail, D)
        else
            fuse!(Di, D)
        end
        return nothing
    end



    n, N = extrema(V)
    if i == n-1
        # Di * Uj = Ujalpha * Di * Dj
        j  =  i  + 1
        Uj = V.x[1]
        c, s = vals(Uj)
        Ujalpha = Rotator(c*conj(alpha), s, j)
        Dj = DiagonalRotator(alpha, j)
        V.x[1] = Ujalpha

        ## Di
        if length(Vs) > 0
            passthrough_phase(Di, Vs[1], Vs[2:end], D)
        else
            fuse!(Di, D)
        end

        ##  Dj
        if length(V.x) > 1
            ## pass Dj through
            VV = V.x[2:end]
            passthrough_phase(Dj,  DescendingChain(VV), Vs, D)
            V.x[2:end] = VV
        else
            if length(Vs) > 0
                passthrough_phase(Dj, Vs[1], Vs[2:end], D)
            else
                fuse!(Dj, D)
            end
        end

    elseif i == n
        # D_i(alpha, i) * R(c,s,i) = R(c alpha/conj(alpha), s, i) * D(conj(alpha),i)
        # fuse with U
        U = V.x[1]
        c,s = vals(U)
        Ualpha = Rotator(c*alpha/conj(alpha), s, i)
        V.x[1] = Ualpha
        Dii = DiagonalRotator(conj(alpha), i)


        if length(V.x) > 1
            VV  = V.x[2:end]
            passthrough_phase(Dii, DescendingChain(VV), Vs, D)
            V.x[2:end] = VV

        else
            ## D = Di*Dj*D
            fuse!(Dii, D)
        end
    elseif n < i <= N
        # turnover
        V.x[i-n], V.x[i-n+1], Di = turnover(Di, V.x[i-n], V.x[i-n+1])
        if length(Vs) > 0
            Vhead, Vtail = Vs[1], Vs[2:end]
            passthrough_phase(Di, Vhead, Vtail, D)
        else
            fuse!(Di, D)
        end
    elseif i == N+1
        U = V.x[end]
        c,s = vals(U)
        Ualpha = Rotator(c*conj(alpha), s, N)
        Dh = DiagonalRotator(alpha, i-1)
        V.x[end] = Ualpha

        if length(Vs) > 0
            Vhead, Vtail = Vs[1], Vs[2:end]
            passthrough_phase(Dh, Vhead, Vtail, D)
            passthrough_phase(Di, Vhead, Vtail, D)
        else
            ## D = Di*Dj*D
            fuse!(Dh, D)
            fuse!(Di, D)
        end
    end


end

function passthrough_phase(Di::DiagonalRotator{T}, V::AscendingChain, Vs::Tuple, D) where {T}

    i = idx(Di)
    alpha, _ = vals(Di)

    if length(V) == 0
        if length(Vs) > 0
            Vhead, Vtail = Vs[1], Vs[2:end]
            passthrough_phase(Di, Vhead, Vtail, D)
        else
            fuse!(Di, D)
        end
        return nothing
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
        if length(Vs) > 0
            passthrough_phase(Di, Vs[1], Vs[2:end], D)
        else
            fuse!(Di, D)
        end

        ##  Dj
        if length(V.x) > 1
            ## pass Dj through
            VV = V.x[2:end]
            passthrough_phase(Dj,  AscendingChain(VV), Vs, D)
            V.x[2:end] = VV
        else
            if length(Vs) > 0
                passthrough_phase(Dj, Vs[1], Vs[2:end], D)
            else
                fuse!(Dj, D)
            end
        end

    elseif i == N
        # D_i(alpha, i) * R(c,s,i) = R(c alpha/conj(alpha), s, i) * D(conj(alpha),i)
        # fuse with U
        U = V.x[1]
        c,s = vals(U)
        Ualpha = Rotator(c*alpha/conj(alpha), s, i)
        V.x[1] = Ualpha
        Dii = DiagonalRotator(conj(alpha), i)


        if length(V.x) > 1
            VV  = V.x[2:end]
            passthrough_phase(Dii, AscendingChain(VV), Vs, D)
            V.x[2:end] = VV

        else
            ## D = Di*Dj*D
            fuse!(Dii, D)
        end
    elseif n <= i < N
        # turnover
        V.x[N+1-i-1], V.x[N+1-i], Di = turnover(Di, V.x[N+1-i-1], V.x[N+1-i])
        if length(Vs) > 0
            Vhead, Vtail = Vs[1], Vs[2:end]
            passthrough_phase(Di, Vhead, Vtail, D)
        else
            fuse!(Di, D)
        end
    elseif i == n-1
        U = V.x[end]
        c,s = vals(U)
        Ualpha = Rotator(c*conj(alpha), s, n)
        Dh = DiagonalRotator(alpha, i+1)
        V.x[end] = Ualpha

        if length(Vs) > 0
            Vhead, Vtail = Vs[1], Vs[2:end]
            passthrough_phase(Dh, Vhead, Vtail, D)
            passthrough_phase(Di, Vhead, Vtail, D)
        else
            ## D = Di*Dj*D
            fuse!(Dh, D)
            fuse!(Di, D)
        end
    end


end


function passthrough_phase(Di::DiagonalRotator{T}, V::TwistedChain, Vs::Tuple, D) where {T}
    i = idx(Di)
    alpha, _ = vals(Di)
    n, N = extrema(V)
    error("implement me")
end

##################################################


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
    # moddify Ms, RF
    bulge_step!(Ms, D, RF, Asc)

end

## Can enter here for more general usage
function bulge_step!(Ms, D, RF, Asc)

    Decoupled = eltype(Ms)[]

    Des = reverse(adjoint.(Asc))
    n = length(Ms) + 1
    m = length(Asc)
    sigma = [idx(R) for R in Ms]
    ps = position_vector(sigma)

    limb, limb_side = step_0!(ps, m, Ms, Des, Asc, D)
    @show :step0, eigvals(limb * (Des * (Ms * (Asc * diagm(0 => D.x)))))
    for k in 1:(n-m-2)

        limb_side = step_k!(k, n, m, ps, limb_side, limb, Des, Asc, Decoupled, Ms, RF)

    end
    @show :stepk, eigvals(Decoupled * (limb * (Des * (Ms * (Asc * diagm(0 => D.x))))))
    @show :stepk, eigvals(Decoupled * (limb * (Des * (Asc * (Ms * diagm(0 => D.x))))))

    # now we choose :left
    limb_side = step_k!(n-m-1, n, m, :left, limb_side, limb, Des, Asc, Decoupled, Ms, RF)


    # now knit in  limb, Des, Asc
    step_knit!(n, m,  limb_side, limb, Des,  Asc, Decoupled, Ms, D, RF)

    ## Decoupled is now Ms, assign somehow?

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
        passthrough_phase(TwistedChain(Ms), Di, (AscendingChain(Asc),), D)
    else
        Asc[1], Di = _fuse(U, Asc[1])      # fuse with  Ascending
        if length(Asc) > 1
            AA = Asc[2:end]
            passthrough_phase(nothing, Di, (AscendingChain(AA),), D) # XXX Descending too?
            Asc[2:end] = AA
        end
    end



    limb, limb_side



end

function step_k!( k, n, m, ps, limb_side, limb, Des, Asc, Decoupled, Ms, RF) where {T}

    phatk = isa(ps, Symbol) ? ps :  ps[k]

    if phatk == :left
        U = iget!(Ms, m + k )
        push!(Des, U)  ## Augement Descending
        passthrough!(DescendingChain(Des), Asc) ## translate ascending

    elseif phatk == :right

        U = iget!(Ms, m + k)
        pushfirst!(Asc, U)
        passthrough!(Des, AscendingChain(Asc))

    end


    if length(limb) > 0
        if limb_side == :left

            passthrough!(limb, AscendingChain(Asc))

        elseif limb_side == :right

            passthrough!(DescendingChain(Des), limb)

        end
    end
    limb_side = :middle

    if phatk == :left
        ## similarity moving Acending from left to right
        passthrough!(RF, Asc) # <--
        U = popfirst!(Des)
        push!(Decoupled, U)  # position consistent with :left?
        limb_side = :left
    elseif phatk == :right
        passthrough!(Asc, RF) #-->
        U = pop!(Asc)
        pushfirst!(Decoupled, U)
        limb_side = :right
    end

    return limb_side
end



function step_knit!( n, m,  limb_side, limb, Des,  Asc, Decoupled, Ms, D, RF) where {T}

    U = pop!(Des)
    V = popfirst!(Asc)
    bottom, Di = _fuse(U,V)

    if length(limb) > 0
        if limb_side == :left
            passthrough_phase(nothing, Di, (AscendingChain(Asc),), D)
        elseif limb_side == :right
            passthrough_phase(nothing, Di, (AscendingChain(Asc),TwistedChain(limb)), D)
        end
    else
        passthrough_phase(nothing, Di, (AscendingChain(Asc),), D)
    end


    ## Final Steps, knit in
    ## setps k=n-m to  n-2


    #bottom = Des[end]
    for k in (n-m):(n-2)


        # set pk = 1
        phatk = :left

        Des = push!(Des, bottom)
        passthrough!(DescendingChain(Des), Asc)

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
                        passthrough_phase(nothing, Di, (AscendingChain(Asc),), D)  # XXX and something prior?
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


        # similarity to move Asc to right side
        passthrough!(RF, Asc)  # <--
        push!(Decoupled,  popfirst!(Des))  # take Des and  move top to decoupled
        U = pop!(Des)
        V = popfirst!(Asc)
        bottom, Di = _fuse(U,V)

        if length(limb) > 0
            if limb_side == :left
                passthrough_phase(nothing, Di, (AscendingChain(Asc),), D)
            elseif limb_side == :right
                passthrough_phase(nothing, Di, (AscendingChain(Asc),TwistedChain(limb)), D)
            end
        else
            passthrough_phase(nothing, Di, (AscendingChain(Asc),), D)
        end

    end

    push!(Decoupled, bottom)

    return nothing

end


## For fun....
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
