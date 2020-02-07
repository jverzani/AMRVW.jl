##################################################

## For the complex-real rotator case we need a diagonal matrix to store the phase
## This is a  means to make this generic with respect to rotator type
## The factorization can have an identity diagonal or a real one.
abstract type AbstractSparseDiagonalMatrix{S} end


struct SparseDiagonal{S} <: AbstractSparseDiagonalMatrix{S}
x::Vector{S}
SparseDiagonal{S}(n) where {S} = new(ones(S, n))
SparseDiagonal(::Type{T}, n) where {T <: Real} = SparseDiagonal{T}(0)
SparseDiagonal(::Type{S}, n) where {S} = SparseDiagonal{S}(n)
end

Base.getindex(D::SparseDiagonal{T}, k) where {T <: Real} = one(T)
@inbounds Base.getindex(D::SparseDiagonal, k)  = D.x[k]
Base.@propagate_inbounds Base.setindex!(D::SparseDiagonal{T}, X, inds...) where {T <: Real} = X
Base.@propagate_inbounds Base.setindex!(D::SparseDiagonal, X, inds...) = setindex!(D.x, X, inds...)

Base.Matrix(D::SparseDiagonal{T}) where {T <: Real} = I
Base.Matrix(D::SparseDiagonal{S}) where {S} = diagm(0 => D.x)

*(D::SparseDiagonal{T}, M::Matrix) where {T <: Real} = M
*(D::SparseDiagonal, M::Matrix) = diagm(0=>D.x) * M
*(M::Matrix, D::SparseDiagonal)  = M * diagm(0=>D.x)




function fuse!(Di::DiagonalRotator, D::SparseDiagonal{S}) where {S <: Complex}
    i = idx(Di)
    alpha, _ = vals(Di)
    D[i] *= alpha
    D[i+1] *= conj(alpha)
    return nothing
end
fuse!(Di::DiagonalRotator, D::SparseDiagonal{T}) where {T} = nothing


## passthrough
## Pass a rotator through a diagonal matrix with phase shifts
## D U -> U' D' (right, as in U starts on right)
## U D -> D' U' (left)


@inline function passthrough!(D::SparseDiagonal{S}, U::Rotator) where {S <: Complex}
    i = idx(U)
    c,s = vals(U)

    alpha, beta = D[i], D[i+1]
    beta1 = alpha * conj(beta)

    D[i] = beta
    D[i+1] = alpha

    return Rotator(beta1 * c, s, i)

end


## U D -> D U
@inline function passthrough!(U::Rotator, D::SparseDiagonal{S}) where {S <: Complex}
    return passthrough!(D, U)
    i = idx(U)
    c,s = vals(U)

    alpha, beta = D[i], D[i+1]
    beta1 = alpha/beta

    D[i] *= conj(beta1)
    D[i+1] *= beta1
    return Rotator(beta1 * c, s, i)
end

passthrough!(D::SparseDiagonal{T}, U::Rotator) where {T <: Real} = U
passthrough!(U::Rotator, D::SparseDiagonal{T}) where {T <: Real}  = U


# Move a chain through a diagonal
function passthrough!(D::SparseDiagonal{S}, Ch::Union{DescendingChain, AscendingChain}) where {S <: Complex}
    for i in 1:length(Ch) # pass  ch through D
        Ch[i] = passthrough!(D, Ch[i])
    end
end

function passthrough!(Ch::Union{DescendingChain, AscendingChain}, D::SparseDiagonal{S}) where {S <: Complex}
    for i in length(Ch):-1:1
        Ch[i] = passthrough!(Ch[i],D)
    end
end

## could do two others...

## noop when Identity Diagonal Matrix
passthrough!(D::SparseDiagonal{T}, C::DescendingChain) where {T <: Real} =  nothing
passthrough!(D::SparseDiagonal{T}, C::AscendingChain) where {T <: Real} =  nothing
passthrough!(D::SparseDiagonal{T}, C::TwistedChain) where {T <: Real} =  nothing
passthrough!(C::DescendingChain, D::SparseDiagonal{T}) where {T <: Real} =  nothing
passthrough!(C::AscendingChain, D::SparseDiagonal{T}) where {T <: Real} =  nothing
passthrough!(C::TwistedChain, D::SparseDiagonal{T}) where {T <: Real} =  nothing

## Di U = U' Di'
passthrough!(D::IdentityRotator, U::AbstractRotator) = (U,D)
function passthrough!(Di::DiagonalRotator, U::AbstractRotator)
    @assert idx(Di) == idx(U)
    i = idx(Di)
    alpha, _ = vals(Di)
    c, s = vals(U)
    Rotator(c * alpha/conj(alpha), s, i), DiagonalRotator(conj(alpha), i)
end

passthrough!(U::AbstractRotator, D::IdentityRotator) = (U,D)
function passthrough!(U::AbstractRotator, Di::DiagonalRotator)
    @assert idx(Di) == idx(U)
    i = idx(Di)
    alpha, _ = vals(Di)
    c, s = vals(U)
    DiagonalRotator(conj(alpha), i), Rotator(c * alpha/conj(alpha), s, i)
end

# Handle the case |i-j| =  1
# return Dj, Vj where Di * Uj = Vj * Dj * Di
function tip(Di::DiagonalRotator, Uj::AbstractRotator)
    i, j = idx(Di),  idx(Uj)
    @assert abs(i-j) ==  1
    alpha, _ = vals(Di)
    c,s =  vals(Uj)
    Rotator(c*conj(alpha), s, j), DiagonalRotator(alpha, j)
end

function tip(Uj::AbstractRotator, Di::DiagonalRotator)
    i, j = idx(Di),  idx(Uj)
    @assert abs(i-j) ==  1
    alpha, _ = vals(Di)
    c,s =  vals(Uj)
    Rotator(c*conj(alpha), s, j), DiagonalRotator(alpha, j)
end
tip(Di::IdentityRotator,  Uj::AbstractRotator) = (Uj, Di)
tip(Uj::AbstractRotator, Di::IdentityRotator) = (Di, Uj)

##################################################
##
## Move a diagonanol rotator from a fuse operation to the right to merge in with D.
##
## pass Di thorugh A; fuses with D
##
## In general, a recursive algorithm can handle this task;
## This function  is speficic to this task, and hence is more efficient.
function passthrough_phase!(Di::DiagonalRotator, A::DescendingChain, D::SparseDiagonal)
    ## We have two rules here
    ## D_i(alpha) * R_{i+1}(c,s) = R_{i+1}(c, conj(alpha)s) D_i(alpha)
    ## R_{i+1}(c, conj(alpha) s) = R_{i+1}(conj(alpha) c, s) D_{i+1}(alpha)
    ## so D_i(alpha) * R_{i+1}(c,s) = R_{i+1}(conj(alpha) c, s) D{{i+1}(alpha) D_i(alpha)
    ## Also D_i(alpha) * R_i(c,s) = R_i(c*alpha/conj(alpha), s) D_i(conj(alpha)
    i, n = idx(Di), length(A)
    alpha, _ = vals(Di)

    ## D_i will passthrough and can be merged into D
    ## we will have D_i(alpha) * D_{i+1}(alpha) * ... * D_j(alpha) which will only have alpha at i and conj(alpha) at j+1
    D[i] *= alpha


    while i < n
        j = i + 1
        c, s = vals(A[j])
        #@assert j == idx(A[j])

        iszero(s) && break
        A[j] = Rotator(conj(alpha)*c, s, j)
        i = i + 1
    end

    D[i+1] *= conj(alpha)

    return nothing

end

function passthrough_phase!(Di::DiagonalRotator, A::AscendingChain, D::SparseDiagonal)

    i, n = idx(Di), length(A)
    alpha, _ = vals(Di)

    ## D_i will passthrough and can be merged into D
    ## we will have D_i(alpha) * D_{i+1}(alpha) * ... * D_j(alpha) which will only have alpha at i and conj(alpha) at j+1
    D[i] *= alpha

    i -= 1
    while i >= 1
        c, s = vals(A[i])
        iszero(s) && break
        A[i] = Rotator(conj(alpha)*c, s, i)
        i -= 1
    end

    if i > 1
        D[i-1] *= conj(alpha)
    end

    return nothing
end


# pass a diagonal rotator from a fuse operation through one or more rotator chains
# merge with D a diagonal matrix
passthrough_phase!(Di, Vs::Tuple, D::SparseDiagonal{T}) where {T <: Real} = nothing
passthrough_phase!(Di, C::AbstractRotatorChain, Vs::Tuple, D::SparseDiagonal{T}) where {T <: Real} = nothing

function passthrough_phase!(Di::DiagonalRotator, Vs::Tuple, D)
    if length(Vs) > 0
        Vhead, Vtail = Vs[1], Vs[2:end]
        passthrough_phase!(Di, Vhead, Vtail, D)
    else
        fuse!(Di, D)
    end
    return nothing
end

# peeled one off
function passthrough_phase!(Di::DiagonalRotator, V::DescendingChain, Vs::Tuple, D)

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



function passthrough_phase!(Di::DiagonalRotator, V::AscendingChain, Vs::Tuple, D)

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

## # two special functions useful to simplify the below:
## function passthrough_phase!(Di::DiagonalRotator, V::TwistedChain, Vs::Tuple, D, j, ::Val{:Des})
##     @show :Des
##     inds = iget(V, j, Val(:Des))
##     passthrough_phase!(Di, DescendingChain(view(V.x, inds)), Vs, D)
## #    passthrough_phase!(Di, DescendingChain(Des), Vs, D)
## #    for (i,j) in enumerate(inds)
## #        V[j] = Des[i]
## #    end
## end

## function passthrough_phase!(Di::DiagonalRotator, V::TwistedChain, Vs::Tuple, D, j, ::Val{:Asc})
##     @show :Asc
##     inds = iget(V, j, Val(:Asc))
##     passthrough_phase!(Di, AscendingChain(view(V.x, inds)), Vs, D)
## #    passthrough_phase!(Di, AscendingChain(Asc), Vs, D)
## #    for (i,j) in enumerate(inds)
## #        V[j] = Asc[i]
## #    end
## end

function passthrough_phase!(Di::DiagonalRotator, V::TwistedChain, Vs::Tuple, D)

    # Handle case where limb is small (length 0,1, or 2)
    if length(V.x) == 0
        return passthrough_phase!(Di, Vs, D)
    elseif length(V.x) <= 2
        if isempty(V.pv) || V.pv[1] == :left
            return passthrough_phase!(Di, DescendingChain(V.x), Vs, D)
        else
            VV = reverse(V.x)
            passthrough_phase!(Di, AscendingChain(VV), Vs, D)
            V.x[1] = VV[2]
            V.x[2] = VV[1]
            return nothing
        end
    end


    n, N = extrema(V)
    pv = V.pv

    i = idx(Di)  # pv[i-n+1] informs if Ui+1 is left or right of Ui

    if i < n - 1
        passthrough_phase!(Di, Vs, D) # misses chain
    elseif i == n - 1
        passthrough_phase!(Di, descending_part(V, n-1), Vs, D)#,  n-1, Val(:Des))
    elseif i == n
        if pv[1] == :left
            ind, Ui = iget(V, i)
            Ui, Di = passthrough!(Di, Ui)
            V[ind] = Ui
            passthrough_phase!(Di, descending_part(V,i), Vs, D)#, i, Val(:Des))
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
                ii,  Ui =  iget(V,  i)
                ij,  Uj =  iget(V, i+1)
                Uj, Ui, Dj = turnover(Di, Uj, Ui)
                V[ii] = Ui
                V[ij] = Uj
                passthrough_phase!(Dj, descending_part(V,i+1), Vs,  D)#,  i+1, Val(:Des))
            end
        end
    elseif i == N
        if pv[end] == :right
            # flip, get Ascending
            ind, Ui = iget(V, i)
            Ui, Di = passthrough!(Di, Ui)
            V[ind] = Ui
            passthrough_phase!(Di, ascending_part(V,i), Vs, D)#,   i, Val(:Asc))
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
            passthrough_phase!(Dh, ascending_part(V,i-1), Vs,  D)#,  i-1, Val(:Asc))
        end
    elseif i == N + 1
        passthrough_phase!(Di, ascending_part(V,N+1), Vs, D)#, N+1, Val(:Asc))
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
                passthrough_phase!(Dh, ascending_part(V,i-1), Vs, D)#, i-1, Val(:Asc))
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
                passthrough_phase!(Dj, descending_part(V,i+1), Vs, D)#, i+1, Val(:Des))
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
                passthrough_phase!(Dh, ascending_part(V,i-1), Vs, D)#, i-1, Val(:Asc))
            else
                passthrough_phase!(Dh, Vs, D)
            end

            ## Move Dj passed Des, if present
            if i-n + 2 <= length(pv)  && pv[i-n+2] == :left
                passthrough_phase!(Dj, descending_part(V,i+1), Vs, D)#, i+1, Val(:Des))
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
                passthrough_phase!(Dh, ascending_part(V,i-1), Vs, D)#, i-1, Val(:Asc))
            else
                passthrough_phase!(Dh, Vs, D)
            end

            ## Move Dj passed Des, if present

            if i-n+2 <= length(pv)  && pv[i-n+2] == :left
                passthrough_phase!(Dj, descending_part(V,i+1), Vs, D)#, i+1, Val(:Des))
            else
                passthrough_phase!(Dj, Vs, D)
            end
        end

    end

    return nothing
end
