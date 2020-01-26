## A factorization type for upper diagonal matrices
## RankOne{T,S}
## Z{T,S}, pencil
## UpperTriangular{T,S}
## Identity{T,S}

abstract type AbstractRFactorization{T, S} end

# Basic interface for an RFactorization:
#
# getindex (need just k-2 <= j <= k, but as general as possible)
# passthrough!(RF, U) to pass a rotator through from left to right or right to left
# simple_passthrough!(RF, U) for a shortcut passthrough
# we also might need:
# length
# eltype
# Matrix
# zero(RF{T, Rt}) to give 0
# one(RF{T, Rt}) to give 1
Base.eltype(RF::AbstractRFactorization{T, S}) where {T,S} = S
Base.zero(RF::AbstractRFactorization) = zero(eltype(RF))
Base.one(RF::AbstractRFactorization) = one(eltype(RF))

function passthrough!(RF::AbstractRFactorization, Us::Union{AscendingChain, DescendingChain})
    for i in  eachindex(Us.x)
        Us.x[i] =  passthrough!(RF, Us.x[i])
    end
end


# passthrough Us -> R
function passthrough!(Us::Union{AscendingChain, DescendingChain}, RF::AbstractRFactorization)
    for i in length(Us.x):-1:1
         Us.x[i] =  passthrough!(Us.x[i], RF)
    end
end

##################################################

## Rank one decomposition for amrvw algorithm to find eigenvalues of companion matrix
struct RFactorizationRankOne{T,S, V} <: AbstractRFactorization{T, S}
  Ct::AscendingChain{T,S,V}
  B::DescendingChain{T,S,V}
  D::SparseDiagonal{S}
end

Base.length(RF::RFactorizationRankOne) = length(RF.Ct)
Base.size(RF::RFactorizationRankOne) = (length(RF)+1, length(RF)+1)


## getindex

function _w(B, j, k)

    ck, sk = vals(B[k])
    if j == k
        return -sk
    else
        cj, sj = vals(B[j])
        s = one(sj)
        for i in (j+1):(k-1)
            ci, si = vals(B[i])
            s *= si
        end
        return conj(cj) * s *  ck
    end
end


## We have W = C * R, where W invovles B and D, and C involves Ct; This
## solves for R terms, as in the paper
## To see the W pattern
## function rotm(a::T,b, i, N) where {T}
##      r = diagm(0 => ones(T, N))
##      r[i:i+1, i:i+1] = [a b; -b conj(a)]
##      r
## end
## @vars bgc bhc bic bjc bkc blc
## @vars bgs bhs bis bjs bks bls
## @vars dg dh di dj dk dl dm
## Bg = rotm(bgc, bgs, 1, 7)
## Bh = rotm(bhc, bhs, 2, 7)
## Bi = rotm(bic, bis, 3, 7)
## Bj = rotm(bjc, bjs, 4, 7)
## Bk =  rotm(bkc, bks, 5, 7)
## Bl = rotm(blc, bls, 6, 7)
## D = diagm(0  => [dg, dh, di, dj, dk,dl, dm])
## en = [0,0,0,0, 0, 1,0]
## Ws = Bg * Bh  * Bi * Bj *  Bk *  *Bl * D   *  en
## To see the C pattern:
## @vars bc_i bs_i bc_j bs_j bc_k bs_k
## @vars cc_g cs_g cc_h cs_h cc_i cs_i cc_j cs_j cc_k cs_k
## @vars d_i d_j d_k d_l
## @vars what wg wh wi wj wk wl
## u = rotm(cc_k, cs_k, 1,2) * [what, wl]
## kk = u[1](what => solve(u[2], what)[1]) |> simplify

## u = rotm(cc_k, cs_k, 2,3) *  rotm(cc_j, cs_j, 1,3) * [what, wk, wl]
## kk_1 = u[1](what => solve(u[3], what)[1]) |> simplify

## u = rotm(cc_k, cs_k,3, 4)*rotm(cc_j, cs_j, 2,4) *  rotm(cc_i, cs_i, 1,4) * [what, wj, wk, wl]
## kk_2 = u[1](what => solve(u[4], what)[1]) |> simplify

## u = rotm(cc_k, cs_k, 4,5) * rotm(cc_j, cs_j,3, 5)*rotm(cc_i, cs_i, 2,5) *  rotm(cc_h, cs_h, 1,5) * [what, wi, wj, wk, wl]
## kk_3 = u[1](what => solve(u[5], what)[1]) |> simplify

## u = rotm(cc_k, cs_k, 5,6) * rotm(cc_j, cs_j, 4,6) * rotm(cc_i, cs_i,3, 6)*rotm(cc_h, cs_h, 2,6) *  rotm(cc_g, cs_g, 1,6) * [what, wh, wi, wj, wk, wl]
## kk_4 = u[1](what => solve(u[6], what)[1]) |> simplify
function Base.getindex(RF::RFactorizationRankOne, j, k)
    # neeed to compute Cts and Ws(B,D)
    Ct = RF.Ct
    B = RF.B
    D = RF.D
    N = length(Ct.x)

    ## return 0 if request is out of bounds
    ## or in lower triangular part (k < j)
    (k < j || j == 0 || k > N+1) && return zero(RF)

    cj, sj  = vals(Ct[N+1-j])
    s = sj
    par = -1


    dk = D[k]
    w  =  _w(B, j,  k) * dk


    tot = w/s

    for i in (j+1):k
        w =  _w(B,  i,  k)
        ci, si = vals(Ct[N+1-i])
        s  = s*si
        tot += par *  cj*conj(ci)/s * w * dk
        par *= -1
    end

    tot
end

function Base.Matrix(RF::RFactorizationRankOne{T,S}) where {T,S}

    n = length(RF) + 1

    M = diagm(0 => ones(S, n))
    e1 = vcat(1, zeros(S, n-1))
    en1 = vcat(zeros(S,n-1), 1)
    en = vcat(zeros(S,n-2), 1, 0)

    ## compute yt
    Ct, B = RF.Ct, RF.B
    D = Matrix(RF.D)
    D = D * M
    rho = (en1' * (Ct * M) * e1)
    yt = -(1/rho * en1' * (Ct*(B*D)))

    # Compute R
    R =   (Ct * (B * D)) +  (Ct * (e1 * yt))
end





## Pass a rotator through Rfactorization from left or right
function passthrough!(RF::RFactorizationRankOne, U::AbstractRotator)
    U = passthrough!(RF.D, U)
    U = passthrough!(RF.B, U)
    U = passthrough!(RF.Ct, U)

    U
end

function passthrough!(U::AbstractRotator, RF::RFactorizationRankOne)

    U = passthrough!(U, RF.Ct)
    U = passthrough!(U, RF.B)
    U = passthrough!(U, RF.D)

    U
end




##################################################

# hold upper triangular matrix as full matrix
struct RFactorizationUpperTriangular{T, S, Rt <: AbstractArray{S,2}} <: AbstractRFactorization{T, S}
R::Rt
RFactorizationUpperTriangular{T, S, Rt}(M) where {T, S, Rt}= new(M)
RFactorizationUpperTriangular(M::AbstractArray{S}) where {S} = RFactorizationUpperTriangular{real(S), S, typeof(M)}(M)
end

Base.length(RF::RFactorizationUpperTriangular) = size(RF.R)[1]
Base.size(RF::RFactorizationUpperTriangular) = size(RF.R)
function Base.getindex(RF::RFactorizationUpperTriangular, i, j)
    if i > 0 && j > 0
        RF.R[i,j]
    else
        zero(eltype(RF.R))
    end
end
Base.Matrix(RF::RFactorizationUpperTriangular) = RF.R

function passthrough!(U::Rt, RF::RFactorizationUpperTriangular) where {Rt <: AbstractRotator}

    R = RF.R
    n = size(R)[2]
    i = idx(U); j = i+1

#    R .= U * R
    c,s = vals(U)
    for k in i:n
        rik, rjk =  R[i,k],  R[k,k]
        R[i,k] = c * rik + s * rjk
        R[j,k] = -conj(s) * rik + conj(c) * rjk
    end

    g2 = givens(R[i+1,i], R[i+1,i+1],1,2)[1]
    c, s = conj(g2.s), real(g2.c)
    V = Rt(c,s,i)
#    R .= R*V
    for k in 1:j
        rki, rkj =  R[k,i], R[k,j]
        R[k, i] = c * rki - conj(s) * rkj
        R[k, j] = s * rki + conj(c) * rkj
    end
    R[j,i] = zero(eltype(R))
    V'

end


function passthrough!(RF::RFactorizationUpperTriangular, V::Rt) where {Rt <: AbstractRotator}

    R = RF.R
    n = size(R)[2]


    #R .= R * V
    c, s = vals(V)
    i = idx(V); j = i+1
    for k in 1:j
        rki, rkj =  R[k,i], R[k,j]
        R[k,i] = c * rki - conj(s) * rkj
        R[k,j] = s * rki + conj(c) * rkj
    end


    g2 = givens(R[i+1,i], R[i,i],1,2)[1]
    c, s = g2.s, real(g2.c)
    U = Rt(c,s,i)
    #R .= U * R
    for k in i:n
        rik, rjk = R[i,k], R[j,k]
        R[i,k ]  = c * rik + s * rjk
        R[j, k]  = -conj(s) * rik + conj(c) * rjk
    end
    R[j,i] = zero(eltype(R))
    U'
end


simple_passthrough!(RF::RFactorizationUpperTriangular, U::AbstractRotator) = false
simple_passthrough!(RF::RFactorizationUpperTriangular, U::AbstractRotator, V::AbstractRotator) = false


##################################################
# Wrapper for I
struct RFactorizationIdentity{T, S} <: AbstractRFactorization{T, S}
RFactorizationIdentity{T, S}() where {T, S}= new()
end

Base.size(RF::RFactorizationIdentity) = (-1, -1)
Base.getindex(RF::RFactorizationIdentity{T, S}, i, j) where {T, S} = i == j ? one(S) : zero(S)
Base.Matrix(::RFactorizationIdentity) = I
Base.length(RF::RFactorizationIdentity) = error("No dimension known")
passthrough!(RF::RFactorizationIdentity, U::AbstractRotator) = U
passthrough!(U::AbstractRotator, RF::RFactorizationIdentity) = U
passthrough!(RF::RFactorizationIdentity, C::DescendingChain) = nothing
passthrough!(RF::RFactorizationIdentity, C::AscendingChain) = nothing
passthrough!(RF::RFactorizationIdentity, C::TwistedChain) = nothing
passthrough!(C::DescendingChain,RF::RFactorizationIdentity) = nothing
passthrough!(C::AscendingChain,RF::RFactorizationIdentity) = nothing
passthrough!(C::TwistedChain,RF::RFactorizationIdentity) = nothing

simple_passthrough!(RF::RFactorizationIdentity, U::AbstractRotator) = true
simple_passthrough!(RF::RFactorizationIdentity, U::AbstractRotator, V::AbstractRotator) = true
