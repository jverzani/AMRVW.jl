## A factorization type for upper diagonal matrices
## T is floating point type
## S is T or Complex{T}
## Rt is rotator type, e.g. RealRotator{T}
## Pt is pencil type: Val(:no_pencil) or Val(:pencil)
abstract type AbstractRFactorization{T, Rt, Pt} end


# Basic interface for an RFactorization:
#
# getindex (need just k-2 <= j <= k, but as general as possible)
# passthrough(RF, U, dir) to pass a rotator through from left to right or right to left
# simple_passthrough(RF, U, dir) for a shortcut passthrough
# we also might need:
# length
# eltype
# Matrix (in diagonostics)
# zero(RF{T, Rt}) to give 0
# one(RF{T, Rt}) to give 1

Base.length(RF::AbstractRFactorization) = length(RF.Ct)

Base.eltype(RF::AbstractRFactorization{T, RealRotator{T},Pt}) where {T,Pt} = T
Base.eltype(RF::AbstractRFactorization{T, ComplexRealRotator{T},Pt}) where {T,Pt} = Complex{T}

Base.zero(RF::AbstractRFactorization) = zero(eltype(RF))
Base.one(RF::AbstractRFactorization) = one(eltype(RF))


## R factorization can include just R or VW' [sp]
## depending if the factorization is pencil type or not

## A factorization of a non pencil type
struct RFactorization{T, Rt} <: AbstractRFactorization{T, Rt, Val{:no_pencil}}
  Ct::AscendingChain{Rt}
  B::DescendingChain{Rt}
  D::AbstractSparseDiagonalMatrix{T,Rt}
end




## xs are decompose(ps)

## Constructor for Real coefficients; no pencil
function r_factorization(xs::Vector{S}) where {S}
    N = length(xs) - 1
    Ct = AscendingChain(Vector{RotatorType(S)}(undef, N))
    B =  DescendingChain(Vector{RotatorType(S)}(undef, N))
    D = sparse_diagonal(S, N+1)

    _r_factorization(xs, Ct, B, D)

    RFactorization(Ct, B, D)
end



# constructor populating Ct, B, D for Real Case
function _r_factorization(xs::Vector{T}, Ct, B, D) where {T <: Real}
    N = length(xs) - 1
    c, s, tmp = givensrot(xs[N], -one(T))
    C = Rotator(c, s, N)
    Ct[1] = C' #Rotator(c, -s, N)
    B[N] = Rotator(s, -c, N)  # not the adjoint, so that Ct_1 * B_n = [0 -1; 1 0]

    @inbounds for i in (N-1):-1:1
        c,s,tmp = givensrot(xs[i], tmp)
        C = Rotator(c, s, i)
        Ct[N + 1 - i], B[i] = C', C
    end

    nothing
end


# constructor populating Ct, B, D for Complex Case (phase is impt)
function _r_factorization(xs::Vector{Complex{T}}, Ct, B, D) where {T <: Real}
    N = length(xs) - 1
    S = Complex{T}

    c, s, tmp = givensrot(xs[N], -one(S))

    nrm = norm(c)
    alpha = c/nrm

    C = Rotator(c, s, N)
    Ct[1] = C'
    B[N] = Rotator(s * alpha, -nrm, N)
    D[N] = conj(alpha)  # Ct*B*D gives [0 -1; 1 0] block at tail
    D[N+1] = alpha


    @inbounds for i in (N-1):-1:1
        c,s,tmp = givensrot(xs[i], tmp)
        C = Rotator(c, s, i)
        Ct[N + 1 - i], B[i] = C', C
    end


    nothing
end


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
function Base.getindex(RF::RFactorization, j, k)
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






## Pass a rotator through Rfactorization from left or right
## function passthrough(RF::RFactorization{T, St}, U::AbstractRotator, ::Val{:right}) where {T, St}
##     passthrough!(RF, U)
## end

function passthrough!(RF::RFactorization{T, St}, U::AbstractRotator) where {T, St}
    U = passthrough!(RF.D, U)
    U = passthrough!(RF.B, U)
    U = passthrough!(RF.Ct, U)

    U
end

#function passthrough(RF::RFactorization{T, St}, U::AbstractRotator, ::Val{:left}) where {T, St}
#    passthrough!(U, RF)
#end

function passthrough!(U::AbstractRotator, RF::RFactorization{T, St}) where {T, St}

    U = passthrough!(U, RF.Ct)
    U = passthrough!(U, RF.B)
    U = passthrough!(U, RF.D)

    U
end

# pass thorugh  R <- Us
function passthrough!(RF::AbstractRFactorization, Us::Vector)
    for i in  eachindex(Us)
        Us[i] =  passthrough!(RF, Us[i])
    end
end

# passthrough Us -> R
function passthrough!(Us::Vector, RF::AbstractRFactorization)
    for i in length(Us):-1:1
         Us[i] =  passthrough!(Us[i], RF)
    end
end



## An Identity R factorization
## used for generic purposes
struct IdentityRFactorization{T, Rt} <: AbstractRFactorization{T, Rt, Val{:no_pencil}}
IdentityRFactorization{T, Rt}() where {T, Rt} = new()
end



Base.getindex(RF::IdentityRFactorization, i, j) = i==j ? one(RF) : zero(RF)

#passthrough(RF::IdentityRFactorization, U::AbstractRotator, dir) = U
passthrough!(RF::IdentityRFactorization, U::AbstractRotator) = U
passthrough!(U::AbstractRotator, RF::IdentityRFactorization) = U
simple_passthrough(RF::IdentityRFactorization, U::AbstractRotator, dir) = true
simple_passthrough(RF::IdentityRFactorization, U::AbstractRotator, V::AbstractRotator, dir) = true
