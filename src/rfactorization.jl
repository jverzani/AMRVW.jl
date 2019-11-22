## A factorization type for upper diagonal matrices
## T is floating point type
## S is T or Complex{T}
## Rt is rotator type, e.g. RealRotator{T}
## Pt is pencil type: Val(:no_pencil) or Val(:pencil)
abstract type AbstractRFactorization{T, Rt, Pt} end


# Basic interface for an RFactorization:
# getindex (k-2 <= j <= k)
# passthrough(RF, U, dir) to pass a rotator through from left to right or right to left
# simple_passthrough(RF, U, dir) for a shortcut passthrough
# we also might need:
# length
# eltype
# zero(RF{T, Rt}) to give 0
# one(RF{T, Rt}) to give 1

Base.length(RF::AbstractRFactorization) = length(RF.Ct)

Base.eltype(RF::AbstractRFactorization{T, RealRotator{T},Pt}) where {T,Pt} = T
Base.eltype(RF::AbstractRFactorization{T, ComplexRealRotator{T},Pt}) where {T,Pt} = Complex{T}

Base.zero(RF::AbstractRFactorization{T, RealRotator{T},Pt}) where {T,Pt} = zero(T)
Base.zero(RF::AbstractRFactorization{T, ComplexRealRotator{T},Pt}) where {T,Pt} = zero(Complex{T})

Base.one(RF::AbstractRFactorization{T, RealRotator{T},Pt}) where {T,Pt} = one(T)
Base.one(RF::AbstractRFactorization{T, ComplexRealRotator{T},Pt}) where {T,Pt} = one(Complex{T})



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
##
## Find R using Ct*B*D decompostion
## writes R as [what, wj, wk,wl] and solves
##
## function rotm(a::T,b, i, N) where {T}
##      r = diagm(0 => ones(T, N))#Matrix{T}(I, N,  N)
##      r[i:i+1, i:i+1] = [a b; -conj(b) conj(a)]
##      r
## end
## @vars bc_i bs_i bc_j bs_j bc_k bs_k
## @vars cc_i cs_i cc_j cs_j cc_k cs_k
## @vars d_i d_j d_k d_l
## D = diagm(0 => [d_i,d_j,d_k,d_l])
## Bi = rotm(bc_i, bs_i, 1, 4)
## Bj = rotm(bc_j, bs_j, 2, 4)
## Bk = rotm(bc_k, bs_k, 3, 4)
## what, wj, wk, wl = Bi * Bj * Bk * D * [0,0,1,0]
## @show wj wk wl
## @vars what wj wk wl
## u = rotm(cc_k, cs_k, 1,2) * [what, wl]
## kk = u[1](what => solve(u[2], what)[1]) |> simplify
## @show kk

## u =  rotm(cc_j, cs_j, 2,3) * rotm(cc_k, cs_k, 1,3) * [what, wk, wl]
## kk_1 = u[1](what => solve(u[2], what)[1]) |> simplify
## @show kk_1

## u =  rotm(cc_k, cs_k, 3,4) * rotm(cc_j, cs_j, 2,4) * rotm(cc_i, cs_i, 1,4) * [what, wi, wk, wl]
## kk_2 = u[1](what => solve(u[2], what)[1]) |> simplify
## @show kk_2
##
## wj = bc_k*bs_j*d_k*conjugate(bc_i)
## wk = bc_k*d_k*conjugate(bc_j)
## wl = -d_k*conjugate(bs_k)
## kk = cc_k*wl*conjugate(cc_k)/conjugate(cs_k) + cs_k*wl
## kk_1 = cc_k*wk*conjugate(cc_k)/conjugate(cs_k) + cs_k*wk + cc_k*cs_j*wl/(cc_j*conjugate(cs_k))
## kk_2 = bc_k*bs_j*cc_i*d_k*conjugate(bc_i)*conjugate(cc_i)/conjugate(cs_i) + bc_k*bs_j*cs_i*d_k*conjugate(bc_i) + cc_i*cs_j*wk/(cc_j*conjugate(cs_i))

function Rjk(Ct, B, dk, j, k)

    N = length(Ct)
    # @assert  0 <= k - j <= 2

    if k - j == 0 # (k,k) case is easiest
        bkc, bks = vals(B[k])
        ckc, cks = vals(Ct[N+1-k])

        wl = -conj(bks) * dk
        return wl / conj(cks)

    elseif k - j == 1

        j = k-1
        bkc, bks = vals(B[k])
        bjc, bjs = vals(B[j])
        ckc, cks = vals(Ct[N+1-k])
        cjc, cjs = vals(Ct[N+1-j])

        wl = -conj(bks) * dk
        wk = bkc * conj(bjc) * dk

        return wk/conj(cjs) - (wl * cjc * conj(ckc))/conj(cjs*cks)


    else# k - j == 2

        i, j = k-2, k-1

        bkc, bks = vals(B[k])
        bjc, bjs = vals(B[j])
        bic, bis = vals(B[i])
        ckc, cks = vals(Ct[N+1-k])
        cjc, cjs = vals(Ct[N+1-j])
        cic, cis = vals(Ct[N+1-i])

        wl = -conj(bks) * dk
        wk = bkc * conj(bjc) * dk
        wj = bjs * bkc * conj(bic) * dk

        return wj/conj(cis) - (wk * cic * conj(cjc)) /conj(cis*cjs) + (wl*cic * conj(ckc))/conj(cis*cjs*cks) #+ wj*cis

    end
end

## RFactorization is upper triangular, so i>j is  zero
function Base.getindex(RF::RFactorization, i, j)
    D = RF.D

    if j < i
        r_ij = zero(RF)
    elseif i + 2 >= j
        r_ij = Rjk(RF.Ct, RF.B, D[j], i, j)
    else
        ## we don't compute this, as it isn't needed
        ## when Q is Hessenberg
        r_ij = zero(RF)*NaN
    end
    r_ij
end



## Pass a rotator through Rfactorization from left or right
function passthrough(RF::RFactorization{T, St}, U::AbstractRotator, ::Val{:right}) where {T, St}

    U = passthrough(RF.D, U, Val(:right))
    U = passthrough(RF.B, U, Val(:right))
    U = passthrough(RF.Ct, U, Val(:right))

    U
end

function passthrough(RF::RFactorization{T, St}, U::AbstractRotator, ::Val{:left}) where {T, St}

    U = passthrough(RF.Ct, U, Val(:left))
    U = passthrough(RF.B, U, Val(:left))
    U = passthrough(RF.D, U, Val(:left))

    U
end



## An Identity R factorization
## used for generic purposes
struct IdentityRFactorization{T, Rt}
end


Base.getindex(RF::IdentityRFactorization, i, j) = i==j ? one(RF) : zero(RF)

passthrough(RF::IdentityRFactorization, U::AbstractRotator, dir) = U
simple_passthrough(RF::IdentityRFactorization, U::AbstractRotator, dir) = true
