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
# zero(RF{T, Rt}) to give 0
# one(RF{T, Rt}) to give 1

Base.length(RF::AbstractRFactorization) = length(RF.Ct)

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

function Base.getindex(RF::RFactorization, i, j)
    D = RF.D
    @assert j-2 <= i
    if j < i
        r_ij = zero(RF)
    else
        r_ij = Rjk(RF.Ct, RF.B, D[j], i, j)
    end
    r_ij
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
