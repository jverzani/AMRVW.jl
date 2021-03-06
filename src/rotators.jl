
## Rotators Our rotators have field names c, s where c nad s are
## either T or Complex{T} It proved to be faster to have immutable
## rotators, rather than mutable ones so there are no "setter" methods
## anymore, a new instance (via `Rotator(c,s,i)`) should be used

abstract type CoreTransform{T,S} end
abstract type AbstractRotator{T,S} <: CoreTransform{T,S} end

## our rotators are [c s; -s conj(c)]
## We *could* just directly use `LinearAlgebra.Givens` here, but
## instead implement Rotator below to ensure the types Real/Real and Complex/Real for c and s
## XXX: for twisted, it would have been easier to have Complex/Complex rotators, as then
## `passthrough_phase` could have been dramatically simplified.
##
## Interface
## vals: return c,s
## idx: return i
## adjoint: return adjoint

## get values
@inline vals(r::AbstractRotator{T}) where {T} = (r.c, r.s)
@inline idx(r::AbstractRotator) = r.i


Base.copy(a::AbstractRotator) = AbstractRotator(a.c, a.s, a.i)
is_diagonal(r::AbstractRotator{T,S}) where {T,S} = norm(r.s) <= eps(T)

## Multiplication of a rotator and a matrix
## this uses [c s; -conj(s) conj(c)] for rotator!
function *(a::AbstractRotator, M::AbstractArray)
    N = copy(M)
    lmul!(a, N)
    N
end

function LinearAlgebra.lmul!(a::AbstractRotator, M::AbstractArray)
    c, s = vals(a)
    i = idx(a); j = i+1
    n = size(M)[2]
    for k in 1:n
        mik, mjk = M[i,k], M[j,k]
        M[i,k]  =  c * mik + s * mjk
        M[j,k]  =  -conj(s) * mik + conj(c) * mjk
    end
    M
end

function *(M::AbstractArray, a::AbstractRotator)
    N = copy(M)
    rmul!(N, a)
    N
end

function LinearAlgebra.rmul!(M::AbstractArray, a::AbstractRotator)
    c, s = vals(a)
    i = idx(a); j = i+1
    n = size(M)[1]
    for k in 1:n
        mki, mkj = M[k,i], M[k,j]
        M[k, i] = c * mki - conj(s) * mkj
        M[k, j] = s * mki + conj(c) * mkj
    end
    M
end

*(Qs::Vector{R}, M::Matrix) where {R <: CoreTransform} = foldr(*, Qs, init=M)
*(M::Matrix, Qs::Vector{R}) where {R <: CoreTransform} = foldl(*, Qs, init=M)
function LinearAlgebra.lmul!(Qs::Vector{R}, M::Matrix) where {R <: CoreTransform}
    foldr(lmul!, Qs, init=M)
end
function LinearAlgebra.rmul!(M::Matrix, Qs::Vector{R}) where {R <: CoreTransform}
    foldl(rmul!, Qs, init=M)
end

# A rotator
# c is real or complex
# s is always real. (The Complex/Complex case was dropped)
# i the index, adds storage space, but simplifies some algorithms
struct Rotator{T,S} <: AbstractRotator{T,S}
c::S
s::T
i::Int
end

function LinearAlgebra.adjoint(U::Rotator)
    c,s = vals(U)
    Rotator(conj(c), -conj(s), idx(U))
end

Base.eltype(::Type{Rotator{T,S}}) where {T, S} = S
Base.one(::Type{Rotator{T,S}}) where {T, S} = Rotator(one(S), zero(T), 0)
Base.copy(U::Rotator) = Rotator(vals(U)..., idx(U))

## A diagonal rotator has iszero(s) being true
struct DiagonalRotator{T,S} <: AbstractRotator{S, S}
c::S
i::Int
DiagonalRotator{T,S}(c,i) where {T,S} = new(convert(S,c),i)
DiagonalRotator{T,S}(c,s,i) where {T,S} = new(convert(S,c),i)
DiagonalRotator{S}(c,i) where {S} = DiagonalRotator{real(S),S}(c,i)
DiagonalRotator(c::S, i) where {S} = DiagonalRotator{real(S),S}(c,i)
end

vals(D::DiagonalRotator{T,S}) where {T,S} = D.c, zero(T)
LinearAlgebra.adjoint(U::DiagonalRotator) = DiagaonalRotator(conj(U.c), U.i)

## for real case, we have identity diagonal
struct IdentityRotator{T,S} <: AbstractRotator{T,S}
i::Int
IdentityRotator{T,S}(i) where {T,S} = new(i)
end

vals(D::IdentityRotator{T,S}) where {T,S} = (one(S), zero(T))
idx(D::IdentityRotator{T,S}) where {T,S} = D.i
LinearAlgebra.adjoint(U::IdentityRotator) = U
