
## Rotators Our rotators have field names c, s where c nad s are
## either T or Complex{T} It proved to be faster to have immutable
## rotators, rather than mutable ones so there are no "setter" methods
## anymore, a new instance (via `Rotator(c,s,i)`) should be used

abstract type CoreTransform{T,S} end
abstract type AbstractRotator{T,S} <: CoreTransform{T,S} end

## our rotators are [conj(c) s; -s c]
## Interface
## vals: return c,s
## idx: return i
## adjoint: return adjoing

## get values
@inline vals(r::AbstractRotator{T}) where {T} = (r.c, r.s)
@inline idx(r::AbstractRotator) = r.i


Base.copy(a::AbstractRotator) = AbstractRotator(a.c, a.s, a.i)
is_diagonal(r::AbstractRotator{T,S}) where {T,S} = norm(r.s) <= eps(T)

## XXX
## this uses [c s; -conj(s) conj(c)] for rotator!
function *(a::AbstractRotator, M::AbstractArray)
    c, s = vals(a)
    i = idx(a); j = i+1
    N = copy(M)
    N[i,:]  =  c * M[i,:] + s * M[j,:]
    N[j,:]  =  -conj(s) * M[i,:] + conj(c) * M[j,:]
    N
end

function *(M::AbstractArray, a::AbstractRotator)
    c, s = vals(a)
    i = idx(a); j = i+1
     N = copy(M)
    N[:, i] = c * M[:,i] - conj(s) * M[:,j]
    N[:, j] = s * M[:,i] + conj(c) * M[:,j]
    N
end
*(Qs::Vector{R}, M::Matrix) where {R <: CoreTransform} = foldr(*, Qs, init=M)
*(M::Matrix, Qs::Vector{R}) where {R <: CoreTransform} = foldl(*, Qs, init=M)

#the index is superflous for now, and a bit of a hassle to keep immutable
#but might be of help later if twisting is approached. Shouldn't effect speed, but does mean 3N storage (Q, Ct, B)
#so may be
#

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

struct DiagonalRotator{T,S} <: AbstractRotator{S, S}
c::S
i::Int
DiagonalRotator{T,S}(c,i) where {T,S} = new(convert(S,c),i)
DiagonalRotator{S}(c,i) where {S} = DiagonalRotator{real(S),S}(c,i)
DiagonalRotator(c::S, i) where {S} = DiagonalRotator{real(S),S}(c,i)
end

vals(D::DiagonalRotator{T,S}) where {T,S} = D.c, zero(real(S))
LinearAlgebra.adjoint(U::DiagonalRotator) = DiagaonalRotator(conj(U.c), U.i)

## for real case, we have identity diagonal
struct IdentityRotator{T,S} <: AbstractRotator{T,S}
i::Int
IdentityRotator{T,S}(i) where {T,S} = new(i)
end

vals(D::IdentityRotator{T,S}) where {T,S} = (one(S), zero(T))
idx(D::IdentityRotator{T,S}) where {T,S} = D.i
LinearAlgebra.adjoint(U::IdentityRotator) = U
