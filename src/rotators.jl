
## Rotators Our rotators have field names c, s where c nad s are
## either T or Complex{T} It proved to be faster to have immutable
## rotators, rather than mutable ones so there are no "setter" methods
## anymore, a new instance (via `Rotator(c,s,i)`) should be used

abstract type CoreTransform{T} end
abstract type AbstractRotator{T} <: CoreTransform{T} end

## our rotators are [conj(c) s; -s c]
## get values
@inline vals(r::AbstractRotator{T}) where {T} = (r.c, r.s)
@inline idx(r::AbstractRotator) = r.i


Base.copy(a::AbstractRotator) = AbstractRotator(a.c, a.s, a.i)
is_diagonal(r::AbstractRotator{T}) where {T} = norm(r.s) <= eps(T)

## XXX
## this uses [c s; -conj(s) conj(c)] for rotator!
function *(a::AbstractRotator, M::AbstractArray)
    c, s = vals(a)
    i = idx(a); j = i+1
    N = copy(M)
    N[i,:]  =  round.(c * M[i,:] + s * M[j,:], digits=16)
    N[j,:]  =   round.(-conj(s) * M[i,:] + conj(c) * M[j,:], digits=16)
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
struct RealRotator{T} <: AbstractRotator{T}
c::T
s::T
i::Int
RealRotator(c::T, s::T, i::Int) where {T} = new{T}(c,s,i)
RealRotator{T}(c,s,i) where {T} = new(convert(T,c),convert(T,s),i)
end

function adjoint(r::RealRotator)
    RealRotator(r.c, -r.s, r.i)
end

Base.eltype(::Type{RealRotator{T}}) where {T} = T

Base.one(::Type{RealRotator{T}}) where {T} = RealRotator(one(T), zero(T), 0)# end

##################################################
### Okay, now try with complex c, real s

struct ComplexRealRotator{T} <: AbstractRotator{T}
c::Complex{T}
s::T
i::Int
ComplexRealRotator(c::Complex{T}, s::T, i::Int) where {T} = new{T}(c,s,i)
ComplexRealRotator{T}(c,s,i) where {T} = new(convert(Complex{T},c),convert(T,s),i)
end

function adjoint(r::ComplexRealRotator)
    ComplexRealRotator(conj(r.c), -r.s, r.i)
end

Base.eltype(::Type{ComplexRealRotator{T}}) where {T} = complex(T)

Base.one(::Type{ComplexRealRotator{T}}) where {T} = ComplexRealRotator(complex(one(T), zero(T)), zero(T), 0)

Base.copy(a::ComplexRealRotator) = ComplexRealRotator(a.c, a.s, a.i)




## Easier constructor

Rotator(c::Complex{T}, s::Complex{T}, i::Int) where {T <: Real} = ComplexComplexRotator(c, s, i)
Rotator(c::Complex{T}, s::T, i::Int) where {T <: Real} = ComplexRealRotator(c,s,i)
Rotator(c::T, s::T, i::Int) where {T <: Real} = RealRotator(c,s,i)

## rotator type from a numeric type
RotatorType(::Type{S}) where {S <: Complex} = ComplexRealRotator{real(S)}
RotatorType(::Type{T}) where {T <: Real} = RealRotator{T}

## A diagonal rotator [c 0; 0 c] witht |c| = 1
struct DiagonalRotator{T} <: AbstractRotator{T}
c::Complex{T}
i::Int
end

vals(D::DiagonalRotator{T}) where {T} = D.c, zero(real(T))

struct IdentityDiagonalRotator{T} <: AbstractRotator{T}
i::Int
IdentityDiagonalRotator{T}(i) where {T} = new(i)
end

vals(D::IdentityDiagonalRotator{T}) where {T} = (one(T), zero(T))
idx(D::IdentityDiagonalRotator{T}) where {T} = D.i
