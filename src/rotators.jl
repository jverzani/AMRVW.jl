
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


#the index is superflous for now, and a bit of a hassle to keep immutable
#but might be of help later if twisting is approached. Shouldn't effect speed, but does mean 3N storage (Q, Ct, B)
#so may be
#
struct RealRotator{T} <: AbstractRotator{T}
c::T
s::T
i::Int
RealRotator(c::T, s::T, i::Int) where {T} = new{T}(c,s,i)
RealRotator{T}() where {T} = new()
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

##################################################
struct DiagonalRotator{T} <: AbstractRotator{T}
c::T
i::Int
end

vals(D::DiagonalRotator) = D.c, zero(D.c)


#### Group rotators

abstract type AbstractRotatorChain{T} end

Base.length(A::AbstractRotatorChain) = length(A.x)

Base.@propagate_inbounds Base.getindex(A::AbstractRotatorChain, i::Int) = getindex(A.x, i)
Base.@propagate_inbounds Base.setindex!(A::AbstractRotatorChain, X, inds...) = setindex!(A.x, X, inds...)


Base.iterate(A::AbstractRotatorChain) = iterate(A.x)
Base.iterate(A::AbstractRotatorChain, st) = iterate(A.x, st)


struct DescendingChain{T} <: AbstractRotatorChain{T}
  x::Vector{T}
end

struct AscendingChain{T} <: AbstractRotatorChain{T}
  x::Vector{T}
end

struct TwistedChain{T} <: AbstractRotatorChain{T}
   x::Vector{T}
end

function Base.size(C::AbstractRotatorChain)
    N = length(C.x)
    (N+1, N+1)
end

Base.adjoint(A::AscendingChain) = DescendingChain(reverse(adjoint.(A.x)))
Base.adjoint(A::DescendingChain) = AscendingChain(reverse(adjoint.(A.x)))
Base.adjoint(A::TwistedChain) = TwistedChain(reverse(adjoint.(A.x)))



function Base.getindex(A::AscendingChain{T}, i, j) where {T}
    N = length(A.x)
    S = eltype(first(A.x))

    if i == N + 1
        i -= 1; j -= 1
        ci, si = vals(A.x[N+1-i])

        j == N && return conj(ci)

        s = -S(si)
        for l in (j+1):(i-1)
            cj, sjj = vals(A.x[N+1-l])
            s *= -sjj
        end
        cj, sj = j > 0 ? vals(A.x[N+1-j]) : (one(S), zero(S))
        return s * conj(cj)

    else


        j > i+1 && return zero(S) # ascending chain's aree lower Hessian

        ci, si = vals(A.x[N+1-i])

        if j == i + 1
            return S(si)
        elseif j == i
            cj, sj  = j > 1 ? vals(A.x[N+1-(j-1)]) : (one(S), zeros(S))
            return conj(cj)*ci
        end

        s = one(S)
        for l in j:(i-1)
            cl, sl = vals(A.x[N+1-l])
            s *= -sl
        end

        cl, sl = j > 1 ? vals(A.x[N+1 - (j-1)]) : (one(S), zero(S))
        return ci * s * conj(cl)
    end

end

function Base.getindex(A::DescendingChain, i, j)

    N = length(A.x)
    S = eltype(first(A.x))

    if i > j + 1
        return zero(S)
    end

    if j == N + 1
        j -= 1
        cj, sj = vals(A.x[j])

        i == N+1 && return conj(cj)

        i -= 1
        s::S = sj
        for l in (i+1):j-1
            cj, sj = vals(A.x[l])
            s *= S(sj)
        end
        ci, si = i > 0 ? vals(A.x[i]) : (one(S), zero(real(S)))
        return s * conj(ci)

    else


        cj, sj = vals(A.x[j])

        if i == j + 1

            return -S(sj)

        elseif j == i

            ci, si  = i >= 2 ? vals(A.x[i-1]) : (one(S), zero(real(S)))
            return conj(ci)*cj

        end

        s = one(S)
        for l in i:(j-1)
            cl, sl = vals(A.x[l])
            s *= sl
        end

        cl, sl = i >= 2 ? vals(A.x[i-1]) : (one(S), zero(real(S)))
        return cj * s * conj(cl)
    end
end


## For the ComplexRealRotator case we need a diagonal matrix to store the phase
## This is a  means to make this generic with respect to rotator type
## The factorization can have an identity diagonal or a real one.
abstract type AbstractSparseDiagonalMatrix{T, Rt} end

struct IdentityDiagonal{T} <: AbstractSparseDiagonalMatrix{T, RealRotator{T}}
IdentityDiagonal{T}() where {T} = new()
end

struct SparseDiagonal{T} <:  AbstractSparseDiagonalMatrix{T, ComplexRealRotator{T}}
x::Vector{Complex{T}}
SparseDiagonal{T}(x::Vector) where {T} = new(x)
SparseDiagonal{T}(N::Int) where {T} = new(ones(Complex{T}, N))
end

sparse_diagonal(::Type{T}, N) where {T <: Real} = IdentityDiagonal{T}()
sparse_diagonal(::Type{S}, N) where {S} = SparseDiagonal{real(S)}(N)

# complex case
@inbounds Base.getindex(D::SparseDiagonal, k) = D.x[k]

Base.@propagate_inbounds Base.setindex!(D::SparseDiagonal, X, inds...) = setindex!(D.x, X, inds...)


# real case is a series of noops
@inbounds Base.getindex(D::IdentityDiagonal{T}, k) where {T}= one(T)

Base.@propagate_inbounds Base.setindex!(D::IdentityDiagonal, X, inds...) where {T} = X
