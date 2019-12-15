##################################################

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

function fuse!(Di::DiagonalRotator, D::SparseDiagonal)
    i = idx(Di)
    alpha, _ = vals(Di)
    D[i] *= alpha
    D[i+1] *= conj(alpha)
    return nothing
end
fuse!(Di::DiagonalRotator, D::IdentityDiagonal) = nothing


# real case is a series of noops
@inbounds Base.getindex(D::IdentityDiagonal{T}, k) where {T}= one(T)

Base.@propagate_inbounds Base.setindex!(D::IdentityDiagonal, X, inds...) where {T} = X

## passthrough
## Pass a rotator through a diagonal matrix with phase shifts
## D U -> U' D' (right, as in U starts on right)
## U D -> D' U' (left)
#function passthrough(D::SparseDiagonal, U::ComplexRealRotator{T}, ::Val{:right}) where {T}
#    passthrough!(D, U)
#end
@inline function passthrough!(D::SparseDiagonal, U::ComplexRealRotator{T}) where {T}
    i = idx(U)
    c,s = vals(U)

    alpha, beta = D[i], D[i+1]
    beta1 = alpha * conj(beta)

    D[i] = beta
    D[i+1] = alpha
    return Rotator(beta1 * c, s, i)
end

## U D -> D U
@inline function passthrough!(U::ComplexRealRotator{T}, D::SparseDiagonal) where {T}
    return passthrough!(D, U)
    i = idx(U)
    c,s = vals(U)

    alpha, beta = D[i], D[i+1]
    beta1 = alpha/beta

    D[i] *= conj(beta1)
    D[i+1] *= beta1
    return Rotator(beta1 * c, s, i)
end
