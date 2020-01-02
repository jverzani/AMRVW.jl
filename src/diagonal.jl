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

passthrough!(D::IdentityDiagonal, U::AbstractRotator{T}) where {T} = U
passthrough!(U::AbstractRotator{T}, D::IdentityDiagonal) where {T} = U

# Move a chain through a diagonal
function passthrough!(D::SparseDiagonal, Asc::AscendingChain)
    for i in 1:length(Asc) # pass  Asc through D
        Asc[i] = passthrough!(D, Asc[i])
    end
end

function passthrough!(Des::DescendingChain, D::SparseDiagonal)
    for i in length(Des):-1:1
        Des[i] = passthrough!(Des[i],D)
    end
end

## could do two others...

## noop when Identity Diagonal
passthrough!(D::IdentityDiagonal, C::AbstractRotatorChain) =  nothing
passthrough!(C::AbstractRotatorChain, D::IdentityDiagonal) =  nothing

## Di U = U' Di'
passthrough!(D::IdentityDiagonal, U::AbstractRotator) = (U,D)
function passthrough!(Di::DiagonalRotator, U::AbstractRotator)
    @assert idx(Di) == idx(U)
    i = idx(Di)
    alpha, _ = vals(Di)
    c, s = vals(U)
    Rotator(c * alpha/conj(alpha), s, i), DiagonalRotator(conj(alpha), i)
end

passthrough!(U::AbstractRotator, D::IdentityDiagonal) = (U,D)
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
tip(Di::IdentityDiagonal,  Uj::AbstractRotator) = (Uj, Di)
tip(Uj::AbstractRotator, Di::IdentityDiagonal) = (Di, Uj)

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
