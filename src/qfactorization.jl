##################################################

## Q Factorization
abstract type AbstractQFactorization{T, S} end
# implement Array interface
Base.length(QF::AbstractQFactorization) = length(QF.Q)
Base.eltype(QF::AbstractQFactorization{T, S}) where {T,S} = S

Base.zero(QF::AbstractQFactorization) = zero(eltype(QF))
Base.one(QF::AbstractQFactorization) = one(eltype(QF))


struct QFactorization{T, S} <: AbstractQFactorization{T, S}
  Q::DescendingChain{T,S}
  D::SparseDiagonal{S}
end


# constructor
function q_factorization(Qs::AbstractArray{T,S}) where {T, S}
    N = length(Qs)
    D = SparseDiagonal(S, N+1)
    QFactorization{T,S}(DescendingChain(Qs), D)
end

function q_factorization(Qs::AbstractArray{T,S}, pv::Vector) where {T, S}
    N = length(Qs)
    D = SparseDiagonal(S, N+1)
    QFactorization{T,S}(TwistedChain(Qs, pv), D)
end



## Find Q[j,k] from its factored form
## QFactorization is Hessenber
## We will only need near  diagonal elements, as we multiply by
## an  upper triangular matrix

function Base.getindex(QF::QFactorization, j, k)
    (j <= 0 || k <= 0) && return zero(QF)
    QF.Q[j,k] * QF.D[k]
end


function Base.Matrix(QF::QFactorization{T, S}) where {T, S}
    n = length(QF) + 1
    M = diagm(0 => ones(S, n))
    QF.Q * (Matrix(QF.D) * M)
end




##################################################
##
## Transformations
##
function passthrough!(QF::QFactorization, U::AbstractRotator)

    UU = passthrough!(QF.D, U)
    UUU = passthrough!(QF.Q, UU)
    UUU

end

function passthrough!(U::AbstractRotator, QF::QFactorization)

    U = passthrough!(U, QF.Q)
    passthrough!(U, QF.D)

end

# pass diagonal rotator through and merge into D
function passthrough_phase!(Di::DiagonalRotator, QF::QFactorization)
    passthrough_phase!(Di, QF.Q, QF.D)
end

## fuse! modifies QF and returns Di is needed.
function fuse!(QF::QFactorization{T, S}, U::AbstractRotator) where {T, S <: Real}
    i = idx(U)
    QF.Q[i], Di = fuse(QF.Q[i], U)
end

function fuse!(U::AbstractRotator, QF::QFactorization{T,S}) where {T, S <: Real}
    i = idx(U)
    QF.Q[i], Di = fuse(U, QF.Q[i])
end

function fuse!(QF::QFactorization{T, S}, U::AbstractRotator) where {T, S <: Complex}
    i = idx(U)
    QF.Q[i], Di = fuse(QF.Q[i], U)
    Di

    alpha, _ = vals(Di)
    @inbounds QF.D[i] = QF.D[i] * alpha
    @inbounds QF.D[i+1] = QF.D[i+1] * conj(alpha)

    return nothing
end

function fuse!(U::AbstractRotator, QF::QFactorization{T,S}) where {T, S <: Complex}
    i = idx(U)
    QF.Q[i], Di = fuse(U, QF.Q[i])
    passthrough_phase!(Di, QF)

    return nothing
end
