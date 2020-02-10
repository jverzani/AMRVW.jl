##################################################

## Q Factorization
## We have two  different Q factorizations: one with a Descending chaing,  one with a  Twisted Chain
abstract type AbstractQFactorization{T, S} end
# implement Array interface
Base.length(QF::AbstractQFactorization) = length(QF.Q)
Base.size(QF::AbstractQFactorization) = (length(QF)+1, length(QF)+1)
Base.eltype(QF::AbstractQFactorization{T, S}) where {T,S} = S

Base.zero(QF::AbstractQFactorization) = zero(eltype(QF))
Base.one(QF::AbstractQFactorization) = one(eltype(QF))


# Factorization  for a   Descending Chain
struct QFactorization{T, S, V} <: AbstractQFactorization{T, S}
  Q::DescendingChain{T, S, V}
  D::SparseDiagonal{S}
end

### QFactorization  Twisted
struct QFactorizationTwisted{T, S, Vt, PVt} <: AbstractQFactorization{T, S}
  Q::TwistedChain{T,S, Vt, PVt} ## Twisted
  D::SparseDiagonal{S}
end

Base.eltype(::QFactorization{T, S, V}) where {T, S, V} = S
Base.eltype(::QFactorizationTwisted{T, S, Vt, PVt}) where {T, S, Vt, PVt} = S

# constructor
#  for descending chain
function q_factorization(Q::DescendingChain{T, S, V}) where {T, S, V}

    N = length(Q)
    D = SparseDiagonal(S, N+1)

    QFactorization(Q, D)

end

function q_factorization(Q::TwistedChain{T, S, Vt, PVt}) where {T, S, Vt, PVt}

    N = length(Q)
    D = SparseDiagonal{S}(N+1)

    QFactorizationTwisted(Q, D)

end

## function q_factorization(Qs::AbstractArray{T,S}) where {T, S}
##     N = length(Qs)
##     D = SparseDiagonal(S, N+1)
##     QFactorization{T,S}(DescendingChain(Qs), D)
## end

## #  for twisted chain
## function q_factorization(Qs::AbstractArray{T,S}, pv::Vector) where {T, S}
##     N = length(Qs)
##     D = SparseDiagonal(S, N+1)
##     QFactorization{T,S}(TwistedChain(Qs, pv), D)
## end



## Find Q[j,k] from its factored form
## QFactorization is Hessenberg
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

    U = passthrough!(QF.D, U)
    passthrough!(QF.Q, U)

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
