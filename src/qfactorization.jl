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


# Factorization for a Descending Chain
struct QFactorization{T, S, V} <: AbstractQFactorization{T, S}
  Q::DescendingChain{T, S, V}
  D::SparseDiagonal{S}
  QFactorization(Q::DescendingChain{T, S, V}, D::SparseDiagonal{S}) where {T, S, V} = new{T, S, V}(Q,D)
  QFactorization(Q::DescendingChain{T, S, V}) where {T, S, V} = QFactorization(Q, SparseDiagonal(S, length(Q)+1))
end

### QFactorization  Twisted
struct QFactorizationTwisted{T, S, Vt, PVt} <: AbstractQFactorization{T, S}
  Q::TwistedChain{T,S, Vt, PVt} ## Twisted
  D::SparseDiagonal{S}
  QFactorizationTwisted(Q::TwistedChain{T, S, Vt, PVt}, D::SparseDiagonal{S}) where  {T, S, Vt, PVt} = new{T, S, Vt, PVt}(Q,D)
  QFactorizationTwisted(Q::TwistedChain{T, S, Vt, PVt}) where {T, S, Vt, PVt} = QFactorizationTwisted(Q, SparseDiagonal{S}(length(Q)+1))

end

Base.eltype(::QFactorization{T, S, V}) where {T, S, V} = S
Base.eltype(::QFactorizationTwisted{T, S, Vt, PVt}) where {T, S, Vt, PVt} = S

# constructor
# for descending chain
function q_factorization(Q::DescendingChain{T, S, V}) where {T, S, V}

    N = length(Q)
    D = SparseDiagonal(S, N+1)

    QFactorization(Q, D)

end

function q_factorization(Q::TwistedChain{T, S, Vt, PVt}) where {T, S, Vt, PVt}

    N = length(Q)
    D = SparseDiagonal{S}(N+1) # note subtlety here, for real+twisted case we have a diagonal

    QFactorizationTwisted(Q, D)

end




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

## For Twisted, one would need to do the inefficient: Matrix(QF)[i,j]


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


function passthrough_phase!(Di::DiagonalRotator, QF::QFactorizationTwisted)
    i = idx(Di)

    asc = ascending_part(QF.Q, i)
    des = descending_part(QF.Q, i)

    passthrough_phase!(Di, (des, asc), QF.D)

end

## fuse! modifies QF and handles phase, if needed
function fuse!(QF::QFactorization{T, S}, U::AbstractRotator) where {T, S <: Real}
    i = idx(U)
    QF.Q[i], Di = fuse(QF.Q[i], U)

    return nothing
end

function fuse!(U::AbstractRotator, QF::QFactorization{T,S}) where {T, S <: Real}
    i = idx(U)
    QF.Q[i], Di = fuse(U, QF.Q[i])

    return nothing
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
