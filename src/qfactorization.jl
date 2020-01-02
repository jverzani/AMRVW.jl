##################################################

## Q Factorization
abstract type AbstractQFactorization{T, Rt, Twt} end

# implement Array interface
Base.length(QF::AbstractQFactorization) = length(QF.Q)

Base.eltype(QF::AbstractQFactorization{T, RealRotator{T},Twt}) where {T,Twt} = T
Base.eltype(QF::AbstractQFactorization{T, ComplexRealRotator{T},Twt}) where {T,Twt} = Complex{T}

Base.zero(QF::AbstractQFactorization) = zero(eltype(QF))
Base.one(QF::AbstractQFactorization) = one(eltype(QF))


struct QFactorization{T, Rt} <: AbstractQFactorization{T, Rt, Val{:not_twisted}}
  Q::DescendingChain{Rt}
  D::AbstractSparseDiagonalMatrix{T, Rt}
  W::Vector{Rt}
end

## Find Q[j,k] from its factored form
## QFactorization is Hessenber
## We will only need near  diagonal elements, as we multiply by
## an  upper triangular matrix

function Base.getindex(QF::AbstractQFactorization, j, k)
    (j <= 0 || k <= 0) && return zero(QF)
    QF.Q[j,k] * QF.D[k]
end

## function Base_getindex(QF::QFactorization, j, k)

##     Q = QF.Q

##     if k == 0
##         return zero(QF)
##     end

##     ## We need to compute QR[j:k, j:k]
##     ## for this we use Q is Hessenberg, R if triangular
##     ## so we only need
##     ##
##     ## [ qji qjj qjk
##     ##    0  qkj qkk]

##     Δ = k - j
##     i, j = k-2, k-1

##     if k <= length(QF)
##         ck, sk =  vals(Q[k])
##     else
##         ck, sk = one(QF), real(zero(QF))
##     end

##     if Δ < -1
##         # Hessenberg
##         return zero(QF)

##     elseif Δ == -1 # e,g, qjk this count is off, as k < j

##         dk = QF.D[k]
##         return -sk * dk

##     elseif Δ == 0 # gkk case, gjj

##         cj, sj = j >= 1 ? vals(Q[j]) : (one(QF), real(zero(QF)))
##         dk = QF.D[k]
##         return ck  * conj(cj) * dk

##     elseif Δ == 1 # qjk case

##         cj, sj = vals(Q[j])
##         ci, si =  i >= 1 ? vals(Q[i]) : (one(QF), real(zero(QF)))
##         dk = QF.D[k]

##         return ck * dk * sj * conj(ci)

##     else # Δ > 1

##         ## we don't need this as we multiply by a triangular matrix
##         return zero(QF) * NaN
##     end


## end




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
function passthrough!(Di::DiagonalRotator, QF::QFactorization)
    passthrough_phase!(Di, QF.Q, QF.D)
end

## fuse! modifies QF and returns Di is needed.
function fuse!(QF::QFactorization{T, RealRotator{T}}, U::AbstractRotator) where {T}
    i = idx(U)
    QF.Q[i] = fuse(QF.Q[i], U)
end

function fuse!(U::AbstractRotator, QF::QFactorization{T,RealRotator{T}}) where {T}
    i = idx(U)
    QF.Q[i] = fuse(U, QF.Q[i])
end

function fuse!(QF::QFactorization{T, ComplexRealRotator{T}}, U::AbstractRotator) where {T}
    i = idx(U)
    QF.Q[i], Di = fuse(QF.Q[i], U)
    Di

    alpha, _ = vals(Di)
    @inbounds QF.D[i] = QF.D[i] * alpha
    @inbounds QF.D[i+1] = QF.D[i+1] * conj(alpha)

    return nothing
end

function fuse!(U::AbstractRotator, QF::QFactorization{T, ComplexRealRotator{T}}) where {T}
    i = idx(U)
    QF.Q[i], Di = fuse(U, QF.Q[i])
    passthrough!(Di, QF)

    return nothing
end
