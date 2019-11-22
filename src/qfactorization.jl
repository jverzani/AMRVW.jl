##################################################

## Q Factorization
abstract type AbstractQFactorization{T, Rt, Twt} end
Base.length(QF::AbstractQFactorization) = length(QF.Q)

Base.zero(RF::AbstractQFactorization{T, RealRotator{T},Twt}) where {T,Twt} = zero(T)
Base.zero(RF::AbstractQFactorization{T, ComplexRealRotator{T},Twt}) where {T,Twt} = zero(Complex{T})
Base.one(RF::AbstractQFactorization{T, RealRotator{T},Twt}) where {T,Twt} = one(T)
Base.one(RF::AbstractQFactorization{T, ComplexRealRotator{T},Twt}) where {T,Twt} = one(Complex{T})



struct QFactorization{T, Rt} <: AbstractQFactorization{T, Rt, Val{:not_twisted}}
  Q::DescendingChain{Rt}
  D::AbstractSparseDiagonalMatrix{T, Rt}
  W::Vector{Rt}
end

## return q[i:k, j:k]
function getQ(QF::QFactorization, k)

    Q = QF.Q

    i, j = k-2, k-1
    cj, sj = vals(Q[j])
    ck, sk = vals(Q[k])
    ci, si = (k > 2) ? vals(Q[i]) : (one(cj), zero(sj))

    dj, dk = QF.D[j], QF.D[k]
    di = (k > 2) ? QF.D[i] : one(dj)

    qji = -si * di
    qkj = -sj * dj
    qjj =  cj  * conj(ci) * dj
    qjk = ck * dk * sj * conj(ci)
    qkk = ck  * conj(cj) * dk

    (qji, qjj, qjk, qkj, qkk)
end

function Base.getindex(QF::QFactorization, j,k)

    Q = QF.Q

    if k == 0
        return zero(QF)
    end

    ck, sk = vals(Q[k])

    ## We need to compute QR[j:k, j:k]
    ## for this we use Q is Hessenberg, R if triangular
    ## so we only need
    ##
    ## [ qji qjj qjk
    ##    0  qkj qkk]

    Δ = k - j
    i, j = k-2, k-1
    ## Q is Hessenberg, so we only need
    if Δ > 1  # qik case

        return zero(ck)

    elseif Δ == 1 # qjk case

        cj, sj = vals(Q[j])
        ci, si =  i >= 1 ? vals(Q[i]) : (one(QF), real(zero(QF)))
        dk = QF.D[k]
        return ck * dk * sj * conj(ci)

    elseif Δ == 0 # gkk case, gjj

        cj, sj = j >= 1 ? vals(Q[j]) : (one(QF), real(zero(QF)))
        dk = QF.D[k]
        return ck  * conj(cj) * dk

    elseif Δ == -1 # e,g, qjk this count is off, as k < j

        dk = QF.D[k]
        return -sk * dk

    elseif Δ < -1

        return zero(ck)

    end





    ## i, j = k-2, k-1
    ## cj, sj = vals(Q[j])
    ## ck, sk = vals(Q[k])
    ## ci, si = (k > 2) ? vals(Q[i]) : (one(cj), zero(sj))

    ## dj, dk = QF.D[j], QF.D[k]
    ## di = (k > 2) ? QF.D[i] : one(dj)

    ## qji = -si * di
    ## qkj = -sj * dj
    ## qjj =  cj  * conj(ci) * dj
    ## qjk = ck * dk * sj * conj(ci)
    ## qkk = ck  * conj(cj) * dk

    ## ##
    ## # [ qji qjj qjk
    ## #    0  qkj qkk]

    ## cj, sj = vals(Q[k-1])
    ## if j == k # kk case
    ##     ck, sk = vals(Q[k])
    ##     dk = QF.D[k]
    ##     return ck * conj(cj) * dk
    ## elseif j - 1 == k # qkj case
    ##     dj = QF.D[k-1]
    ##     return -sj * dj
    ## elseif j + 1 == k # qjk case
    ##     if k == 2
    ##         return zero(cj)
    ##     else
    ##         ck, sk = vals(Q[k])
    ##         ci, si = vals(Q[k-2])
    ##         dk = QF.D[k]
    ##         return ck * dk * sj * conj(ci)
    ##     end
    ## end


end


##################################################
##
## factor
##
## xs are decompose(ps)

## We need zero and one to match T,Complex{T} or just T,T depending
function _zero_one(xs::Vector{S}) where {S}
    T = real(S)
    zero(T), one(T), zero(S), one(S)
end

function q_factorization(xs::Vector{S}) where {S}
    N = length(xs) - 1

    Q =  DescendingChain(Vector{RotatorType(S)}(undef, N))
    zt,ot,zs,os = _zero_one(xs)

    @inbounds for ii = 1:(N-1)
        Q[ii] = Rotator(zs, ot, ii)
    end
    Q[N] = Rotator(os, zt, N)   ## for convenience

    D = sparse_diagonal(S, N+1)
    W = Rotator(zs, ot, 1)  #
    QFactorization(Q, D, [W])
end



##################################################
##
## Transformations
##
function passthrough(QF::QFactorization, U::AbstractRotator, dir::Val{:right})

    U = passthrough(QF.D, U, Val(:right))
    passthrough(QF.Q, U, Val(:right))

end

function passthrough(QF::QFactorization, U::AbstractRotator, dir::Val{:left})

    U = passthrough(QF.Q, U, Val(:left))
    passthrough(QF.D, U, Val(:left))

end

# pass diagonal rotator through and merge into D
function passthrough(QF::QFactorization, U::DiagonalRotator, dir::Val{:left})
    passthrough(QF.Q, QF.D, U, Val(:left))
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
    passthrough(QF, Di, Val(:left))

    return nothing
end
