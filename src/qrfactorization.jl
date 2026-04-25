## In LinearAlgebra/src/qr.jl
## the QR factorization stores
## * Q in a compressed form using Householder transformations
## * R as an upper triangular matrix
##
## QRFactorization stores
## * Q in a compressed form using Givens rotations, either descending or twisted
## * R in a compressed form (RankOne, Pencil, or Identity) or as a full matrix that is upper triangular
## * this factorization gets updated by the AMRVW_algorithm!, and hence `eigvals`
##
abstract type AbstractQRRotatorFactorization  end

##  We  have  two types, though could consolidate
##
struct QRFactorization{QF, RF} <: AbstractQRRotatorFactorization
  QF::QF
  RF::RF
  end


# use identity R when not specified
function QRFactorization(QF::AbstractQFactorization{T,S}) where {T,S}
    RF = RFactorizationIdentity{T,S}()
    QRFactorization(QF, RF)
end


function Base.size(state::AbstractQRRotatorFactorization)
    m1, n1 = size(state.QF)
    m2, n2 = size(state.RF)
    return (max(m1,m2), max(n1, n2))
end

Base.copy(F::QRFactorization) = QRFactorization(copy(F.QF), copy(F.RF))
Base.length(state::QRFactorization) = length(state.QF)+1
Base.eltype(state::QRFactorization) = eltype(state.QF)

# return A = QR
function Base.Matrix(state::QRFactorization)

    Q = Matrix(state.QF)
    R = Matrix(state.RF)

    ## May have Q smaller than  R, so will pad in that case
    if R != I
        m,n =  size(Q)[1], size(R)[1]
        n  < m && error("R is  too  small")
        return Q * R[1:m, 1:m]
    else
        return Q
    end

end

"""
    qr_factorization(H::Matrix; unitary=false)

For a Hessenberg matrix `H` return a factorization,`Q⋅R`, where `Q` is a `QFactorization` object of descending rotators and `R` is either a `AMRVW.RFactorizationUnitaryDiagonal` object (if `unitary=true`) or a `RFactorizationUpperTriangular` object.

`H` may be a full matrix or the `.H` component  of a `hessenberg` factorization.

"""
function qr_factorization(H::AbstractMatrix; unitary=false)
    R = copy(H)
    qr_factorization!(R; unitary=unitary)
end

function qr_factorization!(R::AbstractMatrix{S}; unitary=false) where {S}

    n = size(R, 1)
    Qs = Vector{Rotator{real(S),S}}(undef, n-1)
    for i in 1:n-1
        c,s,r = givensrot(R[i,i], R[i+1,i])
        Ui =  Rotator(c,s,i)
        Qs[i] = Ui'
        j = i + 1
        R[i,i] = r
        R[j,i] = zero(S)
        for k in j:n
            rik, rjk =  R[i,k],  R[j,k]
            R[i,k] = c * rik + s * rjk
            R[j,k] = -conj(s) * rik + conj(c) * rjk
        end
    end

    _qr_factorization(DescendingChain(Qs), R, Val(unitary))
end

function _qr_factorization(Des::DescendingChain{T,S,V}, R, unitary::Val{true}) where {T, S, V}

    QF = QFactorization(Des)
    RF = RFactorizationUnitaryDiagonal(sign.(diag(R)))

    QRFactorization(QF, RF)
end

function _qr_factorization(Des, R, unitary::Val{false})

    QF = QFactorization(Des)
    RF = RFactorizationUpperTriangular(R)

    F = QRFactorization(QF, RF)
end
