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
