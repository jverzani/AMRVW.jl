## Types for factorization objects



abstract type AbstractQRFactorizationState end

function Base.size(state::AbstractQRFactorizationState)
    m1, n1 = size(state.QF)
    m2, n2 = size(state.RF)
    return (max(m1,m2), max(n1, n2))
end



##  We  have  two types, though could consolidate
##
struct QRFactorization{QF, RF} <: AbstractQRFactorizationState
  QF::QF
  RF::RF
end

# use identity R when not specified
function QRFactorization(QF::AbstractQFactorization{T,S}) where {T,S}
    RF = RFactorizationIdentity{T,S}()
    QRFactorization(QF, RF)
end


## XXX  remove me
## struct QRFactorizationTwisted{T, S, Vt, PVt, RF} <: AbstractQRFactorizationState
##   QF::QFactorizationTwisted{T, S, Vt, PVt}
##   RF::RF
## end

Base.length(state::AbstractQRFactorizationState) = length(state.QF)+1
Base.eltype(state::AbstractQRFactorizationState) = eltype(state.QF)

# return A
function Base.Matrix(state::AbstractQRFactorizationState)

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
