## Types for factorization objects

## A container for our counters
mutable struct AMRVW_Counter
    zero_index::Int
    start_index::Int
    stop_index::Int
    it_count::Int
    tr::Int
end


##################################################
## Factorization Types ##
##################################################
## A type to hold our factorization state (Q, R, temp storage)
## T is floating point type
## S is T or Complex{T}
## Rt is rotator type, e.g. RealRotator{T}
## Pt is pencil type: Val(:no_pencil) or Val(:pencil)
## Twt is twist type: Val(:twisted) or Val(:not_twisted)
abstract type AbstractFactorizationState{T, S, Rt, QFt, RFt, Twt} end

# Factorization for No_Pencil and Not Twisted
# N size of poly
# QF, RF factorizations
# UV reusable storage for U and [V] in real cose
# A reusable storage for pieces of the full matrix
# REIGS, IEIGS storage for eigen values
struct QRFactorization{T, S, Rt, QFt, RFt} <: AbstractFactorizationState{T, S, Rt, QFt, RFt, Val{:not_twisted}}
  N::Int
  QF::QFt
  RF::RFt
  UV::Vector{Rt}    # temp storage for U[V] Vector{Rt}, ,
  A::Matrix{S}
  REIGS::Vector{T}
  IEIGS::Vector{T}
  ctrs::AMRVW_Counter
end



# constructor for either case
function qrfactorization(N,
                         QF::AbstractQFactorization{T, Rt, Twt},
                         RF::AbstractRFactorization) where {T, Rt, Twt}


    if Rt == ComplexRealRotator{T}
        S = Complex{T}
        UV = Vector{ComplexRealRotator{T}}(undef, 1)
    else
        S = T
        UV = Vector{RealRotator{T}}(undef, 2)
    end

    A = zeros(S, 2, 2)
    reigs = zeros(T, N)
    ieigs = zeros(T, N)
    ctr = AMRVW_Counter(0, 1, N-1, 0, N-2)

    if Twt == Val{:not_twisted}
        state = QRFactorization(N, QF, RF, UV, A, reigs, ieigs, ctr)
    else
        state = QRFactorizationTwisted(N, QF, RF, UV, A, reigs, ieigs, ctr)
    end

    return state
end

Base.length(state::AbstractFactorizationState) = state.N


# return A
function Base.Matrix(state::AbstractFactorizationState)

    Q = Matrix(state.QF)
    R = Matrix(state.RF)

    ## May have Q smaller than  R, so will pad in that case
    if R != I
        m,n =  size(Q)[1], size(R)[1]
        if m < n
            QQ = diagm(0 => ones(eltype(Q), n))
            QQ[1:m, 1:m] .= Q
            return  QQ  *  R
        else
            return Q * R
        end
    else
        return Q
    end

end
