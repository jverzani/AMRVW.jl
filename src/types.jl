## Types for factorization objects

## A container for our counters
mutable struct AMRVW_Counter
    zero_index::Int  # this should be dropped, it is just start_index-1
    start_index::Int
    stop_index::Int
    it_count::Int
    tr::Int
end


##################################################


abstract type AbstractQRFactorizationState{T, S, Twt} end
function Base.size(state::AbstractQRFactorizationState)
    m1, n1 = size(state.QF)
    m2, n2 = size(state.RF)
    return (max(m1,m2), max(n1, n2))
end


struct QRFactorization{T, S, V, Rt<:AbstractRFactorization{T,S}} <: AbstractQRFactorizationState{T, S, Val{:not_twisted}}
  N::Int
  QF::QFactorization{T, S, V}
  RF::Rt
  UV::Vector{Rotator{T,S}}    # Ascending chain for cfreating bulge
  W::Vector{Rotator{T,S}}     # the limb when m > 1
  A::Matrix{S}
  REIGS::Vector{T}
  IEIGS::Vector{T}
  ctrs::AMRVW_Counter
  QRFactorization{T,S,V, Rt}(N, QF, RF, UV, W, A, REIGS, IEIGS, ctrs) where {T, S, V, Rt} = new(N, QF, RF, UV, W, A, REIGS, IEIGS, ctrs)
  QRFactorization(N::Int, QF::QFactorization{T,S, V}, RF::Rt, UV, W, A, REIGS, IEIGS, ctrs) where {T, S, V, Rt} = QRFactorization{T, S, V, Rt}(N,QF,RF, UV, W, A, REIGS, IEIGS, ctrs)
end



# constructor for either case
function QRFactorization(QF::QFactorization{T, S},
                         RF::AbstractRFactorization) where {T, S}

    N = length(QF) + 1
    A = zeros(S, 2, 2)
    reigs = zeros(T, N)
    ieigs = zeros(T, N)
    ctr = AMRVW_Counter(0, 1, N-1, 0, N-2)
    m = S == T ? 2 : 1
    UV = Vector{Rotator{T, S}}(undef, m)
    W = Vector{Rotator{T, S}}(undef, m-1)
    QRFactorization(N, QF, RF, UV, W, A, reigs, ieigs, ctr)
end

Base.length(state::AbstractQRFactorizationState) = state.N


# return A
function Base.Matrix(state::AbstractQRFactorizationState)

    Q = Matrix(state.QF)
    R = Matrix(state.RF)

    ## May have Q smaller than  R, so will pad in that case
    if R != I
        m,n =  size(Q)[1], size(R)[1]
        if m < n
            return Q * R[1:m, 1:m]
        else
            return Q * R
        end
    else
        return Q
    end

end
