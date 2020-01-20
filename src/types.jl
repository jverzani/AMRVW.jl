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


abstract type AbstractFactorizationState{T, S, Twt} end

struct QRFactorization{T, S,  Qt<:QFactorization{T,S}, Rt<:AbstractRFactorization{T,S}} <: AbstractFactorizationState{T, S, Val{:not_twisted}}
  N::Int
  QF::Qt
  RF::Rt
  UV::Vector{Rotator{T,S}}    # Ascending chain for cfreating bulge
  W::Vector{Rotator{T,S}}     # the limb when m > 1
  A::Matrix{S}
  REIGS::Vector{T}
  IEIGS::Vector{T}
  ctrs::AMRVW_Counter
  QRFactorization{T,S,Qt, Rt}(N, QF, RF, UV, W, A, REIGS, IEIGS, ctrs) where {T, S, Qt, Rt} = new(N, QF, RF, UV, W, A, REIGS, IEIGS, ctrs)
  QRFactorization(N::Int, QF::QFactorization{T,S}, RF::Rt, UV, W, A, REIGS, IEIGS, ctrs) where {T, S, Rt} = QRFactorization{T, S, QFactorization{T,S}, Rt}(N,QF,RF, UV, W, A, REIGS, IEIGS, ctrs)
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
