## Twisted refers to the QFactorization
##
struct QFactorizationTwisted{T, Rt} <: AbstractQFactorization{T, Rt, Val{:twisted}}
  Q::TwistedChain{Rt}
  D::AbstractSparseDiagonalMatrix{T, Rt}
W::Vector{Rt}
end


## return q[i:k, j:k]
function Base.getindex(QF::QFactorizationTwisted, j, k)

    ## XXX
    Q = QF.Q

    i, j = k-2, k-1

    ### XXX  Fill this in ....
    ### This isn't going to be correct
    (qji, qjj, qjk, qkj, qkk)
end
