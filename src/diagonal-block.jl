## A diagonal block is used to identify the shifts and to
## find the eigenvalues for the deflated matrix

## For the QFactorization (Descending) case provide a 2 x 2 view of the full matrix, A[j:j+1,j:j+1]
## For the Twisted case, provide a view of A[j:J+1, j:J+1] that is
## approximate unless the rotators at j-1 and J+1 are both identity rotators (or implicitly so)


##
## For QRFactorization, we exploit fact that QF is Hessenberg and
## R  is upper  triangular to simplify the matrix multiplication:
function diagonal_block!(A, QF::QFactorization, RF, j, J)


    i, k = j-1, j+1
    qji, qjj, qjk, qkj, qkk = QF[j,i], QF[j,j], QF[j,k], QF[k,j], QF[k,k]

    rjj, rjk = RF[j,j], RF[j,k]
    rkk = RF[k,k]
    if i >= 1
        rij, rik = RF[i,j], RF[i,k]
    else
        rij,  rik = 0*rkk, 0*rkk
    end

    ## This is matrix multiplication of a Hessenberg matrix times a upper triangular
    ## (triu(Q,-1) * triu(R))[k-1:k, k-1:k]
    ## could be just sum(QF[i,l] * RF[l,j] for l in i-1:j)

    A[1,1] = qji * rij + qjj * rjj
    A[1,2] = qji * rik + qjj * rjk + qjk * rkk
    A[2,1] = qkj * rjj
    A[2,2] = qkj * rjk + qkk * rkk

    return nothing
end

## Twisted case
function diagonal_block!(A, QF::QFactorizationTwisted, RF, j, k)
    approx_diagonal_block!(A, QF, RF, j, k)
end

## For twisted we use thee approximate diagonal block
## this is exact when the the j-1 and  k+1  rotators are deflated or j=1 or k=end
function approx_diagonal_block!(A, QF, RF, j, k) where {T, S}
    R = RF
    D = QF.D
    Q = QF.Q
    k = min(k, length(Q))
    # make 1:(k-j+1) identity
    A .= zero(eltype(A))
    if j == 1
        _approx_diagonal_block!(A, Q, D, R, j, k)
    else
        _approx_diagonal_block!(A, Q, D, R, j-1, k)
        # shift up and left by 1
        m,n = size(A)
        for j in 1:n-1  # across columns, then down rows
            for i in 1:m-1
                A[i,j] = A[i+1, j+1]
            end
        end
    end

    nothing
end

function _approx_diagonal_block!(A, Q, D, R, j, k)

    S = eltype(A)
    for i in j:k+1
        for ii in j:(i-1)
            A[i-j+1,ii-j+1] = zero(S)
        end
        for ii in i:k+1
            A[i-j+1,ii-j+1] = D[i] *  R[i,ii]
        end
    end

    # we get first one
    Tw = view(Q, j:k)
    for i in reverse(position_vector_indices(Tw.pv))  # allocates here
        U = Tw.x[i]
        V = Rotator(vals(U)..., idx(U) - j + 1)
        lmul!(V, A)
    end

    return  nothing
end


##################################################

# [a11 - l a12; a21 a22] -> l^2 -2 * (tr(A)/2) l + det(A)
# so we use b = tr(A)/2 for qdrtc routing
## function eigen_values(state::AbstractQRFactorizationState)

##     # this allocates:
##     # e1, e2 = eigvals(state.A)
##     # return real(e1), imag(e1), real(e2), imag(e2)

##     a11, a12 = state.A[1,1], state.A[1,2]
##     a21, a22 = state.A[2,1], state.A[2,2]

##     eigen_values(a11, a12, a21, a22)
## end

# return tuple
function eigen_values(A::AbstractArray{T, N}) where {T, N}
    eigen_values(A[1,1], A[1,2], A[2,1], A[2,2])
end

## function eigen_values(A::AbstractArray{S, N}) where {S <: Complex, N}
##     eigen_values(A[1,1], A[1,2], A[2,1], A[2,2])
## end

function  eigen_values(a11::T, a12::T, a21::T, a22::T) where {T <: Real}
    b = (a11 + a22) / 2
    c = a11 * a22 - a12 * a21

    e1r, e1i, e2r, e2i = qdrtc(b,c) #qdrtc(one(b), b, c)
    return  complex(e1r, e1i), complex(e2r, e2i)

end

# from `modified_quadratic.f90`
function  eigen_values(a11::S, a12::S, a21::S, a22::S) where {S <: Complex}

    tr = a11 + a22
    detm = a11 * a22 - a21 * a12
    disc::S = sqrt(tr * tr - 4.0 * detm)

    u::S = abs(tr + disc) > abs(tr - disc) ? tr + disc : tr - disc
    if iszero(u)
        return zero(S), zero(S)
    else
        e1::S = u / 2.0
        e2::S = detm / e1
        return e1, e2
    end

end

#  this  makes some  algorithm type stable.
_eigvals(A::AbstractMatrix{T}) where {T <: Real} = eigvals(A)::Vector{Complex{T}}
_eigvals(A::AbstractMatrix{S}) where {S <: Complex} =  eigvals(A)
