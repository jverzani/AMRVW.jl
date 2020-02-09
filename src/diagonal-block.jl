
## Update state.A temporary storage
## with 2x2 view of full matrix; A[k-1:k, k-1:k]
## For QRFactorization, we exploit fact that QF is Hessenberg and
## R  is upper  triangular to simplify the matrix multiplication:
function diagonal_block(state::QRFactorization{T, S}, k) where {T,  S}
    diagonal_block!(state.A, state.QF, state.RF,  k-1, k-1)
end

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




##################################################

# [a11 - l a12; a21 a22] -> l^2 -2 * (tr(A)/2) l + det(A)
# so we use b = tr(A)/2 for qdrtc routing
function eigen_values(state::AbstractQRFactorizationState{T,S, Val{:not_twisted}}) where {T,S}

    # this allocates:
    # e1, e2 = eigvals(state.A)
    # return real(e1), imag(e1), real(e2), imag(e2)

    a11, a12 = state.A[1,1], state.A[1,2]
    a21, a22 = state.A[2,1], state.A[2,2]

    eigen_values(a11, a12, a21, a22)
end

# return tuple
function eigen_values(A::AbstractArray)
    e1r, e1i, e2r, e2i = eigen_values(A[1,1], A[1,2], A[2,1], A[2,2])
    (complex(e1r, e1i), complex(e2r, e2i))
end

function  eigen_values(a11::T, a12::T, a21::T, a22::T) where {T <: Real}
    b = (a11 + a22) / 2
    c = a11 * a22 - a12 * a21

    e1r, e1i, e2r, e2i = qdrtc(b,c) #qdrtc(one(b), b, c)
    return  e1r, e1i, e2r, e2i

end

# from `modified_quadratic.f90`
function  eigen_values(a11::S, a12::S, a21::S, a22::S) where {S <: Complex}

    tr = a11 + a22
    detm = a11 * a22 - a21 * a12
    disc = sqrt(tr * tr - 4.0 * detm)

    u = abs(tr + disc) > abs(tr - disc) ? tr + disc : tr - disc
    if iszero(u)
        e1r, e1i = zero(T), zero(T)
        e2r, e2i = zero(T), zero(T)
    else
        e1 = u / 2.0
        e2 = detm / e1
        e1r, e1i = real(e1), imag(e1)
        e2r, e2i = real(e2), imag(e2)
    end

    return e1r, e1i, e2r, e2i

end

_eigvals(A::AbstractMatrix{T}) where {T <: Real} = eigvals(A)::Vector{Complex{T}}
_eigvals(A::AbstractMatrix{S}) where {S <: Complex} =  eigvals(A)
