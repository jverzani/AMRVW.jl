
## Update state.A temporary storage
## with 2x2 view of full matrix; A[k-1:k, k-1:k]
## For QRFactorization, we exploit fact that QF is Hessenberg and
## R  is upper  triangular to simplify the matrix multiplication:
function diagonal_block(state::QRFactorization{T, S, Rt, QFt, RFt}, k) where {T,  S, Rt, QFt,  RFt}

    A  = state.A
    QF,  RF =  state.QF, state.RF

    i, j = k-2, k-1

    qji, qjj, qjk, qkj, qkk = QF[j,i], QF[j,j], QF[j,k], QF[k,j], QF[k,k]

    rjj, rjk = RF[j,j], RF[j,k]
    rkk = RF[k,k]
    rij, rik = RF[i,j], RF[i,k]

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
function eigen_values(state::AbstractFactorizationState{T,S,RealRotator{T}, QFt, RFt, Twt}) where {T,S,QFt, RFt, Twt}

    # this allocates:
    # e1, e2 = eigvals(state.A)
    # return real(e1), imag(e1), real(e2), imag(e2)

    a11, a12 = state.A[1,1], state.A[1,2]
    a21, a22 = state.A[2,1], state.A[2,2]

    b = (a11 + a22) / 2
    c = a11 * a22 - a12 * a21

    e1r, e1i, e2r, e2i = qdrtc(b,c) #qdrtc(one(b), b, c)
    return  e1r, e1i, e2r, e2i

end

# from `modified_quadratic.f90`
function eigen_values(state::AbstractFactorizationState{T,S,ComplexRealRotator{T}, QFt, RFt, Twt}) where {T,S,QFt, RFt, Twt}

    a11, a12 = state.A[1,1], state.A[1,2]
    a21, a22 = state.A[2,1], state.A[2,2]

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
