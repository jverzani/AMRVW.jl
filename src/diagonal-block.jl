### Related to decompostion QR into Q*D*Ct*(B*DR + e_1 y^t)


## Update state.A temporary storage
## with 2x2 view of full matrix; A[k-1:k, k-1:k]
function diagonal_block(state::AbstractFactorizationState, k)
    compute_QR(state.A,  state.QF, state.RF, k)
    return nothing
end



## Find R using Ct*B*D decompostion
## writes R as [what, wj, wk,wl] and solves
##
## function rotm(a::T,b, i, N) where {T}
##      r = diagm(0 => ones(T, N))#Matrix{T}(I, N,  N)
##      r[i:i+1, i:i+1] = [a b; -conj(b) conj(a)]
##      r
## end
## @vars bc_i bs_i bc_j bs_j bc_k bs_k
## @vars cc_i cs_i cc_j cs_j cc_k cs_k
## @vars d_i d_j d_k d_l
## D = diagm(0 => [d_i,d_j,d_k,d_l])
## Bi = rotm(bc_i, bs_i, 1, 4)
## Bj = rotm(bc_j, bs_j, 2, 4)
## Bk = rotm(bc_k, bs_k, 3, 4)
## what, wj, wk, wl = Bi * Bj * Bk * D * [0,0,1,0]
## @show wj wk wl
## @vars what wj wk wl
## u = rotm(cc_k, cs_k, 1,2) * [what, wl]
## kk = u[1](what => solve(u[2], what)[1]) |> simplify
## @show kk

## u =  rotm(cc_j, cs_j, 2,3) * rotm(cc_k, cs_k, 1,3) * [what, wk, wl]
## kk_1 = u[1](what => solve(u[2], what)[1]) |> simplify
## @show kk_1

## u =  rotm(cc_k, cs_k, 3,4) * rotm(cc_j, cs_j, 2,4) * rotm(cc_i, cs_i, 1,4) * [what, wi, wk, wl]
## kk_2 = u[1](what => solve(u[2], what)[1]) |> simplify
## @show kk_2
##
## wj = bc_k*bs_j*d_k*conjugate(bc_i)
## wk = bc_k*d_k*conjugate(bc_j)
## wl = -d_k*conjugate(bs_k)
## kk = cc_k*wl*conjugate(cc_k)/conjugate(cs_k) + cs_k*wl
## kk_1 = cc_k*wk*conjugate(cc_k)/conjugate(cs_k) + cs_k*wk + cc_k*cs_j*wl/(cc_j*conjugate(cs_k))
## kk_2 = bc_k*bs_j*cc_i*d_k*conjugate(bc_i)*conjugate(cc_i)/conjugate(cs_i) + bc_k*bs_j*cs_i*d_k*conjugate(bc_i) + cc_i*cs_j*wk/(cc_j*conjugate(cs_i))

function Rjk(Ct, B, dk, j, k)

    N = length(Ct)
    # @assert  0 <= k - j <= 2

    if k - j == 0 # (k,k) case is easiest
        bkc, bks = vals(B[k])
        ckc, cks = vals(Ct[N+1-k])

        wl = -conj(bks) * dk
        return wl / conj(cks)

    elseif k - j == 1

        j = k-1
        bkc, bks = vals(B[k])
        bjc, bjs = vals(B[j])
        ckc, cks = vals(Ct[N+1-k])
        cjc, cjs = vals(Ct[N+1-j])

        wl = -conj(bks) * dk
        wk = bkc * conj(bjc) * dk

        return wk/conj(cjs) - (wl * cjc * conj(ckc))/conj(cjs*cks)


    else# k - j == 2

        i, j = k-2, k-1

        bkc, bks = vals(B[k])
        bjc, bjs = vals(B[j])
        bic, bis = vals(B[i])
        ckc, cks = vals(Ct[N+1-k])
        cjc, cjs = vals(Ct[N+1-j])
        cic, cis = vals(Ct[N+1-i])

        wl = -conj(bks) * dk
        wk = bkc * conj(bjc) * dk
        wj = bjs * bkc * conj(bic) * dk

        return wj/conj(cis) - (wk * cic * conj(cjc)) /conj(cis*cjs) + (wl*cic * conj(ckc))/conj(cis*cjs*cks) #+ wj*cis

    end
end


# Compute [k-1:k, k-1:k] part of QR
# We use
## @vars ci si cj sj ck sk rij rik rjj rjk rkj rkk di dj dk dl
## Q1 = [ci si 0 0; -si conj(ci) 0 0; 0 0 1 0; 0 0 0 1]
## Q2 = [1 0 0 0; 0 cj sj 0; 0 -sj conj(cj) 0; 0 0 0 1]
## Q3 = [1 0 0 0; 0 1 0 0; 0 0 ck sk; 0 0 -sk conj(ck)]
## D = diagm(0 => [di,dj,dk,dl])
## R = [rij rik; rjj rjk; 0 rkk; 0 0]
## ((Q1*Q2 * Q3 * D ) * R)[2:3, :]
function compute_QRX(A, QF::QFactorization, RF, k)

    Q = QF.Q

    i, j = k-2, k-1
    cj, sj = vals(Q[j])
    ck, sk = vals(Q[k])
    ci, si = (k > 2) ? vals(Q[i]) : (one(cj), zero(sj))

    rjj, rjk = RF[j,j], RF[j,k]
    rkk = RF[k,k]
    rij, rik = (k > 2) ? (RF[i,j], RF[i,k]) : (zero(cj), zero(cj))

    dj, dk = QF.D[j], QF.D[k]
    di = (k > 2) ? QF.D[i] : one(dj)


    A[1,1] = cj * dj * rjj * conj(ci) - di*rij*si
    A[1,2] = cj * dj * rjk * conj(ci) + ck * dk * rkk * sj * conj(ci) - di * rik * si
    A[2,1] = -dj * rjj * sj
    A[2,2] = ck * dk * rkk * conj(cj) - dj * rjk * sj

    return nothing

end

function compute_QR(A, QF::QFactorization, RF, k)

    Q = QF.Q

    i, j = k-2, k-1

    #(qji, qjj, qjk, qkj, qkk) = getQ(QF, k)
    qji, qjj, qjk, qkj, qkk = QF[j,i], QF[j,j], QF[j,k], QF[k,j], QF[k,k]

    rjj, rjk = RF[j,j], RF[j,k]
    rkk = RF[k,k]
    rij, rik = (k > 2) ? (RF[i,j], RF[i,k]) : (zero(rkk), zero(rkk))

#    @show qji, qjj, qjk
    ## This is matrix multiplication of a Hessenberg matrix times a upper triangular
    ## (triu(Q,-1) * triu(R))[k-1:k, k-1:k]
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
