## complex-real case


##################################################

# after a quick fuse
# we must passthrough the diagonal rotator through the rest of Q
function absorb_Ut(state::QRFactorization{T, S}, storage, ctr) where {T, S <: Complex}

    QF = state.QF
    Ut = storage.VU[1]'
    fuse!(Ut, QF)

    return nothing
end


function passthrough_triu(state::QRFactorization{T, S}, storage, ctr, dir::Val{:right}) where {T, S <: Complex}

    U = storage.VU[1]
    i = idx(U)

    RF = state.RF

    flag = false

    if  i < ctr.tr

        flag = simple_passthrough!(RF, U)

    end
    if i > ctr.tr || !flag

        U = passthrough!(RF, U)
        storage.VU[1] = U

    end

     ctr.tr -= 1

    return nothing
end


## When Ct and B are identical, we can update just one and leave U,V alone
function simple_passthrough!(RF::RFactorizationRankOne{T, S}, U) where {T, S <: Complex}

    i = idx(U)
    N = length(RF)

    _ = passthrough!(RF.B, U)

    for k in 0:1

        a,b = vals(RF.B[i+k])
        ii = N + 1 - (i + k)
        RF.Ct[ii] = Rotator(conj(a), -b, i + k) # adjoint

    end

    true
end


function passthrough_Q(state::QRFactorization{T, S}, storage, ctr,  dir::Val{:right}) where {T, S <: Complex}

    QF = state.QF
    U = storage.VU[1]
    i = idx(U)


    if i < ctr.stop_index

        U = passthrough!(QF, U)

        storage.VU[1] = U

        return false

    else

        # handle details of knitting in
        D = QF.D
        U = passthrough!(D, U)
        Di = fuse!(QF, U) # handles D bit

        return true

    end
end



##################################################
## Deflation
##
# deflate a term
# deflation for ComplexReal is different, as
# we replace Qi with I and move diagonal part into D
function deflate(QF::AbstractQFactorization{T, S}, k) where {T, S <: Complex}

    # when we deflate here we want to leave Q[k] = I and

    # move Dk matrix over to merge with D
    # Qi       Qi            Qi           Qi
    #   Qj   ->   I Dj  -->    I Dj   -->   I Dj
    #     Qk        QK           Qk Dj        Qk Dj    fuse Dj's with D
    #       QL        QL           QL           Ql Dj


    # then the Dk's are collected into D


    alpha, s = vals(QF.Q[k])
    i = idx(QF.Q[k])

    ## Make Q[k] an identity rotator
    QF.Q[k] = Rotator(one(Complex{T}), zero(T), i)

    # absorb Di into D
    Di =  DiagonalRotator(alpha, i)
    passthrough_phase!(Di, QF)

end
