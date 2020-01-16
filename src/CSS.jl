## ComplexReal case


##################################################

# after a quick fuse
# we must passthrough the diagonal rotator through the rest of Q
function absorb_Ut(state::AbstractFactorizationState{T, S,ComplexRealRotator{T}, QFt, RFt, Val{:not_twisted}}) where {T, S, QFt, RFt}

    QF = state.QF
    Ut = state.UV[1]'
    fuse!(Ut, QF)

    return nothing
end


function passthrough_triu(state::AbstractFactorizationState{T, S,ComplexRealRotator{T}, QFt, RFt, Twt}, dir::Val{:right}) where {T, S, QFt, RFt,Twt}

    U = state.UV[1]
    i = idx(U)

    RF = state.RF

    flag = false

    if  i < state.ctrs.tr

        flag = simple_passthrough!(RF, U)

    end
    if i > state.ctrs.tr || !flag

        U = passthrough!(RF, U)
        state.UV[1] = U

    end

     state.ctrs.tr -= 1

    return nothing
end


## When Ct and B are identical, we can update just one and leave U,V alone
function simple_passthrough!(RF::RFactorization{T, ComplexRealRotator{T}}, U) where {T}

    i = idx(U)
    N = length(RF)

    _ = passthrough!(RF.B, U)

    for k in 0:1

        a,b = vals(RF.B[i+k])
        ii = N + 1 - (i + k)
        RF.Ct[ii] = ComplexRealRotator(conj(a), -b, i + k) # adjoint

    end

    true
end


function passthrough_Q(state::AbstractFactorizationState{T, S,ComplexRealRotator{T}, QFt, RFt, Twt}, dir::Val{:right}) where {T, S, QFt, RFt, Twt}

    QF = state.QF
    U = state.UV[1]
    i = idx(U)


    if i < state.ctrs.stop_index

        U = passthrough!(QF, U)

        state.UV[1] = U

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
function deflate(QF::QFactorization{T, ComplexRealRotator{T}}, k) where {T}

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
    passthrough_phase!(Di,QF)

end
