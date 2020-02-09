
## Main algorithm of AMRVW
## This follows that given in the paper very closely
## should trace algorithm better with `verbose`
function AMRVW_algorithm(state::QRFactorization)

    it_max = 20 * length(state)
    kk = 0

    while kk <= it_max

        kk += 1
        state.ctrs.stop_index <= 0 && return     ## finished up!
        state.ctrs.it_count += 1

        ## show_status(state)


        check_deflation(state)

        k = state.ctrs.stop_index
        delta =  state.ctrs.stop_index - state.ctrs.zero_index

        if delta >= 2

            bulge_step(state)

        elseif delta == 1

            diagonal_block(state,  k + 1)

            e1r,e1i, e2r,e2i = eigen_values(state)
            state.REIGS[k], state.IEIGS[k] = e1r, e1i
            state.REIGS[k+1], state.IEIGS[k+1] = e2r, e2i

            # can finish up if near end
            if state.ctrs.stop_index == 2
                diagonal_block(state, 2)
                e1 = state.A[1,1]
                state.REIGS[1] = real(e1)
                state.IEIGS[1] = imag(e1)
            end

            state.ctrs.zero_index = 0
            state.ctrs.start_index = 1
            state.ctrs.stop_index -= 2

        elseif delta == 0

            diagonal_block(state, state.ctrs.stop_index + 1)
            e1, e2 = state.A[1,1], state.A[2,2]
            if state.ctrs.stop_index == 1
                state.REIGS[1], state.IEIGS[1] = real(e1), imag(e1)
                state.REIGS[2], state.IEIGS[2] = real(e2), imag(e2)
                state.ctrs.stop_index = 0

            else

                state.REIGS[k+1], state.IEIGS[k+1] = real(e2), imag(e2)
                state.ctrs.zero_index = 0
                state.ctrs.start_index = 1
                state.ctrs.stop_index = k - 1

            end

        end

    end


    @warn "Not all eigenvalues were found. The first $(state.ctrs.stop_index-1) are missing."
end

## Find eigenvalues from factorization
function LinearAlgebra.eigvals(state::QRFactorization, args...)
    AMRVW_algorithm(state)
    es = complex.(state.REIGS, state.IEIGS)
#    try
#        sorteig!(es)
#    catch err
#    end
    es
end


##################################################
## Deflation
## when a Q[k] matrix becomes a diagonal matrix, we deflate.
## This is checked by the sine term being basically 0.
function check_deflation(state::AbstractQRFactorizationState{T,S,Twt}, ctr = state.ctrs, tol=eps(T)) where {T, S,Twt}
    QF = state.QF
    for k in ctr.stop_index:-1:ctr.start_index
        c, s = vals(QF.Q[k])
        if abs(s) <= tol

            deflate(QF, k)

            ctr.zero_index = k      # points to a matrix Q[k] either
            ctr.start_index = k + 1
            ctr.it_count = 1        # reset counter

            return
        end
    end
end
