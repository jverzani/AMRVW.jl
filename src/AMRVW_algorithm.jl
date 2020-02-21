## Main algorithm of AMRVW
## This follows that given in the paper very closely
## should trace algorithm better with `verbose`
function AMRVW_algorithm!(state::AbstractQRRotatorFactorization)
    m = eltype(state.QF) <: Real ? 2 : 1
    AMRVW_algorithm!(state, m)
end


function AMRVW_algorithm!(state::AbstractQRRotatorFactorization, m)

    QF, RF = state.QF, state.RF
    S = eltype(QF)
    T = real(S)

    # set counters
    n = length(QF.Q)
    ctr = AMRVW_Counter(0, 1, n, 0, n-1)

    # set up storage
    storage = make_storage(QF, m)
    A = storage.A  # cleanup
    EIGS = zeros(Complex{T}, n+1)

    it_max = 20*n
    kk = 0

    check_deflation(state, ctr)

    while kk <= it_max

        if ctr.stop_index <= 0
            break
        end

        check_deflation(state, ctr)

        ##show_status(state, ctr)

        kk += 1

        ## Delta is number of rotators in piece being considered
        ## m must be greater than delta-1 to run a bulge step
        k = ctr.stop_index
        delta = ctr.stop_index - ctr.zero_index  ## number of rotators in Q factorization

        if delta == 1

            diagonal_block!(A, QF, RF, k, k+1)
            e1::Complex{T}, e2::Complex{T} = eigen_values(A)

            EIGS[k] = e1
            EIGS[k+1] = e2

            # can finish up if near end
            if ctr.stop_index == 2
                diagonal_block!(A, QF, RF, 1, 2)
                e1 = A[1,1]
                EIGS[1] = e1

                break
            end

            ctr.stop_index = ctr.zero_index - 1
            ctr.zero_index = 0
            ctr.start_index = 1
            check_deflation(state, ctr)

            ctr.stop_index <= 0 && break

        elseif delta <= 0


            diagonal_block!(A, QF, RF, k, k+1)
            e1, e2 = A[1,1], A[2,2]  # 0s off diagonal

            if ctr.stop_index == 1
                # finish up
                EIGS[1], EIGS[2] = e1, e2
                ctr.stop_index = 0

                break
            end


            EIGS[k+1] = e2

            ctr.stop_index = ctr.zero_index - 1
            ctr.stop_index <= 0 && break

            ctr.zero_index = 0
            ctr.start_index = 1
            check_deflation(state, ctr)

        elseif m > 2 && m >= delta - 1
            ## When finding shifts we find m eigenvalues from some other means
            ## so here we find the eigvals

            δ, Δ = ctr.start_index, ctr.stop_index
            diagonal_block!(A, QF, RF, δ, Δ+1)

            es = _eigvals(view(A,1:Δ-δ+2,1:Δ-δ+2))
            for (ind, ev) in enumerate(es)
                EIGS[δ+ind-1] = ev
            end

            if  δ == 2
                # finish up
                diagonal_block!(A, QF, RF, 1, 2)

                e1 = A[1,1] # δ == 2 means first rotator a diagonal, so product will be diagonal in 1:2x1:2 block
                EIGS[1] = e1

                break
            end

            ctr.zero_index = 0
            ctr.start_index = 1
            ctr.stop_index = δ - 2

            ctr.stop_index <= 0 && break
            check_deflation(state, ctr)

        else

            bulge_step(QF, RF, storage, ctr, m)
            ctr.it_count -= 1

        end


    end

    LinearAlgebra.sorteig!(EIGS)

    EIGS

end

##################################################
## A container for our counters
mutable struct AMRVW_Counter
    zero_index::Int  # this should be dropped, it is just start_index-1
    start_index::Int
    stop_index::Int
    it_count::Int
    tr::Int
end

function make_counter(state)
    n = length(state.QF.Q)
    AMRVW_Counter(0, 1, n, 0, n-1)
end

##################################################

## Used to  create storage based  on Qfactorization  type
make_storage(state::QRFactorization) = make_storage(state.QF, eltype(state) <: Complex ? 1 : 2)

function make_storage(QF::QFactorization{T,S,V}, m) where {T, S, V}
    A = zeros(S, 4, 4)
    VU = Vector{Rotator{T, S}}(undef, m)   # Ascending Chain
    limb = DescendingChain(copy(QF.Q.x[1:(m-1)]))
    (A=A, VU=VU, limb=limb)
end

function make_storage(QF::QFactorizationTwisted{T,S, Vt,  PVt}, m) where {T, S, Vt, PVt}
    A = zeros(S, max(4, m+2), max(4, m+2))
    VU = Vector{Rotator{T, S}}(undef, m)   # Ascending Chain
    limb = TwistedChain(copy(QF.Q.x[1:(m-1)]), copy(QF.Q.pv[1:(m-2)]))
    (A=A, VU=VU, limb=limb)
end




##################################################
## Deflation
## when a Q[k] matrix becomes a diagonal matrix, we deflate.
## This is checked by the sine term being basically 0.
function check_deflation(state::AbstractQRRotatorFactorization, ctr)

    QF = state.QF
    tol = eps(real(eltype(QF)))

    for k in ctr.stop_index:-1:ctr.start_index

        c, s = vals(QF.Q[k])

        if abs(s) <= tol

            deflate(QF, k, ctr)

            ctr.zero_index = k      # points to a matrix Q[k] either
            ctr.start_index = k + 1
            ctr.it_count = 15        # reset counter

            return
        end
    end

end



## Find eigenvalues from factorization
function eigvals(state::AbstractQRRotatorFactorization)

    es = AMRVW_algorithm!(copy(state))
    LinearAlgebra.sorteig!(es)
    es

end
