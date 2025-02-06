##################################################
## Create Bulge
## 4  cases DoubleShift/SingleShift x Descending/Twisted Q Factorizations
## RFactorization type doesn't enter here

# ## The bulge is created by  (A-rho1) * (A - rho2) * e_1 where rho1 and rho2 are eigenvalue or random
# ## for real case, we take the real part of this result

const IT_COUNT = 15  # if in 15 steps deflation does not occur, this kicks in a random rotator

#=
   create_bulge(state)

Finds m=1 or 2 shifts (m=1 in the CSS case, m=2 in the RDS case) based on the eigenvalues of the lower 2x2 block (using `stop_index`). The
vector `x = alpha (A - rho_1 I) e_1` or `x = alpha (A-rho_1 I) (A-rho_2 I) e1` is found. From this, one or two core transforms
are found so that `U_1' x = gamma e_1` or `U_1' U_2' x = gamma e_1`. The values `U_1` or `U_1`, `U_2` are store in `state`.


=#
function create_bulge(QF::QFactorization{T,S,VV}, RF, storage,  ctr) where {T, S <: Real,VV}

    A = storage.A


    if iszero(ctr.it_count)
        ctr.it_count = IT_COUNT

        t = rand(T) * pi

        re1, ie1 = sincos(t)
        re2, ie2 = re1, -ie1

        i = ctr.start_index
        j = i + 1

        U =  Rotator(re1, ie1, i)
        V = Rotator(re2, ie2, j)

    else


        Δ = ctr.stop_index
        diagonal_block!(A, QF, RF, Δ, Δ) #+1)
        e1,  e2 = eigen_values(A)

        l1r, l1i  =  real(e1), imag(e1)
        l2r, l2i =  real(e2), imag(e2)

        delta = ctr.start_index

        diagonal_block!(A, QF, RF,  delta, delta)
        bk11, bk12 = A[1,1], A[1,2]
        bk21, bk22 = A[2,1], A[2,2]

        diagonal_block!(A, QF, RF, delta+1, delta+1)
        bk32 = A[2,1]


        # compute `x`
        c1 = -l1i * l2i + l1r*l2r -l1r*bk11 -l2r * bk11 + bk11^2 + bk12 * bk21
        c2 = -l1r * bk21 - l2r * bk21 + bk11* bk21 + bk21 * bk22
        c3 = bk21 * bk32

        c, s, nrm = givensrot(c2, c3)
        i = ctr.start_index #1
        j = i + 1
        V = Rotator(c, -s, j)

        c, s, tmp = givensrot(c1, nrm)
        U = Rotator(c, -s, i)
    end

    storage.VU[2] = U
    storage.VU[1] = V

    return nothing
end

## CSS case
function create_bulge(QF::QFactorization{T, S, VV}, RF, storage, ctr) where {T, S <: Complex, VV}
    ray = true # true seems to take  fewer steps than false

    A = storage.A

    if iszero(ctr.it_count)
        ctr.it_count = IT_COUNT

        t = rand(T) * pi
        if ray
            shift = complex(cos(t), sin(t))
        else
            shift = complex(cos(t), zero(T))
        end

    else

        Delta = ctr.stop_index
        diagonal_block!(A, QF, RF, Delta, Delta)

        if ray
            # Wilkinson
            e1, e2 = eigen_values(A[1,1], A[1,2], A[2,1], A[2,2])
            shift = norm(A[2,2] - e1) < norm(A[2,2] - e2) ? e1 : e2
        else
            shift = A[2,2]
        end

    end

    delta =  ctr.start_index
    diagonal_block!(A, QF, RF, delta,  delta)
    c,s,nrm = givensrot(A[1,1] - shift, A[2,1])
    R = Rotator(c, s, delta)
    storage.VU[1] = R'

    nothing

end

##################################################
##
##  Twisted case



function create_bulge(QF::QFactorizationTwisted{T, S}, RF, storage, ctr, m) where {T, S}

    if iszero(ctr.it_count)
        ctr.it_count = IT_COUNT

        for i in 1:m
            storage.VU[i] = random_rotator.(S, ctr.start_index + m - i)
        end
    else
        _create_bulge(QF, RF, storage, ctr, m)
    end

    return nothing

end

## Real
function create_bulge(QF::QFactorizationTwisted{T, S}, RF, storage, ctr, m) where   {T, S <: Real}

    @assert m >= 2  # for even we need atleast a quadratic poly

    if iszero(ctr.it_count)
        ctr.it_count = IT_COUNT


        t = rand(T) * pi

        re1, ie1 = sincos(t)
        re2, ie2 = re1, -ie1

        i = ctr.start_index
        j = i + 1

        U =  Rotator(re1, ie1, i)
        V = Rotator(re2, ie2, j)

        storage.VU[1], storage.VU[2] = V, U
        return nothing
    end

    delta, Delta = ctr.start_index, ctr.stop_index
    n = Delta - delta
    A = storage.A


    offset = delta - 1

    diagonal_block!(A, QF, RF, Delta-m+1, Delta)

    Am = view(A, 2:m+1, 2:m+1)

    if m == 2

        rhos = eigen_values(Am)

    else

        rhos = _eigvals(Am)  ## sorted by LinearAlgebra.sorteig!

    end

    x = zeros(S, m+1)
    x[1] = one(S)

    diagonal_block!(A, QF, RF, delta, delta + m - 1)
    Am = view(A, 1:m+1, 1:m+1)

    i = 1
    while i <= m
        rhoi = rhos[i]
        pi = (i==1) ? :left : QF.Q.pv[i-1]
        pj = QF.Q.pv[i]

        if i < m && isapprox(rhoi, conj(rhos[i+1]))
            ## XXX deal with left right...
            x = (Am^2 - 2real(rhos[i])*Am + norm(rhos[i])^2*I) * x
            if pi == :right
                x = Am \ x
            end
            if pj == :right
                x = Am \ x
            end
            i += 1
        else
            x = (Am - real(rhos[i])*I) * x
            if pi == :right
                x = Am \ x
            end
        end

        i += 1
    end

    # now roll up
    _rollup!(storage.VU, x, delta-1)
    return nothing

end

## complex
function _create_bulge(QF::QFactorizationTwisted{T, S}, RF, storage, ctr, m) where   {T, S <:  Complex}

    delta, Delta = ctr.start_index, ctr.stop_index
    n = Delta -  delta
    A = storage.A

    diagonal_block!(A,QF, RF, Delta-m+1, Delta)

    if m == 1

        Am = view(A, 1:2, 1:2)
        e1, e2 = eigen_values(Am)

        if norm(Am[2,2] - e1) < norm(Am[2,2] - e2)
            rhos = [e1]
        else
            rhos = [e2]
        end

    elseif m == 2

        Am = view(A, 2:3, 2:3)
        rhos = eigen_values(Am)

    else

        Am = view(A, 2:m+1, 2:m+1)
        rhos = _eigvals(Am)

    end

    x = zeros(S, m+1)
    x[1] = one(S)

    diagonal_block!(A, QF, RF, delta, delta + m - 1)
    Am = view(A, 1:m+1, 1:m+1)

    for i in 1:m
        x = (Am - rhos[i] * I) * x
        pi = (i==1) ? :left : QF.Q.pv[i-1]
        if  pi == :right
            x = Am \ x
        end
    end

    # now roll up
    _rollup!(storage.VU, x, delta - 1)

    return nothing
end

function _rollup!(VU, x, offset=0)
    r = x[end]
    m = length(x)
    for i in (m-1):-1:1
        c, s, r = givensrot(x[i], r)
        U = Rotator(c, s, i+offset)
        VU[m-i] = U'
    end
end
