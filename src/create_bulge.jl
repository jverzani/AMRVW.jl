##################################################
## Create Bulge
## One for DoubleShift/SingleShift
## Pencil or Twisted don't matter as they come out in diagonalblock.

# ## The bulge is created by  (A-rho1) * (A - rho2) * e_1 where rho1 and rho2 are eigenvalue or random
# ## for real case, we take the real part of this result

"""
   create_bulge(state)

Finds m=1 or 2 shifts (m=1 in the CSS case, m=2 in the RDS case) based on the eigenvalues of the lower 2x2 block (using `stop_index`). The
vector `x = alpha (A - rho_1 I) e_1` or `x = alpha (A-rho_1 I) (A-rho_2 I) e1` is found. From this, one or two core transforms
are found so that `U_1' x = gamma e_1` or `U_1' U_2' x = gamma e_1`. The values `U_1` or `U_1`, `U_2` are store in `state`.


"""
function create_bulge(QF::QFactorization{T,S,VV}, RF, storage,  ctr) where {T, S <: Real,VV}

    A = storage.A


    if mod(ctr.it_count, 15) == 0

        t = rand(T) * pi

        re1, ie1 = sincos(t)
        re2, ie2 = re1, -ie1

        i = ctr.start_index
        j = i + 1

        U =  Rotator(re1, ie1, i)
        V = Rotator(re2, ie2, j)

    else


        Δ = ctr.stop_index
        #diagonal_block(state, Δ+1) ??
        diagonal_block!(A, QF, RF, Δ, Δ) #+1)
        #l1r, l1i, l2r, l2i = eigen_values(state)
        e1,  e2 = eigen_values(A)
        l1r, l1i  =  real(e1), imag(e1)
        l2r, l2i =  real(e2), imag(e2)

        delta = ctr.start_index

        #diagonal_block(state,  k+1)
        diagonal_block!(A, QF, RF,  delta, delta)
        bk11, bk12 = A[1,1], A[1,2]
        bk21, bk22 = A[2,1], A[2,2]

        #diagonal_block(state, k+2)
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

    if mod(ctr.it_count, 15) == 0

        t = rand(T) * pi
        if ray
            shift = complex(cos(t), sin(t))
        else
            shift = complex(cos(t), zero(T))
        end

    else

        #diagonal_block(state, state.ctrs.stop_index+1)
        Delta = ctr.stop_index
        diagonal_block!(A, QF, RF, Delta, Delta)

        if ray
            # Wilkinson
            e1, e2 = eigen_values(A[1,1], A[1,2], A[2,1], A[2,2])
            #e1, e2 = complex(e1r, e1i), complex(e2r, e2i)
            #e1, e2 = eigen_values(A)
            shift = norm(A[2,2] - e1) < norm(A[2,2] - e2) ? e1 : e2
        else
            shift = A[2,2]
        end

    end

    delta =  ctr.start_index
    #diagonal_block(state, delta+1)
    diagonal_block!(A, QF, RF, delta,  delta)
    c,s,nrm = givensrot(A[1,1] - shift, A[2,1])
    R = Rotator(c, s, delta)
    storage.VU[1] = R'

    nothing

end
