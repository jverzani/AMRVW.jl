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
function create_bulge(state::QRFactorization{T, S, Vt, Rt}) where {T, S <: Real, Vt, Rt}

    if mod(state.ctrs.it_count, 15) == 0

        t = rand(T) * pi

        re1, ie1 = sincos(t)
        re2, ie2 = re1, -ie1

        i = state.ctrs.start_index
        j = i + 1

        U =  Rotator(re1, ie1, i)
        V = Rotator(re2, ie2, j)

    else

        Δ = state.ctrs.stop_index
        diagonal_block(state, Δ+1)

        l1r, l1i, l2r, l2i = eigen_values(state)

        k = state.ctrs.start_index

        diagonal_block(state,  k+1)
        bk11, bk12 = state.A[1,1], state.A[1,2]
        bk21, bk22 = state.A[2,1], state.A[2,2]

        diagonal_block(state, k+2)
        bk32 = state.A[2,1]


        # compute `x`
        c1 = -l1i * l2i + l1r*l2r -l1r*bk11 -l2r * bk11 + bk11^2 + bk12 * bk21
        c2 = -l1r * bk21 - l2r * bk21 + bk11* bk21 + bk21 * bk22
        c3 = bk21 * bk32

        c, s, nrm = givensrot(c2, c3)
        i = state.ctrs.start_index #1
        j = i + 1
        V = Rotator(c, -s, j)

        c, s, tmp = givensrot(c1, nrm)
        U = Rotator(c, -s, i)
    end

    state.UV[1] = U
    state.UV[2] = V

    return nothing
end

## CSS case
function create_bulge(state::QRFactorization{T, S, V, Rt}) where {T, S <: Complex, V, Rt}
    ray = true  # state.ray?

    if mod(state.ctrs.it_count, 15) == 0

        t = rand(T) * pi
        if ray
            shift = complex(cos(t), sin(t))
        else
            shift = complex(cos(t), zero(T))
        end

    else

        diagonal_block(state, state.ctrs.stop_index+1)

        if ray
            # Wilkinson
            #e1, e2 = eigen_values(state)
            e1r, e1i, e2r, e2i = eigen_values(state)
            e1 = complex(e1r,e1i)
            e2 = complex(e2r, e2i)
            shift = norm(state.A[2,2] - e1) < norm(state.A[2,2] - e2) ? e1 : e2
        else
            shift = state.A[2,2]
        end

    end

    k =  state.ctrs.start_index
    diagonal_block(state, k+1)
    c,s,nrm = givensrot(state.A[1,1] - shift, state.A[2,1])
    R = Rotator(c, s, k)
    state.UV[1] = R'

    nothing

end
