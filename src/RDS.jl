## real double shift code


##################################################
## Bulge things


## Absorb Ut
## We have the bulge is generated by
## Ut * A * U in the singleshift case and
## Ut * Vt * A * V * U in the DoubleShift Case
## This combines the Ut with A (or Ut, Vt, A into W,A)

function absorb_Ut(QF::QFactorization{T,S, VV}, RF, storage, ctr) where {T,S <: Real, VV}


    U, V = storage.VU[2], storage.VU[1]
    Ut, Vt = U', V'

    i = idx(U)
    j  = i + 1

    W = QF.Q[i]

    p = i == 1 ? one(S) : get_parity(QF,i-1)

    W = dflip(W, p)
    W, Ut, Vt = turnover(Ut, Vt, W)
    fuse!(Vt, QF)
    Ut  = dflip(Ut, p)

    QF.Q[i] = Ut
    storage.limb[1] = W

    return nothing
end




## For double shift we have V then U
function passthrough_triu(QF::QFactorization{T, S, Vt}, RF, storage, ctr,  dir::Val{:right}) where {T, S <: Real, Vt}

    U, V = storage.VU[2], storage.VU[1]
    j = idx(V)


    flag = false

    if  j <= ctr.tr
        # this leverages Ct*B*U = I*U when Ct, B not touched
        # so V, U are not touched, but Ct, B are
        # flag is true if done, false if not
        flag = simple_passthrough!(RF, U, V)
    end

    if j > ctr.tr || !flag

        V = passthrough!(RF, V)
        U = passthrough!(RF, U)
        storage.VU[1], storage.VU[2]  = V, U

    end

    ctr.tr -= 2


    return nothing

end

## When Ct and B are identical, we can update just one and leave U,V alone
function simple_passthrough!(RF::RFactorizationRankOne{T, S}, U, V) where {T, S <: Real}

    j = idx(V)
    N = length(RF)

    _ = passthrough!(RF.B, V)
    _ = passthrough!(RF.B, U)

    for k in -1:1
        a,b = vals(RF.B[j+k])
        jj = N+1-(j+k)
        RF.Ct[jj] = Rotator(a, -b, j+k)  ## adjoint
    end

    true
end

## passthrough Q
##
## we have Q U or (W,Q) V U to pass through
## This depends on Twisted but not pencil
##
## If we update indices and use a unitary transform, return false (not absorbed)
## else return true (was absorved)
function passthrough_Q(QF::QFactorization{T, S, Vt}, RF, storage, ctr, dir::Val{:right}) where {T, S <: Real, Vt}

    U, V, W = storage.VU[2], storage.VU[1],  storage.limb[1]
    i = idx(U); j = i + 1

    if j < ctr.stop_index
        V = passthrough!(QF, V)
        U = passthrough!(QF, U)
        V, U, W = turnover(W, V, U)
        storage.VU[2], storage.VU[1],  storage.limb[1] = U, V, W

        return false

    else

        p = get_parity(QF,j+1) # 1 or -1 possibly ## XXX
        V = dflip(V, p)
        fuse!(QF, V)

        # now turnover U, merge with W, unitary over, pass through, and fuse
        U = passthrough!(QF.Q, U)
        i += 1 # after turnover, U moves down
        U, Di = fuse(W, U)

        # pass U through triangle then fuse
        U = passthrough!(RF, U)

        p = get_parity(QF,i+1)
        U = dflip(U, p)
        fuse!(QF, U)

        return true

    end
end

##################################################
function deflate(QF::AbstractQFactorization{T, S}, k, ctr) where {T,S <: Real}

    c,s = vals(QF.Q[k])
    i = idx(QF.Q[k])
    QF.Q[k] = Rotator(sign(c), zero(T), i) # ± 1, not just 1


end

## The zero_index and stop_index+1 point at diagonal rotators
## [1 0; 0 1] or [-1 0; 0 -1]
## this recovers 1 or -1
function get_parity(QF::QFactorization{T, S}, k) where {T, S <: Real}
    if k > length(QF)
        return one(S)
    else
        c::S, s::T = vals(QF.Q[k])
        ## @assert iszero(s)
        sign(c)
    end
end

##################################################
## Transformations

## Let a D matrix be one of [1 0; 0 1] or [-1 0; 0 1] (D^2 = I). Then we have this move
## D    --->   D  (we update the rotator)
##   [       [
##
## this is `dflip`
## @vars c s
## D=[-1 0 0; 0 -1 0; 0 0 1]
## U = [1 0 0; 0 c s; 0 -s c]
## U1 = [1 0 0; 0 c -s; 0 s c]
## D*U - U1*D
function dflip(a::Rotator{T,S}, d=one(T)) where {T, S<:Real}
    c::S, s::T = vals(a)
    p = sign(d)
    i = idx(a)
    return Rotator(c, p*s, i)
end
