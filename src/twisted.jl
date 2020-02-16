##
## Implementation of the algorithm from
## A GENERALIZATION OF THE MULTISHIFT QR ALGORITHM by RAF VANDEBRIL AND DAVID S. WATKINS.
## https://doi.org/10.1137/11085219X
## pdf file found at http://www.sci.wsu.edu/math/faculty/watkins/pdfiles/vw12_SIMAX.pdf
##
## Q is neither descending or ascending, but rather twisted
## the unitary transformation is not limited to 1 or 2 rotators, but rather m
##
## TODO:
##
## * Clean up allocations
##
## This algorithm encompasses those in CSS and RDS, but those are more efficient space wise
## This is more general and perhaps of interest as building blocks for experimentation

##################################################


##
## Twisted refers to the QFactorization
##

## unlike faster case, here instead of checking for parity in the diagonal rotators, we move -1 terms into  D
function deflate(QF::QFactorizationTwisted{T, S, Vt, Pvt}, k, ctr) where {T,S <: Real, Vt, Pvt}
    c,s = vals(QF.Q[k])
    i = idx(QF.Q[k])
    #@show :deflate, k
    QF.Q[k] = Rotator(one(T), zero(T), i) # ± 1, not just 1
    if sign(c) < 0
        ## move into D term so no dflip concerns, as with RDS
        Asc = ascending_part(QF.Q, i)
        Des = descending_part(QF.Q, i)

        m, M = i, i+1
        for (i,U) in enumerate(Asc)
            c,s = vals(U)
            idx(U) < ctr.start_index && break
            m -= 1
            j = idx(U)
            Asc[i] = Rotator(-c, s, j)
        end
        for (i, U) in enumerate(Des)
            c,s = vals(U)
            idx(U) > ctr.stop_index && break
            M += 1
            j = idx(U)
            Des[i] = Rotator(-c, s, j)
        end

        QF.D.x[m] *= -one(T)
        QF.D.x[M] *= -one(T)
    end
end

function Base.Matrix(QF::QFactorizationTwisted{T, S}) where {T, S}

    n = length(QF) + 1
    M = diagm(0 => ones(S, n))
    D = Matrix(QF.D) * M
    lmul!(Vector(QF.Q), D)
    D

end









##################################################



## We use `bulge_step!` for more general usage;
function bulge_step(QF::QFactorizationTwisted{T,S}, RF, storage, ctr, m)  where {T,S}

    create_bulge(QF, RF, storage, ctr, m)

    δ, Δ = ctr.start_index, ctr.stop_index
    n = ctr.stop_index - ctr.zero_index + 1

    Asc = copy(storage.VU[1:m])  # stored as V,U ## XXX copy! Could eliminate this
    Ms = QF.Q
    D = QF.D
    RF = RF


    bulge_step!(n, m, view(Ms, δ:Δ), D, RF, Asc)

    return nothing
end

## Can enter here for more general usage
## n is the size of matrix, not the number of rotators
function bulge_step!(n, m, Ms::TwistedChain, D, RF::AbstractRFactorization, Asc, choice = :left)


    ## XXX
    ## We *can* save memory if we
    ## * didn't copy VU to Asc
    ## * put Asc, Des in storage and didn't adjust sizes in step_knit
    ## * used the limb defined in create storage
    ## This would require some reworking of step_knit!

    Des = reverse(adjoint.(Asc))

    # psd[k] = ps[m+k] with padding soecufued by choice
    ps = copy(Ms.pv)
    psd = ps[m+1:end]
    if isa(choice, Symbol)
        append!(psd, repeat([choice], min(m, length(Ms.pv))))
    else
        append!(psd, choice)
    end

    ## Must put AscendingChain into place
    ## py passing through RF and D <--
    passthrough!(RF, AscendingChain(Asc))
    passthrough!(D, AscendingChain(Asc))

    M1 = Des* (Ms * (Asc * (Matrix(D) * Matrix(RF))))
    @show :eigvals
    printtp(eigvals(M1))





    limb, limb_side = step_0!(m, ps, psd, Ms, Des, Asc, D, RF)




    for k in 1:(n-m-2)
        limb_side = step_k!(k, n, m, psd, limb_side, limb, Des, Asc, Ms, D, RF)
    end

    # use new choice
    if n-m-1 > 0
        limb_side = step_k!(n-m-1, n, m, psd, limb_side, limb, Des, Asc, Ms, D, RF)
    end

    @show  limb_side,  psd[end-m+1], n-m-1
    M1 = limb * (Des * (Asc * (view(Ms, 1:(n-m-1)) *  (Matrix(D) * Matrix(RF)))))
    @show :eigvals
    printtp(eigvals(M1))

    # now knit in  limb, Des, Asc
    step_knit!(n, m, psd, limb_side, limb, Des,  Asc, Ms, D, RF)

    Ms.pv[:] = psd

    return nothing
end


############################################################
##
## Implement thee different steps for one pass of the bulge chasing algorithm

function step_0!(m,  ps, psd, Ms::TwistedChain, Des, Asc, D, RF) where {T}

    has_limb = ifelse(m > 1, true, false)

    # grab limb, but keep order
    # limb is a copy, not a view
    limb = TwistedChain(copy(Ms.x[1:(m-1)]), copy(Ms.pv[1:(m-2)]))
    limb_side = :nothing

    if has_limb
        if ps[m-1] == :left
            passthrough!(DescendingChain(Des), limb)  # pass limb through left
            limb_side = :left
        else
            passthrough!(limb, AscendingChain(Asc))   # pass limb to right
            limb_side = :right
        end
    end



    U = Ms.x[m]
    psm = m <= length(ps) ? ps[m] : :right # doesn't matter in this case

    if psm == :left
        Des[end], Di = fuse(Des[end], U)  # fuse with descending; aka Ms[m]

        ## need to passthrough Ms too now
        i = idx(Di)
        dMs = descending_part(Ms, i) # stop_index

        if limb_side == :right
            passthrough_phase!(Di, dMs,  (AscendingChain(Asc), limb), D)
        else
            passthrough_phase!(Di, dMs,  (AscendingChain(Asc), ), D)
        end

    else

        AA = popfirst!(Asc)
        AA, Di = fuse(U, AA)

        if limb_side == :right
            passthrough_phase!(Di, (AscendingChain(Asc), limb), D)
        else
            passthrough_phase!(Di, (AscendingChain(Asc), ), D)
        end
        pushfirst!(Asc, AA)

    end

    limb, limb_side

end


function step_k!(k, n, m, psd, limb_side, limb, Des, Asc, Ms, D, RF) where {T}

    phatk = psd[k]

    ## 4 cases based on:
    ## phatk/limb_side
    ## A bit redundandant, but easier to verify

    U = Ms.x[m+k]

    if phatk == :left && limb_side == :left

        push!(Des, U)  ## Augment Descending

        passthrough!(DescendingChain(Des), AscendingChain(Asc)) ## translate ascending

        passthrough!(limb, AscendingChain(Asc))

        ## similarity transform Asc to right side
        passthrough!(RF, AscendingChain(Asc)) # <-- pass Asc through RF
        passthrough!(D, AscendingChain(Asc))

        U = popfirst!(Des)

        if m > 1
            limb_side = :left
        end


    elseif phatk == :left && (limb_side == :right || limb_side == :nothing)

        push!(Des, U)  ## Augment Descending
        passthrough!(DescendingChain(Des), AscendingChain(Asc)) ## translate ascending
        passthrough!(DescendingChain(Des), limb)

        ## similarity transform Asc to right side
        passthrough!(RF, AscendingChain(Asc)) # <-- pass Asc through RF
        passthrough!(D, AscendingChain(Asc))

        U = popfirst!(Des)

        if m > 1
            limb_side = :left
        end

    elseif phatk == :right && limb_side == :left

        pushfirst!(Asc, U)

        passthrough!(DescendingChain(Des), AscendingChain(Asc))

        passthrough!(limb, AscendingChain(Asc))


        ## similarity transform descending to left side
        passthrough!(DescendingChain(Des), D)
        passthrough!(DescendingChain(Des), RF) #--> pass Des through Rf;  similarity transform

        U = pop!(Asc)

        if m > 1
            limb_side = :right
        end

    elseif phatk == :right && (limb_side == :right || limb_side == :nothing)

        pushfirst!(Asc, U)
        passthrough!(DescendingChain(Des), AscendingChain(Asc))

        passthrough!(DescendingChain(Des), limb)

        ## similarity transform descending to left side
        passthrough!(DescendingChain(Des), D)
        inds = idx.(Des)
        passthrough!(DescendingChain(Des), RF) #--> pass Des through Rf;  similarity transform
        @assert all(idx.(Des) .== inds)
        U = pop!(Asc) # k = idx(U)

        if m > 1
            limb_side = :right
        end

    end
@show k
    Ms.x[k] = U

    return limb_side

end

## steps k=n-m to  n-2
function step_knit!(n, m, psd, limb_side, limb, Des,  Asc, Ms, D, RF) where {T}

    # make bottom
    U = pop!(Des)
    V = popfirst!(Asc)
    bottom, Di = fuse(U,V)

    ## Decoupled too! Decoupled on right? as we augmented Ascending
    ## have a complicted check, as m can be big
    psdm = length(psd)-m+1 > 0 ? psd[end-m+1] : :left

    if psdm == :right
        Decoupled = view(Ms, 1:n-m)
        if limb_side == :right
            passthrough_phase!(Di, (AscendingChain(Asc), Decoupled, limb), D)
        else
            passthrough_phase!(Di, (AscendingChain(Asc),  Decoupled), D)
        end
    else
        if limb_side == :right
            passthrough_phase!(Di, (AscendingChain(Asc), limb), D)
        else
            passthrough_phase!(Di, (AscendingChain(Asc),), D)
        end
    end


    ## Final Steps, knit in
    ## setps k=n-m to  n-2

    for k in (n-m):(n-2)
        phatk = k > 1 ? psd[k-1] : :left
        if k > (n-m) && k > 2 && psd[k-2] != phatk
            ## need to reposition V structure to other side
            if phatk == :left
                passthrough!(RF, DescendingChain(Des))
                passthrough!(D, DescendingChain(Des))
                bottom = passthrough!(RF, bottom)
                bottom = passthrough!(D, bottom)
                passthrough!(RF, AscendingChain(Asc))
                passthrough!(D, AscendingChain(Asc))
            else
                passthrough!(AscendingChain(Asc),D)
                passthrough!(AscendingChain(Asc),RF)
                bottom = passthrough!(bottom, D)
                ind = idx(bottom)
                bottom = passthrough!(bottom, RF)
                @assert idx(bottom) == ind
                passthrough!(DescendingChain(Des),D)
                passthrough!(DescendingChain(Des),RF)
            end
            # flip limb side
            if limb_side == :right
                limb_side = :left
            elseif limb_side == :left
                limb_side = :right
            end
        end

        ## lengthen one side of V
        if phatk == :left
            Des = push!(Des, bottom)
        else
            Asc = pushfirst!(Asc, bottom)
        end
        passthrough!(DescendingChain(Des), AscendingChain(Asc))

        # remove bottom rotator from limb; keep track of if it is right or left
        lpk = length(limb.pv) > 0 ? limb.pv[end] : :nothing
        L, lpk =  pop!(limb)
        ## we have
        ## 4 cases in terms of limb_side and lpk:
        ## :l, :l -> fuse (Asc, Des), translate
        ## :l, :r -> translate, fuse (Asc, limb, Des)
        ## :r, :l -> translate, fuse (limb)
        ## :r, :r -> fuse (.), translate
        if limb_side == :left || limb_side == :nothing
            if lpk == :left || lpk == :nothing

                # fuse then translate
                Asc[1], Di = fuse(L, Asc[1])
                if length(Asc) > 1
                    AA = Asc[2:end]
                    passthrough_phase!(Di,  (AscendingChain(AA), DescendingChain(Des)), D)
                    Asc[2:end] = AA
                else
                    passthrough_phase!(Di, (DescendingChain(Des), ), D)
                end

                ## translate
                length(limb) > 0 && passthrough!(limb, AscendingChain(Asc))

            elseif lpk == :right

                # translate then fuse
                length(limb) > 0 && passthrough!(limb, AscendingChain(Asc))

                # then fuse
                Asc[1], Di = fuse(L, Asc[1])
                if length(Asc) > 1
                    AA = Asc[2:end]
                    passthrough_phase!(Di, (AscendingChain(AA), limb, DescendingChain(Des)), D)
                    Asc[2:end] = AA
                else
                    passthrough_phase!(Di, (limb, DescendingChain(Des)), D)
                end

            end

        else # limb on right

            if  lpk == :left

                # translate then fuse
                length(limb) > 0 && passthrough!(DescendingChain(Des), limb)

                Des[end], Di = fuse(Des[end], L)
                passthrough_phase!(Di, (), D)
            else
                # fuse  then translate
                Des[end], Di = fuse(Des[end], L)
                passthrough_phase!(Di, limb, (), D)

                length(limb) > 0 && passthrough!(DescendingChain(Des), limb)
            end
        end

        if phatk == :left
            # similarity to move Asc to right side
            passthrough!(RF, AscendingChain(Asc))  # <--
            passthrough!(D, AscendingChain(Asc))

            limb_side = :left
        else
            # similarity to move Des to left side
            passthrough!(DescendingChain(Des), D)
            inds = idx.(Des)
            passthrough!(DescendingChain(Des), RF)
            @assert all(idx.(Des) .== inds)

            limb_side = :right

        end

        ## pop off top, and add to Decopuled
        ## we have idx(U) ==   k; so psd[k-1] = phatk indicates direction to  add
        if phatk == :left
            U  = popfirst!(Des)
            Ms.x[k] = U
        else
            U = pop!(Asc)
            Ms.x[k] = U
        end

        U = pop!(Des)
        V = popfirst!(Asc)
        bottom, Di = fuse(U,V)

        if phatk == :right
            Decoupled = view(Ms, 1:k)
            if limb_side == :right
                passthrough_phase!(Di, (AscendingChain(Asc), Decoupled, limb), D)
            else
                passthrough_phase!(Di, (AscendingChain(Asc), Decoupled), D)
            end
        else
            if limb_side == :right
                passthrough_phase!(Di, (AscendingChain(Asc), limb), D)
            else
                passthrough_phase!(Di, (AscendingChain(Asc),), D)
            end
        end

    end

    # just the bottom left, but may be on wrong side if m > 1
    if psd[end] == :left
        if m > 1 && (length(psd) > 1 && psd[end-1] != psd[end])
            bottom = passthrough!(RF, bottom)  # move to other side
            bottom = passthrough!(D, bottom)
        end

        Ms.x[end] = bottom

    else
        if m > 1 && (length(psd) > 1 && psd[end-1] != psd[end])
            bottom = passthrough!(bottom, D)   # move bottom to other side
            bottom = passthrough!(bottom, RF)
        end
        @assert length(Ms.x) == idx(bottom)
        Mx.x[end] = bottom
    end

    return nothing

end
