## Bulge chasing

## A bulge is introduced by a unitary transformation. In the single shift case this is just U^t A U; for
## the double shift case this is (VU)^T A (VU).
## In either case, we have three steps:
##
## * create the unitary transformations (U or (V,U)) -- create_bulge
## * absorb the left side -- absorb_Ut
## * chase the right side down until it is absorbed -- absorb_U
##
function bulge_step(state::QRFactorization{T, S}) where {T, S}

    create_bulge(state)
    absorb_Ut(state) # see CSS or RDS
    absorb_U(state)


    return nothing
end


##################################################
## Absorb U
## After U absorbtion we have either
## A * U or (W,A) * V * U
## THis chases U or (V*U) through until the value can be absorbed
##
## This just pushes work to passthrough_triu and passthrough_Q
##
## * passthrough_triu is basically passthrough(RF) with some efficiencies for first
## time through
##
## * passthrough_Q passes through chain and return false, until fusion when
## it returns true
##
function absorb_U(state::QRFactorization{T, S}) where {T, S}

    flag = false
    while !flag

        passthrough_triu(state, Val(:right))
        flag = passthrough_Q(state, Val(:right))
        ## implicit unitary operation to move U from left to right

    end

    return nothing

end
