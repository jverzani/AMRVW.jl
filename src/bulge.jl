## Bulge chasing for the RDS and CSS case

## A bulge is introduced by a unitary transformation. In the single shift case this is just U^t A U; for
## the double shift case this is (VU)^T A (VU).
## In either case, we have three steps:
##
## * create the unitary transformations (U or (V,U)) -- create_bulge
## * absorb the left side -- absorb_Ut
## * chase the right side down until it is absorbed -- absorb_U
##
function bulge_step(QF::QFactorization{T, S, Vt}, RF, storage, ctr, m=0) where {T, S, Vt}

    create_bulge(QF, RF, storage, ctr)
    absorb_Ut(QF, RF, storage, ctr) # see CSS or RDS
    absorb_U(QF, RF, storage, ctr)


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
function absorb_U(QF::QFactorization{T, S, Vt}, RF, storage, ctr) where {T, S, Vt}

    flag = false
    while !flag

        passthrough_triu(QF, RF, storage, ctr, Val(:right))
        flag = passthrough_Q(QF, RF, storage,  ctr, Val(:right))
        ## implicit unitary operation to move U from left to right

    end

    return nothing

end
