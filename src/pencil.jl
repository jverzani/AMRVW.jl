
##################################################
##
## Factorization

# Pencil is made up of two factorizations
# algorithm is QZ algorithm, so we use Z factorization fir bane
struct ZFactorization{T, Rt} <: AbstractRFactorization{T, Rt, Val{:pencil}}
V::RFactorization{T, Rt}
W::RFactorization{T, Rt}
end

## VV = [Sym("v$i$j") for i in ("i","j","k"), j in ("i","j","k")] |> triu
## WW = [Sym("w$i$j") for i in ("i","j","k"), j in ("i","j","k")] |> triu
## VV * inv(WW)[:,3]
function Base.getindex(RF::ZFactorization, l, k)
    Δ = k - l
    #@assert 0 <= Δ <= 2
    V, W = RF.V, RF.W
    l <= 0 &&  return zero(RF)
    if iszero(Δ)
        return V[k,k] / W[k,k]
    elseif Δ == 1
        j = l
        return - V[j,j] * W[j,k] / (W[j,j]*W[k,k]) + V[j,k]/W[k,k]
    else
        i, j = l, l+1
        return V[i,i] * (W[i,j] * W[j,k] - W[i,k] * W[j,j]) / (W[i,i] * W[j,j] * W[k,k])  - (V[i,j] *  W[j,k]) / ( W[j,j] * W[k,k])  + V[i,k]/W[k,k]


    end
end





## Constructor for Real coefficients; with pencil

function z_factorization(vs, ws)

    V = r_factorization(vs)
    W = r_factorization(ws)
    ZFactorization(V,W)

end


## Pass a rotator through Rfactorization with Pencil from left or right
function passthrough(RF::ZFactorization{T, St}, U::AbstractRotator, ::Val{:right}) where {T, St}

    ## Pass Ut -> W (not W^{-1} <- U
    ## Then  passthrough V
    Ut = passthrough(RF.W, U', Val(:left))
    U = passthrough(RF.V, Ut', Val(:right))

    U
end

function passthrough(RF::ZFactorization{T, St}, U::AbstractRotator, ::Val{:left}) where {T, St}

    U = passthrough(RF.V, U, Val(:left))
    Ut = passthrough(RF.W, U', Val(:right))
    U = Ut'

    U
end



## When Ct and B are identical, we can update just one and leave U,V alone
function simple_passthrough(RF::ZFactorization{T, Rt}, U, V, ::Val{:right}) where {T, Rt}
    false  # might have room for improvement here; o/w consolidate methods
end

function simple_passthrough(RF::ZFactorization{T, Rt}, U, ::Val{:right}) where {T, Rt}
    false  # might have room for improvement here
end
