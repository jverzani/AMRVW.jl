##################################################
##
## Factorization

# Pencil is made up of two factorizations
# notation in paper is
# solve V - \lambda W
# or V*inv(W) - lambda I
# but V = QR, so consider
# Q*R*inv(W)
# to avoid (introduce?) confusion,we let R be V: and have
# Q * V * inv(W) with the R Factorization store both V and W
struct RFactorizationPencil{T, S} <: AbstractRFactorization{T, S}
V::RFactorizationRankOne{T, S}
W::RFactorizationRankOne{T, S}
end

Base.copy(RF::RFactorizationPencil) = RFactorizationPencil(copy(RF.V), copy(RF.W))
Base.size(RF::RFactorizationPencil) = size(RF.V)
function Base.Matrix(RF::RFactorizationPencil)
    V = Matrix(RF.V)
    W = Matrix(RF.W)
    V[1:end-1, 1:end-1] * inv(W[1:end-1, 1:end-1])
end

## VV = [Sym("v$i$j") for i in ("i","j","k"), j in ("i","j","k")] |> triu
## WW = [Sym("w$i$j") for i in ("i","j","k"), j in ("i","j","k")] |> triu
## VV * inv(WW)[:,3]
function Base.getindex(RF::RFactorizationPencil, l, k)
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

function pencil_factorization(vs, ws)

    V = r_factorization(vs)
    W = r_factorization(ws)
    RFactorizationPencil(V,W)

end


## Pass a rotator through Rfactorization with Pencil from left or right
function passthrough!(RF::RFactorizationPencil, U::AbstractRotator)

    ## Pass Ut -> W (not W^{-1} <- U
    ## Then  passthrough V
    Ut = passthrough!(U', RF.W)
    U = passthrough!(RF.V, Ut')

    U
end

function passthrough!(U::AbstractRotator, RF::RFactorizationPencil)

    U = passthrough!(U, RF.V)
    Ut = passthrough!(RF.W, U')
    U = Ut'

    U
end



## When Ct and B are identical, we can update just one and leave U,V alone
function simple_passthrough!(RF::RFactorizationPencil, U, V)
    false  # might have room for improvement here; o/w consolidate methods
end

function simple_passthrough!(RF::RFactorizationPencil, U)
    false  # might have room for improvement here
end
