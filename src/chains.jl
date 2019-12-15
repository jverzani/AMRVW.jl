## A chain is an arrangement of rotators between n and N
## Three types are specialized: Descending, Ascending, Twisted

#### Group rotators

abstract type AbstractRotatorChain{T} end

Base.eltype(A::AbstractRotatorChain) = eltype(A.x)

Base.length(A::AbstractRotatorChain) = length(A.x)

Base.@propagate_inbounds Base.getindex(A::AbstractRotatorChain, i::Int) = getindex(A.x, i)
Base.@propagate_inbounds Base.setindex!(A::AbstractRotatorChain, X, inds...) = setindex!(A.x, X, inds...)


Base.iterate(A::AbstractRotatorChain) = iterate(A.x)
Base.iterate(A::AbstractRotatorChain, st) = iterate(A.x, st)

struct DescendingChain{T} <: AbstractRotatorChain{T}
  x::Vector{T}
end

struct AscendingChain{T} <: AbstractRotatorChain{T}
  x::Vector{T}
end


function Base.size(C::AbstractRotatorChain)
    N = length(C.x)
    (N+1, N+1)
end

Base.extrema(C::AbstractRotatorChain) = extrema(idx.(C.x))

Base.adjoint(A::AscendingChain) = DescendingChain(reverse(adjoint.(A.x)))
Base.adjoint(A::DescendingChain) = AscendingChain(reverse(adjoint.(A.x)))



function Base.getindex(A::AscendingChain{T}, i, j) where {T}
    N = length(A.x)
    S = eltype(first(A.x))

    if i == N + 1
        i -= 1; j -= 1
        ci, si = vals(A.x[N+1-i])

        j == N && return conj(ci)

        s = -S(si)
        for l in (j+1):(i-1)
            cj, sjj = vals(A.x[N+1-l])
            s *= -sjj
        end
        cj, sj = j > 0 ? vals(A.x[N+1-j]) : (one(S), zero(S))
        return s * conj(cj)

    else


        j > i+1 && return zero(S) # ascending chain's aree lower Hessian

        ci, si = vals(A.x[N+1-i])

        if j == i + 1
            return S(si)
        elseif j == i
            cj, sj  = j > 1 ? vals(A.x[N+1-(j-1)]) : (one(S), zeros(S))
            return conj(cj)*ci
        end

        s = one(S)
        for l in j:(i-1)
            cl, sl = vals(A.x[N+1-l])
            s *= -sl
        end

        cl, sl = j > 1 ? vals(A.x[N+1 - (j-1)]) : (one(S), zero(S))
        return ci * s * conj(cl)
    end

end

function Base.getindex(A::DescendingChain, i, j)

    N = length(A.x)
    S = eltype(first(A.x))

    if i > j + 1 || i <= 0 || j <= 0
        return zero(S)
    end

    if j == N + 1
        j -= 1
        cj, sj = vals(A.x[j])

        i == N+1 && return conj(cj)

        i -= 1
        s::S = sj
        for l in (i+1):j-1
            cj, sj = vals(A.x[l])
            s *= S(sj)
        end
        ci, si = i > 0 ? vals(A.x[i]) : (one(S), zero(real(S)))
        return s * conj(ci)

    else


        cj, sj = vals(A.x[j])

        if i == j + 1

            return -S(sj)

        elseif j == i

            ci, si  = i >= 2 ? vals(A.x[i-1]) : (one(S), zero(real(S)))
            return conj(ci)*cj

        end

        s = one(S)
        for l in i:(j-1)
            cl, sl = vals(A.x[l])
            s *= sl
        end

        cl, sl = i >= 2 ? vals(A.x[i-1]) : (one(S), zero(real(S)))
        return cj * s * conj(cl)
    end
end


##################################################

## Various "passthrough" functions
## we have passthrough!(A,B) mutates A,B to return AB = BA
## If A or B is immutatble, return new one
# do not assume DescendingChain is 1 ... n
# left to  right; return U
function passthrough!(A::DescendingChain, U::AbstractRotator)
    i = idx(U)
    n = idx(A.x[1])
    N = idx(A.x[end])
    @assert n <= i < N
    l = (i-n) + 1
    U, A[l], A[l+1] = turnover(A[l], A[l+1], U)
    U
end

function passthrough!(A::AscendingChain, U::AbstractRotator)
    i = idx(U)
    n = idx(A.x[end])
    N = idx(A.x[1])
    @assert n < i <= N
    l = length(A.x)  - (i-n)
    U, A[l], A[l+1] = turnover(A[l], A[l+1], U)
    U

end

## right to left; return U
function passthrough!(U::AbstractRotator, A::DescendingChain)
    i = idx(U)
    n = idx(A.x[1])
    N = idx(A.x[end])
    @assert n < i <= N
    l = (i-n) + 1
    A[l-1], A[l], U = turnover(U, A[l-1], A[l])
    U
end

function passthrough!(U::AbstractRotator, A::AscendingChain)
    i = idx(U)
    n = idx(A.x[end])
    N = idx(A.x[1])
    @assert n <= i < N
    l = length(A.x)  - (i-n)
    A[l-1], A[l], U = turnover(U, A[l-1], A[l])
    U
end

# Need to check  bounds to ensure possible
function passthrough!(A::DescendingChain, B::AscendingChain)
    m, M = idx(A.x[1]), idx(A.x[end])
    n, N = idx(B.x[end]), idx(B.x[1])

    if M < N && n <= m
        for (i, U) in enumerate(B.x)
            B[i] = passthrough!(A, U)
        end
    else
        lA = length(A.x)
        for i in 1:lA
            A[lA+1-1] = passthrough!(A[lA+1-i], B)
        end
    end
    return nothing
end


function passthrough!(A::AscendingChain, B::DescendingChain)
    m, M = idx(A.x[1]), idx(A.x[end])
    n, N = idx(B.x[end]), idx(B.x[1])

    if n <= m-1 && N <= M
        for (i, U) in enumerate(B.x)
            B[i] = passthrough!(A, U)
        end
    else
        lA = length(A.x)
        for i in 1:lA
            j = lA - i + 1
            A[j] = passthrough!(A[j], B)
        end
    end
    return nothing
end


function passthrough!(D::SparseDiagonal, Asc::AscendingChain)
    for i in 1:length(Asc) # pass  Asc through D
        Asc[i] = passthrough(D, Asc[i], Val(:right))
    end
end

function passthrough!(Des::DescendingChain, D::SparseDiagonal)
    for i in length(Des):-1:1
        Des[i] = passthrough(D, Des[i], Val(:left))
    end
end

## could do two others...

## noop when Identity Diagonal
passthrough!(D::IdentityDiagonal, C::AbstractRotatorChain) =  nothing
passthrough!(C::AbstractRotatorChain, D::IdentityDiagonal) =  nothing

passthrough!(D::IdentityDiagonal, U::AbstractRotator{T}) where {T} = U
passthrough!(U::AbstractRotator{T}, D::IdentityDiagonal) where {T} = U



##################################################

## Pass a diagonal rotator through a chain
## [  D --> D [
##  [           [
## function passthrough(A::DescendingChain, D::DiagonalRotator, ::Val{:right})
##     XXX
## end

##  [ D --> D   [
## [          [
## function passthrough(A::AscendingChain, D::DiagonalRotator, ::Val{:right})
##     ## XXX write me, but not needed...
## end


## pass Di thorugh A; fuses with D
## in Twisted we have passthrough_phase with a tuple
## this is speficic to this task
function passthrough_phase!(Di::DiagonalRotator, A::DescendingChain, D::SparseDiagonal)
    ## We have two rules here
    ## D_i(alpha) * R_{i+1}(c,s) = R_{i+1}(c, conj(alpha)s) D_i(alpha)
    ## R_{i+1}(c, conj(alpha) s) = R_{i+1}(conj(alpha) c, s) D_{i+1}(alpha)
    ## so D_i(alpha) * R_{i+1}(c,s) = R_{i+1}(conj(alpha) c, s) D{{i+1}(alpha) D_i(alpha)
    ## Also D_i(alpha) * R_i(c,s) = R_i(c*alpha/conj(alpha), s) D_i(conj(alpha)
    i, n = idx(Di), length(A)
    alpha, _ = vals(Di)

    ## D_i will passthrough and can be merged into D
    ## we will have D_i(alpha) * D_{i+1}(alpha) * ... * D_j(alpha) which will only have alpha at i and conj(alpha) at j+1
    D[i] *= alpha


    while i < n
        j = i + 1
        c, s = vals(A[j])
        @assert j == idx(A[j])

        iszero(s) && break
        A[j] = Rotator(conj(alpha)*c, s, j)
        i = i + 1
    end

    D[i+1] *= conj(alpha)

    return nothing

end
