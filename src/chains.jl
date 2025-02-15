## A chain is an arrangement of rotators between n and N
## Three types are specialized: Descending, Ascending, Twisted

#### Group rotators

abstract type AbstractRotatorChain{T,S} end

Base.eltype(A::AbstractRotatorChain) = eltype(A.x)
Base.length(A::AbstractRotatorChain) = length(A.x)

Base.@propagate_inbounds Base.getindex(A::AbstractRotatorChain, i::Int) = getindex(A.x, i)
Base.@propagate_inbounds Base.setindex!(A::AbstractRotatorChain, X, inds...) = setindex!(A.x, X, inds...)
Base.deleteat!(A::AbstractRotatorChain, j) = deleteat!(A.x, j)

Base.iterate(A::AbstractRotatorChain) = iterate(A.x)
Base.iterate(A::AbstractRotatorChain, st) = iterate(A.x, st)
*(A::AbstractRotatorChain, M::Array) = Vector(A) * M
*(M::Array, A::AbstractRotatorChain) = M * Vector(A)
LinearAlgebra.lmul!(A::AbstractRotatorChain, M::Array) = lmul!(Vector(A), M)
LinearAlgebra.rmul!(M::Array, A::AbstractRotatorChain) = rmul!(M, Vector(A))


function Base.size(C::AbstractRotatorChain)
    n, N = extrema(C)
    (N+1, N+1)
end

Base.extrema(C::AbstractRotatorChain) = extrema(idx.((C.x[1], C.x[end])))

function Base.Matrix(C::AbstractRotatorChain)
    S =  eltype(first(C.x).c)
    M = diagm(0 => ones(S, size(C)[2]))
    lmul!(C, M)
end

## Types of chains
struct DescendingChain{T, S, V <: AbstractVector{Rotator{T,S}}} <: AbstractRotatorChain{T,S}
  x::V
end

## We store ascending chain in order n, n-1, n-2, ...
struct AscendingChain{T,S, V <: AbstractVector{Rotator{T,S}}} <: AbstractRotatorChain{T,S}
  x::V
end

const MonotonicChain = Union{DescendingChain, AscendingChain}

Base.view(D::MonotonicChain, inds...) = DescendingChain(view(D.x, inds...))
Base.copy(D::DescendingChain) = DescendingChain(copy(D.x))
Base.copy(D::AscendingChain) = AscendingChain(copy(D.x))
Base.pop!(D::MonotonicChain) = pop!(D.x)
Base.popfirst!(D::MonotonicChain) = popfirst!(D.x)
Base.push!(D::MonotonicChain) = push!(D.x)
Base.pushfirst!(D::MonotonicChain) = pushfirst!(D.x)

Base.Vector(D::MonotonicChain) = D.x


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


## Structure to hold a twisted representation
## pv is the position vector of length  n-1; we leave parametric so views can be used
struct TwistedChain{T,S, V <: AbstractVector{Rotator{T,S}}, PVt <: AbstractVector{Symbol}} <: AbstractRotatorChain{T,S}
  x::V
  pv::PVt
end

function Base.view(Tw::TwistedChain, inds::UnitRange)
    start, stop = inds.start, inds.stop
    TwistedChain(view(Tw.x, inds), view(Tw.pv, start:stop-1))
end

Base.copy(Tw::TwistedChain) = TwistedChain(copy(Tw.x), copy(Tw.pv))


## Twisted Chains requires a bit more effort
## Constructor from a twisted vector of rotators
function TwistedChain(xs::Vector{Rotator{T,S}}) where {T,S}
    if length(xs) > 0
        sigma = idx.(xs)
        ps = position_vector(sigma)
        TwistedChain(xs[sortperm(sigma)], ps)
    else
        TwistedChain(xs, Symbol[])
    end
end


function position_vector_indices(pv)

    inds = [1]

    for (j,d) in enumerate(pv)
        i = j + 1
        d == :left ? push!(inds, i) : pushfirst!(inds, i)
    end

    inds
end



function Base.Vector(Tw::TwistedChain)
    isempty(Tw.x) && return Tw.x

    return Tw.x[position_vector_indices(Tw.pv)]
end

## Find position vector from a permutation of 1...n
## If i is to the left of i+1, set  ps[i] = :left
## if i is to the right if i+1, set ps[i] = :right
## e.g. 4,3,2,1 -> r,r,r (descending)
##      1,2,3,4 -> l,l,l (ascending)
##      1,3,2,4 -> rlr  (CMV)
function position_vector(sigma)
    ## sigma a perm of 1....n
    sigma = sigma .- minimum(sigma) .+ 1
    n = length(sigma)
    n <= 1 && return Symbol[]
    ps =  repeat([:nothing],n-1)
    for (i,v) in enumerate(sigma)
       if v == 1

            if ps[1] == :nothing
              ps[1] = :left
            end
       elseif  v == n
            if ps[end] == :nothing
                ps[end] = :right
            end
       else
          if ps[v-1] == :nothing
                ps[v-1] =  :right
            end
            if ps[v]  == :nothing
                ps[v] = :left
            end
        end
    end
    ps
end

## note this is expensive
Base.adjoint(A::TwistedChain) = TwistedChain(reverse(adjoint.(Vector(A))))


# In a Twisted chain from n to N, pop off the N rotator
# and modify the position vector
function Base.pop!(A::TwistedChain)
    pop!(A.x), (isempty(A.pv) ? :nothing : pop!(A.pv))
end

## This should be an iterator
##   Indices of Tw.x  to  pass twisted chain from  right side to left side
## grab right terms from bottom
function iterate_rl(pv)
    # indices in Tw.x to move object from l to r
    N = length(pv)+1
    inds =  Int[]
    for (i,pi)  in enumerate(reverse(pv))
        if  pi  ==  :right
            j  = N  - i
            push!(inds, j+1)
        end
    end
    for i  in 1:N
        if  !(i in inds)
            push!(inds, i)
        end
    end
    inds
end


function iterate_lr(pv)
    # indices to move  object from left side to right
      # indices in Tw.x to move object from l to r
    N = length(pv)+1
    inds =  Int[]
    for (i,pi)  in enumerate(reverse(pv))
        if  pi  ==  :left
            j = N - i
            push!(inds, j+1)
        end
    end
    for i  in 1:N
        if  !(i in inds)
            push!(inds, i)
        end
    end
    inds
end


## Fish out of M the rotator with idx i
## return index in Ms.x and the rotator
function iget(Ms::TwistedChain, i)

    m,M = extrema(Ms)
    j = i - (m-1)
    U = Ms[j]

    j, U
end



## From a twisted chain and index i, return ascending part of i-1, i-2, ....
function ascending_part(Ms::TwistedChain, i, Δ=1)
    n, N = extrema(Ms)
    pv = Ms.pv


    start = i
    stop = i-1

    (i <  n || i > N + 1) && return  AscendingChain(view(Ms.x, 1:0)) #start:stop # empyt inds #(inds, Asc)
    while i > Δ
        i -= 1
        if i >= n && (i == N || pv[i-n+1] == :right) ## i-n+1 seems right when n != 1
            start = i
        else
            break
        end
    end

    AscendingChain(view(Ms.x, (stop-n+1):-1:(start-n+1)))

end


## Return (i+1):j where Ui, Ui+1, ... Uj is descending chain
function descending_part(Ms::TwistedChain, i, Δ=extrema(Ms)[2])
    n, N = extrema(Ms)
    pv = Ms.pv

    start = i+1
    stop = i
    (i <  n - 1 || i >= N) && return  DescendingChain(view(Ms.x, 1:0))
    while i < Δ
        i += 1
        if i <= N && (i <= n || pv[i-1-(n-1)] == :left)
            stop = i
        else
            break
        end
    end
    DescendingChain(view(Ms.x, (start-n+1):(stop-n+1)))

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
    #@assert n <= i < N
    l = (i-n) + 1
    U, A[l], A[l+1] = turnover(A[l], A[l+1], U)
    U
end

function passthrough!(A::AscendingChain, U::AbstractRotator)
    i = idx(U)
    n = idx(A.x[end])
    N = idx(A.x[1])
    #@assert n < i <= N
    l = length(A.x)  - (i-n)
    U, A[l], A[l+1] = turnover(A[l], A[l+1], U)
    U

end

## right to left; return U
function passthrough!(U::AbstractRotator, A::DescendingChain)
    i = idx(U)
    n = idx(A.x[1])
    N = idx(A.x[end])
    #@assert n < i <= N
    l = (i-n) + 1
    A[l-1], A[l], U = turnover(U, A[l-1], A[l])
    U
end

function passthrough!(U::AbstractRotator, A::AscendingChain)
    i = idx(U)
    isempty(A.x) && return U
    n, N = extrema(A)
    #@assert n <= i < N
    l = length(A.x)  - (i-n)
    A[l-1], A[l], U = turnover(U, A[l-1], A[l])
    U
end

# Need to check  bounds to ensure possible
function passthrough!(A::DescendingChain, B::AscendingChain)
    m, M = extrema(A)
    n, N = extrema(B)
    if n >= m && N < M
        for (i, U) in enumerate(B.x)
            B[i] = passthrough!(A, U)
        end
    else
        lA = length(A.x)
        for i in 1:lA
            j = lA + 1 - i
            A[j] = passthrough!(A[j], B)
        end
    end
    return nothing
end


function passthrough!(A::AscendingChain, B::DescendingChain)
    m, M = idx(A.x[end]), idx(A.x[1])
    n, N = idx(B.x[1]), idx(B.x[end])

    if m+1 <= n && M >= N
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

## Pass chain (Des, Asc, Twisted) through a chain (Des, Asc)
## Move twisted through down R <- L
## Pass chain (Des, Asc, Twisted) through a chain (Des, Asc)
## Move twisted through down R <- L
function passthrough!(A::DescendingChain, B::TwistedChain)
    length(B) == 0 && return nothing
    n,N = extrema(A)
    m,M = extrema(B)
    #@assert (m >= n && M < N)
    for i in iterate_rl(B.pv)
        B.x[i] = passthrough!(A, B.x[i])
    end
end

function passthrough!(A::AscendingChain, B::TwistedChain)
    length(B) == 0 && return nothing
    n,N = extrema(A)
    m,M = extrema(B)
    #@assert (m > n && M <= N)
    for i in iterate_rl(B.pv)
        B.x[i] = passthrough!(A, B.x[i])
    end
end

function passthrough!(B::TwistedChain, A::DescendingChain)
    length(B) == 0 && return nothing
    n,N = extrema(A)
    m,M = extrema(B)
    #@assert (m > n && M <= N)
    for i in iterate_lr(B.pv)
        B.x[i] = passthrough!(B.x[i], A)
    end
end

function passthrough!(B::TwistedChain, A::AscendingChain)
    length(B) == 0 && return nothing
    n,N = extrema(A)
    m,M = extrema(B)
    #@assert (m >= n && M < N)
    for i in iterate_lr(B.pv)
        B.x[i] = passthrough!(B.x[i], A)
    end
end
