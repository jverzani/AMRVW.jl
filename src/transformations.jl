##
## Transformations are
## * givens rotations  Find [c s; -conj(s) conj(c)][a,b] = [x 0];
## * fuse U,V -> UV or UV,D (when complex real)
## * turnover: [   [ -->  [
##               [      [   [
## to this we add abstract "passthrough" functions

##################################################

"""
    c,s,r  = givensrot(a,b)

where `[c s; -s conj(c)] * [a,b] = [r, 0]
"""
@inline function givensrot(a::T, b::T) where {T <: Real}
    G, r = givens(a,b,1,2)
    s = sign(r)
    s*G.c, s*G.s, s*r
end

## This leaves s real
## a, b= rand(Complex{Float64},2)
## G,r = givens(b,a, 1,2)
## c, s = G.s, real(G.c)
## [c s; -s conj(c)] * [a,b]
function givensrot(a::Complex{T},b::Complex{T}) where {T <: Real}
    G, r = givens(b, a, 1, 2)
    G.s, real(G.c), r
end

givensrot(a::Complex{T},b::T) where {T <: Real} =  givensrot(a, complex(b, zero(T)))

## This is suggested in the book. A bit slower, but only one function needed...
## function givensrotX(a::S, b::T) where {T <: Real, S <: Union{T, Complex{T}}}
##     m = max(norm(a), norm(b))
##     x,y = a/m, b/m
##     r = sqrt(norm(x)^2 + y^2)
##     c = x/r
##     s = y/r
##     r = m*r
##     return conj(c),real(s),r
## end

# avoid square root when norm is near 1
#
# 5.1. Backward Stability of the Turnover. Suggests
# this modification to c,s to ensure unitarity
# Assuming that c and s have been produced by valid formulas, the
# magnitude of the correction will be extremely tiny, on the order of
# the unit roundoff. This seemingly trivial correction is crucial. If it
# is neglected, the core transformations will gradually drift away from
# unitarity over the course of many operations, and the algorithm will
# fail. If it is not neglected, the programmer has a good chance of
# producing an accurate turnover.
function approx_givensrot(a::S,b::T) where {S,T}
    delta::T = abs(a*conj(a)) + b^2 - one(T)
    r::T = one(T) - delta/(one(T)+one(T))  # taylor approx for 1/sqrt(1 + delta) when delta ~ 0
    c = a*r
    s = b*r
    return conj(c),s

end

## extra polish; not used
function polish_givens(c::S, s::T) where {S, T}

    delta = real(conj(c)*c) + s*s - one(T)
    r = one(T) -  delta/(one(T) + one(T))
    r * c, r * s

end


##################################################
## Fuse
## fuse combines two rotations, a and b, into one,
## For general rotation, we have no phase consideration
@inline function fuse(a::AbstractRotator{T}, b::AbstractRotator{T}) where {T}
    #    idx(a) == idx(b) || error("can't fuse")
    i = idx(a)

    c1, s1 = vals(a)
    c2, s2 = vals(b)

    u = c1 * c2 - s1 * conj(s2)
    v = c1 * s2 + s1 * conj(c2)

    return Rotator(u,v,i)

end

## U * V -> UV [ alpha 0; 0 conj(alpha)]

## For ComplexRealRotator, the result of a*b will not have a real sign
## we output by rotating by alpha.
## We have U*V = (UV) * D
## To get D * (UV) a flip is needed
@inline function fuse(a::ComplexRealRotator{T}, b::ComplexRealRotator{T}) where {T}

    i = idx(a)
    #    idx(a) == idx(b) || error("can't fuse")
    c1, s1 = vals(a)
    c2, s2 = vals(b)

    u = c1 * c2 - s1 * s2
    v = c1*s2 + s1 * conj(c2)

    s = norm(v)
    alpha =  iszero(v) ? one(Complex{T}) : v/s
    c = u * alpha


    Rotator(c, s, i), DiagonalRotator(conj(alpha), i)
end

## Fuse for diagonal matrices
fuse(D::IdentityDiagonalRotator, U) = U
function fuse(D::DiagonalRotator{T}, U) where {T}
    alpha,_ = vals(D)
    i = idx(D)
    @assert i == idx(U)
    c, s = vals(U)
    Rotator(c*alpha, s, i)
end

##################################################
# Turnover: Q1    Q3   | x x x |      Q1
#              Q2    = | x x x | = Q3    Q2
#                      | x x x |
#
# This is the key computation once matrices are written as rotators
# This function returns
# c4,s4,c5,s5,c6,s6

# following pp15-17 in book
# and the paper
# M = U    W
#       V
# we find U1 by computing a givens rotation to clear [3,1]
#
#     * M =| a x x |
# U1'      | b x x |
#          | 0 x x |
# Then use fact that a^2 + b^2 ~ 1 to use approx givens to find V1 to clear [2,1]
# V1' *     * M  = | 1 0 0 |
#       U1'        | 0 c s |
#                  | 0 -s cbar|
# The 1 in [1,1] must be 1 and not -1;  so the  lower part is a rotator
# from here [2:3,2:3] is unitary, we pick [3,3] for c and use
# the invariant s1*s2 = s5*s6 to find s6 unless s5=0
#
# question: polishing slows things down, and doesn't seem necessary,
# as approx_givensrot does this; but does make Wilkinson case a bit better
@inline function _turnover(c1::S, s1::T, c2::S, s2::T, c3::S, s3::T) where {S,T}


    ## UVW has first row [UVW11, UVW21, UVW31]
    UVW21 = -c2 * s3 * conj(c1) - c3*s1
    UVW31 = s2 * s3

    c4, s4, r4 = givensrot(UVW21, UVW31)
    c4, s4 = conj(c4), -s4  # conjugate, as [c4 s4; -s4 conj(s4)]*[a,b] = [r,0]
#    c4, s4 =  polish_givens(c4, s4)


    UVW11 = c1*c3 - c2*s1*s3
    c5, s5::T = approx_givensrot(UVW11, real(r4))
    c5, s5 = conj(c5), -s5
#    c5, s5 =  polish_givens(c5, s5)



    ##  use s1*s2 = s5 * s6 if we can
    ## get a from  M  = V1' * U1' * U *  V  *  W; a = M[3,3]
    a = conj(c1) * s2 * s4 + conj(c2) * c4
    if !iszero(s5)
        c6, s6 = approx_givensrot(a, s1*s2/s5)
    else
        b = -c5*s4*conj(c2) + s2 * c5 * conj(c1)*conj(c4) + s2*s1*s5#; M[3,2] when s5=>0
        c6, s6 = approx_givensrot(a, b)
    end
#    c6, s6 = approx_givensrot(a, b)
#    c6, s6 =  polish_givens(c6, s6)

    return (c4, s4, c5, s5, c6, s6)

end


##
##  Turnover interface for rotators
##
function turnover(Q1::AbstractRotator,
                  Q2::AbstractRotator,
                  Q3::AbstractRotator)

    c1, s1 = vals(Q1); c2, s2 = vals(Q2); c3,s3 = vals(Q3)
    i,j,k = idx(Q1), idx(Q2), idx(Q3)
    # @assert i == k && (abs(j-i) == 1)

    c4,s4,c5,s5,c6,s6 = _turnover(c1,s1, c2,s2, c3,s3)
    R1 = Rotator(c4, s4, j)
    R2 = Rotator(c5, s5, i)
    R3 = Rotator(c6, s6, j)

    # we have Q1*Q2*Q3 = R1*R2*R3
    R1, R2, R3

end

## function turnoverXXX(Q1::Rt,
##                   Q2::Rt,
##                   Q3::Rt) where {Rt}

##     c1, s1 = vals(Q1); c2, s2 = vals(Q2); c3,s3 = vals(Q3)
##     i,j,k = idx(Q1), idx(Q2), idx(Q3)
##     # @assert i == k && (abs(j-i) == 1)

##     c4,s4,c5,s5,c6,s6 = _turnover(c1,s1, c2,s2, c3,s3)
##     R1::Rt = Rotator(c4, s4, j)
##     R2::Rt = Rotator(c5, s5, i)
##     R3::Rt = Rotator(c6, s6, j)

##     # we have Q1*Q2*Q3 = R1*R2*R3
##     R1, R2, R3

## end

##################################################

## Various "passthrough" functions
## These passthrough the chains
## XXX use passthrough!; remove these calls
function passthrough(A::DescendingChain, U::AbstractRotator, ::Val{:right})
    return passthrough!(A, U)
end


function passthrough(A::DescendingChain, U::AbstractRotator, ::Val{:left})
    return passthrough!(U, A)
end


function passthrough(A::AscendingChain, U::AbstractRotator, ::Val{:right})
    return passthrough!(A, U)
end


function passthrough(A::AscendingChain, U::AbstractRotator, ::Val{:left})
    return passthrough!(U, A)
end

## XX these to replace the 4 above
# left to  right; return U
# do not assume DescendingChain is 1 ... n
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


## passthrough
## Pass a rotator through a diagonal matrix with phase shifts
## D U -> U' D' (right, as in U starts on right)
## U D -> D' U' (left)
@inline function passthrough(D::SparseDiagonal, U::ComplexRealRotator{T}, ::Val{:right}) where {T}
    i = idx(U)
    c,s = vals(U)

    alpha, beta = D[i], D[i+1]
    beta1 = alpha * conj(beta)

    D[i] = beta
    D[i+1] = alpha
    return Rotator(beta1 * c, s, i)
end

## U D -> D U
@inline function passthrough(D::SparseDiagonal, U::ComplexRealRotator{T}, ::Val{:left}) where {T}
    return passthrough(D, U, Val(:right))
    i = idx(U)
    c,s = vals(U)

    alpha, beta = D[i], D[i+1]
    beta1 = alpha/beta

    D[i] *= conj(beta1)
    D[i+1] *= beta1
    return Rotator(beta1 * c, s, i)
end

passthrough(D::IdentityDiagonal, U::AbstractRotator{T}, args...) where {T} = U

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


## U D -> D U
## Identity
@inline function passthrough(D::IdentityDiagonal, U::AbstractRotator, dir)
    U
end



## Pass a diagonal rotator through a chain
## [  D --> D [
##  [           [
function passthrough(A::DescendingChain, D::DiagonalRotator, ::Val{:right})
end

##  [ D --> D   [
## [          [
function passthrough(A::AscendingChain, D::DiagonalRotator, ::Val{:right})
    ## XXX write me, but not needed...
end


## [ D          [      D     [        D     [
##     [   -->    [ D    -->   [    D   -->   [     * [alpha, I..., conj(alpha]
##       [            [          [ D             [
function passthrough(A::DescendingChain, D::SparseDiagonal, Di::DiagonalRotator, ::Val{:left})
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


function passthrough(A::AscendingChain, D::DiagonalRotator, ::Val{:left})
    ## XXX write me
end
