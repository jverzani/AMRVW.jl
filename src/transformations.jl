##
## Transformations are
## * givens rotations  Find [c s; -conj(s) conj(c)][a,b] = [x 0];
## * fuse U,V -> UV or UV,D (when complex real)
## * turnover: [   [ -->  [
##               [      [   [
## to this we add abstract "passthrough!" functions

##################################################

"""
    c,s,r  = givensrot(a,b)

where `[c s; -s conj(c)] * [a,b] = [r, 0]
"""
@inline function givensrot(a::T, b::T) where {T <: Real}
    G, r = givens(a,b,1,2)
    s = r >= 0 ? one(T) : -one(T)   #sign(r)
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
## For general rotation, we have no phase consideration, so return an identity rotator
@inline function fuse(a::AbstractRotator{T,S}, b::AbstractRotator{T,S}) where {T,S <: Real}
    #    idx(a) == idx(b) || error("can't fuse")
    i = idx(a)

    c1, s1 = vals(a)
    c2, s2 = vals(b)

    u = c1 * c2 - s1 * conj(s2)
    v = c1 * s2 + s1 * conj(c2)

    u,v = polish_givens(u,v)
    return Rotator(u,v,i), IdentityRotator{T, S}(i)

end

## U * V -> UV [ alpha 0; 0 conj(alpha)]

## For ComplexRealRotator, the result of a*b will not have a real sign
## we output by rotating by alpha.
## We have U*V = (UV) * D
## To get D * (UV) a flip is needed
@inline function fuse(a::Rotator{T,S}, b::Rotator{T,S}) where {T, S <: Complex{T}}

    i = idx(a)
    #    idx(a) == idx(b) || error("can't fuse")
    c1, s1 = vals(a)
    c2, s2 = vals(b)

    u = c1 * c2 - s1 * s2
    v = c1*s2 + s1 * conj(c2)

    s = norm(v)
    alpha =  iszero(v) ? one(Complex{T}) : v/s
    c = u * alpha

    c, s = polish_givens(c, s)

    Rotator(c, s, i), DiagonalRotator(conj(alpha), i)
end

## Fuse for diagonal matrices
fuse(D::IdentityRotator, U) = U, D
function fuse(D::DiagonalRotator, U)
    alpha,_ = vals(D)
    i = idx(D)
    @assert i == idx(U)
    c, s = vals(U)
    Rotator(c*alpha, s, i)
end

##################################################
# Turnover: Q1    Q3   | x x x |      U1
#              Q2    = | x x x | = U3    U2
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
# V1' *     * M  = | 1  0 0   |
#       U1'        | 0  c s   |
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

    if iszero(UVW21) && iszero(UVW31)
        if iszero(s2*conj(c3)) ## this is the case of diagonal matrix Q1*Q2*Q3
            c4 = one(S)
            c5 = c1*c3 - c2*s1*s3
            c6 = c2
            return c4, zero(T), c5, zero(T), c6, zero(T)
        end
    end

    c4, s4, r4 = givensrot(UVW21, UVW31)

    c4, s4 = conj(c4), -s4  # conjugate, as [c4 s4; -s4 conj(s4)]*[a,b] = [r,0]
#    c4, s4 =  polish_givens(c4, s4)

    UVW11 = c1*c3 - c2*s1*s3
    c5, s5::T = approx_givensrot(UVW11, real(r4))
    c5, s5 = conj(c5), -s5
#    c5, s5 =  polish_givens(c5, s5)



    ##  use s1*s2 = s5 * s6 if we can
    ## get a from  M  = V1' * U1' * U *  V  *  W; a = M[3,3]
    if !iszero(s5)
        a = conj(c1) * s2 * s4 + conj(c2) * c4
        b = s1*s2/s5
        c6, s6 = approx_givensrot(a, b)
    else
        a = c2*conj(c1)
        b = -s2
        c6, s6 = givensrot(a, b)
    end
#    c6, s6 = approx_givensrot(a, b)
#    c6, s6 =  polish_givens(c6, s6)

    return (c4, s4, c5, s5, c6, s6)

end


##
##  Turnover interface for rotators
##
function turnover(Q1::R1t,
                  Q2::R2t,
                  Q3::R3t) where {R1t, R2t, R3t}

    c1, s1 = vals(Q1); c2, s2 = vals(Q2); c3,s3 = vals(Q3)
    i,j,k = idx(Q1), idx(Q2), idx(Q3)

    # @assert i == k && (abs(j-i) == 1)

    c4,s4,c5,s5,c6,s6 = _turnover(c1,s1, c2,s2, c3,s3)




    R1 = R3t(c4, s4, j)
    R2 = R2t(c5, s5, i)
    R3 = R1t(c6, s6, j)


    R1, R2, R3

end
