function deflate_leading_zeros(ps::Vector{S}) where {S}
    ## trim any 0s from the end of ps
    N = findlast(!iszero, ps)
    K = findfirst(!iszero, ps)

    N == 0 && return(zeros(S,0), length(ps))
    # XXX  should  we make  a view?
    ps = ps[K:N]
    ps, K-1
end

## a pencil has
## v[1] = a0
## v[i+1] + w[i] = ai
## w[n] = an
# we want to trim both ends here
function deflate_leading_zeros(vs::Vector{S}, ws::Vector{S}) where {S}

    # get N
    N = length(vs)
    if iszero(ws[end])
        N -= 1
        while N >= 1
            ai = vs[N+1] + ws[N]
            if !iszero(ai)
                ws[N] = ai
                break
            end
            N -= 1
        end
    end
    K = 1
    if iszero(vs[1])
        K = 2
        while K <= N
            ai = vs[K] + ws[K-1]
            if !iszero(ai)
                vs[K] = ai
                break
            end
            K += 1
        end
    end
    ## Should  this be  a  view?
    vs[K:N], ws[K:N], K-1
end

##################################################
##
## solve degree 2 or less case
function solve_simple_cases(ps::Vector{T}) where {T}
    S = T <: Complex ? T : Complex{T}

    N = length(ps)

    if N <= 1
        return S[]
    elseif N == 2
        return S[-ps[1]/ps[2]]
    elseif N == 3
        c,b,a = ps
        r1,i1, r2,i2 = quadratic_equation(a,b,c)
        return S[complex(r1, i1), complex(r2, i2)]
    end

    error("solve_simple_cases: too many coefficients")
end


#
function quadratic_equation(a::T, b::T, c::T) where {T <: Real}
    qdrtc(a, -b/2, c)
end

## make more robust
function quadratic_equation(a::Complex{T}, b::Complex{T}, c::Complex{T}) where {T}

    d = sqrt(b^2 - 4*a*c)
    e1 = (-b + d)/(2a); e2 = (-b-d)/(2a)
    return (real(e1), imag(e1), real(e2), imag(e2))

end

## Kahan quadratic equation with fma
##  https://people.eecs.berkeley.edu/~wkahan/Qdrtcs.pdf

## solve ax^2 - 2bx + c
function qdrtc(a::T, b::T, c::T) where {T <: Real}
    # z1, z2 roots of ax^2 - 2bx + c
    d = discr(a,b,c)  # (b^2 - a*c), as 2 removes 4

    if d <= 0
        r = b/a  # real
        s = sqrt(-d)/a #imag
        return (r,s,r,-s)
    else
        r = sqrt(d) * (sign(b) + iszero(b)) + b
        return (r/a, zero(T), c/r, zero(T))
    end
end

## more work could be done here.
function discr(a::T,b::T,c::T) where {T}
    pie = 3.0 # depends on 53 or 64 bit...
    d = b*b - a*c
    e = b*b + a*c

    pie*abs(d) > e && return d

    p = b*b
    dp = muladd(b,b,-p)
    q = a*c
    dq = muladd(a,c,-q)

    (p-q) + (dp - dq)
end

function qdrtc(b::T, c::T) where {T <: Real}
    d = b*b - c
    if d <= 0
        r = b  # real
        s = sqrt(-d) #imag
        return (r,s,r,-s)
    else
        r = sqrt(d) * (sign(b) + iszero(b)) + b
        return (r, zero(T), c/r, zero(T))
    end
end
