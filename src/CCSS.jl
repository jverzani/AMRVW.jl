## ComplexComplexRotator case
## XXX Deprecate this, it isn't any faster nor easier.
# ##################################################
# ## We use two complex, rather than 3 reals here.
# ## Will be basically the ame storage, as we don't need to include a D, but not quite (12N, not 11N)

struct ComplexComplexRotator{T} <: AbstractRotator{T}
c::Complex{T}
s::Complex{T}
i::Int
end

function adjoint(r::ComplexComplexRotator)
    ComplexComplexRotator(conj(r.c), -r.s, r.i)
end

Base.eltype(::Type{ComplexComplexRotator{T}}) where {T} = complex(T)
Base.one(::Type{ComplexComplexRotator{T}}) where {T} = ComplexComplexRotator(complex(one(T), zero(T)), complex(zero(T), zero(T)), 0)

function random_rotator(::Type{ComplexComplexRotator{T}}, i) where {T}
    a,b = rand(Complex{T}, 2)
    c,s,_ = givensrot(a,b)
    ComplexComplexRotator(c,s,i)
    ##ns = norm(s)
    ## alpha = conj(s)/ns
    ##ComplexComplexRotator(c*alpha,complex(ns),i)
end




###################################################
## The type and convert method

# ComplexReal Single Shift, no pencil, not twisted
mutable struct ComplexComplex_SingleShift_NoPencil_NotTwisted{T} <: FactorizationType{T, Val{:SingleShift}, Val{:NoPencil}, Val{:NotTwisted}}

N::Int
POLY::Vector{Complex{T}}
Q::Vector{ComplexComplexRotator{T}}
D::Vector{Complex{T}}
Ct::Vector{ComplexComplexRotator{T}}  # We use C', not C here
B::Vector{ComplexComplexRotator{T}}
#
REIGS::Vector{T}
IEIGS::Vector{T}
# reusable storage
U::ComplexComplexRotator{T}
Ut::ComplexComplexRotator{T}
A::Matrix{Complex{T}}    # for parts of A = QR
R::Matrix{Complex{T}}    # temp storage, sometimes R part of QR
e1::Vector{T}   # eigen values e1, e2, store as (re,imag)
e2::Vector{T}
ray::Bool
ctrs::AMRVW_Counter
end

function Base.convert(::Type{FactorizationType{T, Val{:SingleShift}, Val{:NoPencil}, Val{:NotTwisted}}}, ps::Vector{Complex{T}}; ray=true) where {T}

    N = length(ps) - 1

    ComplexComplex_SingleShift_NoPencil_NotTwisted(N, ps,
                                                Vector{ComplexComplexRotator{T}}(undef,N), #Q
                                                ones(Complex{T}, N+1), # D
                                                Vector{ComplexComplexRotator{T}}(undef,N), #Ct
                                                Vector{ComplexComplexRotator{T}}(undef,N), #B
                                                zeros(T, N),  zeros(T, N), #EIGS
    one(ComplexComplexRotator{T}), one(ComplexComplexRotator{T}), #U, Ut
    zeros(Complex{T}, 2, 2),zeros(Complex{T}, 3, 2), # A R
    zeros(T,2), zeros(T,2),
    #    true,  # true for Wilkinson, false for Rayleigh.Make adjustable!
    ray,
    AMRVW_Counter(0,1,N-1, 0, N-2)
    )
end

##################################################
### Factorization


function Q_factorization(state::FactorizationType{T, St, P, Tw}) where {T, St, P, Tw}
    N = state.N
    Q = state.Q
    S = St == Val{:SingleShift} ? Complex{T} : T # type of cosine term in rotator
    zS, oS,zT,oT = zero(S), one(S), zero(T), one(T)
    @inbounds for ii = 1:(N-1)
        Q[ii] = Rotator(zS, oS, ii)

    end
    Q[N] = Rotator(oS, zS, N) # not needed, but helps us in diagonal_block
end


function R_factorization(xs::Vector{Complex{T}}, Ct, B, side=Val{:left}) where {T}


    N = length(xs) - 1 # ps has 1 in it?

    c,s,tmp = givensrot(xs[N], -one(Complex{T}))

    C = Rotator(c, s, N)
    Ct[1] = C'
    B[N] = Rotator(s, -c, N)


    @inbounds for i in (N-1):-1:1
        c,s,tmp = givensrot(xs[i], tmp)
        C = Rotator(c, s, i)
        Ct[N+1-i] = C'
        B[i] = C
    end

    return nothing
end

## populate (D,Ct,B)
function init_triu(state::FactorizationType{T, Val{:SingleShift}, Val{:NoPencil}, Tw}, decompose) where {T, Tw}
    xs = decompose(state.POLY)
    R_factorization(xs, state.Ct, state.B, Val(:right))
end

##################################################
## transformations


## givens rotation
##################################################
# Compute Givens rotation zeroing b
#
# G1 [ ar + i*ai ] = [ nrm ]
# G1 [    b      ] = [     ]
#
# all variables real (nrm complex)
# returns (copmlex(cr, ci), s) with
# u=complex(cr,ci), v=s; then [u -v; v conj(u)] * [complex(ar, ai), s] has 0
#
function givensrot(a::Complex{T},b::Complex{T}) where {T <: Real}

    G,r = givens(a, b, 1, 2)
    return G.c, G.s, r

end


@inline function choose_phase(c::Complex{T}, s::Complex{T}, r) where {T <: Real}
    alpha = conj(r)/norm(r)
    alpha .* (c, s, r)
end


@inline function passthrough(D::SparseDiagonal, U::ComplexComplexRotator{T}, ::Val{:right}) where {T}
    i = idx(U)
    c,s = vals(U)

    alpha, beta = D[i], D[i+1]
    beta1 = alpha/beta

    return Rotator(c, s*beta1, i)
end
@inline function passthrough(D::SparseDiagonal, U::ComplexComplexRotator{T}, ::Val{:left}) where {T}
    XXX("write this")
end


## # D U -> U D

###################################################
## Bulge related
# absorb_UT

function absorb_Ut(state::FactorizationType{T, Val{:SingleShift}, P, Val{:NotTwisted}}) where {T, P}
    i = idx(state.U)
    state.Q[i] = fuse(state.U', state.Q[i])

    nothing
end


## XXX pencil case will need modification
function passthrough_triu(state::FactorizationType{T, Val{:SingleShift}, P, Tw}, ::Val{:right}) where {T, P, Tw}

    i, N = idx(state.U), state.N
    if  i <= state.ctrs.tr
        # this leverages Ct*B*U = I*U when Ct, B not touched
        # So we need no modify U, and Ct and B will change together
        _ = passthrough(Val(:descending), state.B, state.U, Val(:right))

        for k in 0:1
            j = i+k
            a,b = vals(state.B[j])
            R =  Rotator(a, b, idx(state.Ct[N+1-(j)]))
            state.Ct[N+1-(j)] = R'
        end
    else
        state.U = passthrough(Val(:descending), state.B, state.U, Val(:right))
        state.U = passthrough(Val(:ascending),  state.Ct, state.U, Val(:right))

    end
    state.ctrs.tr -= 1
    nothing
end


function passthrough_Q(state::FactorizationType{T, Val{:SingleShift}, P, Val{:NotTwisted}}) where {T,P}

    state.U = passthrough(Val(:diagonal), state.D, state.U, Val(:right))

    i = idx(state.U)
    if i < state.ctrs.stop_index
        state.U = passthrough(Val(:descending), state.Q, state.U, Val(:right))
        return false
    else
        state.Q[i] = fuse(state.Q[i], state.U)
        return true
    end
end

##################################################
## Deflation
##
# deflate a term
# deflation for ComplexReal is different, as
# we replace Qi with I and move diagonal part into D
function deflate(state::FactorizationType{T,Val{:SingleShift},P, Val{:NotTwisted}}, k) where {T, P}

    # when we deflate here we want to leave Q[k] = I and

    # move Dk matrix over to merge with D
    # Qi       Qi            Qi           Qi
    #   Qj   ->   I Dj  -->    I Dj   -->   I Dj
    #     Qk        QK           Qk Dj        Qk Dj    fuse Dj's with D
    #       QL        QL           QL           Ql Dj


    # then the Dk's are collected into D

    alpha, s = vals(state.Q[k])
    ## we express Q_k Q_k1 = ID_k Q_k1 = I Q'_k1 D_k'
    ## XXX This should be a "passthrough"

    state.Q[k] = Rotator(one(Complex{T}), zero(Complex{T}), idx(state.Q[k]))
    state.D[k] *= alpha
    state.D[k+1] *= conj(alpha)

    c, s = vals(state.Q[k+1])
    state.Q[k+1] = Rotator(c, conj(alpha)*s, idx(state.Q[k+1]))
    #cascade(state.Q, state.D, alpha, k, state.ctrs.stop_index)

    # shift zero counter
    state.ctrs.zero_index = k      # points to a matrix Q[k] either RealRotator(-1, 0) or RealRotator(1, 0)
    state.ctrs.start_index = k + 1

    # reset counter
    state.ctrs.it_count = 1
end

##################################################

## Diagaonal block
getDR(state::FactorizationType{T, Val{:SingleShift}, P, Tw}, k) where {T, P, Tw} = one(Complex{T})
