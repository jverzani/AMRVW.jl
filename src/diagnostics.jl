## Diagonostic code
##
import LinearAlgebra: diagm

## debugging tools, basically

XXX(st) = println("XXX: $st")
round1(x) = round.(x, digits=1)
round2(x) = round.(x, digits=2)
round3(x) = round.(x, digits=3)
round5(x) = round.(x, digits=5)
export round1, round2, round3, round5
printtp(x) = println(sprint(io -> show(io, "text/plain", x)))

## Multiplication of rotators
## This is useful for diagnostics

*(A::AbstractRotatorChain, M::Array) = A.x * M


*(D::SparseDiagonal, M::Matrix) = diagm(0=>D.x) * M
*(M::Matrix, D::SparseDiagonal)  = M * diagm(0=>D.x)

*(D::IdentityDiagonal, M::Matrix) = M
*(M::Matrix, D::IdentityDiagonal) = M


### Random rotators are useful for testing
random_rotator(R::AbstractRotator, i) = random_rotator(typeof(R), i)
function random_rotator(::Type{RealRotator{T}}, i) where {T}
    a,b = rand(T, 2)
    c,s,_ = givensrot(a,b)
    RealRotator(c,s,i)
end

function random_rotator(::Type{ComplexRealRotator{T}}, i) where {T}
    a,b = rand(Complex{T}, 2)
    c,s,_ = givensrot(a,b)
    ns = norm(s)
    alpha = conj(s)/ns
    ComplexRealRotator(c*alpha,ns,i)
end



## ### Code to create full matrix from factorization

tilde(A) = A[1:end-1, 1:end-1]

Base.Matrix(IdentityRFactorization) = I

function Base.Matrix(RF::RFactorization{T, Rt}) where {T, Rt}
    S = Rt == RealRotator{T} ? T : Complex{T}

    n = length(RF) + 1

    M = diagm(0 => ones(S, n))
    e1 = vcat(1, zeros(S, n-1))
    en1 = vcat(zeros(S,n-1), 1)
    en = vcat(zeros(S,n-2), 1, 0)

    ## compute yt
    Ct, B = RF.Ct, RF.B
    D = isa(RF.D, IdentityDiagonal) ? I : diagm(0=>RF.D.x)
    D = D * M
    rho = (en1' * (Ct * M) * e1)
    yt = -(1/rho * en1' * (Ct*(B*D)))

    # Compute R
    R =  tilde( (Ct * (B * D)) +  (Ct * (e1 * yt)) )
end


function Base.Matrix(RF::ZFactorization)
    V = Matrix(RF.V)
    W = Matrix(RF.W)
    V * inv(W)
end

function Base.Matrix(QF::QFactorization{T, Rt}) where {T, Rt}
    S = Rt == RealRotator{T} ? T : Complex{T}
    n = length(QF) + 2  ## note 2
    @show n
    M = diagm(0 => ones(S, n))
    D = isa(QF.D, IdentityDiagonal) ? I : diagm(0=>QF.D.x)
    D = D*M
    tilde( QF.Q * D )
end


function Base.Matrix(QF::QFactorizationTwisted{T, Rt}) where {T, Rt}

    S = Rt == RealRotator{T} ? T : Complex{T}
    n = length(QF) + 1
    M = diagm(0 => ones(S, n))
    D = isa(QF.D, IdentityDiagonal) ? I : diagm(0=>QF.D.x)
    D = D*M
    return QF.Q * D
    tilde( QF.Q * D)

end



# return A
function Base.Matrix(state::AbstractFactorizationState)

    Q = Matrix(state.QF)
    R = Matrix(state.RF)
    return Q * R

end

# simple graphic to show march of algorithm
function show_status(state)
    qs = [norm(u.s) for u in state.QF.Q.x[[state.ctrs.start_index:state.ctrs.stop_index...]]]

    minq = length(qs) > 0 ?  minimum(qs) : 0.0


    x = fill(".", state.N+2)
    x[state.ctrs.zero_index+1] = "α"
    x[state.ctrs.start_index+1] = "x"
    x[state.ctrs.stop_index+2] = "Δ"
    println(join(x, ""), " ($minq)")
end
