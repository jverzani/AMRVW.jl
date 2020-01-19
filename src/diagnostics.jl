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
