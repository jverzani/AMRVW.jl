## Diagonostic code
##

## debugging tools, basically
round2(x) = round.(x, digits=2)
round5(x) = round.(x, digits=5)
export round2, round5
printtp(x) = println(sprint(io -> show(io, "text/plain", x)))


### Random rotators are useful for testing
random_rotator(R::AbstractRotator, i) = random_rotator(typeof(R), i)
function random_rotator(t::Type{T}, i) where {T <: Real}
    a,b = rand(T, 2)
    c,s,_ = givensrot(a,b)
    Rotator(c,s,i)
end

function random_rotator(s::Type{S}, i) where {S <: Complex}
    a,b = rand(S, 2)
    c,s,_ = givensrot(a,b)
    ns = norm(s)
    alpha = conj(s)/ns
    Rotator(c*alpha,ns,i)
end



# simple graphic to show march of algorithm
function show_status(state, ctr)
    qs = [norm(u.s) for u in state.QF.Q.x[ctr.start_index:ctr.stop_index]]

    minq = length(qs) > 0 ?  minimum(qs) : 0.0

    x = fill(".", length(state.QF) + 2)
    x[ctr.zero_index+1] = "α"
    x[ctr.start_index+1] = "δ"
    x[ctr.stop_index+2] = "Δ"
    println(join(x, ""), " ($minq)")

end
