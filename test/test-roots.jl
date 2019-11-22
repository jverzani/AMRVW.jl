using AMRVW
using Test

const  A = AMRVW

T = Float64
S = Complex{T}

## Test roots function
p4 = Float64[24, -50, 35, -10, 1] #a0 a1 a2 a4 a4=1
p5 = Float64[-120.0, 274.0, -225.0, 85.0, -15.0, 1.0] #a0,a1,a2,a4,a4,a5=1
p6 = Float64[720, -1764, 1624, -735, 175, -21, 1]
pc = Complex{Float64}[40.0 + 10.0im, -76.0 + 60.0im, 10.0 - 95.0im, 25.0 + 40.0im, -10.0 - 5.0im, 1.0 + 0.0im]
pc_rts = S[im, 1+im,2+im, 3.0+im, 4+im]


csort(xs) = sort(sort(xs, by=x -> imag(x)), by=x -> real(x))


@testset "simple cases" begin

    ## constants
    T = Float64; S = Complex{T}
    ps = [1.0]
    rts = T[]
    @test A.roots(ps) == rts
    @test csort(A.roots(vcat(zeros(2), ps))) == vcat(rts, zeros(T,2))
    @test A.roots(vcat(ps, zeros(2))) == rts
    @test A.roots(vcat(zeros(2), ps, zeros(2))) == vcat(rts, zeros(T,2))

    ps = [1.0, 1.0]
    rts = [-1.0]
    @test A.roots(ps) == rts
    @test csort(A.roots(vcat(zeros(2), ps))) == vcat(rts, zeros(T,2))
    @test A.roots(vcat(ps, zeros(2))) == rts
    @test A.roots(vcat(zeros(2), ps, zeros(2))) == vcat(rts, zeros(T,2))

    ps = [1.0, 2, 1]
    rts = [-1.0, -1.0]
    @test A.roots(ps) == rts
    @test csort(A.roots(vcat(zeros(2), ps))) == vcat(rts, zeros(T,2))
    @test A.roots(vcat(ps, zeros(2))) == rts
    @test A.roots(vcat(zeros(2), ps, zeros(2))) == vcat(rts, zeros(T,2))

    ###

    ps = S[1.0]
    rts = S[]
    @test A.roots(ps) == rts
    @test csort(A.roots(vcat(zeros(2), ps))) == vcat(rts, zeros(S,2))
    @test A.roots(vcat(ps, zeros(2))) == rts
    @test A.roots(vcat(zeros(2), ps, zeros(2))) == vcat(rts, zeros(S,2))



    ps = S[1.0, 1.0]
    rts = S[-1.0]
    @test A.roots(ps) == rts
    @test csort(A.roots(vcat(zeros(2), ps))) == vcat(rts, zeros(S,2))
    @test A.roots(vcat(ps, zeros(2))) == rts
    @test A.roots(vcat(zeros(2), ps, zeros(2))) == vcat(rts, zeros(S,2))


    ps = S[1.0, 2, 1]
    rts = S[-1.0, -1.0]
    @test A.roots(ps) == rts
    @test csort(A.roots(vcat(zeros(2), ps))) == vcat(rts, zeros(S,2))
    @test A.roots(vcat(ps, zeros(2))) == rts
    @test A.roots(vcat(zeros(2), ps, zeros(2))) == vcat(rts, zeros(S,2))


    ## Pencil
    ps = S[1.0, 1]
    rts = S[-1.0]

    vs, ws = A.basic_pencil(ps)
    @test A.roots(vs, ws) == rts
    vs, ws = A.basic_pencil(vcat(zeros(2), ps))
    @test csort(A.roots(vs, ws)) == vcat(rts, zeros(S,2))
    vs, ws = A.basic_pencil(vcat(ps, zeros(2)))
    @test A.roots(vs, ws) == rts
    vs, ws = A.basic_pencil(vcat(zeros(2), ps, zeros(2)))
    @test csort(A.roots(vs, ws)) == vcat(rts, zeros(T,2))



    ps = S[1.0, 2, 1]
    rts = S[-1.0, -1.0]

    vs, ws = A.basic_pencil(ps)
    @test A.roots(vs, ws) == rts
    vs, ws = A.basic_pencil(vcat(zeros(2), ps))
    @test csort(A.roots(vs, ws)) == vcat(rts, zeros(S,2))
    vs, ws = A.basic_pencil(vcat(ps, zeros(2)))
    @test A.roots(vs, ws) == rts
    vs, ws = A.basic_pencil(vcat(zeros(2), ps, zeros(2)))
    @test csort(A.roots(vs, ws)) == vcat(rts, zeros(T,2))

end


@testset "RDS" begin

    @test all(csort(A.roots(p4)) .≈ 1.0:4)
    @test all(csort(A.roots(p5)) .≈ 1.0:5)
    @test all(csort(A.roots(p6)) .≈ 1.0:6)

    rs = rand(T, 500)
    rts = A.roots(rs)
    @test all(!iszero(rts))


end

@testset "CSS" begin

    @test all(isapprox.(csort(A.roots(pc)), pc_rts, atol=sqrt(eps(T))))

    rs = rand(S, 500)
    rts = A.roots(rs)
    @test all(!iszero(rts))


end

@testset "RDS, QZ" begin

    vs, ws = A.basic_pencil(p4)
    @test all(csort(A.roots(vs, ws)) .≈ 1.0:4)

    vs, ws = A.basic_pencil(p5)
    @test all(csort(A.roots(vs, ws)) .≈ 1.0:5)

    vs, ws = A.basic_pencil(p6)
    @test all(csort(A.roots(vs, ws)) .≈ 1.0:6)

    rs = rand(T, 500)
    vs, ws = A.basic_pencil(rs)
    rts = A.roots(vs, ws)
    @test all(!iszero(rts))


end

@testset "CSS, QZ" begin

    vs, ws = A.basic_pencil(pc)
    @test all(csort(round.(A.roots(vs, ws), digits=12)) .≈  pc_rts)


    rs = rand(S, 500)
    vs, ws = A.basic_pencil(rs)
    rts = A.roots(vs, ws)
    @test all(!iszero(rts))





end
