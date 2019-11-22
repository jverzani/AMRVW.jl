using Test
import AMRVW
const A = AMRVW
using LinearAlgebra


# test polynomials
p4 = [24.0, -50.0, 35.0, -10.0, 1.0]
p5 = [-120.0, 274.0, -225.0, 85.0, -15.0, 1.0]
p6 = [720.0, -1764.0, 1624.0, -735.0, 175.0, -21.0, 1.0]
pc = [0.0 + 1.0im, -1.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im, 0.0 - 1.0im, 1.0 + 0.0im]


T = Float64
S = Complex{T}

csort(xs) = sort(sort(xs, by=x -> imag(x)), by=x -> real(x))

@testset "transformations" begin

    Ts = [BigFloat, Float64, Float32, Float16]


    # Real case first
    ## givensrot
    for T in Ts
        a,b = rand(T, 2)
        c,s,r = A.givensrot(a,b)
        M =  diagm(0 => ones(T, 2))
        Rx = (A.RealRotator(c,s,1)*M) * [a,b]
        @test abs(Rx[2]) <= sqrt(eps(T))
    end

    ## turnover
    for T in Ts
        M = diagm(0 => ones(T, 4))
        for i in 1:10
            U,V,W  = A.random_rotator.(A.RealRotator{T}, [2,3,2])
            U1, V1, W1 = A.turnover(U,V,W)
            @test all(isapprox.([U1,V1,W1]*M , [U,V,W]*M, atol=4eps(T)))
        end
    end


    ## fuse
    for T in Ts
        U,V  = A.random_rotator.(A.RealRotator{T}, [2,2])
        M = diagm(0 => ones(T, 4))
        UV = A.fuse(U,V)
        @test all(isapprox.(UV*M, [U,V]*M, atol=2eps(T)))
    end


    # ComplexReal case
    ## givensrot
    for T in Ts
        a,b = rand(Complex{T}, 2)
        c,s,r = A.givensrot(a,b)
        M =  diagm(0 => ones(Complex{T}, 2))
        Rx = (A.ComplexRealRotator(c,s,1)*M) * [a,b]
        @test abs(Rx[2]) <= sqrt(eps(T))
    end

    ## turnover
    ### a special case that had caused issues
    T = Float64
    S = Complex{T}
    U = A.ComplexRealRotator{T}(-0.9990813749617048 + 0.042846941864777804im,
                                   -0.0007387675316362891, 2)
    V = A.ComplexRealRotator{T}(0.20398474677465908 - 0.08295034712193844im,
                                         -0.9754534653153004, 1)
    W = A.ComplexRealRotator{T}(0.03541646337113308 + 0.28843087042807297im,
                                   -0.9568454980331911, 2)

    M = diagm(0 => ones(Complex{T}, 4))
    U1, V1, W1 = A.turnover(U,V,W)
    @test all(isapprox.([U1,V1,W1]*M , [U,V,W]*M))

    for T in Ts
        M = diagm(0 => ones(Complex{T}, 4))
        for i in 1:10
            U,V,W  = A.random_rotator.(A.ComplexRealRotator{T}, [2,3,2])
            U1, V1, W1 = A.turnover(U,V,W)
            @test norm([U,V,W]*M - [U1,V1,W1]*M) <= sqrt(eps(T))
        end
    end

    ## when terms are 0
    for T in Ts
        M = diagm(0 => ones(Complex{T}, 3))

        V,W = A.random_rotator.(A.ComplexRealRotator{T}, [2,1])
        a,b =  sincos(T(pi/4))
        U = A.Rotator(complex(a,b), zero(T), 1)

        A.turnover(U,V,W)
        U1, V1, W1 = A.turnover(U,V,W)
        @test all(isapprox.([U1,V1,W1]*M , [U,V,W]*M, atol=4eps(T)))

        U,W  = A.random_rotator.(A.ComplexRealRotator{T}, [1,1])
        a,b =  sincos(T(pi/4))
        V = A.Rotator(complex(a,b), zero(T), 2)

        A.turnover(U,V,W)
        U1, V1, W1 = A.turnover(U,V,W)
        @test all(isapprox.([U1,V1,W1]*M , [U,V,W]*M, atol=4eps(T)))

        U,V = A.random_rotator.(A.ComplexRealRotator{T}, [1,2])
        a,b =  sincos(T(pi/4))
        W = A.Rotator(complex(a,b), zero(T), 1)

        A.turnover(U,V,W)
        U1, V1, W1 = A.turnover(U,V,W)
        @test all(isapprox.([U1,V1,W1]*M , [U,V,W]*M, atol=4eps(T)))

    end

    ## fuse
    for T in Ts
        U,V  = A.random_rotator.(A.ComplexRealRotator{T}, [2,2])
        M = diagm(0 => ones(Complex{T}, 4))
        UV,Di = A.fuse(U,V)
        alpha = Di.c
        D = diagm(0 => [1, alpha, conj(alpha), 1])
        @test norm(UV*D*M - [U,V]*M) <= sqrt(eps(T))
    end



end


@testset "factorization" begin
    # RDS case first
    ps = [p4, p5, p6]
    Ts = [Float64]
    for p in ps
        for T in Ts
            d = length(p)-1
            state = A.amrvw(p)

            # Cxs = y
            C = state.RF.Ct'
            xs = A.basic_decompose(p)
            M = diagm(0 => ones(T, d+1))
            y = (C*M)*xs
            @test abs(y[1] - norm(xs)) <= sqrt(eps(T)) # check that y[1] is right
            @test norm(y[2:end]) <= sqrt(eps(T))

            # Matrix gives correct eigenvalues
            F = A.Matrix(state)
            rts = sort(eigvals(F))
            @test norm(sort(eigvals(F))  .- [1:d...] ) <= sqrt(eps(T))

            # state.Ct * state.B is I if we skip some
            M = diagm(0 => ones(T, d+1))
            @test norm(state.RF.Ct.x[2:end] * (state.RF.B.x[1:end-1] * M) - I) <= sqrt(eps(T))

            # state.Ct[1] * state.B[end] is twisted
            @test norm((state.RF.Ct.x[1] * (state.RF.B.x[end]*M))[end-1:end, end-1:end] - [0 -1; 1 0]) <= sqrt(eps(T))



            M = diagm(0=> ones(T,d+1))
            state = A.amrvw(p)
            e_1 = vcat(1, zeros(d))
            e_n = vcat(zeros(d-1),1,0)

            xs = A.basic_decompose(p)
            # check that R = Z + x*en'
            Z = state.RF.Ct * (state.RF.B * M) # tested above
            X = state.QF.Q  * (Z + xs*e_n')
            @test norm(sort(eigvals(X))  .- [0:d...] ) <= sqrt(eps(T))
        end
    end



    # CSS case
    ps = [pc]
    Ts = [Float64]
    for p in ps
        for T in Ts
            S = Complex{T}
            d = length(p)-1
            state = A.amrvw(p)

            # Cxs = y
            C = state.RF.Ct'
            xs = A.basic_decompose(p)
            M = diagm(0 => ones(S, d+1))
            y = (C*M)*xs
            @test abs(y[1]) ≈ norm(xs)  # check that y[1] is right
            @test norm(y[2:end]) <= sqrt(eps(T))

            # Matrix gives correct eigenvalues
            F = Matrix(state)
            rts = csort(eigvals(F))
            @test norm(csort(round2(eigvals(F)))  .- S[-1,-im,im,im,1] ) <= sqrt(eps(T))

        end
    end


end


@testset "bulge step" begin

    # RDS First
    ps = [p4, p5, p6]
    Ts = [Float64]

    # check that bulge step works
    for p in ps
        for T in Ts
            state = A.amrvw(p)
            A.bulge_step(state)
            d = length(p)-1
            F = A.Matrix(state)
            rts = sort(eigvals(F))
            @test norm(sort(eigvals(F))  .- [1:d...] ) <= sqrt(eps(T))

            allocs = @allocated A.bulge_step(state)
            @test allocs == 0
        end
    end
end


@testset "deflation" begin

end


@testset "diagonal block" begin
    # RDS First
    ps = [p4, p5, p6]
    Ts = [Float64]


    for p in ps
        for T in Ts
            state = A.amrvw(p)
            A.bulge_step(state)
            state.ctrs.it_count += 1
            A.bulge_step(state)
            A.create_bulge(state)

            k = length(p) - 1
            @test Matrix(state.RF)[k,k] ≈ state.RF[k,k]
            @test Matrix(state.RF)[k-1,k] ≈ state.RF[k-1,k]
            @test Matrix(state.RF)[k-1,k-1] ≈ state.RF[k-1,k-1]

            A.diagonal_block(state, k)
            @test all(Matrix(state)[k-1:k, k-1:k] .≈ state.A)
        end
    end

    # Css code

end
