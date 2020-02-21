using Test
import AMRVW
const A = AMRVW
using LinearAlgebra
using Random

# test polynomials
p4 = [24.0, -50.0, 35.0, -10.0, 1.0]
p5 = [-120.0, 274.0, -225.0, 85.0, -15.0, 1.0]
p6 = [720.0, -1764.0, 1624.0, -735.0, 175.0, -21.0, 1.0]
pc = [0.0 + 1.0im, -1.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im, 0.0 - 1.0im, 1.0 + 0.0im]

# could compare many ways
compare_eigenvalues(M1, M2) = norm(prod(eigvals(M1)) - prod(eigvals(M2))) <= sqrt(eps())


@testset "twisted-bulge-step-any-m" begin
    # bulge step should not change eigenvalues

    N = 10

    ## Real first
    T = Float64
    S = T;
    D = A.SparseDiagonal(T, 0)

    for m in 1:N-1
        @assert 0 < m < N

        ps = randperm(N)
        Ms = A.TwistedChain(A.random_rotator.(T, ps))
        Asc = A.random_rotator.(T, m:-1:1)

        RF = A.RFactorizationIdentity{T, S}()
        M = diagm(0 => ones(S, N+1))
        M1 = Ms * M
        A.bulge_step!(N,m, Ms, D, RF, Asc)
        M2 = Ms * (D * M)
        @test compare_eigenvalues(M1, M2)
    end



    ## complex
    S = Complex{T}
    D =  A.SparseDiagonal(S,N+1)

    for m in 1:N-1
        ps = randperm(N)
        Ms = A.TwistedChain(A.random_rotator.(S, ps))
        Asc = A.random_rotator.(S, m:-1:1)

        RF = A.RFactorizationIdentity{T, S}()
        M = diagm(0 => ones(S, N+1))
        M1 = Ms * M

        A.bulge_step!(N, m, Ms, D, RF, Asc)
        M2 = Ms * (D * M)
        @test compare_eigenvalues(M1, M2)
    end

end


@testset "eigenvalues twisted Q" begin

    # Basic algorithm
    ## real
    T = Float64
    S = Complex{T}
    Qs = A.random_rotator.(T,  [6, 10, 5, 20, 12, 14, 7, 3, 17, 19, 1, 13, 8, 18, 16, 4, 2, 9, 15, 11])
    QF = A.q_factorization(A.TwistedChain(Qs))
    state = A.QRFactorization(QF)
    e1 = eigvals(Matrix(state))
    @test norm(e1 - eigvals(state)) <= sqrt(eps(T))

    ## Complex
    Qs = A.random_rotator.(S,  [6, 10, 5, 20, 12, 14, 7, 3, 17, 19, 1, 13, 8, 18, 16, 4, 2, 9, 15, 11])
    QF = A.q_factorization(A.TwistedChain(Qs))
    state = A.QRFactorization(QF)
    e1 = eigvals(Matrix(state))
    es = eigvals(state)
    @test norm(e1 - es) <= sqrt(eps(T))

    # Different RF
    ## Real
    state = A.amrvw(p6)
    RF1 = state.RF
    RF2 = A.amrvw(rand(6), rand(6)).RF
    R = Matrix(RF1)[1:end-1, 1:end-1]
    RF3 = A.RFactorizationUpperTriangular(R)
    RF4 = A.RFactorizationUpperTriangular(LinearAlgebra.UpperTriangular(Matrix(R)))
    RF5 = A.RFactorizationUnitaryDiagonal(rand([-1.0, 1.0],6))
    RF6 = A.RFactorizationIdentity{T, S}()

    pv = [:right, :left, :right, :left]
    QF = A.QFactorizationTwisted(A.TwistedChain(state.QF.Q.x, pv))
    for RF in (RF1, RF2, RF3, RF4, RF5, RF6)
        state = A.QRFactorization(QF, RF)
        e1 = eigvals(Matrix(state))
        es = eigvals(state)
        @test norm(prod(e1) - prod(es)) < sqrt(eps(T))
    end


    ## Complex
    state = A.amrvw(pc)
    Q = state.QF.Q
    RF1 = state.RF
    RF2 = A.amrvw(rand(S,5), rand(S,5)).RF
    R = Matrix(RF1)#[1:end-1, 1:end-1]
    RF3 = A.RFactorizationUpperTriangular(R)
    RF4 = A.RFactorizationUpperTriangular(LinearAlgebra.UpperTriangular(Matrix(R)))
    ps = [complex(sincos(x)...) for x in rand(T,5)]
    RF5 = A.RFactorizationUnitaryDiagonal(ps)

    RF6 = A.RFactorizationIdentity{T, S}()

    pv = [:right, :left, :right]
    QF = A.QFactorizationTwisted(A.TwistedChain(state.QF.Q.x, pv))
    for RF in (RF1, RF2, RF3, RF4, RF5)
        state = A.QRFactorization(QF, RF)
        e1 = eigvals(Matrix(state))
        es = eigvals(state)
        @test norm(prod(e1) - prod(es)) < sqrt(eps(T))
    end


end
