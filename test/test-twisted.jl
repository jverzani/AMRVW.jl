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
        A.bulge_step!(Ms, D, RF, Asc)
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

        A.bulge_step!(Ms, D, RF, Asc)
        M2 = Ms * (D * M)
        @test compare_eigenvalues(M1, M2)
    end

end


@testset "eigenvalues twisted Q" begin

    ## Real
    state = A.amrvw(p6)
    D = state.QF.D
    Q = state.QF.Q
    RF = state.RF
    N = length(Q)
    Q.x[randperm(N)] = Q.x[:]
    QF = A.QFactorizationTwisted(A.TwistedChain(Q.x), D)
    state = A.QRFactorization(QF, RF)
    es = eigvals(state)
    LinearAlgebra.sorteig!(es)
    ### maximum(norm.(es - eigvals(Q.x * Matrix(RF)))) < sqrt(eps())

    ## Complex

end
