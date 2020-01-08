using Test
import AMRVW
const A = AMRVW
using LinearAlgebra

@testset "passthrough!-chains"  begin

    T = Float64
    S = Complex{T}
    Rt = A.ComplexRealRotator{T}

    # check passthorugh! with Desc; Asc; and Twisted
    N = 10
    Des = A.random_rotator.(Rt, 1:N) |> A.DescendingChain
    Asc = A.random_rotator.(Rt, N:-1:1) |> A.AscendingChain
    #Tw = A.random_rotator.(Rt, 0 .+ reverse([1,3,2,4,7,6,5 ])) |> A.TwistedChain
    Tw = A.random_rotator.(Rt, 0 .+ reverse([3,4,7,6,5 ])) |> A.TwistedChain

    M = diagm(0 => ones(S, N+1))

    ## We have Ascending, Descending, and Twisted to check
    N = 10
    Des = A.random_rotator.(Rt, 1:N) |> A.DescendingChain
    Asc = A.random_rotator.(Rt, N-1:-1:1) |> A.AscendingChain
    Tw = A.random_rotator.(Rt, 0 .+ reverse([3,4,7,6,5 ])) |> A.TwistedChain

    # <--
    M1 = Des * (Tw * M)
    A.passthrough!(Des, Tw)
    M2 = Tw * (Des * M)
    @test M1 - M2 .|> round5 .|> iszero |> all



    M1 = Asc * (Tw * M)
    A.passthrough!(Asc, Tw)
    M2 = Tw * (Asc * M)
    @test M1 - M2 .|> round5 .|> iszero |> all

    # -->
    M1 = Tw * (Des * M)
    A.passthrough!(Tw, Des)
    M2 = Des * (Tw * M)
    @test M1 - M2 .|> round5 .|> iszero |> all


    M1 = Tw * (Asc * M)
    A.passthrough!(Tw, Asc)
    M2 = Asc * (Tw * M)
    @test M1 - M2 .|> round5 .|> iszero |> all

    ## check passthorugh Des, Asc
    M1 = Des * (Asc * M)
    A.passthrough!(Des, Asc)
    M2 = Asc * (Des * M)
    @test M1 - M2 .|> round5 .|> iszero |> all

    Des = A.random_rotator.(Rt, 1:N-1) |> A.DescendingChain
    Asc = A.random_rotator.(Rt, N:-1:1) |> A.AscendingChain
    M1 = Des * (Asc * M)
    A.passthrough!(Des, Asc)
    M2 = Asc * (Des * M)
    @test M1 - M2 .|> round5 .|> iszero |> all

    Des = A.random_rotator.(Rt, 1:N) |> A.DescendingChain
    Asc = A.random_rotator.(Rt, N:-1:2) |> A.AscendingChain
    M1 = Asc * (Des * M)
    A.passthrough!(Asc, Des)
    M2 = Des * (Asc * M)
    @test M1 - M2 .|> round5 .|> iszero |> all

    Des = A.random_rotator.(Rt, 2:N) |> A.DescendingChain
    Asc = A.random_rotator.(Rt, N:-1:1) |> A.AscendingChain
    M1 = Asc * (Des * M)
    A.passthrough!(Asc, Des)
    M2 = Des * (Asc * M)
    @test M1 - M2 .|> round5 .|> iszero |> all
end

@testset "passthrough_phase!"  begin

    # test passthrough_phase!
    T = Float64; S = Complex{T}; Rt = A.ComplexRealRotator{T}
    n = 4

    ## Descending
    for i in 1:6
        Ms = A.random_rotator.(Rt,  [2,3,4])
        Des = A.random_rotator.(Rt,  [1,2,3,4,5])
        D =  A.sparse_diagonal(S,6+1)
        Di = A.DiagonalRotator(complex(sincos(rand())...), i)
        M1 = Di * (Ms * (Des *diagm(0 => D.x)))

        A.passthrough_phase!(Di, A.DescendingChain(Ms), (A.DescendingChain(Des),), D)

        M2 = Ms * (Des *diagm(0 => D.x))
        @test M1 - M2 |> round5  .|> iszero |> all
    end

    ## Ascending
    for i in 1:6
        Ms = A.random_rotator.(Rt,  [4,3,2])
        Des = A.random_rotator.(Rt,  [1,2,3,4,5])
        D =  A.sparse_diagonal(S,6+1)
        Di = A.DiagonalRotator(complex(sincos(rand())...), i)
        M1 = Di * (Ms * (Des *diagm(0 => D.x)))

        A.passthrough_phase!(Di, A.AscendingChain(Ms), (A.DescendingChain(Des),), D)

        M2 = Ms * (Des *diagm(0 => D.x))
        @test M1 - M2 |> round5  .|> iszero |> all
    end

    ## Twisted
    for i in 1:6
        Ms = A.random_rotator.(Rt,  [2,3,4,5])
        D =  A.sparse_diagonal(S,6+1)
        Di = A.DiagonalRotator(complex(sincos(rand())...), i)
        M1 = Di * (Ms * diagm(0 => D.x))

        A.passthrough_phase!(Di, A.TwistedChain(Ms), (), D)

        M2 = Ms * diagm(0 => D.x)
        @test M1 - M2 |> round5  .|> iszero |> all
    end

    ##
    for i in 1:6
        Ms = A.random_rotator.(Rt,  [5,4,3,2])
        D =  A.sparse_diagonal(S,6+1)
        Di = A.DiagonalRotator(complex(sincos(rand())...), i)
        M1 = Di * (Ms * diagm(0 => D.x))

        A.passthrough_phase!(Di, A.TwistedChain(Ms), (), D)

        M2 = Ms * diagm(0 => D.x)
        @test M1 - M2 |> round5  .|> iszero |> all
    end

    ## CMV
    for i in 1:6
        Ms = A.random_rotator.(Rt,  [2,4,3,5]) # XX
        Des = A.random_rotator.(Rt,  [1,2,3,4,5])
        D =  A.sparse_diagonal(S,6+1)
        Di = A.DiagonalRotator(complex(sincos(rand())...), i)
        M1 = Di * (Ms * (Des *diagm(0 => D.x)))

        A.passthrough_phase!(Di, A.TwistedChain(Ms), (A.DescendingChain(Des),), D)

        M2 = Ms * (Des * diagm(0 => D.x))
        @test M1 - M2 |> round5  .|> iszero |> all
    end

    ##
    for i in 1:6

        Ms = A.random_rotator.(Rt,  [3,2,4,5]) # XXX 5 is issue
        Des = A.random_rotator.(Rt,  [1,2,3,4,5])
        D =  A.sparse_diagonal(S,6+1)
        Di = A.DiagonalRotator(complex(sincos(rand())...), i)
        M1 = Di * (Ms * (Des * diagm(0 => D.x)))

        A.passthrough_phase!(Di, A.TwistedChain(Ms), (A.DescendingChain(Des),), D)

        M2 = Ms * (Des * diagm(0 => D.x))
        @test M1 - M2 |> round5  .|> iszero |> all
    end


end
