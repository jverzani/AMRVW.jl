using Test
using JET
using LinearAlgebra

@testset "JET" begin
    # JET loads as empty stubs on Julia versions it does not yet support (e.g. nightly)
    if JET.JET_AVAILABLE
        JET.test_package(AMRVW, ignored_modules=(AnyFrameModule(LinearAlgebra), AnyFrameModule(Base)))
    else
        @test_skip "JET not available on Julia $VERSION"
    end
end
