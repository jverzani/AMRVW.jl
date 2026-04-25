using Test
using JET
using LinearAlgebra

@testset "JET" begin
    JET.test_package(AMRVW, ignored_modules=(AnyFrameModule(LinearAlgebra), AnyFrameModule(Base)))
end
