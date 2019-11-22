module AMRVW


using LinearAlgebra
import Base: adjoint

# package code goes here
include("rotators.jl")
include("transformations.jl")
include("utils.jl")
include("rfactorization.jl")
include("qfactorization.jl")
include("types.jl")
include("diagonal-block.jl")
include("bulge.jl")
include("create_bulge.jl")
include("RDS.jl")
include("CSS.jl")
include("pencil.jl")
include("AMRVW_algorithm.jl")

include("diagnostics.jl")


end # module
