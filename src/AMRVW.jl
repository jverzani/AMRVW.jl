module AMRVW


using LinearAlgebra
import LinearAlgebra: eigvals
import Base: adjoint, *

# package code goes here
include("utils.jl")
include("rotators.jl")
include("transformations.jl")
include("chains.jl")
include("diagonal.jl")

include("rfactorization.jl")
include("qfactorization.jl")
include("types.jl")
include("diagonal-block.jl")
include("bulge.jl")
include("create_bulge.jl")
include("RDS.jl")
include("CSS.jl")
include("pencil.jl")
include("twisted.jl")
include("AMRVW_algorithm.jl")
include("companion.jl")

include("diagnostics.jl")


end # module
