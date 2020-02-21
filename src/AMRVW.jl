module AMRVW


using LinearAlgebra
import LinearAlgebra: eigvals, diagm
import Base: adjoint, *

# package code goes here

include("rotators.jl")
include("transformations.jl")
include("chains.jl")
include("diagonal.jl")

include("rfactorization.jl")
include("qfactorization.jl")
include("qrfactorization.jl")

include("diagonal-block.jl")
include("bulge-step.jl")
include("create-bulge.jl")

include("RDS.jl")
include("CSS.jl")
include("pencil.jl")
include("twisted.jl")

include("AMRVW_algorithm.jl")
include("companion.jl")


include("utils.jl")
include("diagnostics.jl")


end # module
