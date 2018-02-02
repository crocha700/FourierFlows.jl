__precompile__()

module FourierFlows

export AbstractGrid,
       AbstractParams,
       AbstractVars,
       AbstractEquation,
       AbstractTimeStepper,
       AbstractProblem
       
# Abstract supertypes
abstract type AbstractGrid end
abstract type AbstractParams end
abstract type AbstractVars end
abstract type AbstractTimeStepper end
abstract type AbstractProblem end
abstract type AbstractEquation end
abstract type AbstractState end

"""
    cxeltype(a)

Returns Complex{eltype(a)} if eltype(a) <: Real; eltype(a) otherwise.
"""
cxeltype(a) = eltype(a) <: Real ? Complex{eltype(a)} : eltype(a)

# Include core functionality
include("probtypes.jl")
include("domains.jl")
include("diagnostics.jl")
include("output.jl")
include("utils.jl")
include("timesteppers.jl")

# Include physics modules
include("physics/twodturb.jl")
include("physics/barotropicqg.jl")
include("physics/verticallycosineboussinesq.jl")
include("physics/traceradvdiff.jl")

#include("physics/verticallyfourierboussinesq.jl")
#include("physics/tracerpatcheqn.jl")

end # module
