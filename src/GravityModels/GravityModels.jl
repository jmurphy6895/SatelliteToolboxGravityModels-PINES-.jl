# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Submodule to define the gravity model API.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

module GravityModels

using SpecialFunctions
using Dates
using ReferenceFrameRotations
using SatelliteToolboxBase
using SatelliteToolboxLegendre
using StaticArrays

############################################################################################
#                                          Types
############################################################################################

include("./types.jl")

############################################################################################
#                                         Includes
############################################################################################

include("./api.jl")

include("./classical/accelerations.jl")
include("./classical/gravitational_field_derivative.jl")
include("./classical/gravitational_potential.jl")

include("./pines/accelerations.jl")
include("./pines/gravitational_field.jl")
include("./pines/graviational_potential.jl")

end # module GravityModels
