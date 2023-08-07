# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Submodule to define the gravity model API.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

module GravityModels

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

end # module GravityModels
