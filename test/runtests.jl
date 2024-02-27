using Test

using SatelliteToolboxGravityModels

using Dates
using DelimitedFiles
using LinearAlgebra
using SatelliteToolboxTransformations
using Scratch

# We must clear the scratch space before the tests to avoid errors.
clear_scratchspaces!(SatelliteToolboxGravityModels)

model = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))
r_itrf = geodetic_to_ecef(0.0, 0.0, 0)

r_itrf = [8478312.10019284, 7871539.26514707, 2304751.40932291]

U = [0.0]
F = zeros(3)
Uₜ = [0.0]

for N in 1:5
    println(N)
    g_norm_pines = GravityModels.gravitational_field!(model, r_itrf, Val(:Pines); max_order=N, max_degree=N, U=U, F=F, Uₜ=Uₜ)
    println(F)
end
g_norm_class = GravityModels.gravitational_acceleration(model, r_itrf, Val(:Classical); max_order=1, max_degree=1)


@testset "ICGEM" verbose = true begin
    include("./icgem.jl")
end

@testset "GravityModels API" verbose = true begin
    include("./gravity_models.jl")
end
