# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Function to compute the gravitational field potential using the Pines Formulation.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

"""
    gravitational_potential(model::AbstractGravityModel{T}, r::AbstractVector, time_since_JD2000_TT::Number = 0.0; kwargs...) where T<:Number -> T

Compute the gravitational field potential [SI] using the `model` in the position `r` [m], represented in ITRF,
at instant `time`. If the latter argument is omitted, the J2000.0 epoch is used.

!!! info
    In this case, `ϕ` is the geocentric latitude and `λ` is the longitude.

# Keywords

- `max_degree::Int`: Maximum degree used in the spherical harmonics when computing the
    gravitational field derivative. If it is higher than the available number of
    coefficients in the `model`, it will be clamped. If it is lower than 0, it will be set
    to the maximum degree available. (**Default** = -1)
- `max_order::Int`: Maximum order used in the spherical harmonics when computing the
    gravitational field derivative. If it is higher than `max_degree`, it will be clamped.
    If it is lower than 0, it will be set to the same value as `max_degree`.
    (**Default** = -1)
- `P::Union{Nothing, AbstractMatrix}`: An optional matrix that must contain at least
    `max_degree + 1 × max_degree + 1` real numbers that will be used to store the Legendre
    coefficients, reducing the allocations. If it is `nothing`, the matrix will be created
    when calling the function.
"""
function gravitational_potential(
    model::AbstractGravityModel,
    r::AbstractVector{T},
    formulation::Val{:Pines},
    time::DateTime = DateTime("2000-01-01");
    max_degree::Number = -1,
    max_order::Number = -1) where T<:Number

    return gravitational_potential(
        model,
        r,
        formulation,
        (time |> datetime2julian) - time_since_JD2000;
        max_degree = max_degree,
        max_order = max_order)

end


function gravitational_potential(
    model::AbstractGravityModel,
    r::AbstractVector{T},
    formulation::Val{:Pines},
    time_since_JD2000::Number = 0.0;
    max_degree::Number = -1,
    max_order::Number = -1) where T<:Number

    U = 0.0

    gravitational_field!(
        model,
        r,
        formulation,
        time_since_JD2000;
        max_degree = max_degree,
        max_order = max_order,
        U = U,
        F = nothing,
        Uₜ = nothing)

    
    return U

end

function gravitational_potential_time_derivative(
    model::AbstractGravityModel,
    r::AbstractVector{T},
    formulation::Val{:Pines},
    time::DateTime = DateTime("2000-01-01");
    max_degree::Number = -1,
    max_order::Number = -1) where T<:Number

    return gravitational_potential_time_derivative(
        model,
        r,
        formulation,
        (time |> datetime2julian) - time_since_JD2000;
        max_degree = max_degree,
        max_order = max_order)

end

function gravitational_potential_time_derivative(
    model::AbstractGravityModel,
    r::AbstractVector{T},
    formulation::Val{:Pines},
    time_since_JD2000::Number = 0.0;
    max_degree::Number = -1,
    max_order::Number = -1) where T<:Number

    Uₜ = 0.0

    gravitational_field!(
        model,
        r,
        formulation,
        time_since_JD2000;
        max_degree = max_degree,
        max_order = max_order,
        U = nothing,
        F = nothing,
        Uₜ = Uₜ)

    
    return Uₜ

end