# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Function to compute accelerations.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==========================================================================================
#
#   [1] Barthelmes, F (2013). Definition of Functions of the Geopotential and Their
#       Calculation from Spherical Harmonic Models. Scientific Technical Report STR09/02.
#       GeoForschungsZentrum (GFZ).
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

"""
    gravitational_acceleration(model::AbstractGravityModel{T}, r::AbstractVector, time::DateTime = DateTime("2000-01-01"); kwargs...) where T<:Number -> NTuple{3, T}

Compute the gravitational acceleration [m / s²] represented in ITRF using the `model` in the
position `r` [m], also represented in ITRF, at instant `time`. If the latter argument is
omitted, the J2000.0 epoch is used.

!!! note
    Gravitational acceleration is the acceleration caused by the central body mass only,
    i.e., without considering the centrifugal potential.

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
- `dP::Union{Nothing, AbstractMatrix}`: An optional matrix that must contain at least
    `max_degree + 1 × max_degree + 1` real numbers that will be used to store the Legendre
    derivative coefficients, reducing the allocations. If it is `nothing`, the matrix will
    be created when calling the function.

# Returns

- `T`: The derivative of the gravitational field w.r.t. the radius (`∂U/∂r`).
- `T`: The derivative of the gravitational field w.r.t. the geocentric latitude (`∂U/∂ϕ`).
- `T`: The derivative of the gravitational field w.r.t. the longitude (`∂U/∂λ`).
"""
function gravitational_acceleration(
    model::AbstractGravityModel{T},
    r::AbstractVector,
    time::DateTime = DateTime("2000-01-01");
    max_degree::Number = -1,
    max_order::Number = -1,
    P::Union{Nothing, AbstractMatrix} = nothing,
    dP::Union{Nothing, AbstractMatrix} = nothing
) where T<:Number

    # Compute the partial derivatives of the gravitational field w.r.t. the spherical
    # coordinates.
    ∂U_∂r, ∂U_∂ϕ, ∂U_∂λ = gravitational_field_derivative(
        model,
        r,
        time;
        max_degree = max_degree,
        max_order = max_order,
        P = P,
        dP = dP
    )

    # Auxiliary variables.
    ρ²_gc = r[1]^2 + r[2]^2
    r²_gc = ρ²_gc  + r[3]^2
    r_gc  = √r²_gc
    ρ_gc  = √ρ²_gc
    ϕ_gc  = atan(r[3], ρ_gc)
    λ_gc  = atan(r[2], r[1])

    # Acceleration represented in the ITRF
    # ======================================================================================

    # Compute the partial derivatives in spherical coordinate systems [1, p. 22] (eq. 120):
    #
    #     ∂U        1      ∂U     1   ∂U
    #    ---- , ---------.---- , ---.----
    #     ∂r     r.cos ϕ   ∂λ     r   ∂ϕ
    #
    # Notice that the singularity is not a problem here. When computing `cos(π / 2)` a very
    # small number will be returned and `∂U / ∂λ` is 0. Hence, the 2nd component will be 0.

    a_uen = @SVector [
        ∂U_∂r,
        ∂U_∂λ / (r_gc * cos(ϕ_gc)),
        ∂U_∂ϕ / r_gc
    ]

    # The vector `a_uen` is represented in the local UEN (Up-Earth-North) reference frame.
    # Hence, we need to describe the unitary vectors of this frame in the ECEF reference
    # frame. This can be accomplished by the following rotations matrix.
    D_itrf_uen = angle_to_dcm(ϕ_gc, -λ_gc, :YZ)
    a_itrf = D_itrf_uen * a_uen

    return a_itrf
end