# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Function to compute the gravitational field derivative.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

"""
    gravitational_field_derivative(model::AbstractGravityModel{T}, r::AbstractVector, time::DateTime = DateTime("2000-01-01"); kwargs...) where T<:Number -> NTuple{3, T}

Compute the gravitational field derivative [SI] with respect to the spherical coordinates
(`∂U/∂r`, `∂U/∂ϕ`, `∂U/∂λ`) using the `model` in the position `r` [m], represented in ITRF,
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
- `dP::Union{Nothing, AbstractMatrix}`: An optional matrix that must contain at least
    `max_degree + 1 × max_degree + 1` real numbers that will be used to store the Legendre
    derivative coefficients, reducing the allocations. If it is `nothing`, the matrix will
    be created when calling the function.

# Returns

- `T`: The derivative of the gravitational field w.r.t. the radius (`∂U/∂r`).
- `T`: The derivative of the gravitational field w.r.t. the geocentric latitude (`∂U/∂ϕ`).
- `T`: The derivative of the gravitational field w.r.t. the longitude (`∂U/∂λ`).
"""
function gravitational_field_derivative(
    model::AbstractGravityModel{T},
    r::AbstractVector,
    time::DateTime = DateTime("2000-01-01");
    max_degree::Number = -1,
    max_order::Number = -1,
    P::Union{Nothing, AbstractMatrix} = nothing,
    dP::Union{Nothing, AbstractMatrix} = nothing
) where T<:Number

    # Unpack gravity model data
    # ======================================================================================

    μ  = gravity_constant(model)
    R₀ = radius(model)
    norm_type = coefficient_norm(model)
    model_max_degree = maximum_degree(model)

    # Process the inputs
    # ======================================================================================

    # Check maximum degree value.
    if (max_degree < 0) || (max_degree > model_max_degree)
        max_degree = model_max_degree
    end

    # Check maximum order value.
    if (max_order < 0) || (max_order > max_degree)
        max_order = max_degree
    end

    # Obtain the degree and order used for the computation.
    n_max = max_degree
    m_max = max_order

    # Check if the matrices related to Legendre must be computed.
    if isnothing(P)
        P = Matrix{T}(undef, n_max + 1, n_max + 1)

    else
        # If the user passed a matrix, we must check if there are enough space to store the
        # coefficients.
        rows, cols = size(P)

        if (rows < n_max + 1) || (cols < n_max + 1)
            throw(ArgumentError("Matrix `P` must have at least $(n_max + 1) rows and columns."))
        end
    end

    if isnothing(dP)
        dP = Matrix{T}(undef, n_max + 1, n_max + 1)

    else
        # If the user passed a matrix, we must check if there are enough space to store the
        # coefficients.
        rows, cols = size(dP)

        if (rows < n_max + 1) || (cols < n_max + 1)
            throw(ArgumentError("Matrix `P` must have at least $(n_max + 1) rows and columns."))
        end
    end

    # Geocentric latitude and longitude
    # ======================================================================================

    ρ²_gc = r[1]^2 + r[2]^2
    r²_gc = ρ²_gc  + r[3]^2
    r_gc  = √r²_gc
    ρ_gc  = √ρ²_gc
    ϕ_gc  = atan(r[3], ρ_gc)
    λ_gc  = atan(r[2], r[1])

    # Auxiliary variables
    # ======================================================================================

    # Sine and cosine of the geocentric longitude.
    #
    # This values were be used in the algorithm to decrease the computational burden.

    sin_λ,  cos_λ  = sincos(λ_gc)
    sin_2λ, cos_2λ = sincos(2λ_gc)

    # First derivative of the non-spherical portion of the gravitational field
    # ======================================================================================

    ∂U_∂r = T(1)  # ........................................... Derivative w.r.t. the radius
    ∂U_∂ϕ = T(0)  # .............................. Derivative w.r.t. the geocentric latitude
    ∂U_∂λ = T(0)  # ............................. Derivative w.r.t. the geocentric longitude

    # Compute the associated Legendre functions `P_n,m[sin(ϕ_gc)]` with the required
    # normalization and its time derivative.
    #
    # Notice that cos(ϕ_gc - π / 2) = sin(ϕ_gc).
    legendre!(Val(norm_type), P, ϕ_gc - T(π / 2), n_max, m_max; ph_term = false)
    dlegendre!(Val(norm_type), dP, ϕ_gc - T(π / 2), P, n_max, m_max; ph_term = false)

    # Auxiliary variables.
    ratio = R₀ / r_gc
    fact  = ratio

    # Compute the derivatives.
    @inbounds for n in 2:n_max
        aux_∂U_∂r = T(0)
        aux_∂U_∂ϕ = T(0)
        aux_∂U_∂λ = T(0)

        # Sine and cosine with m = 1
        # ==================================================================================
        #
        # This values will be used to update recursively `sin(m * λ_gc)` and
        # `cos(m * λ_gc)`, reducing the computational burden.
        #
        # TODO: Cache the computation.
        # We tried to compute those values only once using an external vector to store the
        # values. However, it leads to a worst performance. This behavior need further
        # investigation.
        sin_mλ   = T(0)      # sin( 0 * λ_gc)
        sin_m_1λ = -sin_λ    # sin(-1 * λ_gc)
        sin_m_2λ = -sin_2λ   # sin(-2 * λ_gc)
        cos_mλ   = T(1)      # cos( 0 * λ_gc)
        cos_m_1λ = +cos_λ    # cos(-1 * λ_gc)
        cos_m_2λ = +cos_2λ   # cos(-2 * λ_gc)

        # Compute the contributions when `m ∈ [1, min(n, m_max)]`
        # ==================================================================================

        for m in 0:min(n, m_max)
            # Compute recursively `sin(m * λ_gc)` and `cos(m * λ_gc)`.
            sin_mλ = 2cos_λ * sin_m_1λ - sin_m_2λ
            cos_mλ = 2cos_λ * cos_m_1λ - cos_m_2λ

            # Get the spherical harmonics coefficients
            # ==============================================================================

            clm, slm = coefficients(model, n, m, time)

            CcSs_nm = clm * cos_mλ + slm * sin_mλ
            ScCs_nm = slm * cos_mλ - clm * sin_mλ

            # Compute the contributions for `m`
            # ==============================================================================

            aux_∂U_∂r +=     P[n+1, m+1] * CcSs_nm
            aux_∂U_∂ϕ +=    dP[n+1, m+1] * CcSs_nm
            aux_∂U_∂λ += m * P[n+1, m+1] * ScCs_nm

            # Update the values for the next step
            # ==============================================================================

            sin_m_2λ = sin_m_1λ
            sin_m_1λ = sin_mλ
            cos_m_2λ = cos_m_1λ
            cos_m_1λ = cos_mλ
        end

        # fact = (a / r)^(n + 1)
        fact *= ratio

        # aux_<> *= (a / r)^(n + 1)
        aux_∂U_∂r *= fact
        aux_∂U_∂ϕ *= fact
        aux_∂U_∂λ *= fact

        ∂U_∂r += (n + 1) * aux_∂U_∂r
        ∂U_∂ϕ += aux_∂U_∂ϕ
        ∂U_∂λ += aux_∂U_∂λ
    end

    ∂U_∂r *= -μ / r²_gc
    ∂U_∂ϕ *= +μ / r_gc
    ∂U_∂λ *= +μ / r_gc

    return ∂U_∂r, ∂U_∂ϕ, ∂U_∂λ
end
