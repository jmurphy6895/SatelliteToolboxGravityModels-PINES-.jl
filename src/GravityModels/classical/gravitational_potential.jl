# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Function to compute the gravitational field potential using the Classical Formulation.
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
    formulation::Val{:Classical},
    time_since_JD2000::Number = 0.0;
    max_degree::Number = -1,
    max_order::Number = -1,
    P::Union{Nothing, AbstractMatrix} = nothing,
) where T<:Number

    # Unpack gravity model coefficients.
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

    # Obtain the required sizes for the matrices P and dP.
    #
    # Notice that, to compute the derivative if `m_max < n_max`, we need that `P` has an
    # order at least one time higher than `dP`. Otherwise, we will access regions with
    # undefined numbers.
    if n_max == m_max
        n_max_P  = m_max_P  = n_max
    else
        n_max_P  = n_max
        m_max_P  = m_max + 1
    end

    # Check if the matrices related to Legendre must be computed.
    if isnothing(P)
        P = Matrix{T}(undef, n_max_P + 1, m_max_P + 1)

    else
        # If the user passed a matrix, we must check if there are enough space to store the
        # coefficients.
        rows, cols = size(P)

        if (rows < n_max_P + 1) || (cols < m_max_P + 1)
            throw(ArgumentError("Matrix `P` must have at least $(n_max_P + 1) rows and $(m_max_P + 1) columns."))
        end
    end


    # Auxiliary variables.
      # ======================================================================================
    ρ²_gc = r[1]^2 + r[2]^2
    r²_gc = ρ²_gc  + r[3]^2
    r_gc  = √r²_gc
    ρ_gc  = √ρ²_gc
    ϕ_gc  = atan(r[3], ρ_gc)
    λ_gc  = atan(r[2], r[1])

    # Sine and cosine of the geocentric longitude.
    #
    # This values were be used in the algorithm to decrease the computational
    # burden.
    sin_λ,  cos_λ  = sincos(λ_gc)
    sin_2λ, cos_2λ = sincos(2λ_gc)

    # Consider the zero degree term.
    # ======================================================================================
    U = T(1)

    # Compute the associated Legendre functions `P_n,m[sin(ϕ_gc)]` with the
    # required normalization.
    #
    # Notice that cos(ϕ_gc - π / 2) = sin(ϕ_gc).
    legendre!(Val(norm_type), P, ϕ_gc - T(π / 2), n_max_P, m_max_P; ph_term = false)

    # Auxiliary variables.
    ratio = R₀ / r_gc
    fact  = ratio

    @inbounds for n in 2:n_max
        aux_U = T(0)

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

        for m in 0:min(n, m_max)
            # Compute recursively `sin(m*λ_gc)` and `cos(m*λ_gc)`.
            sin_mλ = 2cos_λ * sin_m_1λ - sin_m_2λ
            cos_mλ = 2cos_λ * cos_m_1λ - cos_m_2λ

            # Get the spherical harmonics coefficients
            # ==============================================================================

            clm, slm = coefficients(model, n, m, time_since_JD2000)

            aux_U += P[n+1, m+1] * (clm * cos_mλ + slm * sin_mλ)

            # Update the values for the next step.
            sin_m_2λ = sin_m_1λ
            sin_m_1λ = sin_mλ
            cos_m_2λ = cos_m_1λ
            cos_m_1λ = cos_mλ
        end

        # fact = (a / r)^(n + 1)
        fact *= ratio

        U += fact * aux_U
    end

    U *= μ / r_gc

    return U
end