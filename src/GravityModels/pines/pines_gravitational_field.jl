# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Function to compute the gravitational field derivative.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
"""
    gravitational_field!(model::AbstractGravityModel{T}, r::AbstractVector, formulation::Val{:Pines}, time_since_JD2000_TT::Number = 0.0; kwargs...) where T<:Number -> NTuple{3, T}

Compute the gravitational field derivative [SI] with respect to the spherical coordinates
(`∂U/∂r`, `∂U/∂ϕ`, `∂U/∂λ`) using the `model` in the position `r` [m], represented in ITRF,
at instant `time`. If the latter argument is omitted, the J2000.0 epoch is used. This will
be referred to as the classical model.

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
function gravitational_field_datetime!(
    model::AbstractGravityModel,
    r::AbstractVector{T},
    formulation::Val{:Pines},
    time::DateTime = DateTime("2000-01-01");
    max_degree::Number = -1,
    max_order::Number = -1,
    U::Union{Nothing, Number} = nothing,
    F::Union{Nothing, AbstractVector} = nothing,
    Uₜ::Union{Nothing, Number} = nothing
) where T

    JD = time |> datetime2julian
    time_since_JD2000 = JD - JD_J2000

    return gravitational_field!(
        model,
        r,
        formulation,
        time_since_JD2000;
        max_degree = max_degree,
        max_order = max_order,
        U = U,
        F = F,
        Uₜ= Uₜ)

end


function gravitational_field!(
    model::AbstractGravityModel,
    r::AbstractVector{T},
    formulation::Val{:Pines},
    time_since_JD2000::Number = 0.0;
    max_degree::Number = -1,
    max_order::Number = -1,
    U::Union{Nothing, Number} = nothing,
    F::Union{Nothing, AbstractVector} = nothing,
    Uₜ::Union{Nothing, Number} = nothing
) where T<:Number


    # Unpack gravity model data
    # ======================================================================================
    μ  = gravity_constant(model)
    R₀ = radius(model)
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

    r_gc = norm(r)
    s, t, u = r ./ r_gc

    # Auxiliary variables.
    ratio = R₀ / r_gc

    aux1, aux2, aux3, aux4 = auxilary_variable(n_max)

    # Fill Anm Matrix
    Anm = @MMatrix zeros(n_max + 3, n_max + 3)
    Anm[1, 1] = 1.0
    Anm[2, 2] = 1.0
    Anm[2, 1] = u

    @inbounds for n ∈ 2:n_max + 2
        # Fill the Diagonal
        Anm[n+1, n+1] = aux1[n] * Anm[n, n]
        # Fill the First columns
        Anm[n+1, 1] = aux2[n] * u * Anm[n,1] - aux3[n] * Anm[n-1, 1]
        # Fill the Subdiagonal
        Anm[n+1, n] = u * Anm[n+1, n+1]
    end

    @inbounds for n ∈ 3:n_max + 2
        for m ∈ 1:n-2
            Anm[n+1, m+1] = aux4[n,m] * Anm[n,m] + u * Anm[n, m+1]
        end
    end

    #println(Anm)

    # Fill R, I, and P Vectors
    Rm = zeros(n_max)
    Im = zeros(n_max)
    Pn = zeros(n_max + 1)

    Rm[1] = 1.0
    Pn[1] = μ/r_gc
    Pn[2] = ratio * Pn[1]
    for n ∈ 2:n_max
        Rm[n] = s * Rm[n-1] - t * Im[n-1]
        Im[n] = s * Im[n-1] + t * Rm[n-1]
        Pn[n+1] = ratio * Pn[n]
    end

    Cnm = zeros(n_max, n_max + 1)
    Snm = zeros(n_max, n_max + 1)

    for i ∈ 1:n_max
        for j ∈ 1:min(i, m_max)
            Cnm[i,j], Snm[i,j] = coefficients(model, i, j, time_since_JD2000) ./ norm_fact(i, j)
        end
    end

    # Fill D, E, and F Matrices
    Dnm = zeros(n_max, n_max+1)
    Enm = zeros(n_max, n_max+1)
    Fnm = zeros(n_max, n_max+1)

    skip_EF = F === nothing && Uₜ === nothing

    for m ∈ 2:m_max
        for n ∈ m:n_max
            Dnm[n, m] = Cnm[n, m]*Rm[m - 1] + Snm[n, m]*Im[m-1]
            if !skip_EF
                Enm[n, m] = Cnm[n, m]*Rm[m-1] + Snm[n, m]*Im[m-1]
                Fnm[n, m] = Snm[n, m]*Rm[m-1] - Cnm[n, m]*Im[m-1]
            end
        end
    end
    for n ∈ 1:n_max
        Dnm[n, 1] = Cnm[n, 1]*Rm[1]
    end

    println(Dnm)

    ###############################################################################
    #* Perturbing potential
    ###############################################################################
    if U !== nothing
        U = 0.0
        for m ∈ 1:m_max
            for n ∈ m:n_max
                U += Pn[n] * Anm[n, m] * Dnm[n, m]
            end
        end

        for n ∈ 1:n_max
            U += Pn[n] * Anm[n, 1] * Dnm[n, 1]
        end

        U *= -1.0
    end

    ###############################################################################
    #* Perturbing acceleration
    ###############################################################################
    if F !== nothing
        a = zeros(4)
        for m ∈ 1:m_max
            for n ∈ m:n_max
                a += [Pn[n+1] * Anm[n,m] * m * Enm[n,m];
                      Pn[n+1] * Anm[n,m] * m * Fnm[n,m];
                      Pn[n+1] * Anm[n,m+1]   * Dnm[n,m];
                     -Pn[n+1] * Anm[n+1,m+1] * Dnm[n,m]]
            end
        end

        for n ∈ 1:n_max
            a += [0.0;
                  0.0;
                  Pn[n+1] * Anm[n, 2]   * Dnm[n, 1];
                 -Pn[n+1] * Anm[n+1, 2] * Dnm[n, 1]]
        end

        F .= (@view(a[1:3]) + [s; t; u] * a[4]) / R₀
    end

    ###############################################################################
    #* Time derivative of potential in body-fixed frame
    ###############################################################################
    if Uₜ !== nothing
        Uₜ = 0.0
        err_nd = IAU06_err(0.0, MDJ_TT) / TU

        Gnm = zeros(n_max, n_max)
        for m ∈ 1:m_max
            for n ∈ m:n_max
                Gnm[n,m] = (m * (t * Enm[n,m] - s * Fnm[n,m])) * err_nd
                Uₜ += Pn[n] * Anm[n,m] * Gnm[n,m]
            end
        end
    end

    return nothing

end

function norm_fact(l::Int, m::Int)

    kron = (m == 0) ? 1.0 : 0.0

    numer = gamma(l + m + 1.0)
    denom = (2.0 - kron) * (2.0 * l + 1.0) * gamma(l - m + 1.0)

    return √(numer/denom)

end


function auxilary_variable(n_max::Int)
    aux1 = zeros(n_max)
    aux2 = zeros(n_max)
    aux3 = zeros(n_max)
    aux4 = zeros(n_max, n_max + 1)
    for n ∈ 1:n_max
        aux1[n] = 2.0*n + 1.0
        aux2[n] = aux1[n] / (n + 1.0)
        aux3[n] = n / (n + 1)
        for m ∈ 1:n
            aux4[n, m] = n + m
        end
    end

    return aux1, aux2, aux3, aux4
end

function IAU06_err(MJDA_TT, MJDB_TT)

    c0 =   4612.156534
    c1 =   2.7831634
    c2 =  -0.00000132
    c3 =  -0.000119824
    c4 =  -0.000000184

    arcsec2rad = π / (180.0 * 3600.0)
    sec2cy = 3155760000.0

    JC_TT = ((MJDA_TT - MJD_J2000) + MJDB_TT) / 36524.25

    δ̇  = @evalpoly(JC_TT, c0, c1, c2, c3, c4) * arcsec2rad / sec2cy

    return ωE + δ̇

end