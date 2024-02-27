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
MJD_J2000 = 51545.0
ωE = 7.292_115_146_706_979e-5

function gravitational_field!(
    model::AbstractGravityModel,
    r::AbstractVector,
    formulation::Val{:Pines},
    time::DateTime;
    max_degree::Number = -1,
    max_order::Number = -1,
    U::Union{Nothing, AbstractVector} = nothing,
    F::Union{Nothing, AbstractVector} = nothing,
    Uₜ::Union{Nothing, AbstractVector} = nothing)

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
    U::Union{Nothing, AbstractVector} = nothing,
    F::Union{Nothing, AbstractVector} = nothing,
    Uₜ::Union{Nothing, AbstractVector} = nothing) where T<:Number


    # Unpack gravity model data
    # ======================================================================================
    μ  = gravity_constant(model)
    R₀ = radius(model)
    model_max_degree = maximum_degree(model)

    # Process the inputs
    # ======================================================================================

    DU = R₀
    TU = √(μ/DU^3)

    # Check maximum degree value.
    if (max_degree < 0) || (max_degree > model_max_degree)
        max_degree = model_max_degree
    end

    # Check maximum order value.
    if (max_order < 0) || (max_order > max_degree)
        max_order = max_degree
    end

    # Obtain the degree and order used for the computation.
    n_max = max_degree + 1
    m_max = max_order + 1

    r_gc = norm(r)
    s, t, u = r ./ r_gc

    # Auxiliary variables.
    ratio = R₀ / r_gc

    # Fill Anm Matrix
    Anm = fully_normalized_derived_legendre(u, n_max, m_max)

    # Fill R, I, and P Vectors
    Rm = zeros(n_max+1)
    Im = zeros(n_max+1)
    ρn = zeros(n_max+2)

    Rm[1] = 1.0
    ρn[1] = μ/r_gc
    ρn[2] = ratio * ρn[1]
    for n ∈ 2:n_max+1
        Rm[n] = s * Rm[n-1] - t * Im[n-1]
        Im[n] = s * Im[n-1] + t * Rm[n-1]
        ρn[n+1] = ratio * ρn[n]
    end

    Cnm = zeros(n_max, m_max)
    Snm = zeros(n_max, m_max)

    @inbounds for i ∈ 1:n_max
        for j ∈ 1:min(i, m_max)
            Cnm[i,j], Snm[i,j] = coefficients(model, i, j, time_since_JD2000)
        end
    end

    # Fill D, E, and F Matrices
    Dnm = zero(Anm)
    Enm = zero(Anm)
    Fnm = zero(Anm)

    for n in 1:n_max
        for m in 1:min(n, m_max)
            Dnm[n, m] = Cnm[n, m]*Rm[m] + Snm[n, m]*Im[m]
        end
    end

    skip_EF = F === nothing && Uₜ === nothing
    if !skip_EF
        for n in 1:n_max
            for m in 2:min(n, m_max)
                Enm[n, m] = Cnm[n, m]*Rm[m-1] + Snm[n, m]*Im[m-1]
                Fnm[n, m] = Snm[n, m]*Rm[m-1] - Cnm[n, m]*Im[m-1]
            end
        end
    end

    ###############################################################################
    #* Perturbing potential
    ###############################################################################
    if U !== nothing
        U[1] = 0.0
        for n ∈ 2:n_max
            for m ∈ 1:min(m_max, n)
                U[1] += ρn[n] * Anm[n, m] * Dnm[n, m]
            end
        end

        U[1] *= -1.0
    end

    ###############################################################################
    #* Perturbing acceleration
    ###############################################################################
    if F !== nothing
        a = zeros(4)
        for n ∈ 2:n_max
            for m ∈ 1:min(m_max, n)
                n1q = √(((n - m) * ((m == 1) + 1.0) * (n + m - 1.0)) / 2.0)
                n2q = √(((n + m + 2.0) * (n + m + 1.0) * (2.0*n - 1.0) * ((m == 1) + 1.0)) / ((2.0*n + 1.0) * 2.0))

                println(n, ", ", m, ": ", n1q)
                println(n, ", ", m, ": ", n2q)

                a += [ρn[n+1] * Anm[n,m] * m * Enm[n,m];
                      ρn[n+1] * Anm[n,m] * m * Fnm[n,m];
                      ρn[n+1] * n1q * Anm[n,m+1]   * Dnm[n,m];
                     -ρn[n+1] * n2q * Anm[n+1,m+1] * Dnm[n,m]]
            end
        end

        F[1:3] = (@view(a[1:3]) + [s; t; u] * a[4]) / R₀
    end

    ###############################################################################
    #* Time derivative of potential in body-fixed frame
    ###############################################################################
    if Uₜ !== nothing
        Uₜ[1] = 0.0
        err_nd = IAU06_err(0.0, time_since_JD2000) / TU

        Gnm = zeros(n_max, n_max)
        for m ∈ 1:m_max
            for n ∈ 1:m
                Gnm[n,m] = (m * (t * Enm[n,m] - s * Fnm[n,m])) * err_nd
                Uₜ[1] += ρn[n] * Anm[n,m] * Gnm[n,m]
            end
        end

        Uₜ[1] *= -1
    end

    return nothing
    
end

function IAU06_err(MJDA_TT::Number, MJDB_TT::Number)

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


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Compute the derived Legendre functions with full normalization.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==========================================================================================
#
#   [1] Holmes, S. A. and W. E. Featherstone, 2002. A unified approach to the Clenshaw
#       summation and the recursive computation of very high degree and order normalised
#       associated Legendre functions. Journal of Geodesy, 76(5), pp. 279-299.
#
#       For more info.: http://mitgcm.org/~mlosch/geoidcookbook/node11.html
#
#   [2] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. Microcosm
#       Press, Hawthorn, CA, USA.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

