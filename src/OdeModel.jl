module OdeModel

using SparseArrays
using LinearAlgebra

export MetacommunityParams, metacommunity_ode!, precompute_dispersal_matrix

struct MetacommunityParams{T<:Real, M<:AbstractMatrix{T}, V<:AbstractVector{T}, S<:AbstractMatrix{T}}
    n_sites::Int
    n_species::Int
    interaction_matrix::M
    dispersal_matrix::S
    dispersal_scaling::V
    intrinsic_growth_rates::M
    temperatures::V
    habitat_suitability::V
    thermal_optima::V
    thermal_sigmas::V
    carrying_capacity::V
end

function gaussian_thermal_filter(temp, opt, sigma)
    return exp(-(temp - opt)^2 / (2 * sigma^2))
end

function metacommunity_ode!(du, u, p::MetacommunityParams, t)
    U = reshape(u, p.n_sites, p.n_species)
    dU = reshape(du, p.n_sites, p.n_species)

    for s in 1:p.n_species
        opt = p.thermal_optima[s]
        sigma = p.thermal_sigmas[s]

        for i in 1:p.n_sites
            N_is = max(U[i, s], 0.0)

            total_biomass_i = 0.0
            for j in 1:p.n_species
                total_biomass_i += max(U[i, j], 0.0)
            end

            env_filter = gaussian_thermal_filter(p.temperatures[i], opt, sigma) * p.habitat_suitability[i]
            r_eff = p.intrinsic_growth_rates[i, s] * env_filter

            interaction_term = 0.0
            for j in 1:p.n_species
                interaction_term += p.interaction_matrix[s, j] * max(U[i, j], 0.0)
            end
            interaction_term /= max(p.carrying_capacity[i], 1e-6)

            logistic_term = clamp(1.0 - total_biomass_i / max(p.carrying_capacity[i], 1e-6), -1.0, 2.0)

            dU[i, s] = N_is * (r_eff * logistic_term + interaction_term)
        end
    end

    for s in 1:p.n_species
        species_pop = @view U[:, s]
        dispersal_scale = p.dispersal_scaling[s]

        immigration = p.dispersal_matrix * species_pop

        for i in 1:p.n_sites
            emigration_rate = 0.0
            col_start = p.dispersal_matrix.colptr[i]
            col_end = p.dispersal_matrix.colptr[i+1] - 1
            for idx in col_start:col_end
                emigration_rate += p.dispersal_matrix.nzval[idx]
            end

            dU[i, s] += dispersal_scale * (immigration[i] - emigration_rate * max(U[i, s], 0.0))
        end
    end

    return nothing
end

function precompute_dispersal_matrix(n_sites, distances, elevations, c, dams)
    I = Int[]
    J = Int[]
    V = Float64[]

    for j in 1:n_sites
        for i in 1:n_sites
            d_ij = distances[i, j]
            if i == j || d_ij == 0 || isinf(d_ij)
                continue
            end

            e_j = elevations[j]
            e_i = elevations[i]

            x = 1.0
            if e_i > e_j
                x = 1.0 / (1.0 + c * (e_i - e_j))
            end

            d_km = d_ij / 1000.0

            rate = x * dams[j, i] / max(d_km, 0.001)

            push!(I, i)
            push!(J, j)
            push!(V, rate)
        end
    end

    return sparse(I, J, V, n_sites, n_sites)
end

function precompute_dispersal_matrix(n_sites, distances, elevations, c, dams, species_codes)
    return precompute_dispersal_matrix(n_sites, distances, elevations, c, dams)
end

end # module OdeModel