# ============================================================================
# Sensitivity Analysis: Parameter Sweep
# 
# This script runs simulations with different combinations of:
#   1. UPSTREAM_COST - additional cost factor for upstream dispersal
#   2. Exploitation scenarios (per-subcatchment exploitation factors)
#   3. Passability scenarios (per-subcatchment dam passability multipliers)
#
# SETUP INSTRUCTIONS:
#   1. First activate the local environment:
#      julia> using Pkg
#      julia> Pkg.activate(".")
#      julia> Pkg.instantiate()
#
#   2. Generate the serialized data file (one-time setup):
#      julia> include("serialize_data.jl")
#
#   3. Run this sensitivity analysis:
#      julia> include("run_sensitivity.jl")
#
# Results are saved in the results/sensitivity directory.
# ============================================================================

using Pkg
Pkg.activate(".")

using DifferentialEquations
using DataFrames
using CSV
using SparseArrays
using LinearAlgebra
using JLD2

include("src/OdeModel.jl")
include("src/Visualization.jl")

using .OdeModel
using .Visualization

const DAYS_PER_YEAR = 365
const SIMULATION_YEARS = 1
const T_SPAN = (0.0, Float64(SIMULATION_YEARS * DAYS_PER_YEAR))

all_subcatchments = [1.1, 1.2, 1.3, 1.4, 2.2, 3.0, 6.0, 7.0, 9.0, 10.0, 11.3, 11.5, 11.2, 11.1, 11.4, 12.0, 13.3, 13.2, 13.1, 14.0, 15.2, 15.1, 16.0, 17.0, 18.1, 18.3, 18.2, 19.0, 20.0, 21.0, 22.1, 22.4, 22.2, 23.0, 24.1, 24.2, 25.0, 26.1, 26.2, 26.3, 26.4, 26.5, 26.6, 27.0, 28.1, 28.2, 28.3, 28.4, 29.0, 30.1, 30.2, 30.3, 30.4, 30.5, 30.6, 30.7, 31.0, 32.4, 32.3, 32.1, 32.2, 33.0, 34.1, 34.2, 34.3, 34.4, 35.0, 36.1, 36.2, 36.3, 36.4, 37.0, 38.0, 39.0, 40.1, 40.2, 41.0]

upstream_costs = [0.01, 0.05, 0.1, 0.5]

exploitation_scenarios = Dict{String, Dict{Float64, Float64}}(
    "baseline" => Dict{Float64, Float64}(a => 1.0 for a in all_subcatchments),
    "low_exploitation" => Dict{Float64, Float64}(a => 0.9 for a in all_subcatchments),
    "medium_exploitation" => Dict{Float64, Float64}(a => 0.7 for a in all_subcatchments),
    "high_exploitation" => Dict{Float64, Float64}(a => 0.5 for a in all_subcatchments),
)

passability_scenarios = Dict{String, Dict{Float64, Float64}}(
    "baseline" => Dict{Float64, Float64}(a => 1.0 for a in all_subcatchments),
    "improved_passability" => Dict{Float64, Float64}(a => 1.5 for a in all_subcatchments),
    "reduced_passability" => Dict{Float64, Float64}(a => 0.5 for a in all_subcatchments),
    "blocked" => Dict{Float64, Float64}(a => 0.1 for a in all_subcatchments),
)

function build_exploitation_vector(site_df, sites; exploitation_per_subcatchment::Dict{Float64, Float64}=Dict{Float64, Float64}())
    n_sites = length(sites)
    exploitation_vector = ones(n_sites)
    site_to_sc = Dict(row.CODIGO => row.CODIGO_S for row in eachrow(site_df))
    for (i, site) in enumerate(sites)
        if haskey(site_to_sc, site) && haskey(exploitation_per_subcatchment, site_to_sc[site])
            exploitation_vector[i] = exploitation_per_subcatchment[site_to_sc[site]]
        end
    end
    return exploitation_vector
end

function build_dam_passability_vector(site_df, sites; passability_per_subcatchment::Dict{Float64, Float64}=Dict{Float64, Float64}())
    n_sites = length(sites)
    passability_vector = ones(n_sites)
    site_to_sc = Dict(row.CODIGO => row.CODIGO_S for row in eachrow(site_df))
    for (i, site) in enumerate(sites)
        if haskey(site_to_sc, site) && haskey(passability_per_subcatchment, site_to_sc[site])
            passability_vector[i] = passability_per_subcatchment[site_to_sc[site]]
        end
    end
    return passability_vector
end

function apply_exploitation_to_growth_rates!(intrinsic_growth_rates, exploitation_vector)
    n_sites = length(exploitation_vector)
    for i in 1:n_sites
        intrinsic_growth_rates[i, :] .*= exploitation_vector[i]
    end
    return intrinsic_growth_rates
end

function run_single_simulation(data_base, upstream_cost, exploitation_dict, passability_dict;
    simulation_years=SIMULATION_YEARS)

    t_span = (0.0, Float64(simulation_years * DAYS_PER_YEAR))

    params_copy = MetacommunityParams(
        data_base.params.n_sites,
        data_base.params.n_species,
        copy(data_base.params.interaction_matrix),
        copy(data_base.params.dispersal_matrix),
        copy(data_base.params.dispersal_scaling),
        copy(data_base.params.intrinsic_growth_rates),
        copy(data_base.params.temperatures),
        copy(data_base.params.habitat_suitability),
        copy(data_base.params.thermal_optima),
        copy(data_base.params.thermal_sigmas),
        copy(data_base.params.carrying_capacity)
    )

    exploitation_vector = build_exploitation_vector(data_base.site_df, data_base.sites;
        exploitation_per_subcatchment=exploitation_dict)
    passability_vector = build_dam_passability_vector(data_base.site_df, data_base.sites;
        passability_per_subcatchment=passability_dict)

    apply_exploitation_to_growth_rates!(params_copy.intrinsic_growth_rates, exploitation_vector)

    modified_dams = copy(data_base.dams)
    for j in 1:params_copy.n_sites
        for i in 1:params_copy.n_sites
            if i != j && modified_dams[i, j] < 1.0
                modified_dams[i, j] *= passability_vector[j]
                modified_dams[i, j] = min(1.0, modified_dams[i, j])
            end
        end
    end

    new_dispersal_matrix = precompute_dispersal_matrix(
        params_copy.n_sites,
        Matrix(data_base.distance_matrix),
        data_base.elevations,
        upstream_cost,
        modified_dams,
        data_base.species
    )

    params_final = MetacommunityParams(
        params_copy.n_sites,
        params_copy.n_species,
        params_copy.interaction_matrix,
        new_dispersal_matrix,
        params_copy.dispersal_scaling,
        params_copy.intrinsic_growth_rates,
        params_copy.temperatures,
        params_copy.habitat_suitability,
        params_copy.thermal_optima,
        params_copy.thermal_sigmas,
        params_copy.carrying_capacity
    )

    u0 = data_base.u0

    prob = ODEProblem(metacommunity_ode!, u0, t_span, params_final)

    function positivity_condition(u, t, integrator)
        any(x -> x < 0, u)
    end
    function positivity_affect!(integrator)
        integrator.u .= max.(integrator.u, 0.0)
    end
    positivity_cb = DiscreteCallback(positivity_condition, positivity_affect!; save_positions=(false, true))

    sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6, saveat=0:1.0:t_span[2], callback=positivity_cb)

    return sol, exploitation_vector, passability_vector
end

println("="^60)
println("Sensitivity Analysis: Parameter Sweep")
println("="^60)

base_output_dir = "results/sensitivity"
mkpath(base_output_dir)

println("\nLoading serialized data from data/data.jld2...")
data_file = joinpath("data", "data.jld2")
if !isfile(data_file)
    error("Data file not found. Please run 'serialize_data.jl' first to generate the serialized data.")
end

loaded = load(data_file)
data_base = (
    params = loaded["params"],
    sites = loaded["sites"],
    species = loaded["species"],
    distance_matrix = loaded["distance_matrix"],
    elevations = loaded["elevations"],
    dams = loaded["dams"],
    site_df = loaded["site_df"],
    u0 = loaded["u0"]
)

println("Loaded data for $(data_base.params.n_sites) sites and $(data_base.params.n_species) species")

total_runs = length(upstream_costs) * length(exploitation_scenarios) * length(passability_scenarios)
current_run = 0

for uc in upstream_costs
    for (exp_name, exp_dict) in exploitation_scenarios
        for (pass_name, pass_dict) in passability_scenarios
            global current_run += 1
            run_label = "uc=$(uc)_exp=$(exp_name)_pass=$(pass_name)"
            println("\n[$current_run/$total_runs] Running: $run_label")

            run_dir = joinpath(base_output_dir, "upstream_$(uc)", exp_name, pass_name)
            mkpath(run_dir)

            sol, exp_vec, pass_vec = run_single_simulation(
                data_base, uc, exp_dict, pass_dict
            )

            output_jld2 = joinpath(run_dir, "simulation_output.jld2")
            jldsave(output_jld2;
                sol_t = sol.t,
                sol_u = sol.u,
                sites = data_base.sites,
                species = data_base.species,
                upstream_cost = uc,
                simulation_years = SIMULATION_YEARS,
                exploitation_scenario = exp_name,
                passability_scenario = pass_name,
                exploitation_vector = exp_vec,
                passability_vector = pass_vec
            )

            fig_biomass = plot_avg_total_biomass(sol, data_base.sites, data_base.species)
            save_figure(fig_biomass, joinpath(run_dir, "avg_total_biomass.png"))

            fig_richness = plot_avg_species_richness(sol, data_base.sites, data_base.species)
            save_figure(fig_richness, joinpath(run_dir, "avg_species_richness.png"))

            fig_combined = plot_combined_analysis(sol, data_base.site_df, data_base.sites, data_base.species, data_base.distance_matrix)
            save_figure(fig_combined, joinpath(run_dir, "combined_analysis.png"))

            println("  Saved to: $run_dir")
        end
    end
end

println("\n" * "="^60)
println("Sensitivity analysis complete. $total_runs runs saved to '$base_output_dir'")
println("="^60)