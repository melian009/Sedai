module Visualization

using CairoMakie
using Makie
using GeoMakie
using SparseArrays

export plot_avg_total_biomass, plot_avg_species_richness, plot_combined_analysis, save_figure

function utm_to_latlon(easting, northing, zone::Int=30)
    a = 6378137.0
    f = 1/298.257223563
    k0 = 0.9996

    e = sqrt(2*f - f^2)
    e2 = e^2 / (1 - e^2)

    FE = 500000.0
    FN = (zone >= 33) ? 0.0 : 0.0

    x = easting - FE
    y = northing - FN

    M = y / k0
    mu = M / (a * (1 - e^2/4 - 3*e^4/64 - 5*e^6/256))

    e1 = (1 - sqrt(1 - e^2)) / (1 + sqrt(1 - e^2))

    phi1 = mu +
           (3*e1/2 - 27*e1^3/32) * sin(2*mu) +
           (21*e1^2/16 - 55*e1^4/32) * sin(4*mu) +
           (151*e1^3/96) * sin(6*mu) +
           (1097*e1^4/512) * sin(8*mu)

    sin_phi1 = sin(phi1)
    cos_phi1 = cos(phi1)
    tan_phi1 = tan(phi1)

    N1 = a / sqrt(1 - e^2 * sin_phi1^2)
    T1 = tan_phi1^2
    C1 = e2 * cos_phi1^2
    R1 = a * (1 - e^2) / ((1 - e^2 * sin_phi1^2)^1.5)
    D = x / (N1 * k0)

    lat = phi1 - (N1 * tan_phi1 / R1) * (
        D^2/2 -
        (5 + 3*T1 + 10*C1 - 4*C1^2 - 9*e2) * D^4/24 +
        (61 + 90*T1 + 298*C1 + 45*T1^2 - 252*e2 - 3*C1^2) * D^6/720
    )

    lon = (D -
        (1 + 2*T1 + C1) * D^3/6 +
        (5 - 2*C1 + 28*T1 - 3*C1^2 + 8*e2 + 24*T1^2) * D^5/120
    ) / cos_phi1

    lon = lon + (zone - 1) * 6 - 180 + 3

    lat_deg = lat * 180 / π
    lon_deg = lon * 180 / π

    return (lat_deg, lon_deg)
end

function plot_avg_total_biomass(sol, sites, species; figure_size::Tuple = (900, 600))
    n_sites = length(sites)
    n_species = length(species)
    n_timepoints = length(sol.t)

    total_biomass = zeros(n_timepoints, n_sites)

    for t in 1:n_timepoints
        u_matrix = reshape(sol.u[t], n_sites, n_species)
        total_biomass[t, :] .= sum(u_matrix, dims = 2)[:]
    end

    avg_biomass = vec(mean(total_biomass, dims = 2))
    std_biomass = vec(std(total_biomass, dims = 2))
    upper = avg_biomass .+ std_biomass
    lower = avg_biomass .- std_biomass
    lower = max.(lower, 0.0)

    time = sol.t

    fig = Figure(figure_size = figure_size)
    ax = Axis(fig[1, 1],
        title = "Average Total Biomass Over Time (± 1 SD)",
        xlabel = "Time",
        ylabel = "Average Total Biomass")

    band!(ax, time, lower, upper, color = (:blue, 0.2))
    lines!(ax, time, avg_biomass, color = :blue, linewidth = 2)

    return fig
end

function plot_avg_species_richness(sol, sites, species; threshold::Float64 = 0.1, figure_size::Tuple = (900, 600))
    n_sites = length(sites)
    n_species = length(species)
    n_timepoints = length(sol.t)

    richness = zeros(n_timepoints, n_sites)

    for t in 1:n_timepoints
        u_matrix = reshape(sol.u[t], n_sites, n_species)
        for site_idx in 1:n_sites
            richness[t, site_idx] = sum(u_matrix[site_idx, :] .> threshold)
        end
    end

    avg_richness = vec(mean(richness, dims = 2))
    std_richness = vec(std(richness, dims = 2))
    upper = avg_richness .+ std_richness
    lower = max.(avg_richness .- std_richness, 0.0)

    time = sol.t

    fig = Figure(figure_size = figure_size)
    ax = Axis(fig[1, 1],
        title = "Average Species Richness Over Time (± 1 SD, threshold = $threshold)",
        xlabel = "Time",
        ylabel = "Average Number of Species")

    band!(ax, time, lower, upper, color = (:green, 0.2))
    lines!(ax, time, avg_richness, color = :green, linewidth = 2)

    return fig
end

function plot_combined_analysis(sol, site_df, sites, species, distance_matrix; figure_size::Tuple = (1400, 1000))
    n_sites = length(sites)
    n_species = length(species)

    fig = Figure(figure_size = figure_size)

    ax1 = Axis(fig[1, 1], title = "Total Biomass Over Time", xlabel = "Time", ylabel = "Total Biomass")
    time = sol.t
    total_biomass = [sum(reshape(sol.u[t], n_sites, n_species)) for t in 1:length(sol.t)]
    lines!(ax1, time, total_biomass, linewidth = 2, color = :blue)

    ax2 = Axis(fig[1, 2], title = "Total Species Richness Over Time", xlabel = "Time", ylabel = "Richness")
    richness = [sum(reshape(sol.u[t], n_sites, n_species) .> 0.1) for t in 1:length(sol.t)]
    lines!(ax2, time, richness, linewidth = 2, color = :green)

    filtered_df = filter(row -> row.CODIGO in sites, site_df)

    lats = Float64[]
    lons = Float64[]
    for (ex, ey) in zip(Float64.(filtered_df.UTMX), Float64.(filtered_df.UTMY))
        lat, lon = utm_to_latlon(ex, ey, 30)
        push!(lats, lat)
        push!(lons, lon)
    end

    elevations = Float64.(filtered_df.ALTITUD)

    ax3 = GeoAxis(fig[2, 1:2];
        title = "Site Network",
        xlabel = "Longitude",
        ylabel = "Latitude",
        dest = "+proj=latlong",
        limits = (extrema(lons) .+ (-0.5, 0.5), extrema(lats) .+ (-0.5, 0.5)))

    land = GeoMakie.land()
    poly!(ax3, land; color = (:lightgray, 0.3), strokecolor = :gray, strokewidth = 0.5)
    lines!(ax3, GeoMakie.coastlines(); color = :darkgray, linewidth = 1)

    sc = scatter!(ax3, lons, lats;
        color = elevations,
        markersize = 15,
        colormap = :terrain,
        strokecolor = :black,
        strokewidth = 1)

    Colorbar(fig[2, 3], sc, label = "Elevation (m)")

    site_to_pos = Dict{String, Int}()
    for (pos, row) in enumerate(eachrow(filtered_df))
        site_to_pos[row.CODIGO] = pos
    end

    I, J, V = findnz(distance_matrix)
    valid_connections = []
    for idx in 1:length(I)
        i, j = I[idx], J[idx]
        if i > length(sites) || j > length(sites)
            continue
        end
        site_i = sites[i]
        site_j = sites[j]
        if haskey(site_to_pos, site_i) && haskey(site_to_pos, site_j)
            pos_i = site_to_pos[site_i]
            pos_j = site_to_pos[site_j]
            push!(valid_connections, (pos_i, pos_j, V[idx]))
        end
    end

    sort!(valid_connections, by = c -> c[3])
    max_shown = min(100, length(valid_connections))

    for idx in 1:max_shown
        pos_i, pos_j, dist = valid_connections[idx]
        strength = min(1.0, 5000.0 / (dist + 1))
        lines!(ax3, [lons[pos_j], lons[pos_i]], [lats[pos_j], lats[pos_i]],
            color = (:blue, 0.6 * strength),
            linewidth = 1.0 + strength * 2)
    end

    return fig
end

function save_figure(fig, filename::String; size::Tuple = (1200, 900))
    Makie.save(filename, fig, size = size)
    println("Figure saved to: $filename")
end

end # module Visualization