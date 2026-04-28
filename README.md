# Sedai
Fish Metacommunity simulation on a 3D dendritic network to show the exploitation-to-recovery dynamics

This readme provides an introductory description of the spatially-explicit modeling framework used to simulate fish metacommunity dynamics in an hypothetical River basin. The model is implemented as a system of Ordinary Differential Equations (ODEs) that integrate local population dynamics, environmental filtering, and 3D dendritic dispersal.

## 1. General Model Structure

The rate of change for the population size $N_{i,s}$ of species $s$ at site $i$ is defined by the sum of local dynamics and net spatial flux:

$$\frac{dN_{i,s}}{dt} = \underbrace{f_{local}(N_{i,s}, \mathbf{N}_i, E_i)}_{\text{Local Dynamics}} + \underbrace{\sum_{j \in \text{neighbors}} (m_{ji} N_{j,s} - m_{ij} N_{i,s})}_{\text{Spatial Flux}}$$

where:
- $N_{i,s}$ is the population density of species $s$ at site $i$.
- $\mathbf{N}_i$ is the vector of all species populations at site $i$.
- $E_i$ represents the local environmental conditions.
- $m_{ij}$ is the dispersal rate from site $i$ to site $j$.

## 2. Local Dynamics and Environmental Filtering

The local component follows a Generalized Lotka-Volterra (GLV) structure, modified by environmental suitability filters.

### 2.1. Population Growth and Biotic Interactions

The local dynamics are expressed as:

$$f_{local} = N_{i,s} \left( r_{i,s}^{eff} + \sum_{j=1}^{S} \alpha_{sj} N_{i,j} \right)$$

where:
- $r_{i,s}^{eff}$ is the effective intrinsic growth rate.
- $\alpha_{sj}$ is the interaction coefficient representing the effect of species $j$ on species $s$ (asymmetric interaction matrix).

### 2.2. Environmental Filtering

The effective growth rate $r_{i,s}^{eff}$ is determined by the base growth rate $r_{s}$ and a site-specific suitability factor $W_{i,s}$:

$$r_{i,s}^{eff} = r_{s} \cdot W_{i,s}$$

The suitability factor $W_{i,s}$ incorporates thermal tolerance and habitat quality:

$$W_{i,s} = \text{ThermalFilter}(T_i, \text{opt}_s, \sigma_s) \cdot h_i$$

#### Thermal Filter
We use a Gaussian-shaped "thermal window" to model species-specific temperature tolerances:

$$\text{ThermalFilter}(T_i, \text{opt}_s, \sigma_s) = \exp\left( -\frac{(T_i - \text{opt}_s)^2}{2\sigma_s^2} \right)$$

where:
- $T_i$ is the temperature at site $i$.
- $\text{opt}_s$ is the optimal temperature for species $s$.
- $\sigma_s$ is the thermal tolerance (standard deviation) of species $s$.

#### Habitat Suitability
- $h_i$ is the habitat suitability index at site $i$, which accounts for anthropogenic stressors (e.g., conductivity), river width, and mesohabitat types.

## 3. Spatial Flux (3D Dendritic Dispersal)

The spatial component models the movement of individuals through the river network, accounting for hydrological distance, elevation changes, and physical barriers.

### 3.1. Dispersal Rate Formulation

The dispersal rate $m_{ij}$ from site $i$ to site $j$ is defined as:

$$m_{ij} = m \cdot \left( \frac{1}{d_{ij}} \right) \cdot x_{ij} \cdot p_{ij}$$

where:
- $m$ is the base dispersal intensity.
- $d_{ij}$ is the hydrological distance between site $i$ and $j$.
- $x_{ij}$ is the 3D connectivity factor (elevation cost).
- $p_{ij}$ is the dam passability factor (ranging from 0.0 to 1.0).

### 3.2. Elevation Cost (Upstream vs. Downstream)

The cost of movement depends on the relative elevation of the sites:

- **Downstream or Equal Elevation ($e_j \le e_i$):**
  $$x_{ij} = 1.0$$
- **Upstream ($e_j > e_i$):**
  $$x_{ij} = \frac{1}{1 + c(e_j - e_i)}$$

where:
- $e_i$ and $e_j$ are the elevations of the source and destination sites, respectively.
- $c$ is the upstream migration cost constant.

## 4. Data Preparation

The model relies on a comprehensive data preparation pipeline, implemented in `src/data_preparation.jl`, to integrate heterogeneous datasets into the ODE framework. The process involves:

- **Species Characteristics (`load_species_characteristics`):** Loads data from `data/ABIOTIC/basinexample.csv`. Thermal optima are calculated as the midpoint of temperature ranges (`TEMPERATURE_C`), and thermal breadth ($\sigma_s$) is approximated as $\sigma_s \approx \text{range}/6$. Max size (`MAX_SIZE_mm`) is used to scale intrinsic growth rates.
- **Site Data (`load_site_data`):** Merges connectivity data (`data/ConnectivityUTM.csv`) and environmental data (`data/ABIOTIC/Matriz_Ambiental_Data.csv`) on site codes (`CODIGO`). Missing values in dam distance columns (`Demb arr.(m)`, `Demb ab.(m)`) are handled, and data is converted to numeric formats.
- **Interaction Matrix (`load_interaction_matrix`):** Parses qualitative interaction descriptions from `data/BIOTIC/FishInteractions.csv` into a quantitative matrix using `parse_interaction_string`. The values are assigned based on the interaction type:
    - "No coexist": -1.0
    - "Displaces": -0.8
    - "Predation": -0.5
    - "Competition" or "Interfere": -0.3
    - "Affects" (general): -0.2
    - "Coexist" or "Neutral": 0.0
- **Spatial Connectivity (`build_distance_matrix`):** Constructs a sparse distance matrix from `data/Matrix_distances_1037puntos_BRUTO_FINAL.csv` by grouping sites by subcatchment (`CODIGO_S`) and sorting them by distance to the river (`Distance.(m)`). Connections are restricted to adjacent sites within the same subcatchment to reflect the dendritic river network.
- **Environmental Filtering:** Site-specific temperatures (`extract_site_temperatures`) and habitat suitability indices (`extract_habitat_suitability`) are extracted. Habitat suitability is normalized based on the Trophic State Index (`IET`), where higher IET values indicate lower suitability.
- **Intrinsic Growth Rates (`build_intrinsic_growth_rates`):** Base growth rates are calculated as proportional to $1/\text{max\_size}$, then scaled by species presence/absence at each site (using columns ending in `_DEN` from `data/BIOTIC/FishDensity_and_Juveniles_Matrix.csv`).
- **Dispersal Matrix (`precompute_dispersal_matrix`):** A sparse dispersal matrix $M$ is pre-computed, incorporating hydrological distance, elevation-dependent migration costs (using an upstream cost constant $c$), and dam passability factors (ranging from 0.0 to 1.0).

## 5. Implementation Details

- **State Vector:** The system state $u$ is a matrix of size $N_{sites} \times N_{species}$.
- **Sparse Matrix Optimization:** To handle large-scale river networks efficiently, the dispersal rates are pre-calculated into a sparse matrix $M$, where $M_{ji}$ represents the rate from $i$ to $j$.
- **Numerical Integration:** Given the potential stiffness arising from the coupling of fast dispersal and slower local dynamics, the system is solved using stiff-aware solvers (e.g., `Rodas5` or `TRBDF2`).

## 6. Model Assumptions

1. **Constant Interaction Matrix:** Biotic interaction coefficients ($\alpha_{sj}$) are assumed to be constant across the entire basin, though the realized interaction depends on local densities.
2. **Passive/Active Dispersal Balance:** Dispersal is modeled as a density-dependent flux proportional to the population at the source site.
3. **Gaussian Thermal Niche:** Species performance is assumed to decay symmetrically as temperatures deviate from the optimum.
4. **Dam Passability:** Dams act as semi-permeable filters that reduce the effective dispersal rate between connected nodes.
