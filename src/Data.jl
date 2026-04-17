# copied manually from jameel-institute/daedalus
# will be moved over to package data eventually

module Data

using ..Constants
using ..DataLoader
using LinearAlgebra

export prepare_contacts, contacts3d, get_settings, total_contacts,
       prepare_community_contacts, worker_contacts, consumer_worker_contacts,
       prepare_demog, initial_state

"""
    worker_contacts(workers; scaled=true)::Vector{Float64}

Get per-capita social contacts within each economic sector.

Data sourced from `sectorcontacts.csv` via `DataLoader`. When `scaled=true`
(default), values are divided element-wise by sector workforce counts so that
the result is contacts per worker (as used in the ODE force-of-infection).

# Arguments
- `workers`: Worker counts for each of the 45 economic sectors
- `scaled::Bool`: If true (default), contacts are divided by worker counts

# Returns
A `Vector{Float64}` of length 45 with per-capita within-sector contact rates
"""
function worker_contacts(workers; scaled = true)::Vector{Float64}
    x = Vector{Float64}(DataLoader.get_economic_contacts().contacts_workplace)

    if scaled
        x = x ./ workers
    end

    return x
end

"""
    consumer_worker_contacts(demography; scaled=true)::Matrix{Float64}

Get the consumer-worker contact matrix (45 sectors × 4 age groups).

Represents contacts between consumer-sector workers and age groups.

# Arguments
- `demography::Vector{Float64}`: Population vector for 4 age groups
- `scaled::Bool`: If true (default), contacts are scaled by age-group demographics

# Returns
A `Matrix{Float64}` of size (45, 4) with consumer-worker contact rates
"""
function consumer_worker_contacts(demography; scaled = true)::Matrix{Float64}
    ccw = repeat([1.0], N_ECON_GROUPS * N_AGE_GROUPS)
    ccw = reshape(ccw, N_ECON_GROUPS, N_AGE_GROUPS)

    # colwise div by size of from groups
    if scaled
        ccw *= Diagonal(1 ./ demography[i_AGE_GROUPS])
    end

    return ccw
end

"""
    prepare_community_contacts(cm; scaled=true)::Matrix{Float64}

Get a 49×49 community-only contact matrix for all age-groups and economic
sectors. Unlike `prepare_contacts`, this function does **not** add
within-sector workplace contacts or consumer-worker contacts to the matrix.
Those routes are kept separate for use in the ODE force-of-infection
calculation (see `plan_ode.md`).

# Arguments
- `cm::Matrix{Float64}`: 4×4 contact matrix for age groups
- `scaled::Bool`: If true (default), contact matrix is scaled by demographics

# Returns
A `Matrix{Float64}` of size (49, 49) with community contacts only
"""
function prepare_community_contacts(cm; scaled = true)::Matrix{Float64}
    cm_x = ones(N_TOTAL_GROUPS, N_TOTAL_GROUPS) .* cm[i_WORKING_AGE, i_WORKING_AGE]
    cm_x[i_AGE_GROUPS, i_AGE_GROUPS] = cm
    cm_x[i_AGE_GROUPS, i_ECON_GROUPS] .= cm[:, i_WORKING_AGE]
    cm_x[i_ECON_GROUPS, i_AGE_GROUPS] .= reshape(cm[i_WORKING_AGE, :], 1, N_AGE_GROUPS)
    # NOTE: workplace and consumer-worker contacts intentionally excluded

    if scaled
        cm_x *= Diagonal(1 ./ prepare_demog())
    end

    return cm_x
end

"""
    prepare_demog(demog, workers)::Vector{Float64}

Get a 49-element population vector for all age-groups and economic sectors.

Concatenates 4 age groups with 45 economic sector worker counts to form the
population vector used in force-of-infection calculations.

# Arguments
- `demog::Vector{Float64}`: Demographics (4 age groups)
- `workers::Vector{Int}`: Worker counts (45 economic sectors)

# Returns
A `Vector{Float64}` of length 49 (4 age groups + 45 workers)
"""
function prepare_demog(demog, workers)::Vector{Float64}
    return [demog; workers]
end

"""
    prepare_demog(cd::CountryData) -> Vector

Get the 49-element population vector (4 age groups + 45 worker sectors) for
any country from a [`DataLoader.CountryData`](@ref) struct.

Worker counts are clamped to a minimum of 1 to avoid division-by-zero when
this vector is used as a denominator (e.g. in contact-matrix scaling and the
Rt callback). This matches the `+1` padding applied in `initial_state`.
"""
function prepare_demog(cd::CountryData)
    return prepare_demog(cd.demography, max.(cd.workers, 1))
end

"""
    worker_contacts(cd::CountryData; scaled=true) -> SVector

Get per-capita social contacts within each economic sector using workforce
counts from `cd`. Sectors with zero workers are treated as having 1 worker
to avoid division by zero.
"""
function worker_contacts(cd::CountryData; scaled = true)
    x = Vector{Float64}(DataLoader.get_economic_contacts().contacts_workplace)
    if scaled
        x = x ./ max.(cd.workers, 1)
    end
    return x
end

"""
    expand_contacts(cm::Matrix{Float64})::Matrix{Float64}

Expand a 4×4 contact matrix to 49×49 covering all age groups and economic sectors.

Maps the 4×4 age-group contact matrix to the full 49×49 population, assuming
economic sectors adopt the contact patterns of the working-age group.

# Arguments
- `cm::Matrix{Float64}`: 4×4 contact matrix for age groups

# Returns
A `Matrix{Float64}` of size (49, 49)
"""
function expand_contacts(cm::Matrix{Float64})::Matrix{Float64}
    cm_x = ones(N_TOTAL_GROUPS, N_TOTAL_GROUPS) .* cm[i_WORKING_AGE, i_WORKING_AGE]
    cm_x[i_AGE_GROUPS, i_AGE_GROUPS] = cm
    cm_x[i_AGE_GROUPS, i_ECON_GROUPS] .= cm[:, i_WORKING_AGE]
    cm_x[i_ECON_GROUPS, i_AGE_GROUPS] .= reshape(cm[i_WORKING_AGE, :], 1, N_AGE_GROUPS)

    # cm_x[i_ECON_GROUPS, i_AGE_GROUPS] += consumer_worker_contacts(
    #     cd.demography; scaled = false)
    # cm_x[i_ECON_GROUPS, i_ECON_GROUPS] += Diagonal(worker_contacts(
    #     cd.workers; scaled = false))

    return cm_x
end

"""
    consumer_worker_contacts(cd::CountryData; scaled=true)

Get the 45×4 consumer-worker contact matrix scaled by the age-group
demography in `cd`.
"""
function consumer_worker_contacts(cd::CountryData; scaled = true)
    return consumer_worker_contacts(cd.demography; scaled = scaled)
end

"""
    prepare_contacts(cd::CountryData; scaled=true)

Get a 49×49 contact matrix for all age-groups and economic sectors using
the demographic and contact data in `cd`.
"""
function prepare_contacts(cd::CountryData; scaled = true)
    cm = cd.contact_matrix

    if isa(cm, Vector)
        cm_x = [expand_contacts(i) for i in cm]
        if scaled
            demog_diag = Diagonal(1 ./ prepare_demog(cd))
            cm_x = [i * demog_diag for i in cm_x]
        end

        return cm_x
    else
        cm_x = expand_contacts(cm)
        if scaled
            cm_x *= Diagonal(1 ./ prepare_demog(cd))
        end

        return cm_x
    end
end

"""
    total_contacts(contacts) -> Matrix{Float64}

Reduce a contact matrix representation to a single aggregated 2D matrix.

When `contacts` is a `Vector{Matrix{Float64}}` (one matrix per closure setting),
the matrices are summed element-wise to produce the total contact matrix. When
`contacts` is already a single `Matrix{Float64}`, it is returned unchanged.

# Arguments
- `contacts`: Either a `Vector{Matrix{Float64}}` or a `Matrix{Float64}`

# Returns
A single `Matrix{Float64}` representing total contacts across all settings.
"""
function total_contacts(contacts::Union{
        Vector{Matrix{Float64}}, Matrix{Float64}})::Matrix{Float64}
    isa(contacts, Vector) ? sum(contacts) : contacts
end

"""
    contacts3d(cd::CountryData) -> Array{Float64, 3}

Return the scaled contact matrices for a country as a single 3D array of size
`(N_TOTAL_GROUPS, N_TOTAL_GROUPS, K)`, where `K` is the number of closure
settings (economic contact scenarios).

Calls `prepare_contacts(cd)` and normalises the result into 3D form:
- If the country has multiple contact matrices (one per closure setting),
  they are stacked along the third dimension to give shape `(49, 49, K)`.
- If the country has a single contact matrix, it is reshaped to `(49, 49, 1)`
  so the return type is always a 3D array regardless of `K`.

The third dimension of the returned array is the one contracted by
`weighted_slice_sum!` in the ODE force-of-infection calculation.

# Arguments
- `cd::CountryData`: Country data struct from `DataLoader`

# Returns
An `Array{Float64, 3}` of size `(N_TOTAL_GROUPS, N_TOTAL_GROUPS, K)`.
"""
function contacts3d(cd::CountryData)::Array{Float64, 3}
    cm = prepare_contacts(cd)

    if isa(cm, Vector)
        settings = length(cm)
        return stack(cm) # type is Array{Float64, 3} with dims 49, 49, K
    else
        dim_main = first(size(cm))
        return reshape(cm, dim_main, dim_main, 1) # K = 1
    end
end

"""
    get_settings(cd::CountryData)::Int

Return the number of contact matrix settings (closure strategies) for a country.

# Arguments
- `cd::CountryData`: Country data struct

# Returns
The number of contact matrix settings: typically 1 (single matrix) or >1 (multiple scenarios)
"""
function get_settings(cd::CountryData)::Int
    isa(cd.contact_matrix, Vector) ? length(cd.contact_matrix) : 1
end

"""
    prepare_community_contacts(cd::CountryData; scaled=true)

Get a 49×49 community-only contact matrix for all age-groups and economic
sectors using `cd`. Workplace and consumer-worker contacts are excluded.
"""
function prepare_community_contacts(cd::CountryData; scaled = true)
    cm = cd.contact_matrix
    cm_x = ones(N_TOTAL_GROUPS, N_TOTAL_GROUPS) .* cm[i_WORKING_AGE, i_WORKING_AGE]
    cm_x[i_AGE_GROUPS, i_AGE_GROUPS] = cm
    cm_x[i_AGE_GROUPS, i_ECON_GROUPS] .= cm[:, i_WORKING_AGE]
    cm_x[i_ECON_GROUPS, i_AGE_GROUPS] .= reshape(cm[i_WORKING_AGE, :], 1, N_AGE_GROUPS)

    if scaled
        cm_x *= Diagonal(1 ./ prepare_demog(cd))
    end

    return cm_x
end

"""
    initial_state(cd::CountryData) -> Array{Float64,3}

Construct the initial epidemic state (N_TOTAL_GROUPS × N_COMPARTMENTS ×
N_VACCINE_STRATA) for a country given as a [`DataLoader.CountryData`](@ref)
struct. A fraction `1e-6` of each group is seeded as symptomatic infectious.

The state has 9 compartments: 7 primary epidemiological (S, E, Is, Ia, H, R, D)
and 2 data tracking (newI, newH for cumulative infections and hospitalisations).
"""
function initial_state(cd::CountryData)
    demography = cd.demography
    workers = cd.workers .+ 1.0

    p_infected = 1e-6
    p_susc = 1.0 - p_infected

    zero_compartments = zeros(N_COMPARTMENTS - 3)
    init = [p_susc, 0.0, p_infected]
    init = [init; zero_compartments]
    init = reshape(init, 1, N_COMPARTMENTS)
    init = repeat(init, N_TOTAL_GROUPS)

    dummy = zeros(N_TOTAL_GROUPS, N_COMPARTMENTS, N_VACCINE_STRATA)
    dummy[:, :, i_UNVAX_STRATUM] = init

    demog = copy(demography)
    inactive_workers = demog[i_WORKING_AGE] - sum(workers)
    demog[i_WORKING_AGE] = inactive_workers
    demog = [demog; workers]

    return dummy .* demog
end

"""
    initial_state(country::String) -> Array{Float64,3}

Construct the initial epidemic state for a named country. Calls
[`DataLoader.get_country`](@ref) and delegates to `initial_state(::CountryData)`.

# Example
```julia
state = initial_state("United Kingdom")
```
"""
function initial_state(country::String)
    return initial_state(DataLoader.get_country(country))
end

"""
    prepare_contacts(country::String; scaled=true)

Get the 49×49 contact matrix for a named country.
"""
function prepare_contacts(country::String; scaled = true)
    return prepare_contacts(DataLoader.get_country(country); scaled = scaled)
end

"""
    prepare_community_contacts(country::String; scaled=true)

Get the 49×49 community-only contact matrix for a named country.
"""
function prepare_community_contacts(country::String; scaled = true)
    return prepare_community_contacts(DataLoader.get_country(country); scaled = scaled)
end

"""
    worker_contacts(country::String; scaled=true) -> SVector

Get per-capita within-sector contact rates for a named country.
"""
function worker_contacts(country::String; scaled = true)
    return worker_contacts(DataLoader.get_country(country); scaled = scaled)
end

"""
    prepare_demog(country::String) -> Vector

Get the 49-element population vector for a named country.
"""
function prepare_demog(country::String)
    return prepare_demog(DataLoader.get_country(country))
end

end
