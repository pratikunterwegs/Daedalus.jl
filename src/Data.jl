# copied manually from jameel-institute/daedalus
# will be moved over to package data eventually

module Data

using ..Constants
using ..DataLoader
using LinearAlgebra

export prepare_contacts, contacts3d, get_settings, total_contacts,
       worker_contacts, consumer_worker_contacts,
       prepare_demog, initial_state

"""
    worker_contacts(cd::CountryData; scaled=true) -> SVector

Get per-capita social contacts within each economic sector from `cd`.
Sectors with zero workers are treated as having 1 worker
to avoid division by zero.

Data sourced from `sectorcontacts.csv` via `DataLoader`. When `scaled=true`
(default), values are divided element-wise by sector workforce counts so that
the result is contacts per worker (as used in the ODE force-of-infection).

# Arguments
- `cd`: Used to access worker counts for each of the 45 economic sectors
- `scaled::Bool`: If true (default), contacts are divided by worker counts

# Returns
A `Vector{Float64}` of length 45 with per-capita within-sector contact rates
"""
function worker_contacts(cd::CountryData; scaled = true)
    x = Vector{Float64}(DataLoader.get_economic_contacts().contacts_workplace)
    if scaled
        x = x ./ max.(cd.workers, 1)
    end
    return x
end

"""
    prepare_demog(cd::CountryData) -> Vector

Get a 49-element population vector for all age-groups and economic sectors.

Concatenates 4 age groups with 45 economic sector worker counts to form the
population vector used in force-of-infection calculations.

# Arguments
- `cd::CountryData` - gives the country data

# Returns
A `Vector{Float64}` of length 49 (4 age groups + 45 workers)
"""
function prepare_demog(cd::CountryData)
    return [cd.demography; max.(cd.workers, 1.0)] # prevent zero division
end

"""
    consumer_worker_contacts(demography; scaled=true)::Matrix{Float64}

Get the consumer-worker contact matrix (45 sectors × 4 age groups).

Represents contacts between consumer-sector workers and age groups.

# Arguments
- `cd::CountryData`: Population vector for 4 age groups
- `scaled::Bool`: If true (default), contacts are scaled by age-group demographics

# Returns
A `Matrix{Float64}` of size (45, 4) with consumer-worker contact rates
"""
function consumer_worker_contacts(cd::CountryData; scaled = true)::Matrix{Float64}
    ccw = repeat([1.0], N_ECON_GROUPS * N_AGE_GROUPS)
    ccw = reshape(ccw, N_ECON_GROUPS, N_AGE_GROUPS)

    demog = prepare_demog(cd)

    # TODO: NEEDS TO READ REAL DATA, CURRENTLY USES DUMMY DATA
    # colwise div by size of from groups
    if scaled
        ccw *= Diagonal(1 ./ demog[i_AGE_GROUPS])
    end

    return ccw
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
    prepare_contacts(cd::CountryData; scaled=true)

Get contact matrices for all age-groups and economic sectors using demographic
and contact data in `cd`.

Delegates to `contacts3d()` which includes both community and workplace layers.
For use with `total_contacts()` in beta calibration and NGM calculations.
"""
function prepare_contacts(cd::CountryData; scaled = true)
    contacts3d(cd; scaled = scaled)
end

"""
    total_contacts(contacts) -> Matrix{Float64}

Reduce a contact matrix representation to a single aggregated 2D matrix.

When `contacts` is a `Vector{Matrix{Float64}}` (one matrix per closure setting),
the matrices are summed element-wise to produce the total contact matrix. When
`contacts` is a `Matrix{Float64}`, it is returned unchanged. When `contacts` is
a 3D `Array{Float64,3}` (e.g., from `contacts3d`), the slices are summed.

# Arguments
- `contacts`: Either a `Vector{Matrix{Float64}}`, a `Matrix{Float64}`, or `Array{Float64,3}`

# Returns
A single `Matrix{Float64}` representing total contacts across all settings/layers.
"""
function total_contacts(contacts::Union{
        Vector{Matrix{Float64}}, Matrix{Float64}, Array{Float64, 3}})::Matrix{Float64}
    if isa(contacts, Array{Float64, 3})
        dropdims(sum(contacts; dims = 3); dims = 3)
    elseif isa(contacts, Vector)
        sum(contacts)
    else
        contacts
    end
end

"""
    contacts3d(cd::CountryData; scaled=true) -> Array{Float64, 3}

Return contact matrices for a country as a 3D array of size
`(N_TOTAL_GROUPS, N_TOTAL_GROUPS, K+1)`, where K is the number of community
contact settings and the final slice is the workplace contact matrix.

Structure:
- Slices 1..K: Community contact matrices (expanded from 4×4 age-specific matrices)
- Slice K+1: Workplace contact matrix (worker-worker diagonal + consumer-worker block)

For countries with a single contact matrix (common case), returns `(49, 49, 2)`:
- `[:, :, 1]` = community contacts
- `[:, :, 2]` = workplace contacts

The third dimension is contracted by `weighted_slice_sum!` in the ODE
force-of-infection calculation.

# Arguments
- `cd::CountryData`: Country data struct from `DataLoader`
- `scaled::Bool`: If true (default), scales contacts by demographics for FOI calculation

# Returns
An `Array{Float64, 3}` of size `(N_TOTAL_GROUPS, N_TOTAL_GROUPS, K+1)`.
"""
function contacts3d(cd::CountryData; scaled::Bool = true)::Array{Float64, 3}
    cm = cd.contact_matrix
    demog = prepare_demog(cd)

    # Community layer(s): expand each 4×4 matrix to 49×49
    cms = isa(cm, Vector) ? [expand_contacts(c) for c in cm] : [expand_contacts(cm)]
    if scaled
        demog_diag = Diagonal(1.0 ./ demog)
        cms = [c * demog_diag for c in cms]
    end

    # Workplace layer: worker-worker diagonal + consumer-worker block
    cm_work = zeros(N_TOTAL_GROUPS, N_TOTAL_GROUPS)
    cm_work[i_ECON_GROUPS, i_ECON_GROUPS] .= Diagonal(worker_contacts(cd; scaled = scaled))
    cm_work[i_ECON_GROUPS, i_AGE_GROUPS] .= consumer_worker_contacts(cd; scaled = scaled)

    return stack([cms..., cm_work])
end

"""
    get_settings(cd::CountryData)::Int

Return the number of contact matrix layers (community settings + workplace) for a country.

# Arguments
- `cd::CountryData`: Country data struct

# Returns
The number of contact layers: K (community settings) + 1 (workplace).
Typically 2 for countries with one community contact matrix, or K+1 for multiple scenarios.
"""
function get_settings(cd::CountryData)::Int
    isa(cd.contact_matrix, Vector) ? length(cd.contact_matrix) + 1 : 2
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

    u = zeros(N_TOTAL_GROUPS, N_COMPARTMENTS, N_VACCINE_STRATA)
    u[:, :, i_UNVAX_STRATUM] = init

    demog = copy(demography)
    inactive_workers = demog[i_WORKING_AGE] - sum(workers)
    demog[i_WORKING_AGE] = inactive_workers
    demog = [demog; workers]

    u = u .* demog # N_TOT_GRP x N_COMP x N_VAX
    u = reshape(u, length(u))

    # add new vax and Reff
    u = [u; zeros(N_TOTAL_GROUPS); 0.0]

    return u
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

Get contact matrices for a named country. Delegates to `prepare_contacts(cd)`.
"""
function prepare_contacts(country::String; scaled = true)
    prepare_contacts(DataLoader.get_country(country); scaled = scaled)
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
