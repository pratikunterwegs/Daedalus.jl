# copied manually from jameel-institute/daedalus
# will be moved over to package data eventually

module Data

using ..Constants
using ..DataLoader
using LinearAlgebra
using StaticArrays

export australia_initial_state, australia_demography, prepare_contacts,
       prepare_community_contacts, worker_contacts, consumer_worker_contacts,
       prepare_demog, australia_contacts, initial_state

"""
    australia_demography()

Get number of individuals in each DAEDALUS demography group.
"""
function australia_demography()
    return [1669707, 4769747, 14926119, 4134308]
end

"""
    australia_contacts()

Get contacts between DAEDALUS demography groups as a matrix.
"""
function australia_contacts()
    cm = [[3.7187500 2.5982168 5.739112 0.2728101]
          [0.9095369 13.0623374 5.741992 0.5229291]
          [0.6420045 1.8348941 11.256655 1.0003495]
          [0.1347582 0.6540519 3.760931 2.5421895]]

    return cm
end

"""
    aus_workers()

Get number of workers in DAEDALUS economic sectors.

All sectors have one additional worker to avoid zero division.
"""
function aus_workers()
    return [
        331714, 10730, 110740, 113170, 65189, 224517, 35631, 47879, 45956,
        8604, 31481, 18322, 30266, 30367, 52435, 69643, 38263, 21966,
        59115, 33101, 32696, 195972, 79259, 73591, 1199314, 1783483, 319770,
        66707, 62051, 151331, 89280, 899587, 91912, 90495, 270878, 418059,
        166920, 823060, 567062, 859198, 1059016, 1686004, 277154, 268246, 0
    ] .+ 1
end

"""
    worker_contacts()

Get per-capita social contacts within each economic sector.

Data sourced from `sectorcontacts.csv` via `DataLoader`. When `scaled=true`
(default), values are divided element-wise by sector workforce counts so that
the result is contacts per worker (as used in the ODE force-of-infection).
"""
function worker_contacts(workers = aus_workers(); scaled = true)
    x = Vector{Float64}(DataLoader.get_economic_contacts().contacts_workplace)

    if scaled
        x = x ./ workers
    end

    return SVector{N_ECON_GROUPS}(x)
end

function consumer_worker_contacts(demography = australia_demography(); scaled = true)
    ccw = repeat([1.0], N_ECON_GROUPS * N_AGE_GROUPS)
    ccw = reshape(ccw, N_ECON_GROUPS, N_AGE_GROUPS)

    # colwise div by size of from groups
    if scaled
        ccw *= Diagonal(1 ./ demography[i_AGE_GROUPS])
    end

    return SMatrix{N_ECON_GROUPS, N_AGE_GROUPS}(ccw)
end

"""
    prepare_contacts()

Get a contact matrix for all age-groups and economic sectors as a StaticArrays
    `SMatrix`. Operates on default data.
"""
# Function to prepare 49x49 community contacts matrix
function prepare_contacts(cm = australia_contacts(); scaled = true)
    # community contacts, colwise div by size of from-group
    cm_x = ones(N_TOTAL_GROUPS, N_TOTAL_GROUPS) .* cm[i_WORKING_AGE, i_WORKING_AGE]
    cm_x[i_AGE_GROUPS, i_AGE_GROUPS] = cm
    cm_x[i_AGE_GROUPS, i_ECON_GROUPS] .= cm[:, i_WORKING_AGE]
    cm_x[i_ECON_GROUPS, i_AGE_GROUPS] .= reshape(cm[i_WORKING_AGE, :], 1, N_AGE_GROUPS)

    cm_x[i_ECON_GROUPS, i_AGE_GROUPS] += consumer_worker_contacts(scaled = false)
    cm_x[i_ECON_GROUPS, i_ECON_GROUPS] += Diagonal(worker_contacts(scaled = false))

    if scaled
        cm_x *= Diagonal(1 ./ prepare_demog())
    end

    return SMatrix{N_TOTAL_GROUPS, N_TOTAL_GROUPS}(cm_x)
end

"""
    prepare_community_contacts()

Get a 49×49 community-only contact matrix for all age-groups and economic
sectors. Unlike `prepare_contacts`, this function does **not** add
within-sector workplace contacts or consumer-worker contacts to the matrix.
Those routes are kept separate for use in the ODE force-of-infection
calculation (see `plan_ode.md`).
"""
function prepare_community_contacts(cm = australia_contacts(); scaled = true)
    cm_x = ones(N_TOTAL_GROUPS, N_TOTAL_GROUPS) .* cm[i_WORKING_AGE, i_WORKING_AGE]
    cm_x[i_AGE_GROUPS, i_AGE_GROUPS] = cm
    cm_x[i_AGE_GROUPS, i_ECON_GROUPS] .= cm[:, i_WORKING_AGE]
    cm_x[i_ECON_GROUPS, i_AGE_GROUPS] .= reshape(cm[i_WORKING_AGE, :], 1, N_AGE_GROUPS)
    # NOTE: workplace and consumer-worker contacts intentionally excluded

    if scaled
        cm_x *= Diagonal(1 ./ prepare_demog())
    end

    return SMatrix{N_TOTAL_GROUPS, N_TOTAL_GROUPS}(cm_x)
end

"""
    prepare_demog()

Get a population vector for all age-groups and economic sectors for the force of
    infection calculation for community infections.
"""
# Function to prepare 49 element demography vector for I/N in FOI calculation
function prepare_demog(demog = australia_demography(), workers = aus_workers())
    return [demog; workers]
end

"""
    australia_initial_state()

Initial state with seven compartments: S, E, Is, Ia, R, D, and H, and four
age groups and 49 economic sectors (comprised of working age individuals), and
two vaccination strata for unvaccinated and vaccinated individuals.
"""
function australia_initial_state(
        demography = australia_demography(), workers = aus_workers())

    # assume all infectious go into Is
    p_infected = 1e-6
    p_susc = 1.0 - p_infected

    zero_compartments = zeros(N_COMPARTMENTS - 3)
    init = [p_susc, 0.0, p_infected]
    init = [init; zero_compartments]
    init = reshape(init, 1, N_COMPARTMENTS)
    init = repeat(init, N_TOTAL_GROUPS)

    # assign to first layer
    dummy = zeros(N_TOTAL_GROUPS, N_COMPARTMENTS, N_VACCINE_STRATA)
    dummy[:, :, i_UNVAX_STRATUM] = init

    # process demography to calculate non-working working age, and concat
    # demography and worker counts
    inactive_workers = demography[i_WORKING_AGE] - sum(workers)
    demography[i_WORKING_AGE] = inactive_workers

    demography = [demography; workers]

    # multiply by demography - row i times element i
    init = dummy .* demography

    return init
end

"""
    prepare_demog(cd::CountryData) -> Vector

Get the 49-element population vector (4 age groups + 45 worker sectors) for
any country from a [`DataLoader.CountryData`](@ref) struct.
"""
function prepare_demog(cd::CountryData)
    return prepare_demog(cd.demography, cd.workers)
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
    return SVector{N_ECON_GROUPS}(x)
end

"""
    consumer_worker_contacts(cd::CountryData; scaled=true) -> SMatrix

Get the 45×4 consumer-worker contact matrix scaled by the age-group
demography in `cd`.
"""
function consumer_worker_contacts(cd::CountryData; scaled = true)
    return consumer_worker_contacts(cd.demography; scaled = scaled)
end

"""
    prepare_contacts(cd::CountryData; scaled=true) -> SMatrix

Get a 49×49 contact matrix for all age-groups and economic sectors using
the demographic and contact data in `cd`.
"""
function prepare_contacts(cd::CountryData; scaled = true)
    cm = cd.contact_matrix
    cm_x = ones(N_TOTAL_GROUPS, N_TOTAL_GROUPS) .* cm[i_WORKING_AGE, i_WORKING_AGE]
    cm_x[i_AGE_GROUPS, i_AGE_GROUPS] = cm
    cm_x[i_AGE_GROUPS, i_ECON_GROUPS] .= cm[:, i_WORKING_AGE]
    cm_x[i_ECON_GROUPS, i_AGE_GROUPS] .= reshape(cm[i_WORKING_AGE, :], 1, N_AGE_GROUPS)

    cm_x[i_ECON_GROUPS, i_AGE_GROUPS] += consumer_worker_contacts(
        cd.demography; scaled = false)
    cm_x[i_ECON_GROUPS, i_ECON_GROUPS] += Diagonal(worker_contacts(
        cd.workers; scaled = false))

    if scaled
        cm_x *= Diagonal(1 ./ prepare_demog(cd))
    end

    return SMatrix{N_TOTAL_GROUPS, N_TOTAL_GROUPS}(cm_x)
end

"""
    prepare_community_contacts(cd::CountryData; scaled=true) -> SMatrix

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

    return SMatrix{N_TOTAL_GROUPS, N_TOTAL_GROUPS}(cm_x)
end

"""
    initial_state(cd::CountryData) -> Array{Float64,3}

Construct the initial epidemic state (N_TOTAL_GROUPS × N_COMPARTMENTS ×
N_VACCINE_STRATA) for a country given as a [`DataLoader.CountryData`](@ref)
struct. A fraction `1e-6` of each group is seeded as symptomatic infectious.
"""
function initial_state(cd::CountryData)
    demography = cd.demography
    workers = cd.workers

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
result = daedalus(
    initial_state = state,
    contacts      = prepare_contacts("United Kingdom"),
    cw            = worker_contacts("United Kingdom"),
)
```
"""
function initial_state(country::String)
    return initial_state(DataLoader.get_country(country))
end

"""
    prepare_contacts(country::String; scaled=true) -> SMatrix

Get the 49×49 contact matrix for a named country.
"""
function prepare_contacts(country::String; scaled = true)
    return prepare_contacts(DataLoader.get_country(country); scaled = scaled)
end

"""
    prepare_community_contacts(country::String; scaled=true) -> SMatrix

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
