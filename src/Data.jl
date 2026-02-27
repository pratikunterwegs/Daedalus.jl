# copied manually from jameel-institute/daedalus
# will be moved over to package data eventually

module Data

using ..Constants
using ..DataLoader
using LinearAlgebra
using StaticArrays

export prepare_contacts, contacts3d, get_settings, total_contacts,
       prepare_community_contacts, worker_contacts, consumer_worker_contacts,
       prepare_demog, initial_state

"""
    worker_contacts()

Get per-capita social contacts within each economic sector.

Data sourced from `sectorcontacts.csv` via `DataLoader`. When `scaled=true`
(default), values are divided element-wise by sector workforce counts so that
the result is contacts per worker (as used in the ODE force-of-infection).
"""
function worker_contacts(workers; scaled = true)
    x = Vector{Float64}(DataLoader.get_economic_contacts().contacts_workplace)

    if scaled
        x = x ./ workers
    end

    return SVector{N_ECON_GROUPS}(x)
end

function consumer_worker_contacts(demography; scaled = true)
    ccw = repeat([1.0], N_ECON_GROUPS * N_AGE_GROUPS)
    ccw = reshape(ccw, N_ECON_GROUPS, N_AGE_GROUPS)

    # colwise div by size of from groups
    if scaled
        ccw *= Diagonal(1 ./ demography[i_AGE_GROUPS])
    end

    return SMatrix{N_ECON_GROUPS, N_AGE_GROUPS}(ccw)
end

"""
    prepare_community_contacts()

Get a 49×49 community-only contact matrix for all age-groups and economic
sectors. Unlike `prepare_contacts`, this function does **not** add
within-sector workplace contacts or consumer-worker contacts to the matrix.
Those routes are kept separate for use in the ODE force-of-infection
calculation (see `plan_ode.md`).
"""
function prepare_community_contacts(cm; scaled = true)
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
function prepare_demog(demog, workers)
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
    return SVector{N_ECON_GROUPS}(x)
end

function expand_contacts(cm::Matrix{Float64})
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

function total_contacts(contacts::Union{
        Vector{Matrix{Float64}}, Matrix{Float64}})::Matrix{Float64}
    isa(contacts, Vector) ? sum(contacts) : contacts
end

"""
    contacts3d(cd::CountryData) -> SArray
"""
function contacts3d(cd::CountryData)::SArray
    cm = prepare_contacts(cd)

    if isa(cm, Vector)
        settings = length(cm)
        return SArray{Tuple{N_TOTAL_GROUPS, N_TOTAL_GROUPS, settings}}(stack(cm))
    else
        return SArray{Tuple{N_TOTAL_GROUPS, N_TOTAL_GROUPS, 1}}(cm)
    end
end

function get_settings(cd::CountryData)::Int
    isa(cd.contact_matrix, Vector) ? length(cd.contact_matrix) : 1
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
