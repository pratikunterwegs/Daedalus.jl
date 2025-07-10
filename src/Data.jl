# copied manually from jameel-institute/daedalus
# will be moved over to package data eventually

module Data

using ..Constants
using StaticArrays

export australia_initial_state, australia_demography, prepare_contacts,
    worker_contacts, prepare_demog

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

Get a dummy value of worker contacts within economic sectors.

This data is synthetic and not generated from the R package {daedalus}.
"""
function worker_contacts(workers=aus_workers())
    # in proportion to workforce and scaled by workforce
    return SVector{N_ECON_GROUPS}(2 .+ workers / sum(workers)) ./ workers
end

"""
    prepare_contacts()

Get a contact matrix for all age-groups and economic sectors as a StaticArrays
    `SMatrix`.
"""
# Function to prepare 49x49 community contacts matrix
function prepare_contacts(cm=australia_contacts())
    cm_x = ones(N_TOTAL_GROUPS, N_TOTAL_GROUPS) .* cm[i_WORKING_AGE, i_WORKING_AGE]
    cm_x[i_AGE_GROUPS, i_AGE_GROUPS] = cm
    cm_x[i_AGE_GROUPS, i_ECON_GROUPS] .= cm[:, i_WORKING_AGE]
    cm_x[i_ECON_GROUPS, i_AGE_GROUPS] .= reshape(cm[i_WORKING_AGE, :], 1, N_AGE_GROUPS)

    return SMatrix{N_TOTAL_GROUPS,N_TOTAL_GROUPS}(cm_x)
end

"""
    prepare_demog()

Get a population vector for all age-groups and economic sectors for the force of
    infection calculation for community infections.
"""
# Function to prepare 49 element demography vector for I/N in FOI calculation
function prepare_demog(demog=australia_demography(), workers=aus_workers())
    demog_x = ones(N_TOTAL_GROUPS) * demog[i_WORKING_AGE]
    demog_x[i_AGE_GROUPS] = demog

    return demog_x
end

"""
    australia_initial_state()

Initial state with seven compartments: S, E, Is, Ia, R, D, and H, and four
age groups and 49 economic sectors (comprised of working age individuals), and
two vaccination strata for unvaccinated and vaccinated individuals.
"""
function australia_initial_state(
    demography=australia_demography(), workers=aus_workers())

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

end
