# copied manually from jameel-institute/daedalus
# will be moved over to package data eventually

module Data

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
    return (2 .+ workers / sum(workers)) ./ workers
end

"""
    prepare_contacts()

Get a contact matrix for all age-groups and economic sectors as a StaticArrays
    `SMatrix`.
"""
# Function to prepare 49x49 community contacts matrix
function prepare_contacts(cm=australia_contacts())
    cm_x = ones(49, 49) .* cm[3, 3]
    cm_x[1:4, 1:4] = cm
    cm_x[1:4, 5:49] .= cm[:, 3]
    cm_x[5:49, 1:4] .= reshape(cm[3, :], 1, 4)

    return SMatrix{49,49}(cm_x)
end

"""
    prepare_demog()

Get a population vector for all age-groups and economic sectors for the force of
    infection calculation for community infections.
"""
# Function to prepare 49 element demography vector for I/N in FOI calculation
function prepare_demog(demog=australia_demography(), workers=aus_workers())
    demog_x = ones(49) * demog[3]
    demog_x[1:4] = demog

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
    p_infected = 1e-6
    p_susc = 1 - p_infected
    compartments = 8
    age_groups = 4
    econ_groups = 45

    init = [p_susc, 0.0, p_infected, 0.0, 0.0, 0.0, 0.0, 0.0]
    init = reshape(init, 1, compartments)
    init = repeat(init, age_groups + econ_groups)

    # process demography to calculate non-working working age, and concat
    # demography and worker counts
    i_working_age = 3
    inactive_workers = demography[i_working_age] - sum(workers)
    demography[i_working_age] = inactive_workers

    demography = [demography; workers]

    # multiply by demography - row i times element i
    init = init .* demography

    return init
end

end
