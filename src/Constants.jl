# DAEDALUS related constants copied over from {daedalus}

module Constants

export iS, iE, iIs, iIa, iH, iR, iD, inewI, inewH, N_COMPARTMENTS, COMPARTMENTS,
       i_AGE_GROUPS, i_WORKING_AGE, i_ECON_GROUPS,
       N_AGE_GROUPS, N_ECON_GROUPS, N_TOTAL_GROUPS, N_VACCINE_STRATA, N_SEXES,
       i_UNVAX_STRATUM, i_VAX_STRATUM, i_rel_Rt, i_rel_Rt_cont,
       get_indices, get_age_group_indices

const iS = 1
const iE = 2
const iIs = 3
const iIa = 4
const iH = 5
const iR = 6
const iD = 7
const inewI = 8
const inewH = 9

const N_COMPARTMENTS = 9
const N_DATA_COMPARTMENTS = 2
const N_PRIMARY_COMPARTMENTS = N_COMPARTMENTS - N_DATA_COMPARTMENTS

const COMPARTMENTS = ["S", "E", "Is", "Ia", "H", "R", "D", "new_inf",
    "new_hosp", "new_vax"]
const PRIMARY_COMPARTMENTS = COMPARTMENTS[1:N_PRIMARY_COMPARTMENTS]

const i_AGE_GROUPS = 1:4
const i_WORKING_AGE = 3
const i_ECON_GROUPS = 5:49
const N_AGE_GROUPS = 4
const N_ECON_GROUPS = 45
const N_TOTAL_GROUPS = N_AGE_GROUPS + N_ECON_GROUPS
const N_VACCINE_STRATA = 2
const i_UNVAX_STRATUM = 1
const i_VAX_STRATUM = 2

const i_rel_Rt = 1
const i_rel_Rt_cont = i_rel_Rt + 1

const N_SEXES = 2  # Male and female for life expectancy data averaging

# Maps compartment name → its 1-based block index in the state vector.
# Avoids eval() in get_indices for fast, type-stable lookups.
const COMPARTMENT_INDICES = Dict{String, Int}(
    "S" => iS,
    "E" => iE,
    "Is" => iIs,
    "Ia" => iIa,
    "H" => iH,
    "R" => iR,
    "D" => iD,
    "new_inf" => inewI,
    "new_hosp" => inewH
)

"""
    get_indices(compartment, groups=nothing)

Return the linear indices into the ODE state vector for `compartment` and
`groups`.

The state vector is laid out as `N_VACCINE_STRATA` repetitions of
`N_COMPARTMENTS` blocks, each block containing `N_TOTAL_GROUPS` consecutive
entries (one per demographic/economic group). The Rt scalar parameters are
appended after all epi-compartment entries.

# Arguments
- `compartment`: compartment name — one of $(join(COMPARTMENTS, ", ")) or `"Rt"`.
- `groups`: groups to select within the compartment block:
  - `nothing` (default) — all groups (`1:N_TOTAL_GROUPS`)
  - `Int` — a single group
  - `AbstractVector{Int}` or `UnitRange{Int}` — an arbitrary subset of groups

# Returns
A `UnitRange{Int}` when all groups are requested or a contiguous range is
passed, otherwise a `Vector{Int}`.
"""
function get_indices(
        compartment::String,
        groups::Union{Nothing, Int, AbstractVector{Int}, UnitRange{Int}} = nothing
)
    if compartment == "Rt"
        return N_TOTAL_GROUPS * N_COMPARTMENTS * N_VACCINE_STRATA +
               N_TOTAL_GROUPS + i_rel_Rt
    end

    if compartment == "vax"
        start = N_TOTAL_GROUPS * N_COMPARTMENTS + 1
        vax_range = start:((start - 1) * 2)
        isnothing(groups) && return vax_range
        return vax_range[groups]
    end

    if compartment == "new_vax"
        start = N_TOTAL_GROUPS * N_COMPARTMENTS * N_VACCINE_STRATA + 1
        new_vax_range = start:(start + N_TOTAL_GROUPS - 1)
        isnothing(groups) && return new_vax_range
        return new_vax_range[groups]
    end

    if compartment == "I"
        # Total infectious: Ia + Is
        return vcat(get_indices("Ia", groups), get_indices("Is", groups))
    end

    idx_comp = get(COMPARTMENT_INDICES, compartment, nothing)
    isnothing(idx_comp) && error("Compartment \"", compartment, "\" not available!")

    first_idx = (idx_comp - 1) * N_TOTAL_GROUPS + 1
    comp_range = first_idx:(first_idx + N_TOTAL_GROUPS - 1)

    isnothing(groups) && return comp_range
    return comp_range[groups]
end

"""
    get_age_group_indices(age_group::Int)

Return the indices within the N_TOTAL_GROUPS vector for a given age group (1-4).

The 49-group vector is partitioned into 4 age groups:
- Age group 1 (0-4 years): indices 1-12
- Age group 2 (5-19 years): indices 13-24
- Age group 3 (20-64 years): indices 25-36
- Age group 4 (65+ years): indices 37-49

# Arguments
- `age_group`: Age group number (1-4)

# Returns
A `UnitRange{Int}` for the indices of that age group
"""
function get_age_group_indices(age_group::Int)
    age_group < 1 ||
        age_group > N_AGE_GROUPS && error("Age group must be 1-$(N_AGE_GROUPS)")
    starts = [1, 13, 25, 37]
    ends = [12, 24, 36, 49]
    return starts[age_group]:ends[age_group]
end

end
