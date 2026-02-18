# DAEDALUS related constants copied over from {daedalus}

module Constants

export iS, iE, iIs, iIa, iH, iR, iD, N_COMPARTMENTS, COMPARTMENTS,
       i_AGE_GROUPS, i_WORKING_AGE, i_ECON_GROUPS,
       N_AGE_GROUPS, N_ECON_GROUPS, N_TOTAL_GROUPS, N_VACCINE_STRATA,
       i_UNVAX_STRATUM, i_VAX_STRATUM, i_rel_Rt, i_rel_Rt_cont,
       get_indices

const iS = 1
const iE = 2
const iIs = 3
const iIa = 4
const iH = 5
const iR = 6
const iD = 7

const N_COMPARTMENTS = 7

const COMPARTMENTS = ["S", "E", "Is", "Ia", "H", "R", "D"]

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

# Maps compartment name → its 1-based block index in the state vector.
# Avoids eval() in get_indices for fast, type-stable lookups.
const COMPARTMENT_INDICES = Dict{String, Int}(
    "S" => iS,
    "E" => iE,
    "Is" => iIs,
    "Ia" => iIa,
    "H" => iH,
    "R" => iR,
    "D" => iD
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
        return N_TOTAL_GROUPS * N_COMPARTMENTS * N_VACCINE_STRATA + i_rel_Rt
    end

    idx_comp = get(COMPARTMENT_INDICES, compartment, nothing)
    isnothing(idx_comp) && error("Compartment \"", compartment, "\" not available!")

    first_idx = (idx_comp - 1) * N_TOTAL_GROUPS + 1
    comp_range = first_idx:(first_idx + N_TOTAL_GROUPS - 1)

    isnothing(groups) && return comp_range
    return comp_range[groups]
end

end
