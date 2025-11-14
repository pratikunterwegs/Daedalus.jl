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

"""
    get_indices(compartment::String)

Get compartment indices for state tests.
    """
get_indices(compartment::String) = begin
    if compartment in COMPARTMENTS
        idx_comp = eval(Symbol("i" * compartment))
        final_idx = idx_comp * N_TOTAL_GROUPS
        first_idx = final_idx - N_TOTAL_GROUPS + 1

        return first_idx:final_idx
    elseif compartment == "Rt"
        return N_TOTAL_GROUPS * N_COMPARTMENTS * N_VACCINE_STRATA + i_rel_Rt
    else
        error("Compartment ", compartment, " not available!")
    end
end

end
