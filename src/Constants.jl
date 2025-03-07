# DAEDALUS related constants copied over from {daedalus}

module Constants

export iS, iE, iIs, iIa, iH, iR, iD, N_COMPARTMENTS,
    i_AGE_GROUPS, i_WORKING_AGE, i_ECON_GROUPS,
    N_AGE_GROUPS, N_ECON_GROUPS, N_TOTAL_GROUPS, N_VACCINE_STRATA,
    i_UNVAX_STRATUM, i_VAX_STRATUM

const iS = 1
const iE = 2
const iIs = 3
const iIa = 4
const iH = 5
const iR = 6
const iD = 7

const N_COMPARTMENTS = 7

const i_AGE_GROUPS = 1:4
const i_WORKING_AGE = 3
const i_ECON_GROUPS = 5:49
const N_AGE_GROUPS = 4
const N_ECON_GROUPS = 45
const N_TOTAL_GROUPS = N_AGE_GROUPS + N_ECON_GROUPS
const N_VACCINE_STRATA = 2
const i_UNVAX_STRATUM = 1
const i_VAX_STRATUM = 2

end
