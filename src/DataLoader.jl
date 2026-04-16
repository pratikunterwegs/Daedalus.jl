
module DataLoader

export CountryData, InfectionData, VaccinationScenario, EconomicContacts,
       get_country, list_countries,
       get_pathogen, list_pathogens,
       get_economic_contacts, get_sector_names,
       get_vaccination_scenario, list_vaccination_scenarios,
       get_closure_strategy, list_closure_strategies

using CSV
using DataFrames
using Statistics

const DATA_DIR = joinpath(@__DIR__, "data")

"""
    CountryData

Demographic and economic data for a single country, matching the structure of
the `country_data` list in the R package `{daedalus.data}`.

This struct is mutable to allow reassigning values to examine different
    scenarios.
"""
mutable struct CountryData
    demography::Vector{Float64} # 4 age groups: 0-4, 5-19, 20-64, 65+
    contact_matrix::Union{Matrix{Float64}, Vector{Matrix{Float64}}} # matrix or list of matrices
    workers::Vector{Int} # 45 economic sectors
    gva::Vector{Float64} # 45 sectors, million USD/day
    hospital_capacity::Float64 # spare pandemic hospital beds
    gni::Float64 # GNI per capita PPP USD
    life_expectancy::Vector{Float64} # 4 age groups
    vsl::Vector{Float64} # value of statistical life, 4 age groups
end

"""
    InfectionData

Epidemiological parameters for a single pathogen, matching the structure of
`infection_data` in the R package `{daedalus.data}`.

This struct is mutable to allow reassigning values to examine different
    scenarios.
"""
mutable struct InfectionData
    r0::Float64
    sigma::Float64
    p_sigma::Float64
    epsilon::Float64
    gamma_Is::Float64
    gamma_Ia::Float64
    gamma_H_recovery::Float64
    gamma_H_death::Float64
    rho::Float64
    eta::Vector{Float64}    # 4 age groups
    hfr::Vector{Float64}    # 4 age groups
    ifr::Vector{Float64}    # 4 age groups
end

"""
    VaccinationScenario

Parameters for a single pre-pandemic vaccine investment scenario.
"""
struct VaccinationScenario
    start_time::Float64       # days from epidemic start
    rate::Float64             # fraction of population vaccinated per day
    uptake_limit::Float64     # maximum fraction willing to be vaccinated
    efficacy::Float64         # vaccine efficacy (fraction)
    waning_period::Float64    # days until immunity wanes
end

"""
    EconomicContacts

Per-sector social contact rates for economic sectors.
"""
struct EconomicContacts
    contacts_workplace::Vector{Float64}           # 45 sectors
    contacts_between_sectors::Matrix{Float64}     # 45×45
end

const _country_cache = Ref{Union{Nothing, Dict{String, CountryData}}}(nothing)
const _pathogen_cache = Ref{Union{Nothing, Dict{String, InfectionData}}}(nothing)
const _econ_cache = Ref{Union{Nothing, EconomicContacts}}(nothing)
const _sectors_cache = Ref{Union{Nothing, Vector{String}}}(nothing)
const _closure_cache = Ref{Union{Nothing, Dict{String, Vector{Float64}}}}(nothing)

# Map 5-year bin index (1-based, where bin k covers ages (k-1)*5 to k*5)
# to one of the 4 DAEDALUS age groups.
# Used for both the demography bins (1-21) and the CM bins (1-16).
function _bin_to_group(k::Int)::Int
    k == 1 && return 1   # 0-4
    1 < k <= 4 && return 2   # 5-19
    4 < k <= 13 && return 3   # 20-64
    return 4                          # 65+
end

function _load_economic_contacts()
    df = CSV.read(joinpath(DATA_DIR, "sectorcontacts.csv"), DataFrame)
    contacts_wp = Float64.(df[!, :n_cnt])
    n = length(contacts_wp)
    contacts_bw = zeros(Float64, n, n)  # off-diagonal zero per daedalus.data
    return EconomicContacts(contacts_wp, contacts_bw)
end

"""
    get_economic_contacts() -> EconomicContacts

Return within-sector and between-sector social contact rates for all 45
economic sectors. Data sourced from `sectorcontacts.csv`.
"""
function get_economic_contacts()::EconomicContacts
    if isnothing(_econ_cache[])
        _econ_cache[] = _load_economic_contacts()
    end
    return _econ_cache[]
end

function _load_sector_names()
    df = CSV.read(joinpath(DATA_DIR, "sectorcontacts.csv"), DataFrame)
    return String.(df[!, :sector])
end

"""
    get_sector_names() -> Vector{String}

Return the names of all 45 economic sectors.
"""
function get_sector_names()::Vector{String}
    if isnothing(_sectors_cache[])
        _sectors_cache[] = _load_sector_names()
    end
    return _sectors_cache[]
end

# Age-bin groups for sevenpathogens.csv: 17 five-year bins (0-4, ..., 80+)
# matching the IHR/IFR row indices ihr1-ihr17, ifr1-ifr17
function _aggregate_age_bins(vals::Vector{Float64})::Vector{Float64}
    # vals[1..17] correspond to 5-year bins: 0-4, 5-9, ..., 75-79, 80+
    # DAEDALUS groups: [0-4], [5-19 = bins 2-4], [20-64 = bins 5-13], [65+ = bins 14-17]
    g1 = vals[1]
    g2 = mean(vals[2:4])
    g3 = mean(vals[5:13])
    g4 = mean(vals[14:17])
    return [g1, g2, g3, g4]
end

function _load_pathogens()
    df = CSV.read(joinpath(DATA_DIR, "sevenpathogens.csv"), DataFrame)

    # Column 1 is "code label", columns 2-8 are pathogen values, col 9 is definition
    pathogen_names = String.(names(df)[2:8])
    code_col = df[!, Symbol("code label")]

    result = Dict{String, InfectionData}()

    for (pi, pname) in enumerate(pathogen_names)
        col = df[!, pi + 1]  # +1 because col 1 is code label

        # Helper: get scalar value for a given code
        getval(code) = Float64(col[findfirst(==(code), code_col)])

        # IHR and IFR: 17 five-year bins
        ihr = [Float64(col[findfirst(==(string("ihr", i)), code_col)]) for i in 1:17]
        ifr = [Float64(col[findfirst(==(string("ifr", i)), code_col)]) for i in 1:17]

        ihr_4 = _aggregate_age_bins(ihr)
        ifr_4 = _aggregate_age_bins(ifr)

        # Scalar parameters
        ps = getval("ps")
        Tlat = getval("Tlat")
        Tay = getval("Tay")
        Tsr = getval("Tsr")
        Tsh = getval("Tsh")
        Threc = getval("Threc")
        Thd = getval("Thd")
        Ti = getval("Ti")
        red = getval("red")
        r0 = getval("R0")

        # Derived rates
        sigma = 1.0 / Tlat
        gamma_Ia = 1.0 / Tay
        gamma_Is = 1.0 / Tsr
        gamma_H_rec = 1.0 / Threc
        gamma_H_death = 1.0 / Thd
        rho = 1.0 / Ti
        epsilon = red

        # eta[age] = (IHR[age] / p_sigma) / Tsh
        eta = (ihr_4 ./ ps) ./ Tsh

        # hfr[age] = IFR[age] / IHR[age]  (clamp to avoid div-by-zero)
        hfr = ifr_4 ./ max.(ihr_4, 1e-12)

        result[lowercase(pname)] = InfectionData(
            r0, sigma, ps, epsilon,
            gamma_Is, gamma_Ia, gamma_H_rec, gamma_H_death, rho,
            eta, hfr, ifr_4
        )
    end

    return result
end

"""
    get_pathogen(name::String) -> InfectionData

Return infection parameters for the named pathogen. Call [`list_pathogens`](@ref)
to see available names.
"""
function get_pathogen(name::String)::InfectionData
    if isnothing(_pathogen_cache[])
        _pathogen_cache[] = _load_pathogens()
    end
    d = _pathogen_cache[]
    name_lower = lowercase(name)
    haskey(d, name_lower) ||
        error("Pathogen not found: $name. Available: " * join(keys(d), ", "))
    return d[name_lower]
end

"""
    list_pathogens() -> Vector{String}

Return the names of all available pathogens.
"""
function list_pathogens()::Vector{String}
    if isnothing(_pathogen_cache[])
        _pathogen_cache[] = _load_pathogens()
    end
    return collect(keys(_pathogen_cache[]))
end

const VACCINATION_SCENARIOS = Dict{String, VaccinationScenario}(
    "none" => VaccinationScenario(365.0, 1.0 / 7 / 100, 0.40, 0.50, 270.0),
    "low" => VaccinationScenario(300.0, 2.0 / 7 / 100, 0.50, 0.50, 270.0),
    "medium" => VaccinationScenario(200.0, 3.0 / 7 / 100, 0.60, 0.50, 270.0),
    "high" => VaccinationScenario(100.0, 3.5 / 7 / 100, 0.80, 0.50, 270.0)
)

"""
    get_vaccination_scenario(name::String) -> VaccinationScenario

Return parameters for the named vaccination scenario. Call
[`list_vaccination_scenarios`](@ref) to see available names.
"""
function get_vaccination_scenario(name::String)::VaccinationScenario
    haskey(VACCINATION_SCENARIOS, name) ||
        error("Scenario not found: $name. Available: " *
              join(keys(VACCINATION_SCENARIOS), ", "))
    return VACCINATION_SCENARIOS[name]
end

"""
    list_vaccination_scenarios() -> Vector{String}

Return names of all vaccination scenarios: `"none"`, `"low"`, `"medium"`, `"high"`.
"""
function list_vaccination_scenarios()::Vector{String}
    return collect(keys(VACCINATION_SCENARIOS))
end

function _load_closure_strategies()
    df = CSV.read(joinpath(DATA_DIR, "economic_closures.csv"), DataFrame)
    n_sectors = nrow(df)

    # Column mapping matching the R closure_data.R script:
    #   "none"              -> all 1.0 (hardcoded)
    #   "elimination"       -> "Elimination" column      (heavy level)
    #   "economic_closures" -> "Economic Closures" column (light level)
    #   "school_closures"   -> "School Closures" column   (light level)
    result = Dict{String, Vector{Float64}}(
        "none" => ones(Float64, n_sectors),
        "elimination" => Float64.(df[!, "Elimination"]),
        "economic_closures" => Float64.(df[!, "Economic Closures"]),
        "school_closures" => Float64.(df[!, "School Closures"])
    )
    return result
end

"""
    get_closure_strategy(name::String) -> Vector{Float64}

Return a 45-element vector of sector openness coefficients (0–1) for the named
pandemic response strategy. Call [`list_closure_strategies`](@ref) for available
names.
"""
function get_closure_strategy(name::String)::Vector{Float64}
    if isnothing(_closure_cache[])
        _closure_cache[] = _load_closure_strategies()
    end
    d = _closure_cache[]
    haskey(d, name) || error("Strategy not found: $name. Available: " * join(keys(d), ", "))
    return d[name]
end

"""
    list_closure_strategies() -> Vector{String}

Return names of available pandemic response strategies.
"""
function list_closure_strategies()::Vector{String}
    if isnothing(_closure_cache[])
        _closure_cache[] = _load_closure_strategies()
    end
    return collect(keys(_closure_cache[]))
end

# Letter position in a-p (a=1, ..., p=16)
const _LETTER_POS = Dict(c => i for (i, c) in enumerate('a':'p'))

"""
Aggregate a 16-element vector of 5-year-bin populations into 4 DAEDALUS groups.
Bins 1-16 correspond to ages 0-4, 5-9, ..., 75-79. Bins 14-16 map to 65+.
"""
function _demog_bins_to_4(pop16::Vector{Float64})::Vector{Float64}
    return [
        pop16[1],
        sum(pop16[2:4]),
        sum(pop16[5:13]),
        sum(pop16[14:16])
    ]
end

"""
Compute the 4×4 population-weighted contact matrix from 16-bin CM values.

Replicates the weighted aggregation in the R `data-raw/country_data.R` script:
  contacts_4x4[i, j] = sum_{to in group_i, from in group_j}(CM[to,from] * pop[to])
                        / sum_{to in group_i} pop[to]
"""
function _aggregate_cm(cm16::Matrix{Float64}, pop16::Vector{Float64})::Matrix{Float64}
    # Group assignments for bins 1-16
    groups = [_bin_to_group(k) for k in 1:16]
    group_pops = [sum(pop16[groups .== g]) for g in 1:4]
    group_pops = max.(group_pops, 1.0)  # guard against zero populations

    cm4 = zeros(Float64, 4, 4)
    for to_bin in 1:16, from_bin in 1:16

        gi = groups[to_bin]
        gj = groups[from_bin]
        cm4[gi, gj] += cm16[to_bin, from_bin] * pop16[to_bin]
    end
    for i in 1:4
        cm4[i, :] ./= group_pops[i]
    end
    return cm4
end

# Load and aggregate life expectancy from CSV: average over sex, aggregate to 4 age groups
function _load_life_expectancy()::Dict{String, Vector{Float64}}
    df = CSV.read(joinpath(DATA_DIR, "life_expectancy.csv"), DataFrame)

    # Extract relevant columns
    countries = String.(df[!, Symbol("Location")])
    sexes = String.(df[!, Symbol("Dim1")])
    age_groups = String.(df[!, Symbol("Dim2")])
    values = Float64.(df[!, Symbol("Value")])

    # Group by country and age group, average over sex
    result_dict = Dict{String, Dict{String, Float64}}()
    for (i, row) in enumerate(eachrow(df))
        country = String(row.Location)
        age_group = String(row.Dim2)
        value = Float64(row.Value)

        if !haskey(result_dict, country)
            result_dict[country] = Dict{String, Float64}()
        end

        if !haskey(result_dict[country], age_group)
            result_dict[country][age_group] = 0.0
        end
        result_dict[country][age_group] += value
    end

    # Average over sex and aggregate to 4 age groups
    agg_result = Dict{String, Vector{Float64}}()
    for (country, age_dict) in result_dict
        # Count how many entries per age group (should be 2: male, female)
        age_counts = Dict(ag => 0 for ag in keys(age_dict))
        for (country_inner, ad) in result_dict
            if country_inner == country
                for ag in keys(ad)
                    age_counts[ag] = get(age_counts, ag, 0) + 1
                end
            end
        end

        # Recompute: average and aggregate to 4 groups
        le_by_age = Dict{String, Vector{Float64}}()
        # Aggregate: 0-4, 5-19, 20-64, 65+
        group1 = Float64[]  # 0-4
        group2 = Float64[]  # 5-19
        group3 = Float64[]  # 20-64
        group4 = Float64[]  # 65+

        for (age_group, val_sum) in age_dict
            # Determine which group this age belongs to
            if lowercase(age_group) in ["<1 year", "0-4 years", "1-4 years"]
                push!(group1, val_sum / 2)  # average over 2 sexes
            elseif occursin(r"^[5-9]", age_group) || occursin(r"^[1][0-9]", age_group) || occursin(r"^20", age_group)
                # 5-9, 10-14, 15-19, 20 (if exists)
                if occursin(r"^20", age_group)
                    # This might be 20-24
                    push!(group3, val_sum / 2)
                else
                    # 5-19
                    push!(group2, val_sum / 2)
                end
            elseif occursin(r"^[2-6][0-9]|^[56][0-9]", age_group) && !occursin(r"^65", age_group)
                push!(group3, val_sum / 2)  # 20-64
            elseif occursin(r"^65", age_group) || occursin(r"^7[0-9]", age_group) || occursin(r"^8[0-9]", age_group) || lowercase(age_group) == "85+ years"
                push!(group4, val_sum / 2)  # 65+
            end
        end

        le_4 = [
            isempty(group1) ? 70.0 : mean(group1),
            isempty(group2) ? 50.0 : mean(group2),
            isempty(group3) ? 45.0 : mean(group3),
            isempty(group4) ? 20.0 : mean(group4)
        ]
        agg_result[country] = le_4
    end

    return agg_result
end

function _load_countries()
    country_df = CSV.read(joinpath(DATA_DIR, "country_data.csv"), DataFrame)
    hosp_df = CSV.read(joinpath(DATA_DIR, "hospital_capacity.csv"), DataFrame)
    # sector_gva_data.csv has a junk first row; row 2 is the real header
    gva_df = CSV.read(
        joinpath(DATA_DIR, "sector_gva_data.csv"), DataFrame;
        header = 2
    )
    # Ensure the country column is named "Country" (strip any leading whitespace)
    rename!(gva_df, first(names(gva_df)) => "Country")

    # Load life expectancy
    le_lookup = _load_life_expectancy()

    # Build hospital capacity lookup: country -> spare_capacity
    hosp_lookup = Dict{String, Float64}(
        String(r.country) => Float64(r.spare_capacity) for r in eachrow(hosp_df)
    )

    # Build GVA lookup: country -> Vector{Float64}(45)
    # Columns 2:46 are the 45 sector GVA values; some cells may be non-numeric
    _parse_gva(x::Number) = abs(Float64(x))
    _parse_gva(x::AbstractString) = abs(something(tryparse(Float64, strip(x)), 0.0))
    _parse_gva(::Missing) = 0.0
    _parse_gva(x) = 0.0

    gva_lookup = Dict{String, Vector{Float64}}()
    for r in eachrow(gva_df)
        cname = String(r.Country)
        vals = [_parse_gva(r[i]) for i in 2:ncol(gva_df)]
        length(vals) == 45 && (gva_lookup[cname] = vals)
    end

    # Identify CM and population columns
    all_cols = names(country_df)
    cm_cols = filter(c -> startswith(c, "CM"), all_cols)
    npop_cols = ["Npop$i" for i in 1:21]
    nns_cols = ["NNs$i" for i in 1:45]

    result = Dict{String, CountryData}()

    for row in eachrow(country_df)
        cname = String(row.country)

        # Skip countries without hospital, GVA, or life expectancy data
        (haskey(hosp_lookup, cname) && haskey(gva_lookup, cname) && haskey(le_lookup, cname)) || continue

        # --- Demography: 21 five-year bins -> 4 DAEDALUS groups ---
        pop21 = Float64[row[Symbol(c)] for c in npop_cols]
        # Fill any missing
        replace!(pop21, NaN => 0.0)

        demog4 = [
            pop21[1],
            sum(pop21[2:4]),
            sum(pop21[5:13]),
            sum(pop21[14:21])
        ]

        # --- Workers: 45 sectors ---
        workers45 = Int[max(0, round(Int, coalesce(row[Symbol(c)], 0))) for c in nns_cols]

        # --- Contact matrix: 16×16 -> 4×4 ---
        # Parse each CMxy column: x = to letter (row), y = from letter (col)
        pop16 = Float64[row[Symbol("Npop$k")] for k in 1:16]
        replace!(pop16, NaN => 0.0)

        cm16 = zeros(Float64, 16, 16)
        for colname in cm_cols
            letters = colname[3:4]   # e.g. "CMab" -> "ab"
            to_pos = get(_LETTER_POS, letters[1], 0)
            frm_pos = get(_LETTER_POS, letters[2], 0)
            (to_pos == 0 || frm_pos == 0) && continue
            v = row[Symbol(colname)]
            cm16[to_pos, frm_pos] = isnan(v) ? 0.0 : Float64(v)
        end

        cm4 = _aggregate_cm(cm16, pop16)

        # Replace any NaN/Inf cells with 1.0 (some countries have missing data)
        replace!(cm4, NaN => 1.0)
        replace!(cm4, Inf => 1.0)

        # --- Hospital capacity ---
        spare = hosp_lookup[cname]
        total_pop = sum(demog4)
        hosp_cap = round(spare * total_pop / 1000.0)

        # --- GVA ---
        gva = gva_lookup[cname]

        # --- GNI and life expectancy ---
        gni = Float64(coalesce(row.gnipc, 0.0))
        life_exp = le_lookup[cname]
        vsl = life_exp .* gni

        result[cname] = CountryData(demog4, cm4, workers45, gva, hosp_cap, gni, life_exp, vsl)
    end

    return result
end

"""
    get_country(name::String) -> CountryData

Return demographic and economic data for the named country. Call
[`list_countries`](@ref) to see available countries.
"""
function get_country(name::String)::CountryData
    if isnothing(_country_cache[])
        _country_cache[] = _load_countries()
    end
    d = _country_cache[]
    haskey(d, name) ||
        error("Country \"$name\" not found. Call list_countries() to see available names.")
    return d[name]
end

"""
    list_countries() -> Vector{String}

Return a sorted list of all available country names.
"""
function list_countries()::Vector{String}
    if isnothing(_country_cache[])
        _country_cache[] = _load_countries()
    end
    return sort(collect(keys(_country_cache[])))
end

end  # module DataLoader
