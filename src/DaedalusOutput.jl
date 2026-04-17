"""
    DaedalusOutput

Container for results from a `daedalus()` simulation, including the ODE solution,
input parameters, and events.

# Fields
- `sol`: ODE solution (`ODESolution` from OrdinaryDiffEq)
- `saves`: Event tracking data (`Union{Nothing, Vector{SavedValues}}`)
- `events`: Event objects (`Vector{Event}` — can contain `Npi` and/or `Vaccination`)
- `country`: Input country data (`CountryData`)
- `infection`: Input infection parameters (`InfectionData`)
- `time_end`: Simulation duration in days (`Float64`)
"""
struct DaedalusOutput
    sol                # ODESolution
    saves              # Union{Nothing, Vector{SavedValues}}
    events             # Vector{Event}
    country            # CountryData
    infection          # InfectionData
    time_end::Float64
end

"""
    get_data(output::DaedalusOutput, field::Symbol)

Access a field of the output struct by name.

# Arguments
- `output`: A `DaedalusOutput` object
- `field`: Field name as a symbol (`:country`, `:infection`, `:sol`, etc.)

# Returns
The requested field value.
"""
function get_data(output::DaedalusOutput, field::Symbol)
    return getfield(output, field)
end

"""
    get_npi(output::DaedalusOutput)

Extract all NPI events from the output.

# Arguments
- `output`: A `DaedalusOutput` object

# Returns
A `Vector{Npi}` containing all NPI events (empty if none)
"""
function get_npi(output::DaedalusOutput)
    return filter(e -> isa(e, Npi), output.events)
end

"""
    get_vaccination(output::DaedalusOutput)

Extract all vaccination events from the output.

# Arguments
- `output`: A `DaedalusOutput` object

# Returns
A `Vector{Vaccination}` containing all vaccination events (empty if none)
"""
function get_vaccination(output::DaedalusOutput)
    return filter(e -> isa(e, Vaccination), output.events)
end

function Base.show(io::IO, output::DaedalusOutput)
    country_name = hasfield(typeof(output.country), :name) ? output.country.name : "Unknown"
    r0 = output.infection.r0
    npi_count = length(get_npi(output))
    vax_count = length(get_vaccination(output))
    duration = output.time_end

    npi_status = npi_count > 0 ? "active ($npi_count)" : "none"
    vax_status = vax_count > 0 ? "active ($vax_count)" : "none"

    print(io, "DaedalusOutput\n")
    print(io, "  Country: $country_name\n")
    print(io, "  R₀: $(round(r0, digits=2))\n")
    print(io, "  NPI: $npi_status\n")
    print(io, "  Vaccination: $vax_status\n")
    print(io, "  Duration: $(Int(duration)) days\n")
end
