"""
    DaedalusOutput

Container for results from a `daedalus()` simulation, including the ODE solution,
input parameters, and interventions.

# Fields
- `sol`: ODE solution (`ODESolution` from OrdinaryDiffEq)
- `saves`: NPI event tracking data (`Union{Nothing, Vector{SavedValues}}`)
- `npi`: Intervention object (`Union{Nothing, Npi}`)
- `country`: Input country data (`CountryData`)
- `infection`: Input infection parameters (`InfectionData`)
- `time_end`: Simulation duration in days (`Float64`)
"""
struct DaedalusOutput
    sol                # ODESolution
    saves              # Union{Nothing, Vector{SavedValues}}
    npi                # Union{Nothing, Npi}
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

function Base.show(io::IO, output::DaedalusOutput)
    country_name = get(output.country, :name, "Unknown")
    r0 = output.infection.r0
    npi_status = isnothing(output.npi) ? "none" : "active"
    duration = output.time_end

    print(io, "DaedalusOutput\n")
    print(io, "  Country: $country_name\n")
    print(io, "  R₀: $(round(r0, digits=2))\n")
    print(io, "  NPI: $npi_status\n")
    print(io, "  Duration: $(Int(duration)) days\n")
end
