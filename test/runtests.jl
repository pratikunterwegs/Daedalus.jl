using Daedalus
using LinearAlgebra
using Random
using Test

@testset "Daedalus.jl" begin
    # Basic model and helper tests
    include("test_basic.jl")

    # Include eigenvalue tests
    include("test_eigenvalue.jl")

    # Include trigger and effect tests
    include("test_triggers_effects.jl")

    # Include time-limited NPI tests
    include("test_timed_npi.jl")

    # Include NPI helper tests
    include("test_npi_helpers.jl")

    # Include vaccination event tests
    include("test_vaccination_events.jl")

    # Include output processing tests
    include("test_outputs.jl")

    include("test_countries.jl")

    include("test_weighted_slice_sum.jl")

    include("test_settings.jl")
end
