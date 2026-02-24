using Daedalus
using LinearAlgebra
using Random
using Test

@testset "Daedalus.jl" begin
    # Basic model and helper tests
    include("test_basic.jl")

    # Include eigenvalue tests
    include("test_eigenvalue.jl")

    # Include time-limited NPI tests
    include("test_timed_npi.jl")

    # Include output processing tests
    include("test_outputs.jl")
end
