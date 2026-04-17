"""
Tests for multi-parameter runs with vector r0 and multi-threading support.
Tests basic functionality and structure without checking numerical correctness.
"""

using Daedalus
using Test

# Helper to create infection with custom r0
function make_infection(pathogen_name::String, r0_value::Float64)
    inf = Daedalus.DataLoader.get_pathogen(pathogen_name)
    inf.r0 = r0_value
    return inf
end

@testset "Multi-run infection parameter support" begin
    @testset "Scalar infection returns single result" begin
        # Scalar infection should return a single NamedTuple, not a vector
        result = Daedalus.daedalus("Canada", make_infection("sars-cov-2 delta", 1.5), time_end = 10.0)

        @test isa(result, NamedTuple)
        @test haskey(result, :sol)
        @test haskey(result, :saves)
        @test haskey(result, :npi)
        @test length(result) == 3  # sol, saves, npi only
    end

    @testset "Vector infection returns vector of results" begin
        # Vector infection should return a Vector of NamedTuples
        infections = [
            make_infection("sars-cov-2 delta", 1.0),
            make_infection("sars-cov-2 delta", 1.5),
            make_infection("sars-cov-2 delta", 2.0)
        ]
        results = Daedalus.daedalus("Canada", infections, time_end = 10.0)

        @test isa(results, Vector)
        @test length(results) == 3

        # Each result should be a NamedTuple with sol, saves, npi, r0 fields
        for (i, result) in enumerate(results)
            @test isa(result, NamedTuple)
            @test haskey(result, :sol)
            @test haskey(result, :saves)
            @test haskey(result, :npi)
            @test haskey(result, :r0)
            @test result.r0 == infections[i].r0
        end
    end

    @testset "Single-element infection vector matches scalar" begin
        # A vector with one element should behave like vector
        infection_scalar = make_infection("sars-cov-2 delta", 1.5)

        infections_vector = [make_infection("sars-cov-2 delta", 1.5)]

        result_scalar = Daedalus.daedalus("Canada", infection_scalar, time_end = 10.0)
        results_vector = Daedalus.daedalus("Canada", infections_vector, time_end = 10.0)

        # Vector should return a vector of length 1
        @test isa(results_vector, Vector)
        @test length(results_vector) == 1

        # The single result should be a NamedTuple with the r0 parameter field
        @test isa(results_vector[1], NamedTuple)
        @test results_vector[1].r0 == infection_scalar.r0
    end

    @testset "Different r0 values produce different results" begin
        # Different r0 values should produce different ODE solutions
        infections = [
            make_infection("sars-cov-2 delta", 1.0),
            make_infection("sars-cov-2 delta", 3.0)
        ]
        results = Daedalus.daedalus("Canada", infections, time_end = 50.0)

        @test length(results) == 2

        # Extract final state from each solution
        final_1 = results[1].sol.u[end]
        final_2 = results[2].sol.u[end]

        # Solutions should differ (higher r0 should produce different dynamics)
        @test final_1 != final_2
        @test length(final_1) == length(final_2)  # Same state dimension
    end

    @testset "Time dimension consistency" begin
        # All runs should have same number of timepoints
        infections = [
            make_infection("sars-cov-2 delta", 1.0),
            make_infection("sars-cov-2 delta", 2.0),
            make_infection("sars-cov-2 delta", 3.0)
        ]
        time_end = 30.0
        increment = 1.0

        results = Daedalus.daedalus(
            "Canada", infections, time_end = time_end, increment = increment
        )

        expected_timepoints = Int(time_end / increment) + 1

        for result in results
            # Check number of saved timepoints
            @test length(result.sol.t) == expected_timepoints
        end
    end
end

@testset "Multi-threading support" begin
    @testset "Serial execution (n_threads=1)" begin
        # With n_threads=1, should execute serially without error
        infections = [
            make_infection("sars-cov-2 delta", 1.0),
            make_infection("sars-cov-2 delta", 1.5),
            make_infection("sars-cov-2 delta", 2.0)
        ]
        results = Daedalus.daedalus(
            "Canada", infections, time_end = 10.0, n_threads = 1
        )

        @test isa(results, Vector)
        @test length(results) == 3

        for result in results
            @test isa(result, NamedTuple)
            @test haskey(result, :sol)
        end
    end

    @testset "Multi-threaded execution (n_threads>1)" begin
        # With n_threads>1, should execute with threads without error
        infections = [
            make_infection("sars-cov-2 delta", 1.0),
            make_infection("sars-cov-2 delta", 1.5),
            make_infection("sars-cov-2 delta", 2.0)
        ]
        n_threads_use = min(4, Threads.nthreads())  # Use up to 4 threads, but not more than available

        results = Daedalus.daedalus(
            "Canada", infections, time_end = 10.0, n_threads = n_threads_use
        )

        @test isa(results, Vector)
        @test length(results) == 3

        for result in results
            @test isa(result, NamedTuple)
            @test haskey(result, :sol)
        end
    end

    @testset "Serial and multi-threaded produce consistent results" begin
        # Results should be identical regardless of threading
        infections_serial = [
            make_infection("sars-cov-2 delta", 1.0),
            make_infection("sars-cov-2 delta", 1.5),
            make_infection("sars-cov-2 delta", 2.0)
        ]

        results_serial = Daedalus.daedalus(
            "Canada", infections_serial, time_end = 10.0, n_threads = 1
        )

        n_threads_use = min(2, Threads.nthreads())
        if n_threads_use > 1
            infections_threaded = [
                make_infection("sars-cov-2 delta", 1.0),
                make_infection("sars-cov-2 delta", 1.5),
                make_infection("sars-cov-2 delta", 2.0)
            ]

            results_threaded = Daedalus.daedalus(
                "Canada", infections_threaded, time_end = 10.0, n_threads = n_threads_use
            )

            @test length(results_serial) == length(results_threaded)

            # Check that r0 values match and are in same order
            for (rs, rt) in zip(results_serial, results_threaded)
                @test rs.r0 == rt.r0
                # Solutions should be very similar (allowing for floating point variations)
                @test length(rs.sol.u) == length(rt.sol.u)
            end
        else
            @test_skip "Multi-threading test skipped (only 1 thread available)"
        end
    end

    @testset "n_threads=1 default for scalar" begin
        # When using scalar infection, n_threads parameter should be ignored (no effect)
        result = Daedalus.daedalus(
            "Canada", make_infection("sars-cov-2 delta", 1.5), time_end = 10.0, n_threads = 4)

        # Should still return single result, not vector
        @test isa(result, NamedTuple)
        @test haskey(result, :sol)
    end
end

@testset "Edge cases" begin
    @testset "Empty infection vector raises error" begin
        # Empty vector should raise an error
        @test_throws BoundsError Daedalus.daedalus(
            "Canada", Daedalus.DataLoader.InfectionData[], time_end = 10.0
        )
    end

    @testset "Large infection vector" begin
        # Should handle moderate-sized parameter sweeps
        r0_vals = collect(1.0:0.5:10.0)  # 19 values
        infections = [make_infection("sars-cov-2 delta", r0) for r0 in r0_vals]
        results = Daedalus.daedalus("Canada", infections, time_end = 5.0)

        @test length(results) == length(r0_vals)

        for (i, result) in enumerate(results)
            @test result.r0 ≈ r0_vals[i]
            @test haskey(result, :sol)
        end
    end

    @testset "Very small time_end" begin
        # Should handle very short simulation times
        infections = [
            make_infection("sars-cov-2 delta", 1.0),
            make_infection("sars-cov-2 delta", 2.0)
        ]
        results = Daedalus.daedalus("Canada", infections, time_end = 1.0)

        @test length(results) == 2
        for result in results
            @test isa(result, NamedTuple)
            @test haskey(result, :sol)
        end
    end

    @testset "Different countries with vector infection" begin
        # Should work with different country data
        countries = ["Australia", "United Kingdom"]

        for country_name in countries
            infections = [
                make_infection("sars-cov-2 delta", 1.0),
                make_infection("sars-cov-2 delta", 2.0),
                make_infection("sars-cov-2 delta", 3.0)
            ]
            results = Daedalus.daedalus(country_name, infections, time_end = 5.0)
            @test length(results) == 3
        end
    end
end

@testset "Parameter preservation" begin
    @testset "r0 values are correctly stored in results" begin
        r0_vals = [1.2, 2.3, 3.4]
        infections = [make_infection("sars-cov-2 delta", r0) for r0 in r0_vals]
        results = Daedalus.daedalus("Canada", infections, time_end = 10.0)

        for (i, result) in enumerate(results)
            @test result.r0 == r0_vals[i]
        end
    end

    @testset "Other parameters are shared across runs" begin
        infections = [
            make_infection("sars-cov-2 delta", 1.0),
            make_infection("sars-cov-2 delta", 2.0)
        ]
        results = Daedalus.daedalus(
            "Canada",
            infections,
            time_end = 10.0,
            npi = nothing,
            log_rt = true
        )

        # All runs should have the same events
        @test results[1].events === results[2].events
        @test isempty(results[1].events)
    end
end

@testset "Solution structure" begin
    @testset "Each solution has valid ODE output" begin
        infections = [
            make_infection("sars-cov-2 delta", 1.0),
            make_infection("sars-cov-2 delta", 2.0)
        ]
        results = Daedalus.daedalus("Canada", infections, time_end = 15.0)

        for result in results
            sol = result.sol

            # Check solution structure
            @test haskey(sol, :t)  # Time vector
            @test haskey(sol, :u)  # State vector
            @test !isempty(sol.t)
            @test !isempty(sol.u)
            @test length(sol.t) == length(sol.u)

            # State should be a vector
            @test isa(sol.u[1], AbstractVector)
        end
    end

    @testset "Solutions have correct state dimension" begin
        # State: N_COMPARTMENTS * N_TOTAL_GROUPS * N_VACCINE_STRATA + N_TOTAL_GROUPS + 1 (Rt)
        infections = [make_infection("sars-cov-2 delta", 1.0)]
        results = Daedalus.daedalus("Canada", infections, time_end = 10.0)

        state = results[1].sol.u[1]
        expected_dim = Daedalus.Constants.N_COMPARTMENTS *
                       Daedalus.Constants.N_TOTAL_GROUPS *
                       Daedalus.Constants.N_VACCINE_STRATA +
                       Daedalus.Constants.N_TOTAL_GROUPS + 1

        @test length(state) == expected_dim
    end
end

@testset "New interface with InfectionData" begin
    @testset "String pathogen name works" begin
        # String pathogen should fetch and work
        result = Daedalus.daedalus(
            "Australia",
            "sars-cov-2 delta",
            time_end = 20.0,
            log_rt = true
        )

        @test isa(result, NamedTuple)
        @test haskey(result, :sol)
        @test haskey(result, :saves)
        @test haskey(result, :npi)
        @test length(result) == 3  # Only sol, saves, npi (no r0 field)
    end

    @testset "InfectionData object works" begin
        # InfectionData object should work directly
        infection = Daedalus.DataLoader.get_pathogen("sars-cov-2 delta")
        infection.r0 = 2.0
        result = Daedalus.daedalus(
            "Australia",
            infection,
            time_end = 20.0,
            log_rt = true
        )

        @test isa(result, NamedTuple)
        @test haskey(result, :sol)
        @test haskey(result, :saves)
        @test haskey(result, :npi)
        @test length(result) == 3  # Only sol, saves, npi (no r0 field)
    end
end
