"""
Tests for multi-parameter runs with vector r0 and multi-threading support.
Tests basic functionality and structure without checking numerical correctness.
"""

using Daedalus
using Test

@testset "Multi-run r0 parameter support" begin
    @testset "Scalar r0 returns single result" begin
        # Scalar r0 should return a single NamedTuple, not a vector
        result = Daedalus.daedalus("Canada", 1.5, time_end = 10.0)

        @test isa(result, NamedTuple)
        @test haskey(result, :sol)
        @test haskey(result, :saves)
        @test haskey(result, :npi)
        @test length(result) == 3  # sol, saves, npi only
    end

    @testset "Vector r0 returns vector of results" begin
        # Vector r0 should return a Vector of NamedTuples
        r0_vals = [1.0, 1.5, 2.0]
        results = Daedalus.daedalus("Canada", r0_vals, time_end = 10.0)

        @test isa(results, Vector)
        @test length(results) == 3

        # Each result should be a NamedTuple with sol, saves, npi, r0 fields
        for (i, result) in enumerate(results)
            @test isa(result, NamedTuple)
            @test haskey(result, :sol)
            @test haskey(result, :saves)
            @test haskey(result, :npi)
            @test haskey(result, :r0)
            @test result.r0 == r0_vals[i]
        end
    end

    @testset "Single-element r0 vector matches scalar" begin
        # A vector with one element should behave like scalar
        r0_scalar = 1.5
        r0_vector = [1.5]

        result_scalar = Daedalus.daedalus(country = "Canada", r0 = r0_scalar, time_end = 10.0)
        results_vector = Daedalus.daedalus(country = "Canada", r0 = r0_vector, time_end = 10.0)

        # Vector should return a vector of length 1
        @test isa(results_vector, Vector)
        @test length(results_vector) == 1

        # The single result should be a NamedTuple with the r0 parameter field
        @test isa(results_vector[1], NamedTuple)
        @test results_vector[1].r0 == r0_scalar
    end

    @testset "Different r0 values produce different results" begin
        # Different r0 values should produce different ODE solutions
        r0_vals = [1.0, 3.0]
        results = Daedalus.daedalus(country = "Canada", r0 = r0_vals, time_end = 50.0)

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
        r0_vals = [1.0, 2.0, 3.0]
        time_end = 30.0
        increment = 1.0

        results = Daedalus.daedalus(
            country = "Canada", r0 = r0_vals, time_end = time_end, increment = increment
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
        r0_vals = [1.0, 1.5, 2.0]
        results = Daedalus.daedalus(
            "Canada", r0_vals, time_end = 10.0, n_threads = 1
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
        r0_vals = [1.0, 1.5, 2.0]
        n_threads_use = min(4, Threads.nthreads())  # Use up to 4 threads, but not more than available

        results = Daedalus.daedalus(
            "Canada", r0_vals, time_end = 10.0, n_threads = n_threads_use
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
        r0_vals = [1.0, 1.5, 2.0]

        results_serial = Daedalus.daedalus(
            "Canada", r0_vals, time_end = 10.0, n_threads = 1
        )

        n_threads_use = min(2, Threads.nthreads())
        if n_threads_use > 1
            results_threaded = Daedalus.daedalus(
                "Canada", r0_vals, time_end = 10.0, n_threads = n_threads_use
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
        # When using scalar r0, n_threads parameter should be ignored (no effect)
        result = Daedalus.daedalus("Canada", 1.5, time_end = 10.0, n_threads = 4)

        # Should still return single result, not vector
        @test isa(result, NamedTuple)
        @test haskey(result, :sol)
    end
end

@testset "Edge cases" begin
    @testset "Empty r0 vector raises error" begin
        # Empty vector should raise an error
        @test_throws BoundsError Daedalus.daedalus(
            "Canada", Float64[], time_end = 10.0
        )
    end

    @testset "Large r0 vector" begin
        # Should handle moderate-sized parameter sweeps
        r0_vals = collect(1.0:0.5:10.0)  # 19 values
        results = Daedalus.daedalus("Canada", r0_vals, time_end = 5.0)

        @test length(results) == length(r0_vals)

        for (i, result) in enumerate(results)
            @test result.r0 ≈ r0_vals[i]
            @test haskey(result, :sol)
        end
    end

    @testset "Very small time_end" begin
        # Should handle very short simulation times
        r0_vals = [1.0, 2.0]
        results = Daedalus.daedalus("Canada", r0_vals, time_end = 1.0)

        @test length(results) == 2
        for result in results
            @test isa(result, NamedTuple)
            @test haskey(result, :sol)
        end
    end

    @testset "Different countries with vector r0" begin
        # Should work with different country data
        countries = ["Australia", "United Kingdom"]
        r0_vals = [1.0, 2.0, 3.0]

        for country_name in countries
            results = Daedalus.daedalus(country_name, r0_vals, time_end = 5.0)
            @test length(results) == 3
        end
    end
end

@testset "Parameter preservation" begin
    @testset "r0 values are correctly stored in results" begin
        r0_vals = [1.2, 2.3, 3.4]
        results = Daedalus.daedalus("Canada", r0_vals, time_end = 10.0)

        for (i, result) in enumerate(results)
            @test result.r0 == r0_vals[i]
        end
    end

    @testset "Other parameters are shared across runs" begin
        r0_vals = [1.0, 2.0]
        results = Daedalus.daedalus(
            "Canada",
            r0_vals,
            time_end = 10.0,
            npi = nothing,
            log_rt = true
        )

        # All runs should share the same npi
        @test results[1].npi === results[2].npi
        @test isnothing(results[1].npi)
    end
end

@testset "Solution structure" begin
    @testset "Each solution has valid ODE output" begin
        r0_vals = [1.0, 2.0]
        results = Daedalus.daedalus("Canada", r0_vals, time_end = 15.0)

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
        # State should have N_COMPARTMENTS * N_TOTAL_GROUPS * N_VACCINE_STRATA + 1 (for Rt)
        r0_vals = [1.0]
        results = Daedalus.daedalus("Canada", r0_vals, time_end = 10.0)

        state = results[1].sol.u[1]
        expected_dim = Daedalus.Constants.N_COMPARTMENTS *
                       Daedalus.Constants.N_TOTAL_GROUPS *
                       Daedalus.Constants.N_VACCINE_STRATA + 1

        @test length(state) == expected_dim
    end
end

@testset "Backward compatibility" begin
    @testset "Scalar r0 behavior unchanged" begin
        # Original scalar interface should work exactly as before
        result = Daedalus.daedalus(
            "Australia",
            2.0,
            sigma = 0.217,
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
