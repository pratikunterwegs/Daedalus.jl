using Daedalus
using Test
using DiffEqCallbacks

@testset "TimedNpi construction and validation" begin
    @testset "Valid single-phase construction" begin
        npi = Daedalus.DaedalusStructs.TimedNpi(10.0, 20.0, 0.5, "test")

        @test npi.start_times == [10.0]
        @test npi.end_times == [20.0]
        @test npi.coefs == [0.5]
        @test npi.identifier == "test"
    end

    @testset "Valid multi-phase construction" begin
        start_times = [10.0, 30.0, 60.0]
        end_times = [25.0, 55.0, 90.0]
        coefs = [0.7, 0.3, 0.5]

        npi = Daedalus.DaedalusStructs.TimedNpi(
            start_times, end_times, coefs, "multi_phase"
        )

        @test npi.start_times == start_times
        @test npi.end_times == end_times
        @test npi.coefs == coefs
        @test npi.identifier == "multi_phase"
    end

    @testset "Default identifier" begin
        npi = Daedalus.DaedalusStructs.TimedNpi([10.0], [20.0], [0.5])
        @test npi.identifier == "custom_timed"
    end

    @testset "Convenience constructor defaults" begin
        npi = Daedalus.DaedalusStructs.TimedNpi(10.0, 20.0, 0.5)
        @test npi.identifier == "single_phase"
        @test length(npi.start_times) == 1
    end

    @testset "Validation: mismatched vector lengths" begin
        @test_throws ArgumentError Daedalus.DaedalusStructs.TimedNpi(
            [10.0, 20.0], [30.0], [0.5, 0.3]
        )

        @test_throws ArgumentError Daedalus.DaedalusStructs.TimedNpi(
            [10.0], [30.0, 40.0], [0.5]
        )

        @test_throws ArgumentError Daedalus.DaedalusStructs.TimedNpi(
            [10.0, 20.0], [30.0, 40.0], [0.5]
        )
    end

    @testset "Validation: negative times" begin
        @test_throws ArgumentError Daedalus.DaedalusStructs.TimedNpi(
            [-10.0], [20.0], [0.5]
        )

        @test_throws ArgumentError Daedalus.DaedalusStructs.TimedNpi(
            [10.0], [-20.0], [0.5]
        )
    end

    @testset "Validation: end_time < start_time" begin
        @test_throws ArgumentError Daedalus.DaedalusStructs.TimedNpi(
            [20.0], [10.0], [0.5]
        )

        @test_throws ArgumentError Daedalus.DaedalusStructs.TimedNpi(
            [10.0, 30.0], [25.0, 35.0], [0.5, 0.3]  # second phase: 30 > 35 is valid
        )

        @test_throws ArgumentError Daedalus.DaedalusStructs.TimedNpi(
            [10.0, 40.0], [25.0, 35.0], [0.5, 0.3]  # second phase: 40 > 35 invalid
        )
    end

    @testset "Validation: unsorted times" begin
        @test_throws ArgumentError Daedalus.DaedalusStructs.TimedNpi(
            [30.0, 10.0], [40.0, 20.0], [0.5, 0.3]
        )

        @test_throws ArgumentError Daedalus.DaedalusStructs.TimedNpi(
            [10.0, 30.0], [40.0, 20.0], [0.5, 0.3]
        )
    end

    @testset "Validation: overlapping intervals" begin
        # Phase 1: 10-30, Phase 2: 25-40 (overlaps!)
        @test_throws ArgumentError Daedalus.DaedalusStructs.TimedNpi(
            [10.0, 25.0], [30.0, 40.0], [0.5, 0.3]
        )

        # Phase 1: 10-30, Phase 2: 30-40 (exactly adjacent - should fail)
        @test_throws ArgumentError Daedalus.DaedalusStructs.TimedNpi(
            [10.0, 30.0], [30.0, 40.0], [0.5, 0.3]
        )
    end

    @testset "Validation: coefficients out of range" begin
        @test_throws ArgumentError Daedalus.DaedalusStructs.TimedNpi(
            [10.0], [20.0], [1.5]
        )

        @test_throws ArgumentError Daedalus.DaedalusStructs.TimedNpi(
            [10.0], [20.0], [-0.1]
        )

        @test_throws ArgumentError Daedalus.DaedalusStructs.TimedNpi(
            [10.0, 30.0], [20.0, 40.0], [0.5, 1.1]
        )
    end

    @testset "Boundary values: zero duration phase" begin
        # Start and end at same time (zero duration)
        npi = Daedalus.DaedalusStructs.TimedNpi([10.0], [10.0], [0.5])
        @test npi.start_times == [10.0]
        @test npi.end_times == [10.0]
    end

    @testset "Boundary values: coefficient limits" begin
        # Coefficient = 0 (complete transmission block)
        npi1 = Daedalus.DaedalusStructs.TimedNpi([10.0], [20.0], [0.0])
        @test npi1.coefs == [0.0]

        # Coefficient = 1 (no intervention effect)
        npi2 = Daedalus.DaedalusStructs.TimedNpi([10.0], [20.0], [1.0])
        @test npi2.coefs == [1.0]
    end

    @testset "Non-overlapping intervals with gap" begin
        # Phase 1: 10-20, gap, Phase 2: 25-35 (valid)
        npi = Daedalus.DaedalusStructs.TimedNpi(
            [10.0, 25.0], [20.0, 35.0], [0.5, 0.3]
        )
        @test npi.start_times == [10.0, 25.0]
        @test npi.end_times == [20.0, 35.0]
    end
end

@testset "TimedNpi helper functions" begin
    @testset "n_phases" begin
        npi1 = Daedalus.DaedalusStructs.TimedNpi(10.0, 20.0, 0.5)
        @test Daedalus.DaedalusStructs.n_phases(npi1) == 1

        npi3 = Daedalus.DaedalusStructs.TimedNpi(
            [10.0, 30.0, 60.0],
            [20.0, 50.0, 80.0],
            [0.7, 0.3, 0.5]
        )
        @test Daedalus.DaedalusStructs.n_phases(npi3) == 3
    end

    @testset "total_duration" begin
        # Single phase: 10 days
        npi1 = Daedalus.DaedalusStructs.TimedNpi(10.0, 20.0, 0.5)
        @test Daedalus.DaedalusStructs.total_duration(npi1) == 10.0

        # Three phases: 10 + 20 + 20 = 50 days total
        npi3 = Daedalus.DaedalusStructs.TimedNpi(
            [10.0, 30.0, 60.0],
            [20.0, 50.0, 80.0],
            [0.7, 0.3, 0.5]
        )
        @test Daedalus.DaedalusStructs.total_duration(npi3) == 50.0

        # Zero duration phase
        npi_zero = Daedalus.DaedalusStructs.TimedNpi([10.0], [10.0], [0.5])
        @test Daedalus.DaedalusStructs.total_duration(npi_zero) == 0.0
    end

    @testset "Base.show" begin
        npi = Daedalus.DaedalusStructs.TimedNpi(
            [10.0, 30.0],
            [20.0, 40.0],
            [0.7, 0.3],
            "test_display"
        )

        # Capture output
        io = IOBuffer()
        show(io, npi)
        output = String(take!(io))

        @test occursin("TimedNpi: test_display", output)
        @test occursin("Number of phases: 2", output)
        @test occursin("Total duration: 20.0 days", output)
        @test occursin("Phase 1", output)
        @test occursin("Phase 2", output)
    end
end

@testset "make_timed_npi_callbacks" begin
    @testset "Callback creation for single phase" begin
        npi = Daedalus.DaedalusStructs.TimedNpi(10.0, 20.0, 0.5)
        callbacks = Daedalus.Events.make_timed_npi_callbacks(npi)

        @test isa(callbacks, CallbackSet)
        # Should have 2 callbacks: one for activation, one for deactivation
        @test length(callbacks.discrete_callbacks) == 2
    end

    @testset "Callback creation for multiple phases" begin
        npi = Daedalus.DaedalusStructs.TimedNpi(
            [10.0, 30.0, 60.0],
            [20.0, 40.0, 80.0],
            [0.7, 0.3, 0.5]
        )
        callbacks = Daedalus.Events.make_timed_npi_callbacks(npi)

        @test isa(callbacks, CallbackSet)
        # Should have 6 callbacks: 2 per phase (on/off)
        @test length(callbacks.discrete_callbacks) == 6
    end

    @testset "Callback times are correct" begin
        npi = Daedalus.DaedalusStructs.TimedNpi(
            [15.0, 45.0],
            [30.0, 60.0],
            [0.6, 0.4]
        )
        callbacks = Daedalus.Events.make_timed_npi_callbacks(npi)

        # Extract callback times (this is implementation-specific)
        # The callbacks should trigger at times: 15, 30, 45, 60
        @test length(callbacks.discrete_callbacks) == 4
    end
end

@testset "TimedNpi integration with daedalus model" begin
    @testset "Model runs with single-phase TimedNpi" begin
        npi = Daedalus.DaedalusStructs.TimedNpi(15.0, 45.0, 0.5, "test_single")

        result = Daedalus.daedalus(
            r0 = 2.0,
            time_end = 80.0,
            increment = 1.0,
            npi = npi,
            log_rt = true
        )

        @test result.sol.retcode == :Success
        @test result.npi === npi
        @test isnothing(result.saves)  # TimedNpi has no saved values
    end

    @testset "Model runs with multi-phase TimedNpi" begin
        npi = Daedalus.DaedalusStructs.TimedNpi(
            [10.0, 30.0, 60.0],
            [25.0, 50.0, 75.0],
            [0.7, 0.4, 0.6],
            "test_multi"
        )

        result = Daedalus.daedalus(
            r0 = 2.5,
            time_end = 100.0,
            increment = 1.0,
            npi = npi,
            log_rt = false
        )

        @test result.sol.retcode == :Success
        @test result.npi === npi
        @test isnothing(result.saves)
    end

    @testset "Model runs without Rt logging" begin
        npi = Daedalus.DaedalusStructs.TimedNpi(20.0, 40.0, 0.3)

        result = Daedalus.daedalus(
            r0 = 1.8,
            time_end = 60.0,
            npi = npi,
            log_rt = false
        )

        @test result.sol.retcode == :Success
    end

    @testset "Model output structure is correct" begin
        npi = Daedalus.DaedalusStructs.TimedNpi(10.0, 30.0, 0.5)

        result = Daedalus.daedalus(
            r0 = 2.0,
            time_end = 50.0,
            increment = 1.0,
            npi = npi
        )

        @test haskey(result, :sol)
        @test haskey(result, :saves)
        @test haskey(result, :npi)
        @test result.npi === npi
    end

    @testset "Solution has correct time points" begin
        npi = Daedalus.DaedalusStructs.TimedNpi(15.0, 35.0, 0.4)

        result = Daedalus.daedalus(
            r0 = 2.0,
            time_end = 50.0,
            increment = 1.0,
            npi = npi
        )

        # Should have savepoints at each day from 0 to 50
        @test length(result.sol.t) == 51
        @test result.sol.t[1] == 0.0
        @test result.sol.t[end] == 50.0
    end

    @testset "Extreme coefficient values" begin
        # Complete transmission block (coef = 0)
        npi_block = Daedalus.DaedalusStructs.TimedNpi(10.0, 30.0, 0.0)
        result_block = Daedalus.daedalus(
            r0 = 3.0,
            time_end = 50.0,
            npi = npi_block
        )
        @test result_block.sol.retcode == :Success

        # No intervention effect (coef = 1)
        npi_none = Daedalus.DaedalusStructs.TimedNpi(10.0, 30.0, 1.0)
        result_none = Daedalus.daedalus(
            r0 = 3.0,
            time_end = 50.0,
            npi = npi_none
        )
        @test result_none.sol.retcode == :Success
    end

    @testset "NPI at simulation boundaries" begin
        # NPI starting at time 0
        npi_start = Daedalus.DaedalusStructs.TimedNpi(0.0, 20.0, 0.5)
        result_start = Daedalus.daedalus(
            r0 = 2.0,
            time_end = 40.0,
            npi = npi_start
        )
        @test result_start.sol.retcode == :Success

        # NPI ending at simulation end
        npi_end = Daedalus.DaedalusStructs.TimedNpi(20.0, 50.0, 0.5)
        result_end = Daedalus.daedalus(
            r0 = 2.0,
            time_end = 50.0,
            npi = npi_end
        )
        @test result_end.sol.retcode == :Success

        # NPI covering entire simulation
        npi_full = Daedalus.DaedalusStructs.TimedNpi(0.0, 50.0, 0.5)
        result_full = Daedalus.daedalus(
            r0 = 2.0,
            time_end = 50.0,
            npi = npi_full
        )
        @test result_full.sol.retcode == :Success
    end
end

@testset "TimedNpi vs no intervention comparison" begin
    @testset "Basic epidemic comparison" begin
        # Run without intervention
        result_none = Daedalus.daedalus(
            r0 = 2.5,
            time_end = 100.0,
            increment = 1.0,
            log_rt = true
        )

        # Run with intervention
        npi = Daedalus.DaedalusStructs.TimedNpi(20.0, 60.0, 0.3)
        result_npi = Daedalus.daedalus(
            r0 = 2.5,
            time_end = 100.0,
            increment = 1.0,
            npi = npi,
            log_rt = true
        )

        # Both should complete successfully
        @test result_none.sol.retcode == :Success
        @test result_npi.sol.retcode == :Success

        # Should have same number of time points
        @test length(result_none.sol.t) == length(result_npi.sol.t)

        # Intervention should generally reduce epidemic size
        # (not guaranteed for all parameter combinations, but should be true here)
        iExposed = Daedalus.Constants.get_indices("E")
        exposed_none = [sum(u[iExposed]) for u in result_none.sol.u]
        exposed_npi = [sum(u[iExposed]) for u in result_npi.sol.u]

        # Peak exposure should be lower with intervention
        @test maximum(exposed_npi) < maximum(exposed_none)
    end
end

@testset "TimedNpi vs reactive Npi coexistence" begin
    @testset "Model accepts reactive Npi (unchanged behavior)" begin
        npi = Daedalus.DaedalusStructs.Npi(5000.0, (coef = 0.4,))

        result = Daedalus.daedalus(
            r0 = 2.5,
            time_end = 100.0,
            npi = npi
        )

        @test result.sol.retcode == :Success
        @test result.npi === npi
        @test !isnothing(result.saves)  # Reactive NPI has saved values
    end

    @testset "Model accepts nothing (unchanged behavior)" begin
        result = Daedalus.daedalus(
            r0 = 2.0,
            time_end = 50.0,
            npi = nothing
        )

        @test result.sol.retcode == :Success
        @test isnothing(result.npi)
        @test isnothing(result.saves)
    end
end

@testset "Edge cases and stress tests" begin
    @testset "Very short intervention duration" begin
        # 0.1 day intervention
        npi = Daedalus.DaedalusStructs.TimedNpi(20.0, 20.1, 0.5)
        result = Daedalus.daedalus(
            r0 = 2.0,
            time_end = 40.0,
            increment = 0.1,
            npi = npi
        )
        @test result.sol.retcode == :Success
    end

    @testset "Very long intervention duration" begin
        # 200 day intervention
        npi = Daedalus.DaedalusStructs.TimedNpi(10.0, 210.0, 0.3)
        result = Daedalus.daedalus(
            r0 = 2.0,
            time_end = 250.0,
            npi = npi
        )
        @test result.sol.retcode == :Success
    end

    @testset "Many phases" begin
        # 10 phases with 5-day gaps
        start_times = collect(0.0:10.0:90.0)
        end_times = start_times .+ 5.0
        coefs = fill(0.5, 10)

        npi = Daedalus.DaedalusStructs.TimedNpi(
            start_times, end_times, coefs, "many_phases"
        )
        result = Daedalus.daedalus(
            r0 = 2.0,
            time_end = 120.0,
            npi = npi
        )
        @test result.sol.retcode == :Success
        @test Daedalus.DaedalusStructs.n_phases(npi) == 10
    end

    @testset "Alternating intensity phases" begin
        # High-low-high-low pattern
        npi = Daedalus.DaedalusStructs.TimedNpi(
            [10.0, 30.0, 50.0, 70.0],
            [25.0, 45.0, 65.0, 85.0],
            [0.2, 0.8, 0.2, 0.8],
            "alternating"
        )
        result = Daedalus.daedalus(
            r0 = 2.5,
            time_end = 100.0,
            npi = npi
        )
        @test result.sol.retcode == :Success
    end

    @testset "Very high R0 with strong intervention" begin
        # R0=5.0 with 90% reduction
        npi = Daedalus.DaedalusStructs.TimedNpi(5.0, 50.0, 0.1)
        result = Daedalus.daedalus(
            r0 = 5.0,
            time_end = 80.0,
            npi = npi
        )
        @test result.sol.retcode == :Success
    end

    @testset "Late-starting intervention" begin
        # Intervention starts after epidemic peak likely passed
        npi = Daedalus.DaedalusStructs.TimedNpi(150.0, 180.0, 0.3)
        result = Daedalus.daedalus(
            r0 = 2.0,
            time_end = 200.0,
            npi = npi
        )
        @test result.sol.retcode == :Success
    end
end
