using Daedalus
using Test
using OrdinaryDiffEq
using DiffEqCallbacks

@testset "make_events with timed effects" begin
    @testset "Callback creation for single timed phase" begin
        effect = Daedalus.DaedalusStructs.TimedEffect(
            :beta, x -> x .* 0.5, x -> x ./ 0.5, 10.0, 20.0
        )
        npi = Daedalus.DaedalusStructs.Npi([effect])
        savepoints = 0.0:1.0:30.0
        callbacks = Daedalus.Events.make_events(npi, savepoints)

        @test isa(callbacks, DiffEqCallbacks.CallbackSet)
        # Should have 2 callbacks: one for activation at start_time, one for deactivation at end_time
        @test length(callbacks.discrete_callbacks) == 2
    end

    @testset "Callback creation for multiple timed phases" begin
        effect1 = Daedalus.DaedalusStructs.TimedEffect(
            :beta, x -> x .* 0.7, x -> x ./ 0.7, 10.0, 20.0
        )
        effect2 = Daedalus.DaedalusStructs.TimedEffect(
            :beta, x -> x .* 0.3, x -> x ./ 0.3, 30.0, 40.0
        )
        effect3 = Daedalus.DaedalusStructs.TimedEffect(
            :beta, x -> x .* 0.5, x -> x ./ 0.5, 60.0, 80.0
        )
        npi = Daedalus.DaedalusStructs.Npi([effect1, effect2, effect3])
        savepoints = 0.0:1.0:100.0
        callbacks = Daedalus.Events.make_events(npi, savepoints)

        @test isa(callbacks, DiffEqCallbacks.CallbackSet)
        # Should have 6 callbacks: 2 per effect (on/off)
        @test length(callbacks.discrete_callbacks) == 6
    end
end

@testset "TimedEffect integration with daedalus model" begin
    @testset "Model runs with single timed effect" begin
        effect = Daedalus.DaedalusStructs.TimedEffect(
            :beta, x -> x .* 0.5, x -> x ./ 0.5, 15.0, 45.0
        )
        npi = Daedalus.DaedalusStructs.Npi([effect])

        infection = Daedalus.DataLoader.get_pathogen("sars-cov-2 delta")
        infection.r0 = 2.0

        result = Daedalus.daedalus(
            "Australia",
            infection,
            time_end = 80.0,
            increment = 1.0,
            npi = npi,
            log_rt = true
        )

        @test result.sol.retcode == OrdinaryDiffEq.ReturnCode.Success
        @test result.npi === npi
        @test isnothing(result.saves)  # TimedEffect has no saved values
    end

    @testset "Model runs with multiple timed effects" begin
        effect1 = Daedalus.DaedalusStructs.TimedEffect(
            :beta, x -> x .* 0.7, x -> x ./ 0.7, 10.0, 25.0
        )
        effect2 = Daedalus.DaedalusStructs.TimedEffect(
            :beta, x -> x .* 0.4, x -> x ./ 0.4, 30.0, 50.0
        )
        effect3 = Daedalus.DaedalusStructs.TimedEffect(
            :beta, x -> x .* 0.6, x -> x ./ 0.6, 60.0, 75.0
        )
        npi = Daedalus.DaedalusStructs.Npi([effect1, effect2, effect3])

        result = Daedalus.daedalus(
            "Australia",
            "sars-cov-1",
            time_end = 100.0,
            increment = 1.0,
            npi = npi,
            log_rt = false
        )

        @test result.sol.retcode == OrdinaryDiffEq.ReturnCode.Success
        @test result.npi === npi
        @test isnothing(result.saves)
    end

    @testset "Model runs without Rt logging" begin
        effect = Daedalus.DaedalusStructs.TimedEffect(
            :beta, x -> x .* 0.3, x -> x ./ 0.3, 20.0, 40.0
        )
        npi = Daedalus.DaedalusStructs.Npi([effect])

        result = Daedalus.daedalus(
            "Australia",
            "sars-cov-1",
            time_end = 60.0,
            npi = npi,
            log_rt = false
        )

        @test result.sol.retcode == OrdinaryDiffEq.ReturnCode.Success
    end
end

@testset "TimedEffect and reactive Npi coexistence" begin
    @testset "Model accepts reactive Npi (unchanged behavior)" begin
        effect = Daedalus.DaedalusStructs.ReactiveEffect(
            :beta, x -> x .* 0.4, x -> x ./ 0.4;
            on = ("H", 5000.0), off = ("Rt", 1.0)
        )
        npi = Daedalus.DaedalusStructs.Npi([effect])

        result = Daedalus.daedalus(
            "Australia",
            "influenza 2009",
            time_end = 100.0,
            npi = npi
        )

        @test result.sol.retcode == OrdinaryDiffEq.ReturnCode.Success
        @test result.npi === npi
        @test !isnothing(result.saves)  # Reactive NPI has saved values
    end

    @testset "Model accepts nothing (unchanged behavior)" begin
        result = Daedalus.daedalus(
            "Australia",
            "influenza 2009",
            time_end = 50.0,
            npi = nothing
        )

        @test result.sol.retcode == OrdinaryDiffEq.ReturnCode.Success
        @test isnothing(result.npi)
        @test isnothing(result.saves)
    end
end

@testset "Npi with flexible parameter effects" begin
    @testset "Single parameter (new style with ReactiveEffect)" begin
        effect = Daedalus.DaedalusStructs.ReactiveEffect(
            :beta, x -> x .* 0.5, x -> x ./ 0.5;
            on = ("H", 5000.0), off = ("Rt", 1.0)
        )
        npi = Daedalus.DaedalusStructs.Npi([effect])
        @test length(npi.effects) == 1
        @test npi.effects[1].target == :beta

        result = Daedalus.daedalus(
            "Australia",
            "influenza 2009",
            time_end = 200.0,
            npi = npi
        )
        @test result.sol.retcode == OrdinaryDiffEq.ReturnCode.Success
    end

    @testset "Multiple parameters" begin
        effect_beta = Daedalus.DaedalusStructs.ReactiveEffect(
            :beta, x -> x .* 0.6, x -> x ./ 0.6;
            on = ("H", 5000.0), off = ("Rt", 1.0)
        )
        effect_omega = Daedalus.DaedalusStructs.ReactiveEffect(
            :omega, x -> x .* 0.6, x -> x ./ 0.6;
            on = ("H", 5000.0), off = ("Rt", 1.0)
        )
        npi = Daedalus.DaedalusStructs.Npi([effect_beta, effect_omega])
        @test length(npi.effects) == 2
        @test npi.effects[1].target == :beta
        @test npi.effects[2].target == :omega

        result = Daedalus.daedalus(
            "Australia",
            "influenza 2009",
            time_end = 200.0,
            npi = npi
        )
        @test result.sol.retcode == OrdinaryDiffEq.ReturnCode.Success
    end

    @testset "Different functions for different parameters" begin
        effect_beta = Daedalus.DaedalusStructs.ReactiveEffect(
            :beta, x -> x .* 0.4, x -> x ./ 0.4;
            on = ("H", 5000.0), off = ("Rt", 1.0)
        )
        effect_omega = Daedalus.DaedalusStructs.ReactiveEffect(
            :omega, x -> x .* 0.7, x -> x ./ 0.7;
            on = ("D", 100.0), off = ("Rt", 1.0)
        )
        npi = Daedalus.DaedalusStructs.Npi([effect_beta, effect_omega])
        @test length(npi.effects) == 2
        @test npi.effects[1].target == :beta
        @test npi.effects[2].target == :omega

        result = Daedalus.daedalus(
            "Australia",
            "influenza 2009",
            time_end = 200.0,
            npi = npi
        )
        @test result.sol.retcode == OrdinaryDiffEq.ReturnCode.Success
    end

    @testset "Per-parameter independent triggers" begin
        effect_beta = Daedalus.DaedalusStructs.ReactiveEffect(
            :beta, x -> x .* 0.4, x -> x ./ 0.4;
            on = ("H", 5000.0), off = ("Rt", 1.0)
        )
        effect_omega = Daedalus.DaedalusStructs.ReactiveEffect(
            :omega, x -> x .* 0.7, x -> x ./ 0.7;
            on = ("D", 100.0), off = ("Rt", 1.0)
        )
        npi = Daedalus.DaedalusStructs.Npi([effect_beta, effect_omega])
        @test length(npi.effects) == 2
        @test npi.effects[1].target == :beta
        @test npi.effects[2].target == :omega
        test_val = 2.0
        @test npi.effects[1].func(test_val) ≈ 0.8
        @test npi.effects[2].func(test_val) ≈ 1.4
        @test npi.effects[1].comp_on.name == "H"
        @test npi.effects[2].comp_on.name == "D"

        result = Daedalus.daedalus(
            "Australia",
            "influenza 2009",
            time_end = 200.0,
            npi = npi
        )
        @test result.sol.retcode == OrdinaryDiffEq.ReturnCode.Success
    end

    @testset "Multiple effects with different on/off conditions" begin
        effect_beta = Daedalus.DaedalusStructs.ReactiveEffect(
            :beta, x -> x .* 0.5, x -> x ./ 0.5;
            on = ("H", 5000.0), off = ("Rt", 1.0)
        )
        effect_omega = Daedalus.DaedalusStructs.ReactiveEffect(
            :omega, x -> x .* 0.9, x -> x ./ 0.9;
            on = ("I", 10000.0), off = ("H", 2000.0)
        )
        npi = Daedalus.DaedalusStructs.Npi([effect_beta, effect_omega])
        @test length(npi.effects) == 2
        param_names = Set(eff.target for eff in npi.effects)
        @test param_names == Set([:beta, :omega])

        beta_eff = [eff for eff in npi.effects if eff.target == :beta][1]
        omega_eff = [eff for eff in npi.effects if eff.target == :omega][1]
        test_val = 10.0
        @test beta_eff.func(test_val) ≈ 5.0
        @test omega_eff.func(test_val) ≈ 9.0

        @test beta_eff.comp_on.name == "H"
        @test beta_eff.comp_on.value == 5000.0
        @test omega_eff.comp_on.name == "I"
        @test omega_eff.comp_on.value == 10000.0
        @test beta_eff.comp_off.name == "Rt"
        @test omega_eff.comp_off.name == "H"

        result = Daedalus.daedalus(
            "Australia",
            "influenza 2009",
            time_end = 200.0,
            npi = npi
        )
        @test result.sol.retcode == OrdinaryDiffEq.ReturnCode.Success
    end
end

@testset "Npi with per-effect trigger conditions" begin
    @testset "Direct ReactiveEffect construction with custom triggers" begin
        effect_beta = Daedalus.DaedalusStructs.ReactiveEffect(
            :beta, x -> x .* 0.4, x -> x ./ 0.4;
            on = ("H", 5000.0), off = ("Rt", 1.0)
        )
        effect_gamma = Daedalus.DaedalusStructs.ReactiveEffect(
            :omega, x -> x .* 0.8, x -> x ./ 0.8;
            on = ("D", 100.0), off = ("Rt", 1.0)
        )

        npi = Daedalus.DaedalusStructs.Npi([effect_beta, effect_gamma])
        @test length(npi.effects) == 2
        @test npi.effects[1].target == :beta
        @test npi.effects[2].target == :omega
        @test npi.effects[1].comp_on.name == "H"
        @test npi.effects[2].comp_on.name == "D"
        @test npi.effects[1].comp_on.value == 5000.0
        @test npi.effects[2].comp_on.value == 100.0

        result = Daedalus.daedalus(
            "Australia",
            "influenza 2009",
            time_end = 200.0,
            npi = npi
        )
        @test result.sol.retcode == OrdinaryDiffEq.ReturnCode.Success
    end

    @testset "Infectious (I) compartment trigger" begin
        effect = Daedalus.DaedalusStructs.ReactiveEffect(
            :beta, x -> x .* 0.5, x -> x ./ 0.5;
            on = ("I", 8000.0), off = ("Rt", 1.0)
        )
        npi = Daedalus.DaedalusStructs.Npi([effect])
        @test npi.effects[1].comp_on.name == "I"

        result = Daedalus.daedalus(
            "Australia",
            "influenza 2009",
            time_end = 200.0,
            npi = npi
        )
        @test result.sol.retcode == OrdinaryDiffEq.ReturnCode.Success
    end

    @testset "Each effect tracks its own activation state" begin
        effect1 = Daedalus.DaedalusStructs.ReactiveEffect(
            :beta, x -> x .* 0.4, x -> x ./ 0.4;
            on = ("H", 5000.0), off = ("Rt", 1.0)
        )
        effect2 = Daedalus.DaedalusStructs.ReactiveEffect(
            :omega, x -> x .* 0.7, x -> x ./ 0.7;
            on = ("D", 100.0), off = ("Rt", 1.0)
        )
        npi = Daedalus.DaedalusStructs.Npi([effect1, effect2])

        @test !npi.effects[1].ison
        @test !npi.effects[2].ison
    end
end
