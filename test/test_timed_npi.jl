using Daedalus
using Test
using OrdinaryDiffEq
using DiffEqCallbacks

@testset "ParamEffect with reactive triggers" begin
    @testset "Construct ParamEffect with ReactiveTrigger" begin
        trigger_on = Daedalus.DaedalusStructs.ReactiveTrigger(5000.0, "H")
        trigger_off = Daedalus.DaedalusStructs.ReactiveTrigger(1.0, "Rt")
        effect = Daedalus.DaedalusStructs.ParamEffect(
            :beta, x -> x .* 0.4, x -> x ./ 0.4,
            trigger_on, trigger_off
        )

        @test effect.target == :beta
        @test effect.trigger_on.value == 5000.0
        @test effect.trigger_on.name == "H"
        @test effect.trigger_off.value == 1.0
        @test effect.trigger_off.name == "Rt"
        @test !effect.ison
    end

    @testset "ParamEffect with timed triggers" begin
        trigger_on = Daedalus.DaedalusStructs.TimeTrigger(10.0, "time")
        trigger_off = Daedalus.DaedalusStructs.TimeTrigger(20.0, "time")
        effect = Daedalus.DaedalusStructs.ParamEffect(
            :beta, x -> x .* 0.5, x -> x ./ 0.5,
            trigger_on, trigger_off
        )

        @test effect.target == :beta
        @test effect.trigger_on.value == 10.0
        @test effect.trigger_off.value == 20.0
        @test !effect.ison
    end
end

@testset "Npi with ParamEffect" begin
    @testset "Npi construction with single reactive effect" begin
        trigger_on = Daedalus.DaedalusStructs.ReactiveTrigger(5000.0, "H")
        trigger_off = Daedalus.DaedalusStructs.ReactiveTrigger(1.0, "Rt")
        effect = Daedalus.DaedalusStructs.ParamEffect(
            :beta, x -> x .* 0.4, x -> x ./ 0.4,
            trigger_on, trigger_off
        )
        npi = Daedalus.DaedalusStructs.Npi([effect])

        @test length(npi.effects) == 1
        @test npi.effects[1].target == :beta
    end

    @testset "Npi construction with multiple effects" begin
        trigger_on1 = Daedalus.DaedalusStructs.ReactiveTrigger(5000.0, "H")
        trigger_off1 = Daedalus.DaedalusStructs.ReactiveTrigger(1.0, "Rt")
        effect1 = Daedalus.DaedalusStructs.ParamEffect(
            :beta, x -> x .* 0.4, x -> x ./ 0.4,
            trigger_on1, trigger_off1
        )

        trigger_on2 = Daedalus.DaedalusStructs.ReactiveTrigger(100.0, "D")
        trigger_off2 = Daedalus.DaedalusStructs.ReactiveTrigger(1.0, "Rt")
        effect2 = Daedalus.DaedalusStructs.ParamEffect(
            :omega, x -> x .* 0.8, x -> x ./ 0.8,
            trigger_on2, trigger_off2
        )

        npi = Daedalus.DaedalusStructs.Npi([effect1, effect2])
        @test length(npi.effects) == 2
        @test npi.effects[1].target == :beta
        @test npi.effects[2].target == :omega
    end
end

@testset "Reactive NPI model integration" begin
    @testset "Model runs with reactive NPI" begin
        trigger_on = Daedalus.DaedalusStructs.ReactiveTrigger(5000.0, "H")
        trigger_off = Daedalus.DaedalusStructs.ReactiveTrigger(1.0, "Rt")
        effect = Daedalus.DaedalusStructs.ParamEffect(
            :beta, x -> x .* 0.4, x -> x ./ 0.4,
            trigger_on, trigger_off
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

    @testset "Model runs without NPI" begin
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

@testset "Npi with multiple reactive effects" begin
    @testset "Multiple parameters with same triggers" begin
        trigger_on = Daedalus.DaedalusStructs.ReactiveTrigger(5000.0, "H")
        trigger_off = Daedalus.DaedalusStructs.ReactiveTrigger(1.0, "Rt")

        effect_beta = Daedalus.DaedalusStructs.ParamEffect(
            :beta, x -> x .* 0.6, x -> x ./ 0.6,
            trigger_on, trigger_off
        )
        effect_omega = Daedalus.DaedalusStructs.ParamEffect(
            :omega, x -> x .* 0.6, x -> x ./ 0.6,
            trigger_on, trigger_off
        )
        npi = Daedalus.DaedalusStructs.Npi([effect_beta, effect_omega])

        @test length(npi.effects) == 2
        @test npi.effects[1].target == :beta
        @test npi.effects[2].target == :omega

        result = Daedalus.daedalus(
            "Australia",
            "influenza 2009",
            time_end = 100.0,
            npi = npi
        )
        @test result.sol.retcode == OrdinaryDiffEq.ReturnCode.Success
    end

    @testset "Different transformation functions per parameter" begin
        trigger_on_beta = Daedalus.DaedalusStructs.ReactiveTrigger(5000.0, "H")
        trigger_off_beta = Daedalus.DaedalusStructs.ReactiveTrigger(1.0, "Rt")
        effect_beta = Daedalus.DaedalusStructs.ParamEffect(
            :beta, x -> x .* 0.4, x -> x ./ 0.4,
            trigger_on_beta, trigger_off_beta
        )

        trigger_on_omega = Daedalus.DaedalusStructs.ReactiveTrigger(100.0, "D")
        trigger_off_omega = Daedalus.DaedalusStructs.ReactiveTrigger(1.0, "Rt")
        effect_omega = Daedalus.DaedalusStructs.ParamEffect(
            :omega, x -> x .* 0.7, x -> x ./ 0.7,
            trigger_on_omega, trigger_off_omega
        )

        npi = Daedalus.DaedalusStructs.Npi([effect_beta, effect_omega])
        @test length(npi.effects) == 2

        test_val = 2.0
        @test npi.effects[1].func(test_val) ≈ 0.8
        @test npi.effects[2].func(test_val) ≈ 1.4

        result = Daedalus.daedalus(
            "Australia",
            "influenza 2009",
            time_end = 100.0,
            npi = npi
        )
        @test result.sol.retcode == OrdinaryDiffEq.ReturnCode.Success
    end

    @testset "Independent triggers for each effect" begin
        trigger_on_beta = Daedalus.DaedalusStructs.ReactiveTrigger(5000.0, "H")
        trigger_off_beta = Daedalus.DaedalusStructs.ReactiveTrigger(1.0, "Rt")
        effect_beta = Daedalus.DaedalusStructs.ParamEffect(
            :beta, x -> x .* 0.5, x -> x ./ 0.5,
            trigger_on_beta, trigger_off_beta
        )

        trigger_on_omega = Daedalus.DaedalusStructs.ReactiveTrigger(10000.0, "I")
        trigger_off_omega = Daedalus.DaedalusStructs.ReactiveTrigger(2000.0, "H")
        effect_omega = Daedalus.DaedalusStructs.ParamEffect(
            :omega, x -> x .* 0.9, x -> x ./ 0.9,
            trigger_on_omega, trigger_off_omega
        )

        npi = Daedalus.DaedalusStructs.Npi([effect_beta, effect_omega])
        @test length(npi.effects) == 2

        beta_eff = npi.effects[1]
        omega_eff = npi.effects[2]

        @test beta_eff.trigger_on.name == "H"
        @test beta_eff.trigger_on.value == 5000.0
        @test beta_eff.trigger_off.name == "Rt"
        @test omega_eff.trigger_on.name == "I"
        @test omega_eff.trigger_on.value == 10000.0
        @test omega_eff.trigger_off.name == "H"

        result = Daedalus.daedalus(
            "Australia",
            "influenza 2009",
            time_end = 100.0,
            npi = npi
        )
        @test result.sol.retcode == OrdinaryDiffEq.ReturnCode.Success
    end

    @testset "Each effect tracks its own activation state" begin
        trigger_on1 = Daedalus.DaedalusStructs.ReactiveTrigger(5000.0, "H")
        trigger_off1 = Daedalus.DaedalusStructs.ReactiveTrigger(1.0, "Rt")
        effect1 = Daedalus.DaedalusStructs.ParamEffect(
            :beta, x -> x .* 0.4, x -> x ./ 0.4,
            trigger_on1, trigger_off1
        )

        trigger_on2 = Daedalus.DaedalusStructs.ReactiveTrigger(100.0, "D")
        trigger_off2 = Daedalus.DaedalusStructs.ReactiveTrigger(1.0, "Rt")
        effect2 = Daedalus.DaedalusStructs.ParamEffect(
            :omega, x -> x .* 0.7, x -> x ./ 0.7,
            trigger_on2, trigger_off2
        )

        npi = Daedalus.DaedalusStructs.Npi([effect1, effect2])
        @test !npi.effects[1].ison
        @test !npi.effects[2].ison
    end
end
