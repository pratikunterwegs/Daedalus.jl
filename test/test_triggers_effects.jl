using Daedalus
using Test

@testset "Trigger types" begin
    @testset "ReactiveTrigger construction and fields" begin
        trigger = Daedalus.DaedalusStructs.ReactiveTrigger(5000.0, "H")
        @test trigger.value == 5000.0
        @test trigger.name == "H"
    end

    @testset "TimeTrigger construction and fields" begin
        trigger = Daedalus.DaedalusStructs.TimeTrigger(10.0, "time")
        @test trigger.value == 10.0
        @test trigger.name == "time"
    end

    @testset "TimeTrigger with default name" begin
        trigger = Daedalus.DaedalusStructs.TimeTrigger(15.0)
        @test trigger.value == 15.0
        @test trigger.name == "time"
    end

    @testset "Trigger is abstract type" begin
        trigger_reactive = Daedalus.DaedalusStructs.ReactiveTrigger(5000.0, "H")
        trigger_time = Daedalus.DaedalusStructs.TimeTrigger(10.0)

        @test trigger_reactive isa Daedalus.DaedalusStructs.Trigger
        @test trigger_time isa Daedalus.DaedalusStructs.Trigger
    end
end

@testset "ParamEffect basic functionality" begin
    @testset "ParamEffect with ReactiveTriggers" begin
        trigger_on = Daedalus.DaedalusStructs.ReactiveTrigger(5000.0, "H")
        trigger_off = Daedalus.DaedalusStructs.ReactiveTrigger(1.0, "Rt")

        effect = Daedalus.DaedalusStructs.ParamEffect(
            :beta, x -> x .* 0.5, x -> x ./ 0.5,
            trigger_on, trigger_off
        )

        @test effect.target == :beta
        @test effect.trigger_on.value == 5000.0
        @test effect.trigger_off.value == 1.0
        @test !effect.ison  # Initially inactive
    end

    @testset "ParamEffect with TimeTriggers" begin
        trigger_on = Daedalus.DaedalusStructs.TimeTrigger(10.0)
        trigger_off = Daedalus.DaedalusStructs.TimeTrigger(20.0)

        effect = Daedalus.DaedalusStructs.ParamEffect(
            :beta, x -> x .* 0.7, x -> x ./ 0.7,
            trigger_on, trigger_off
        )

        @test effect.target == :beta
        @test isa(effect.trigger_on, Daedalus.DaedalusStructs.TimeTrigger)
        @test isa(effect.trigger_off, Daedalus.DaedalusStructs.TimeTrigger)
    end

    @testset "ParamEffect is an Effect" begin
        trigger_on = Daedalus.DaedalusStructs.ReactiveTrigger(5000.0, "H")
        trigger_off = Daedalus.DaedalusStructs.ReactiveTrigger(1.0, "Rt")

        effect = Daedalus.DaedalusStructs.ParamEffect(
            :beta, x -> x .* 0.5, x -> x ./ 0.5,
            trigger_on, trigger_off
        )

        @test effect isa Daedalus.DaedalusStructs.Effect
    end

    @testset "ParamEffect function application" begin
        trigger_on = Daedalus.DaedalusStructs.ReactiveTrigger(5000.0, "H")
        trigger_off = Daedalus.DaedalusStructs.ReactiveTrigger(1.0, "Rt")

        effect = Daedalus.DaedalusStructs.ParamEffect(
            :beta, x -> x .* 0.4, x -> x ./ 0.4,
            trigger_on, trigger_off
        )

        test_val = 10.0
        @test effect.func(test_val) ≈ 4.0
        @test effect.reset_func(4.0) ≈ 10.0
    end
end

@testset "Npi container" begin
    @testset "Npi with single effect" begin
        trigger_on = Daedalus.DaedalusStructs.ReactiveTrigger(5000.0, "H")
        trigger_off = Daedalus.DaedalusStructs.ReactiveTrigger(1.0, "Rt")
        effect = Daedalus.DaedalusStructs.ParamEffect(
            :beta, x -> x .* 0.5, x -> x ./ 0.5,
            trigger_on, trigger_off
        )

        npi = Daedalus.DaedalusStructs.Npi([effect])
        @test length(npi.effects) == 1
        @test npi.effects[1] === effect
    end

    @testset "Npi with multiple effects" begin
        trigger1_on = Daedalus.DaedalusStructs.ReactiveTrigger(5000.0, "H")
        trigger1_off = Daedalus.DaedalusStructs.ReactiveTrigger(1.0, "Rt")
        effect1 = Daedalus.DaedalusStructs.ParamEffect(
            :beta, x -> x .* 0.5, x -> x ./ 0.5,
            trigger1_on, trigger1_off
        )

        trigger2_on = Daedalus.DaedalusStructs.ReactiveTrigger(100.0, "D")
        trigger2_off = Daedalus.DaedalusStructs.ReactiveTrigger(1.0, "Rt")
        effect2 = Daedalus.DaedalusStructs.ParamEffect(
            :omega, x -> x .* 0.8, x -> x ./ 0.8,
            trigger2_on, trigger2_off
        )

        npi = Daedalus.DaedalusStructs.Npi([effect1, effect2])
        @test length(npi.effects) == 2
        @test npi.effects[1].target == :beta
        @test npi.effects[2].target == :omega
    end

    @testset "Npi is an Event" begin
        trigger_on = Daedalus.DaedalusStructs.ReactiveTrigger(5000.0, "H")
        trigger_off = Daedalus.DaedalusStructs.ReactiveTrigger(1.0, "Rt")
        effect = Daedalus.DaedalusStructs.ParamEffect(
            :beta, x -> x .* 0.5, x -> x ./ 0.5,
            trigger_on, trigger_off
        )

        npi = Daedalus.DaedalusStructs.Npi([effect])
        @test npi isa Daedalus.DaedalusStructs.Event
    end
end

@testset "Exported types verification" begin
    # Verify that the correct types are exported
    @test Daedalus.DaedalusStructs.ParamEffect isa DataType
    @test Daedalus.DaedalusStructs.Effect isa DataType
    @test Daedalus.DaedalusStructs.Trigger isa DataType
    @test Daedalus.DaedalusStructs.ReactiveTrigger isa DataType
    @test Daedalus.DaedalusStructs.TimeTrigger isa DataType
    @test Daedalus.DaedalusStructs.Npi isa DataType
end
