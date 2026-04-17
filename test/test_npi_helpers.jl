@testset "NPI Helpers" begin
    import Daedalus: create_time_intervention, create_reactive_intervention, create_multi_intervention

    @testset "Time-limited intervention creation" begin
        npi = create_time_intervention(:beta, x -> x .* 0.5, start=10.0, end_time=40.0)

        @test isa(npi, Daedalus.Npi)
        @test length(npi.effects) == 1
        @test npi.effects[1].target == :beta
        @test isa(npi.effects[1].trigger_on, Daedalus.TimeTrigger)
        @test isa(npi.effects[1].trigger_off, Daedalus.TimeTrigger)
        @test npi.effects[1].trigger_on.value == 10.0
        @test npi.effects[1].trigger_off.value == 40.0
    end

    @testset "Reactive intervention creation" begin
        npi = create_reactive_intervention(:beta, x -> x .* 0.6,
                                          on_threshold=5000.0, off_threshold=2000.0,
                                          compartment="H")

        @test isa(npi, Daedalus.Npi)
        @test length(npi.effects) == 1
        @test npi.effects[1].target == :beta
        @test isa(npi.effects[1].trigger_on, Daedalus.ReactiveTrigger)
        @test isa(npi.effects[1].trigger_off, Daedalus.ReactiveTrigger)
        @test npi.effects[1].trigger_on.value == 5000.0
        @test npi.effects[1].trigger_off.value == 2000.0
        @test npi.effects[1].trigger_on.name == "H"
        @test npi.effects[1].trigger_off.name == "H"
    end

    @testset "Multi-parameter intervention creation" begin
        reductions = Dict(:beta => x -> x .* 0.5, :omega => x -> x .* 1.1)
        npi = create_multi_intervention(reductions, start=20.0, end_time=60.0)

        @test isa(npi, Daedalus.Npi)
        @test length(npi.effects) == 2

        targets = Set(eff.target for eff in npi.effects)
        @test :beta in targets
        @test :omega in targets

        for eff in npi.effects
            @test isa(eff.trigger_on, Daedalus.TimeTrigger)
            @test eff.trigger_on.value == 20.0
            @test eff.trigger_off.value == 60.0
        end
    end

    @testset "Integration with daedalus" begin
        # Create and use a time-limited intervention
        npi = create_time_intervention(:beta, x -> x .* 0.5, start=10.0, end_time=30.0)
        infection = Daedalus.DataLoader.get_pathogen("sars-cov-2 delta")
        result = daedalus("Australia", infection, npi=npi, time_end=50.0, log_rt=false)

        @test isa(result, Daedalus.DaedalusOutput)
        @test !isempty(Daedalus.get_npi(result))
        @test length(result.sol.t) > 0  # ODE solved successfully
        @test all(isfinite, result.sol.u[end])
    end

    @testset "Reactive intervention in simulation" begin
        npi = create_reactive_intervention(:beta, x -> x .* 0.7,
                                          on_threshold=10000.0, off_threshold=5000.0,
                                          compartment="H")
        infection = Daedalus.DataLoader.get_pathogen("sars-cov-2 delta")
        result = daedalus("Australia", infection, npi=npi, time_end=100.0, log_rt=false)

        @test isa(result, Daedalus.DaedalusOutput)
        @test length(result.sol.t) > 0  # ODE solved successfully
        @test all(isfinite, result.sol.u[end])
    end
end
