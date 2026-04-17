@testset "Vaccination Events" begin
    import Daedalus: Vaccination, daedalus, DaedalusOutput, get_vaccination, get_npi

    @testset "Vaccination struct creation" begin
        vax = Vaccination(0.0, 0.01, 0.7, 0.9, 180.0)

        @test isa(vax, Vaccination)
        @test vax.start_time == 0.0
        @test vax.rate == 0.01
        @test vax.uptake_limit == 0.7
        @test vax.efficacy == 0.9
        @test vax.waning_period == 180.0
        @test vax.ison == false
    end

    @testset "Make events for vaccination" begin
        vax = Vaccination(10.0, 0.01, 0.7, 0.9, 180.0)
        savepoints = 0.0:1.0:100.0

        cb_set = Daedalus.make_events(vax, savepoints)

        @test isa(cb_set, DiffEqCallbacks.CallbackSet)
    end

    @testset "Vaccination-only run" begin
        vax = Vaccination(0.0, 0.01, 0.7, 0.9, 180.0)
        infection = Daedalus.DataLoader.get_pathogen("sars-cov-2 delta")
        result = daedalus("Canada", infection, vaccination=vax, time_end=50.0, log_rt=false)

        @test isa(result, DaedalusOutput)
        @test !isempty(get_vaccination(result))
        @test length(result.events) == 1
        @test isa(result.events[1], Vaccination)
    end

    @testset "Combined NPI and vaccination run" begin
        npi = Daedalus.create_time_intervention(:beta, x -> x .* 0.5, start=10.0, end_time=30.0)
        vax = Vaccination(5.0, 0.01, 0.7, 0.9, 180.0)
        infection = Daedalus.DataLoader.get_pathogen("sars-cov-2 delta")
        result = daedalus("Canada", infection, npi=npi, vaccination=vax, time_end=50.0, log_rt=false)

        @test isa(result, DaedalusOutput)
        @test !isempty(get_npi(result))
        @test !isempty(get_vaccination(result))
        @test length(result.events) == 2
    end

    @testset "No events run" begin
        infection = Daedalus.DataLoader.get_pathogen("sars-cov-2 delta")
        result = daedalus("Canada", infection, time_end=50.0, log_rt=false)

        @test isa(result, DaedalusOutput)
        @test isempty(get_npi(result))
        @test isempty(get_vaccination(result))
        @test isempty(result.events)
    end

    @testset "Vaccination affects state vs no vaccination" begin
        vax = Vaccination(0.0, 0.01, 0.7, 0.9, 180.0)
        infection = Daedalus.DataLoader.get_pathogen("sars-cov-2 delta")

        result_vax = daedalus("Canada", infection, vaccination=vax, time_end=100.0, log_rt=false)
        result_no_vax = daedalus("Canada", infection, time_end=100.0, log_rt=false)

        # With vaccination, should have different state at end
        # (vaccinated stratum should have non-zero susceptibles)
        @test result_vax.sol.u[end] != result_no_vax.sol.u[end]
    end

    @testset "DaedalusOutput show method with vaccination" begin
        vax = Vaccination(0.0, 0.01, 0.7, 0.9, 180.0)
        infection = Daedalus.DataLoader.get_pathogen("sars-cov-2 delta")
        result = daedalus("Canada", infection, vaccination=vax, time_end=50.0, log_rt=false)

        # Should not error
        io = IOBuffer()
        show(io, result)
        output = String(take!(io))

        @test contains(output, "DaedalusOutput")
        @test contains(output, "Vaccination: active")
        @test contains(output, "R₀:")
    end
end
