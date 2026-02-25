@testset "Zero-worker sector countries solve without warnings" begin
    # 30 countries have at least one sector with 0 workers in the data.
    # They previously emitted dt_NaN solver warnings because prepare_demog
    # returned zeros, causing Inf in the scaled contact matrix.
    zero_worker_countries = [
        "Australia", "Belgium", "Brunei", "Cambodia", "Chile", "China",
        "Costa Rica", "Cyprus", "Estonia", "Finland", "Hong Kong", "Iceland",
        "Japan", "Kazakhstan", "Laos", "Latvia", "Luxembourg", "Malaysia",
        "Malta", "Mexico", "Morocco", "Myanmar", "New Zealand", "Portugal",
        "Romania", "Rwanda", "Singapore", "Slovenia", "Switzerland", "Tunisia"
    ]
    for country in zero_worker_countries
        # prepare_demog must have no zeros (would cause Inf in contact scaling)
        demog = Daedalus.Data.prepare_demog(country)
        @test all(demog .> 0) broken=false

        # Scaled contact matrix must be finite
        cm = Daedalus.Data.prepare_contacts(country; scaled = true)
        @test all(isfinite, cm)

        # ODE must solve without NaN/Inf in the solution
        result = daedalus(country = country, r0 = 3.0, time_end = 100.0, log_rt = false)
        last_u = result.sol.u[end]

        @test length(result.sol.t) == 101 # hardcoded but oh well
        @test all(isfinite, last_u)
    end
end
