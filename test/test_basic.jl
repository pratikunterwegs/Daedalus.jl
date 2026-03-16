@testset "DAEDALUS model" begin
    try
        infection = Daedalus.DataLoader.get_pathogen("sars-cov-2 delta")
        infection.r0 = 5.0
        daedalus("Australia", infection, time_end = 100.0)
        @test true
    catch e
        @test false
    end
end

@testset "daedalus accepts CountryData struct" begin
    cd = Daedalus.DataLoader.get_country("Australia")
    infection = Daedalus.DataLoader.get_pathogen("sars-cov-2 delta")
    infection.r0 = 3.0
    result = daedalus(cd, infection, time_end = 100.0, log_rt = false)
    @test length(result.sol.t) == 101
    @test all(isfinite, result.sol.u[end])

    # result must be identical to string-based call
    infection_str = Daedalus.DataLoader.get_pathogen("sars-cov-2 delta")
    infection_str.r0 = 3.0
    result_str = daedalus("Australia", infection_str, time_end = 100.0, log_rt = false)
    @test result.sol.u[end] ≈ result_str.sol.u[end]
end

# tests for helper functions
@testset "Daedalus helpers for β and NGM" begin
    r0 = 1.3
    sigma = 0.217
    p_sigma = 0.867
    epsilon = 0.58
    rho = 0.003
    gamma_Ia = 0.476
    gamma_Is = 0.25

    aus = Daedalus.DataLoader.get_country("Australia")
    contacts_unscaled = Daedalus.Data.total_contacts(
        Daedalus.Data.prepare_contacts(aus; scaled = false)
    )

    beta = Daedalus.Helpers.get_beta(
        contacts_unscaled,
        r0, sigma, p_sigma, epsilon, gamma_Ia, gamma_Is
    )
    @test typeof(beta) == Float64

    ngm = Daedalus.Helpers.get_ngm(
        contacts_unscaled,
        beta, sigma, p_sigma, epsilon, gamma_Ia, gamma_Is
    )
    lambda = maximum(real(eigen(ngm).values))

    @test lambda ≈ r0
end
