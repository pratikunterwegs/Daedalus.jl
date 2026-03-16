@testset "Multiple contact settings" begin
    cd = Daedalus.DataLoader.get_country("Australia")
    cm = cd.contact_matrix

    @testset "get_settings returns 1 for a single matrix" begin
        @test Daedalus.Data.get_settings(cd) == 1
    end

    @testset "get_settings returns 2 for a vector of two matrices" begin
        cd2 = deepcopy(cd)
        cd2.contact_matrix = [cm, cm]
        @test Daedalus.Data.get_settings(cd2) == 2
    end

    @testset "contacts3d shape is (49, 49, 1) for a single setting" begin
        c3d = Daedalus.Data.contacts3d(cd)
        @test size(c3d) == (49, 49, 1)
    end

    @testset "contacts3d shape is (49, 49, 2) for two settings" begin
        cd2 = deepcopy(cd)
        cd2.contact_matrix = [cm, cm]
        c3d = Daedalus.Data.contacts3d(cd2)
        @test size(c3d) == (49, 49, 2)
    end

    @testset "total_contacts sums a vector of matrices element-wise" begin
        # A vector of two identical matrices should sum to 2 × the single matrix.
        cd2 = deepcopy(cd)
        cd2.contact_matrix = [cm, cm]
        single = Daedalus.Data.total_contacts(
            Daedalus.Data.prepare_contacts(cd; scaled = false))
        doubled = Daedalus.Data.total_contacts(
            Daedalus.Data.prepare_contacts(cd2; scaled = false))
        @test doubled ≈ 2 .* single
    end

    @testset "daedalus runs with two equal settings" begin
        cd2 = deepcopy(cd)
        cd2.contact_matrix = [cm, cm]
        infection = Daedalus.DataLoader.get_pathogen("sars-cov-2 delta")
        infection.r0 = 2.5
        result = daedalus(cd2, infection, time_end = 100.0, log_rt = false)
        @test length(result.sol.t) == 101
        @test all(isfinite, result.sol.u[end])
    end

    @testset "two equal settings gives same epidemic as one setting" begin
        # beta is calibrated from total_contacts (sum of all settings).
        # With contact_matrix = [cm, cm], total contacts = 2*cm, so beta is
        # halved. In the ODE, weighted_slice_sum! with weights [1, 1] gives
        # 2*scaled_cm; the force of infection (beta/2) * (2*scaled_cm) equals
        # beta * scaled_cm — the same as for a single setting.
        # Therefore, the epidemic trajectory and final size are identical when
        # the same r0 is used.
        cd2 = deepcopy(cd)
        cd2.contact_matrix = [cm, cm]
        infection_1 = Daedalus.DataLoader.get_pathogen("sars-cov-2 delta")
        infection_1.r0 = 2.5
        infection_2 = Daedalus.DataLoader.get_pathogen("sars-cov-2 delta")
        infection_2.r0 = 2.5
        result_1 = daedalus(cd, infection_1, time_end = 300.0, log_rt = false)
        result_2 = daedalus(cd2, infection_2, time_end = 300.0, log_rt = false)
        deaths_1 = last(Daedalus.Outputs.get_values(result_1, "D", 1))
        deaths_2 = last(Daedalus.Outputs.get_values(result_2, "D", 1))
        @test deaths_2 ≈ deaths_1 rtol = 1e-3
    end
end
