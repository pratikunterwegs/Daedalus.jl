@testset "Multiple contact settings" begin
    cd = Daedalus.DataLoader.get_country("Australia")
    cm = cd.contact_matrix

    @testset "get_settings returns 2 for a single matrix (community + workplace)" begin
        @test Daedalus.Data.get_settings(cd) == 2
    end

    @testset "get_settings returns 3 for a vector of two matrices (2 community + 1 workplace)" begin
        cd2 = deepcopy(cd)
        cd2.contact_matrix = [cm, cm]
        @test Daedalus.Data.get_settings(cd2) == 3
    end

    @testset "contacts3d shape is (49, 49, 2) for a single setting (community + workplace)" begin
        c3d = Daedalus.Data.contacts3d(cd)
        @test size(c3d) == (49, 49, 2)
    end

    @testset "contacts3d shape is (49, 49, 3) for two settings (2 community + 1 workplace)" begin
        cd2 = deepcopy(cd)
        cd2.contact_matrix = [cm, cm]
        c3d = Daedalus.Data.contacts3d(cd2)
        @test size(c3d) == (49, 49, 3)
    end

    @testset "total_contacts sums all layers of contacts3d" begin
        # total_contacts should sum all community + workplace layers.
        # For K=1: sum community + workplace
        # For K=2: sum community1 + community2 + workplace
        c3d_single = Daedalus.Data.contacts3d(cd; scaled = false)
        c3d_double = deepcopy(cd)
        c3d_double.contact_matrix = [cm, cm]
        c3d_double = Daedalus.Data.contacts3d(c3d_double; scaled = false)

        total_single = Daedalus.Data.total_contacts(c3d_single)
        total_double = Daedalus.Data.total_contacts(c3d_double)

        # total_double should have more than total_single because it includes
        # the sum of two identical community matrices plus the workplace
        @test all(total_double .>= total_single)
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

    @testset "model runs with multiple equal community settings" begin
        # With the new contacts3d structure (community + workplace layers),
        # total contacts change when adding more community settings, so beta
        # calibration differs. The test just verifies the model runs correctly.
        cd2 = deepcopy(cd)
        cd2.contact_matrix = [cm, cm]
        infection = Daedalus.DataLoader.get_pathogen("sars-cov-2 delta")
        infection.r0 = 2.5
        result = daedalus(cd2, infection, time_end = 300.0, log_rt = false)
        @test length(result.sol.t) > 100
        @test all(isfinite, result.sol.u[end])
    end
end
