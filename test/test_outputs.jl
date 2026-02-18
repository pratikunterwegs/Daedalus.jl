@testset "Outputs.get_values" begin
    output = daedalus(time_end=100.0)
    tmax = 100

    @testset "daily values (timebin=1)" begin
        vals = Daedalus.Outputs.get_values(output, "H", 1)
        @test length(vals) == tmax + 1
        @test eltype(vals) == Float64
        @test all(vals .>= 0)
    end

    @testset "binned values sum matches daily total" begin
        daily = Daedalus.Outputs.get_values(output, "H", 1)
        binned = Daedalus.Outputs.get_values(output, "H", 90)
        # 100 days: bin 1 = days 0–89, bin 2 = days 90–100
        @test length(binned) == 2
        @test all(binned .>= 0)
        @test sum(binned) ≈ sum(daily)
    end

    @testset "bin alignment when tmax divisible by timebin" begin
        output90 = daedalus(time_end=90.0)
        daily = Daedalus.Outputs.get_values(output90, "H", 1)
        binned = Daedalus.Outputs.get_values(output90, "H", 90)
        @test length(binned) == 1
        @test sum(binned) ≈ sum(daily)
    end

    @testset "all standard compartments return correct length" begin
        for comp in ["S", "E", "Is", "Ia", "H", "R", "D"]
            vals = Daedalus.Outputs.get_values(output, comp, 1)
            @test length(vals) == tmax + 1
        end
    end

    @testset "Rt compartment" begin
        rt = Daedalus.Outputs.get_values(output, "Rt", 1)
        @test length(rt) == tmax + 1
        @test eltype(rt) == Float64
        @test all(rt .>= 0)
    end

    @testset "invalid compartment throws" begin
        @test_throws ErrorException Daedalus.Outputs.get_values(output, "X", 1)
        @test_throws ErrorException Daedalus.Outputs.get_values(output, "", 1)
    end
end
