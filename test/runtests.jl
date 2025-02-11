using Daedalus
using Test

@testset "Daedalus.jl" begin
    # Write your tests here.
    @testset "DAEDALUS model" begin
        try
            daedalus()
            @test true
        catch e
            @test false
        end
    end

end
