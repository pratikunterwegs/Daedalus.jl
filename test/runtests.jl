using Daedalus
using LinearAlgebra
using Random
using Test

@testset "Daedalus.jl" begin
    # Write your tests here.
    @testset "DAEDALUS model" begin
        try
            daedalus(r0 = 5.0, time_end = 100.0)
            @test true
        catch e
            @test false
        end
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

        beta = Daedalus.Helpers.get_beta(
            Daedalus.Data.australia_contacts(),
            r0, sigma, p_sigma, epsilon, gamma_Ia, gamma_Is
        )
        @test typeof(beta) == Float64

        ngm = Daedalus.Helpers.get_ngm(
            Daedalus.Data.australia_contacts(),
            r0, sigma, p_sigma, epsilon, gamma_Ia, gamma_Is
        )
        lambda = maximum(eigen(ngm).values)

        @test lambda ≈ r0
    end

    # Include eigenvalue tests
    include("test_eigenvalue.jl")
end
