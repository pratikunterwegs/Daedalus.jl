using Daedalus
using LinearAlgebra
using Random
using Test

@testset "Dominant eigenvalue computation" begin
    # Test 1: Simple 2x2 matrix with known eigenvalues
    @testset "Simple 2x2 matrix" begin
        A = [3.0 1.0; 1.0 3.0]
        # Eigenvalues are 4.0 and 2.0
        λ_dom = Daedalus.Helpers.dominant_eigenvalue(A)
        λ_exact = maximum(real(eigen(A).values))

        @test isapprox(λ_dom, λ_exact, rtol = 1e-5)
        @test isapprox(λ_dom, 4.0, rtol = 1e-5)
    end

    # Test 2: Identity matrix
    @testset "Identity matrix" begin
        n = 10
        A = Matrix{Float64}(I, n, n)
        λ_dom = Daedalus.Helpers.dominant_eigenvalue(A)

        @test isapprox(λ_dom, 1.0, rtol = 1e-5)
    end

    # Test 3: Diagonal matrix
    @testset "Diagonal matrix" begin
        A = Diagonal([5.0, 3.0, 1.0, 2.0])
        λ_dom = Daedalus.Helpers.dominant_eigenvalue(A)
        λ_exact = maximum(real(eigen(A).values))

        @test isapprox(λ_dom, λ_exact, rtol = 1e-5)
        @test isapprox(λ_dom, 5.0, rtol = 1e-5)
    end

    # Test 4: Random matrix
    @testset "Random 10x10 matrix" begin
        Random.seed!(42)
        A = rand(10, 10)
        λ_dom = Daedalus.Helpers.dominant_eigenvalue(A)
        λ_exact = maximum(real(eigen(A).values))

        @test isapprox(λ_dom, λ_exact, rtol = 1e-4)
    end

    # Test 5: Australia NGM (realistic use case)
    @testset "Australia NGM" begin
        r0 = 1.3
        sigma = 0.217
        p_sigma = 0.867
        epsilon = 0.58
        gamma_Ia = 0.476
        gamma_Is = 0.25

        ngm = Daedalus.Helpers.get_ngm(
            Daedalus.Data.australia_contacts(),
            r0, sigma, p_sigma, epsilon, gamma_Ia, gamma_Is
        )

        λ_dom = Daedalus.Helpers.dominant_eigenvalue(ngm)
        λ_exact = maximum(real(eigen(ngm).values))

        # Should match the exact eigenvalue
        @test isapprox(λ_dom, λ_exact, rtol = 1e-4)

        # Should be approximately r0
        @test isapprox(λ_dom, r0, rtol = 1e-2)
    end

    # Test 6: Full contact matrix (49x49) with susceptible scaling
    @testset "Full 49x49 contact matrix with susceptible scaling" begin
        Random.seed!(123)

        # Get the full contact matrix
        contacts = Daedalus.Data.prepare_contacts(scaled = false)

        # Get NGM parameters
        r0 = 2.5
        sigma = 0.217
        p_sigma = 0.867
        epsilon = 0.58
        gamma_Ia = 0.476
        gamma_Is = 0.25

        ngm = Daedalus.Helpers.get_ngm(
            contacts,
            r0, sigma, p_sigma, epsilon, gamma_Ia, gamma_Is
        )

        # Simulate susceptible proportions (as would happen during simulation)
        p_susc = rand(49) .* 0.8 .+ 0.1  # Random proportions between 0.1 and 0.9
        ngm_susc = ngm .* p_susc

        λ_dom = Daedalus.Helpers.dominant_eigenvalue(ngm_susc)
        λ_exact = maximum(real(eigen(ngm_susc).values))

        @test isapprox(λ_dom, λ_exact, rtol = 1e-4)
    end

    # Test 7: Convergence with different tolerances
    @testset "Convergence tolerance" begin
        Random.seed!(456)
        A = rand(20, 20)

        λ_tight = Daedalus.Helpers.dominant_eigenvalue(A, tol = 1e-8)
        λ_loose = Daedalus.Helpers.dominant_eigenvalue(A, tol = 1e-3)
        λ_exact = maximum(real(eigen(A).values))

        # Both should be close to exact
        @test isapprox(λ_tight, λ_exact, rtol = 1e-4)
        @test isapprox(λ_loose, λ_exact, rtol = 1e-2)

        # Tight tolerance should be more accurate
        @test abs(λ_tight - λ_exact) <= abs(λ_loose - λ_exact)
    end

    # Test 8: Warm start initialization
    @testset "Warm start with v_init" begin
        Random.seed!(789)
        A = rand(15, 15)

        # First computation
        λ1 = Daedalus.Helpers.dominant_eigenvalue(A)

        # Use a good initial vector (close to the eigenvector)
        v_init = randn(15)
        v_init = v_init / norm(v_init)

        λ2 = Daedalus.Helpers.dominant_eigenvalue(A, v_init = v_init)

        # Both should give the same result
        @test isapprox(λ1, λ2, rtol = 1e-5)

        # Both should match exact
        λ_exact = maximum(real(eigen(A).values))
        @test isapprox(λ2, λ_exact, rtol = 1e-4)
    end

    # Test 9: Non-negative matrix (Perron-Frobenius property)
    @testset "Non-negative matrix (epidemiological relevance)" begin
        Random.seed!(101)
        # NGM-like non-negative matrix
        A = abs.(rand(20, 20))

        λ_dom = Daedalus.Helpers.dominant_eigenvalue(A)
        λ_exact = maximum(real(eigen(A).values))

        # Should be positive for non-negative matrix
        @test λ_dom > 0

        # Should match exact value
        @test isapprox(λ_dom, λ_exact, rtol = 1e-4)
    end

    # Test 10: Zero matrix edge case
    @testset "Zero matrix edge case" begin
        A = zeros(5, 5)
        λ_dom = Daedalus.Helpers.dominant_eigenvalue(A)

        @test isapprox(λ_dom, 0.0, atol = 1e-10)
    end
end

@testset "Rt calculation in make_rt_logger" begin
    @testset "Rt logging with power iteration" begin
        # Run a short simulation with Rt logging enabled
        result = Daedalus.daedalus(
            r0 = 2.0,
            time_end = 50.0,
            increment = 1.0,
            log_rt = true
        )

        # Extract Rt values
        iRt = Daedalus.Constants.get_indices("Rt")
        rt_values = [u[iRt] for u in result.sol.u]

        # Check that Rt values are reasonable
        @test all(rt_values .> 0)  # Rt should be positive
        @test rt_values[1]≈2.0 atol=0.1  # Initial Rt should be close to r0

        # Check that Rt decreases over time (as susceptibles are depleted)
        # This is expected behavior for an epidemic without interventions
        @test rt_values[end] < rt_values[1]

        # Verify no NaN or Inf values
        @test all(isfinite.(rt_values))
    end

    @testset "Rt with NPI intervention" begin
        # Create a simple NPI for testing
        npi = Daedalus.DaedalusStructs.Npi(20000.0, (coef = 0.7,))

        result = Daedalus.daedalus(
            r0 = 3.0,
            time_end = 100.0,
            increment = 1.0,
            npi = npi,
            log_rt = true
        )

        # Extract Rt values
        iRt = Daedalus.Constants.get_indices("Rt")
        rt_values = [u[iRt] for u in result.sol.u]

        # Check basic properties
        @test all(rt_values .> 0)
        @test all(isfinite.(rt_values))
    end
end
