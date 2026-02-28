@testset "weighted_slice_sum!" begin
    wss! = Daedalus.Helpers.weighted_slice_sum!

    M, K, N = 5, 6, 3
    X = rand(M, K, N)
    result = zeros(M, K)

    @testset "unit weights equal sum over third dimension" begin
        v = ones(N)
        wss!(X, v, result)
        @test result ≈ dropdims(sum(X, dims = 3), dims = 3)
    end

    @testset "single non-zero weight selects one slice" begin
        for i in 1:N
            v = zeros(N)
            v[i] = 1.0
            wss!(X, v, result)
            @test result ≈ X[:, :, i]
        end
    end

    @testset "scalar weight scales the result" begin
        c = 3.7
        v = zeros(N)
        v[1] = c
        wss!(X, v, result)
        @test result ≈ c .* X[:, :, 1]
    end

    @testset "zero weights produce zero result" begin
        v = zeros(N)
        wss!(X, v, result)
        @test all(iszero, result)
    end

    @testset "result is overwritten, not accumulated" begin
        fill!(result, 999.0)
        v = zeros(N)
        v[1] = 1.0
        wss!(X, v, result)
        @test result ≈ X[:, :, 1]
    end

    @testset "matches explicit reference loop" begin
        v = rand(N)
        wss!(X, v, result)
        ref = zeros(M, K)
        for i in 1:N
            ref .+= X[:, :, i] .* v[i]
        end
        @test result ≈ ref
    end
end
