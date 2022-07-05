# Tests of square histogram construction

@testset "Sqhist: Robin Hood alias table" begin
    p = [2/15, 7/15, 6/15]
    K, V = sqhist_robinhood(p)
    @test K == [2,3,3]
    @test V == [2/15, nextfloat(9/15), 15/15]
    p = [2/15, 6/15, 7/15]
    K, V = sqhist_robinhood(p)
    @test K == [3,2,2]
    @test V == [2/15, 2/3, 14/15]
    p = [7/15, 2/15, 6/15]
    K, V = sqhist_robinhood(p)
    @test K == [3,1,3]
    @test V == [nextfloat(4/15), 7/15, 1.0]
    #
    p = [.21, .18, .26, .17, .18]
    K, V = sqhist_robinhood(p)
    @test K == [3, 3, 3, 3, 1]
    @test V == [0.18999999999999997, 0.38, 0.6000000000000001, 0.7700000000000001, 0.98]
    #
    K′, V′ = Vector{Int}(), Vector{Float64}()
    q = Vector{Float64}()
    n = length(p)
    resize!(K′, n); resize!(V′, n); resize!(q, n)
    sqhist_robinhood!(K′, V′, q, p)
    @test K == K′
    @test V == V′
    #
    resize!(V′, 0)
    @test_throws DimensionMismatch sqhist_robinhood!(K′, V′, q, p)
    #
    @testset "numerical stability" begin
        # equal probability
        for i = 1:10
            n = 1 << i
            p = fill(1/n, n)
            K, V = sqhist_robinhood(p)
            @test K == collect(1:n)
            @test V == collect(1/n:1/n:1.0)
        end
        # cases on the verge of instability
        # p₁ = 0.999
        # n = 1 << i
        # p = [p₁; fill((1.0 - p₁) / n, n)]
        # K, V = sqhist_robinhood(p)
    end
end

@testset "Sqhist: non-Robin Hood alias table" begin
    p = [2/15, 7/15, 6/15]
    K, V = sqhist(p)
    @test K == [3,2,2]
    @test V == [2/15, 2/3, 13/15]
    p = [2/15, 6/15, 7/15]
    K, V = sqhist(p)
    @test K == [3,2,2]
    @test V == [2/15, 2/3, 14/15]
    p = [7/15, 2/15, 6/15]
    K, V = sqhist(p)
    @test K == [1,3,1]
    @test V == [1/3, 7/15, 13/15]
    #
    p = [8/15, 3/15, 4/15]
    K, V = sqhist(p)
    @test K == [1,1,1]
    @test V == [1/3, 8/15, 14/15]
    p = fill(1/3, 3)
end
