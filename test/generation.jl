#
Random.seed!(0x123456789abcdef)

# Helper functions
function unsafe_countcategory!(v::AbstractArray, A::AbstractArray)
    @inbounds for i ∈ eachindex(A)
        v[A[i]] += 1
    end
    v
end
unsafe_countcategory(A, n::Int) = unsafe_countcategory!(zeros(Int, n), A)
countcategory(A) = unsafe_countcategory(A, maximum(A))

in13(x) = 1 ≤ x ≤ 3
_normalize1!(A) = A .*= inv(sum(A))

@testset "generate: general case" begin
    p = [2/15, 7/15, 6/15]
    K, V = sqhist(p)
    A = generate(K, V, 1,2,3);
    mn, mx = extrema(A)
    @test mn ≥ 1 && mx ≤ 3
    @test all(1 ≤ generate(K, V) ≤ 3 for _ = 1:10^4)
    A = generate(K, V, 10^7);
    t = unsafe_countcategory(A, 3)
    @test all((≈).(t ./ 10^7, p, atol=1e-3))
    @test all(in13, A)
    u = Vector{Float64}(undef, 10)
    @test_throws DimensionMismatch generate!(A, u, K, V)
    resize!(u, length(A))
    generate!(A, u, K, V)
    t = unsafe_countcategory(A, 3)
    @test all((≈).(t ./ 10^7, p, atol=1e-3))
    @test all(in13, A)
    # multidimensional initialization
    @test all(in13, A)
    for dims ∈ (0, 1, 2, (), (1,), (1,2), (1,1,3), (3,2,4), (1,2,3,4))
        A = generate(K, V, dims)
        @test all(in13, A)
    end
    # views
    M = zeros(Int, 7, 5)
    for i_1 = 1:7
        Mv = view(M, i_1, :)
        resize!(u, length(Mv))
        generate!(Mv, u, K, V)
        @test all(!=(0), Mv)
        @test all(in13, Mv)
        Mv .= 0
    end
    for i_2 = 1:5
        Mv = view(M, :, i_2)
        resize!(u, length(Mv))
        generate!(Mv, u, K, V)
        @test all(!=(0), Mv)
        @test all(in13, Mv)
        Mv .= 0
    end
    for (s_1, s_2) ∈ ((1:3, 1:2), (4:6, 3:5), (2:5, 1:3), (1:2:4, 1:2), (1:3:7, 1:2:5))
        Mv = view(M, s_1, s_2)
        u = similar(Mv, Float64)
        generate!(Mv, u, K, V)
        @test all(!=(0), Mv)
        @test all(in13, Mv)
        Mv .= 0
    end
    # throws based on K, V
    for n = 0:2
        resize!(V, n)
        @test_throws BoundsError generate(K, V)
        @test_throws BoundsError generate(K, V, 1, 2)
        for dims ∈ (0, 1, 2, (), (1,), (1,2))
            @test_throws BoundsError generate(K, V, dims)
        end
        @test_throws BoundsError generate!(A, K, V)
    end
    resize!(V, 3)
    for n = (0, 4)
        resize!(K, n)
        @test_throws BoundsError generate(K, V)
        @test_throws BoundsError generate(K, V, 1, 2)
        for dims ∈ (0, 1, 2, (), (1,), (1,2))
            @test_throws BoundsError generate(K, V, dims)
        end
        @test_throws BoundsError generate!(A, K, V)
    end
end

@testset "generate: SqHist interface" begin
    p = [2/15, 7/15, 6/15]
    K, V = sqhist(p)
    x = SqHist(p)
    A = generate(x, 1,2,3);
    mn, mx = extrema(A)
    @test mn ≥ 1 && mx ≤ 3
    @test all(1 ≤ generate(x) ≤ 3 for _ = 1:10^4)
    A = generate(x, 10^7);
    t = unsafe_countcategory(A, 3)
    @test all((≈).(t ./ 10^7, p, atol=1e-3))
    @test all(in13, A)
    u = Vector{Float64}(undef, 10)
    @test_throws DimensionMismatch generate!(A, u, x)
    resize!(u, length(A))
    generate!(A, u, x)
    t = unsafe_countcategory(A, 3)
    @test all((≈).(t ./ 10^7, p, atol=1e-3))
    @test all(in13, A)
    # multidimensional initialization
    @test all(in13, A)
    for dims ∈ (0, 1, 2, (), (1,), (1,2), (1,1,3), (3,2,4), (1,2,3,4))
        A = generate(x, dims)
        @test all(in13, A)
    end
end

@testset "generate: equal probability" begin
    K, V = [1,2,3], [1/3, 2/3, 3/3]
    p = [1/3, 1/3, 1/3]
    n = 3
    A = generate(n, 1,2,3);
    mn, mx = extrema(A)
    @test mn ≥ 1 && mx ≤ 3
    @test all(1 ≤ generate(n) ≤ 3 for _ = 1:10^4)
    A = generate(n, 10^7);
    t = unsafe_countcategory(A, 3)
    @test all((≈).(t ./ 10^7, p, atol=1e-3))
    @test all(in13, A)
    u = Vector{Float64}(undef, 10)
    @test_throws DimensionMismatch generate!(A, u, n)
    resize!(u, length(A))
    generate!(A, u, n)
    t = unsafe_countcategory(A, 3)
    @test all((≈).(t ./ 10^7, p, atol=1e-3))
    @test all(in13, A)
    # multidimensional initialization
    @test all(in13, A)
    for dims ∈ (0, 1, 2, (), (1,), (1,2), (1,1,3), (3,2,4), (1,2,3,4))
        A = generate(n, dims)
        @test all(in13, A)
    end
    # views
    M = zeros(Int, 7, 5)
    for i_1 = 1:7
        Mv = view(M, i_1, :)
        resize!(u, length(Mv))
        generate!(Mv, u, n)
        @test all(!=(0), Mv)
        @test all(in13, Mv)
        Mv .= 0
    end
    for i_2 = 1:5
        Mv = view(M, :, i_2)
        resize!(u, length(Mv))
        generate!(Mv, u, n)
        @test all(!=(0), Mv)
        @test all(in13, Mv)
        Mv .= 0
    end
    for (s_1, s_2) ∈ ((1:3, 1:2), (4:6, 3:5), (2:5, 1:3), (1:2:4, 1:2), (1:3:7, 1:2:5))
        Mv = view(M, s_1, s_2)
        u = similar(Mv, Float64)
        generate!(Mv, u, n)
        @test all(!=(0), Mv)
        @test all(in13, Mv)
        Mv .= 0
    end
    # throws based on n
    @test_throws ArgumentError generate(0)
    @test_throws ArgumentError generate(0, 1, 2)
    for dims ∈ (0, 1, 2, (), (1,), (1,2))
        @test_throws ArgumentError generate(0, dims)
    end
    @test_throws ArgumentError generate!(A, 0)
end

@testset "generate: SqHistEquiprobable interface" begin
    K, V = [1,2,3], [1/3, 2/3, 3/3]
    p = [1/3, 1/3, 1/3]
    n = 3
    x = SqHistEquiprobable(n)
    A = generate(x, 1,2,3);
    mn, mx = extrema(A)
    @test mn ≥ 1 && mx ≤ 3
    @test all(1 ≤ generate(x) ≤ 3 for _ = 1:10^4)
    A = generate(x, 10^7);
    t = unsafe_countcategory(A, 3)
    @test all((≈).(t ./ 10^7, p, atol=1e-3))
    @test all(in13, A)
    u = Vector{Float64}(undef, 10)
    @test_throws DimensionMismatch generate!(A, u, x)
    resize!(u, length(A))
    generate!(A, u, x)
    t = unsafe_countcategory(A, 3)
    @test all((≈).(t ./ 10^7, p, atol=1e-3))
    @test all(in13, A)
    # multidimensional initialization
    @test all(in13, A)
    for dims ∈ (0, 1, 2, (), (1,), (1,2), (1,1,3), (3,2,4), (1,2,3,4))
        A = generate(x, dims)
        @test all(in13, A)
    end
end

@testset "sqhist, generate: numerical stability" begin
    n_samples = 10^8
    c = inv(n_samples)
    A = Vector{Int}(undef, n_samples)
    u = Vector{Float64}(undef, n_samples)
    for i = 1:15
        n = (1 << i)
        p = fill(1/n, n)
        K, V = sqhist(p)
        @test 1 ≤ generate(K, V) ≤ n
        generate!(A, u, K, V)
        t = unsafe_countcategory(A, n);
        @test all(i -> ≈(t[i] * c, p[i], atol=1e-3), 1:n)
        # @test all((≈).(t .* c, p, atol=1e-3))
        @test all(∈(1:n), A)
        for _ = 1:10
            rand!(p)
            _normalize1!(p)
            K, V = sqhist(p)
            @test 1 ≤ generate(K, V) ≤ n
            generate!(A, u, K, V)
            t = unsafe_countcategory(A, n);
            @test all(i -> ≈(t[i] * c, p[i], atol=1e-3), 1:n)
            # @test all((≈).(t .* c, p, atol=1e-3))
            @test all(∈(1:n), A)
        end
    end
end
