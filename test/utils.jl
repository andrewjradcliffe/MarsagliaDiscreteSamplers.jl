# Basic utils tests
@testset "vfindextrema(f, domain)" begin
    @test vfindextrema(-, 1:10) == (findmin(-, 1:10), findmax(-, 1:10))
    @test vfindextrema(identity, [1, NaN, 3]) == ((1.0, 1), (3.0, 3))
    @test vfindextrema(identity, [1, 3, NaN]) == ((1.0, 1), (3.0, 2))
    @test vfindextrema(cos, 0:π/2:2π) == (findmin(cos, 0:π/2:2π), findmax(cos, 0:π/2:2π))
end

@testset "vfill!" begin
    for T ∈ (Float16, Float32, Float64, Int16, Int32, Int64)
        x = rand(T)
        y = Vector{T}(undef, 10)
        z = Vector{T}(undef, 10)
        @test fill!(y, x) == vfill!(z, x)
    end
end
