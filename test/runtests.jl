using MarsagliaDiscreteSamplers
using Random
using Test

const tests = [
    "squarehistogram.jl",
    "generation.jl"
    "utils.jl"
]
for t in tests
    @testset "Test $t" begin
        include(t)
    end
end
