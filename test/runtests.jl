using MarsagliaDiscreteSamplers
using Random
using Test

const tests = [
    "squarehistogram.jl",
]
for t in tests
    @testset "Test $t" begin
        include(t)
    end
end
