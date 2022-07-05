module MarsagliaDiscreteSamplers

using Random, Test, BenchmarkTools
using LoopVectorization
using VectorizedReduction

export sqhist!, sqhist

include("squarehistogram.jl")

end
