module MarsagliaDiscreteSamplers

using Random, Test, BenchmarkTools
using LoopVectorization

export sqhist!, sqhist

include("squarehistogram.jl")

end
