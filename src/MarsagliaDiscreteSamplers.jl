module MarsagliaDiscreteSamplers

import Base: unsafe_trunc

using Random
using LoopVectorization

export SqHist, SqHistEquiprobable

export sqhist!, sqhist, generate!, generate, vgenerate!, vgenerate

include("utils.jl")
include("squarehistogram.jl")
include("generation.jl")

end
