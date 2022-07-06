#
# Date created: 2022-07-05
# Author: aradclif
#
#
############################################################################################

function generate(K::Vector{Int}, V::Vector{T}) where {T<:AbstractFloat}
    n = length(K)
    u = rand()
    # j = unsafe_trunc(Int, muladd(u, n, 1))
    j = unsafe_trunc(Int, u[i] * n) + 1
    u < V[j] ? j : K[j]
end

function generate!(A::AbstractArray{<:Integer}, u::AbstractArray{Float64}, K::Vector{Int}, V::Vector{T}) where {T<:AbstractFloat}
    n = length(K)
    checkbounds(V, n)
    rand!(u)
    @inbounds for i ∈ eachindex(A, u)
        # j = unsafe_trunc(Int, muladd(u[i], n, 1.0))
        j = unsafe_trunc(Int, u[i] * n) + 1
        # A[i] = u[i] < V[j] ? j : K[j]
        A[i] = ifelse(u[i] < V[j], j, K[j])
    end
    A
end
# Dense arrays: seemingly safe w.r.t. memory deps
function generate!(A::Array{<:Integer}, u::Array{Float64}, K::Vector{Int}, V::Vector{T}) where {T<:AbstractFloat}
    n = length(K)
    checkbounds(V, n)
    rand!(u)
    @inbounds @simd ivdep for i ∈ eachindex(A, u)
        # j = unsafe_trunc(Int, muladd(u[i], n, 1.0))
        j = unsafe_trunc(Int, u[i] * n) + 1
        # A[i] = u[i] < V[j] ? j : K[j]
        A[i] = ifelse(u[i] < V[j], j, K[j])
    end
    A
end

generate!(A::AbstractArray{<:Integer}, K::Vector{Int}, V::Vector{T}) where {T<:AbstractFloat} = generate!(A, similar(A, Float64), K, V)

generate(K::Vector{Int}, V::Vector{T}, dims::Vararg{Int, N}) where {T<:AbstractFloat} where {N} =
    generate!(Array{Int}(undef, dims), K, V)

################
# Why bother with @turbo? Well, unless one can guarantee that @simd ivdep is safe,
# @turbo will exhibit superior performance.
vgenerate(K::Vector{Int}, V::Vector{<:AbstractFloat}) = generate(K, V)

function vgenerate!(A::AbstractArray{<:Integer}, u::AbstractArray{Float64}, K::Vector{Int}, V::Vector{T}) where {T<:AbstractFloat}
    n = length(K)
    checkbounds(V, n)
    rand!(u)
    @turbo for i ∈ eachindex(A, u)
        # j = unsafe_trunc(Int, muladd(u[i], n, 1.0))
        j = unsafe_trunc(Int, u[i] * n) + 1
        A[i] = ifelse(u[i] < V[j], j, K[j])
    end
    A
end
vgenerate!(A::AbstractArray{<:Integer}, K::Vector{Int}, V::Vector{T}) where {T<:AbstractFloat} =
    vgenerate!(A, similar(A, Float64), K, V)

vgenerate(K::Vector{Int}, V::Vector{T}, dims::Vararg{Int, N}) where {T<:AbstractFloat} where {N} =
    vgenerate!(Array{Int}(undef, dims), K, V)


################
# Equal probability case admits an optimized form:
# U ~ Uniform(0,1)
# j = ⌊nU + 1⌋; return j

function generate(n::Int)
    n > 0 || throw(ArgumentError("n must be > 0"))
    unsafe_trunc(Int, rand() * n) + 1
end

function generate!(A::AbstractArray{<:Integer}, u::AbstractArray{Float64}, n::Int)
    n > 0 || throw(ArgumentError("n must be > 0"))
    n == 1 && return fill!(A, 1)
    rand!(u)
    @inbounds for i ∈ eachindex(A, u)
        A[i] = unsafe_trunc(Int, u[i] * n) + 1
    end
    A
end
generate!(A::AbstractArray, n::Int) = generate!(A, similar(A, Float64), n)

generate(n::Int, dims::NTuple{N, Int}) where {N} = generate!(Array{Int}(undef, dims), n)
generate(n::Int, dims::Vararg{Int, N}) where {N} = generate(n, dims)

################
vgenerate(n::Int) = generate(n)
function vfill!(A::AbstractArray, v::Real)
    @turbo for i ∈ eachindex(A)
        A[i] = v
    end
    A
end

function vgenerate!(A::AbstractArray{<:Integer}, u::AbstractArray{Float64}, n::Int)
    n > 0 || throw(ArgumentError("n must be > 0"))
    n == 1 && return vfill!(A, 1)
    rand!(u)
    @turbo for i ∈ eachindex(A, u)
        A[i] = unsafe_trunc(Int, u[i] * n) + 1
    end
    A
end
vgenerate!(A::AbstractArray, n::Int) = vgenerate!(A, similar(A, Float64), n)

vgenerate(n::Int, dims::NTuple{N, Int}) where {N} = vgenerate!(Array{Int}(undef, dims), n)
vgenerate(n::Int, dims::Vararg{Int, N}) where {N} = vgenerate(n, dims)
