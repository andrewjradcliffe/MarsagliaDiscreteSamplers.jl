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
