#
# Date created: 2022-07-10
# Author: aradclif
#
#
############################################################################################
# for use in: Comparison: non-chunked vs. chunked

function option1!(A::AbstractArray, u::Array{Float64})
    rand!(u)
    @inbounds @simd ivdep for i ∈ eachindex(A, u)
        A[i] = u[i]
    end
    A
end
function voption1!(A::AbstractArray, u::Array{Float64})
    rand!(u)
    @turbo for i ∈ eachindex(A, u)
        A[i] = u[i]
    end
    A
end

function option2!(A::AbstractArray, u::Array{Float64})
    n = length(A)
    c = length(u)
    q, r = divrem(n, c)
    i₀ = firstindex(A)
    f = i₀
    l = f - 1
    k = 1
    while k ≤ q
        f = l + 1
        l += c
        j = f:l
        rand!(u)
        @inbounds @simd for i ∈ eachindex(u)
            A[j[i]] = u[i]
        end
        k += 1
    end
    if r != 0
        f = l + 1
        j = f:n
        rand!(u)
        @inbounds @simd for i ∈ eachindex(j)
            A[j[i]] = u[i]
        end
    end
    A
end

function voption2!(A::AbstractArray, u::Array{Float64})
    n = length(A)
    c = length(u)
    q, r = divrem(n, c)
    i₀ = firstindex(A)
    f = i₀
    l = f - 1
    k = 1
    while k ≤ q
        f = l + 1
        l += c
        j = f:l
        rand!(u)
        @turbo for i ∈ eachindex(u)
            A[j[i]] = u[i]
        end
        k += 1
    end
    if r != 0
        f = l + 1
        j = f:n
        rand!(u)
        @turbo for i ∈ eachindex(j)
            A[j[i]] = u[i]
        end
    end
    A
end

function option3!(A::AbstractArray, u::Array{Float64})
    n = length(A)
    c = length(u)
    q, r = divrem(n, c)
    i₀ = firstindex(A)
    f = i₀
    l = f - 1
    k = 1
    while k ≤ q
        f = l + 1
        l += c
        j = f:l
        rand!(u)
        @inbounds for (i, i′) ∈ enumerate(j)
            A[i′] = u[i]
        end
        k += 1
    end
    if r != 0
        f = l + 1
        j = f:n
        rand!(u)
        @inbounds for (i, i′) ∈ enumerate(j)
            A[i′] = u[i]
        end
    end
    A
end

function option4!(A::AbstractArray, u::Array{Float64})
    n = length(A)
    c = length(u)
    q, r = divrem(n, c)
    i₀ = firstindex(A)
    f = i₀
    l = f - 1
    k = 1
    while k ≤ q
        f = l + 1
        l += c
        j = f
        rand!(u)
        @inbounds for i ∈ eachindex(u)
            A[j] = u[i]
            j += 1
        end
        k += 1
    end
    if r != 0
        f = l + 1
        j = f:n
        rand!(u)
        i = 1
        @inbounds for j ∈ f:n
            A[j] = u[i]
            i += 1
        end
    end
    A
end

function generate_option2!(A::AbstractArray{<:Integer}, u::Array{Float64}, K::Vector{Int}, V::Vector{T}) where {T<:AbstractFloat}
    n = length(A)
    c = length(u)
    q, r = divrem(n, c)
    N = length(K)
    checkbounds(V, N)
    N == 1 && return fill!(A, @inbounds K[1])
    i₀ = firstindex(A)
    f = i₀
    l = f - 1
    k = 1
    while k ≤ q
        f = l + 1
        l += c
        j′ = f:l
        rand!(u)
        @inbounds @simd ivdep for i ∈ eachindex(u)
            j = unsafe_trunc(Int, u[i] * N) + 1
            A[j′[i]] = ifelse(u[i] < V[j], j, K[j])
        end
        k += 1
    end
    if r != 0
        f = l + 1
        j′ = f:n
        rand!(u)
        @inbounds @simd ivdep for i ∈ eachindex(j′)
            j = unsafe_trunc(Int, u[i] * N) + 1
            A[j′[i]] = ifelse(u[i] < V[j], j, K[j])
        end
    end
    A
end
