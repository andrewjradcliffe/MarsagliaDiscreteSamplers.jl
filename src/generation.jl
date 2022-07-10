#
# Date created: 2022-07-05
# Author: aradclif
#
#
############################################################################################

"""
    generate(K::Vector{Int}, V::Vector{<:AbstractFloat})

Generate a random category from the discrete distribution which corresponds to the
squared histogram held in the category aliases, `K`, and division points, `V`.
"""
function generate(K::Vector{Int}, V::Vector{T}) where {T<:AbstractFloat}
    n = length(K)
    checkbounds(V, n)
    u = rand()
    j = unsafe_trunc(Int, u * n) + 1
    @inbounds u < V[j] ? j : K[j]
end

"""
    generate!(A, u::AbstractArray{Float64}, K::Vector{Int}, V::Vector{<:AbstractFloat})
    generate!(A, K::Vector{Int}, V::Vector{<:AbstractFloat})

Populate the array `A` with random categories drawn from the discrete distribution which
corresponds to the squared histogram `K` and `V`. If provided, `u` is used as temporary storage
for uniform(0,1) draws; it must be an array which has equal size as `A`.
If `u` is not provided, it will be allocated; repeated callers are encouraged to
pre-allocate this temporary storage once in order to avoid unnecessary allocations.
"""
function generate!(A::AbstractArray{<:Integer}, u::AbstractArray{Float64}, K::Vector{Int}, V::Vector{T}) where {T<:AbstractFloat}
    n = length(K)
    checkbounds(V, n)
    n == 1 && return fill!(A, @inbounds K[1])
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
    n == 1 && return fill!(A, @inbounds K[1])
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

"""
    generate(K, V, dims::Tuple)
    generate(K, V, dims::Integer...)

Generate an array of random categories from the discrete distribution given the
squared histogram `K` and `V`.

See also: [`generate!`](@ref)

# Examples
```julia-repl
julia> p = [2/15, 7/15, 6/15]; K, V = sqhist(p);

julia> generate(K, V, 1,2,3)
1×2×3 Array{Int64, 3}:
[:, :, 1] =
 1  2

[:, :, 2] =
 2  1

[:, :, 3] =
 3  2
```
"""
generate(K::Vector{Int}, V::Vector{<:AbstractFloat}, dims::NTuple{N, Integer}) where {N} =
    generate!(Array{Int}(undef, dims), K, V)
generate(K::Vector{Int}, V::Vector{<:AbstractFloat}, dims::Vararg{Integer, N}) where {N} =
    generate(K, V, dims)

################
# Why bother with @turbo? Well, unless one can guarantee that @simd ivdep is safe,
# @turbo will exhibit superior performance.

"""
    vgenerate(K::Vector{Int}, V::Vector{<:AbstractFloat})

Generate a random category from the discrete distribution which corresponds to the
squared histogram held in the category aliases, `K`, and division points, `V`.
"""
vgenerate(K::Vector{Int}, V::Vector{<:AbstractFloat}) = generate(K, V)

"""
    vgenerate!(A, u::AbstractArray{Float64}, K::Vector{Int}, V::Vector{<:AbstractFloat})
    vgenerate!(A, K::Vector{Int}, V::Vector{<:AbstractFloat})

Populate the array `A` with random categories drawn from the discrete distribution which
corresponds to the squared histogram `K` and `V`. If provided, `u` is used as temporary storage
for uniform(0,1) draws; it must be an array which has equal size as `A`.
If `u` is not provided, it will be allocated; repeated callers are encouraged to
pre-allocate this temporary storage once in order to avoid unnecessary allocations.
"""
function vgenerate!(A::AbstractArray{<:Integer}, u::AbstractArray{Float64}, K::Vector{Int}, V::Vector{T}) where {T<:AbstractFloat}
    n = length(K)
    checkbounds(V, n)
    n == 1 && return vfill!(A, @inbounds K[1])
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


"""
    vgenerate(K, V, dims::Tuple)
    vgenerate(K, V, dims::Integer...)

Generate an array of random categories from the discrete distribution given the
squared histogram `K` and `V`.

See also: [`generate!`](@ref)

# Examples
```julia-repl
julia> p = [2/15, 7/15, 6/15]; K, V = sqhist(p);

julia> vgenerate(K, V, 1,2,3)
1×2×3 Array{Int64, 3}:
[:, :, 1] =
 2  1

[:, :, 2] =
 2  3

[:, :, 3] =
 3  2
```
"""
vgenerate(K::Vector{Int}, V::Vector{<:AbstractFloat}, dims::NTuple{N, Integer}) where {N} =
    vgenerate!(Array{Int}(undef, dims), K, V)
vgenerate(K::Vector{Int}, V::Vector{<:AbstractFloat}, dims::Vararg{Integer, N}) where {N} =
    vgenerate(K, V, dims)

################
# Equal probability case admits an optimized form:
# U ~ Uniform(0,1)
# j = ⌊nU + 1⌋; return j

"""
    generate(n::Int)

Generate a random category from the uniform discrete distribution which has
support on `1` through `n` (inclusive).
"""
function generate(n::Int)
    n > 0 || throw(ArgumentError("n must be > 0"))
    unsafe_trunc(Int, rand() * n) + 1
end

"""
    generate!(A, u::AbstractArray{Float64}, n::Int)
    generate!(A, n::Int)

Populate the array `A` with random categories drawn from the uniform discrete distribution
which has support on `1` through `n` (inclusive). If provided, `u` is used as temporary storage
for uniform(0,1) draws; it must be an array which has equal size as `A`.
If `u` is not provided, it will be allocated; repeated callers are encouraged to
pre-allocate this temporary storage once in order to avoid unnecessary allocations.
"""
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

"""
    generate(n::Int, dims::Tuple)
    generate(n::Int, dims::Integer...)

Generate an array of random categories from the uniform discrete distribution which has
support on `1` through `n` (inclusive).

See also: [`generate!`](@ref)

# Examples
```julia-repl
julia> generate(10, 3, 2)
3×2 Matrix{Int64}:
 10  1
  2  3
  2  1
```
"""
generate(n::Int, dims::NTuple{N, Integer}) where {N} = generate!(Array{Int}(undef, dims), n)
generate(n::Int, dims::Vararg{Integer, N}) where {N} = generate(n, dims)

################
"""
    vgenerate(n::Int)

Generate a random category from the uniform discrete distribution which has
support on `1` through `n` (inclusive).
"""
vgenerate(n::Int) = generate(n)

"""
    vgenerate!(A, u::AbstractArray{Float64}, n::Int)
    vgenerate!(A, n::Int)

Populate the array `A` with random categories drawn from the uniform discrete distribution
which has support on `1` through `n` (inclusive). If provided, `u` is used as temporary storage
for uniform(0,1) draws; it must be an array which has equal size as `A`.
If `u` is not provided, it will be allocated; repeated callers are encouraged to
pre-allocate this temporary storage once in order to avoid unnecessary allocations.
"""
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

"""
    vgenerate(n::Int, dims::Tuple)
    vgenerate(n::Int, dims::Integer...)

Generate an array of random categories from the uniform discrete distribution which has
support on `1` through `n` (inclusive).

See also: [`vgenerate!`](@ref)

# Examples
```julia-repl
julia> vgenerate(10, 3, 2)
3×2 Matrix{Int64}:
 9  9
 4  8
 6  1
```
"""
vgenerate(n::Int, dims::NTuple{N, Integer}) where {N} = vgenerate!(Array{Int}(undef, dims), n)
vgenerate(n::Int, dims::Vararg{Integer, N}) where {N} = vgenerate(n, dims)

####
# ugly, but makes @turbo generation case behave like base
vgenerate(K::Vector{Int}, V::Vector{<:AbstractFloat}, dims::Tuple{}) = generate(K, V, ())
vgenerate(n::Int, dims::Tuple{}) = generate(n, ())

################################################################
# Interface -- for convenience of sampling
"""
    generate(x::SqHist)
    generate(x::SqHistEquiprobable)

Generate a random category from the discrete distribution which corresponds to the
squared histogram `x`.
"""
generate(x::SqHist) = generate(x.K, x.V)

"""
    generate(x::SqHist, dims::Tuple)
    generate(x::SqHist, dims::Integer...)
    generate(x::SqHistEquiprobable, dims::Tuple)
    generate(x::SqHistEquiprobable, dims::Integer...)

Generate an array of random categories from the discrete distribution given the
squared histogram `x`.

See also: [`generate!`](@ref)

# Examples
```julia-repl
julia> p = [2/15, 7/15, 6/15]; x = SqHist(p);

julia> generate(x, 1,2,3)
1×2×3 Array{Int64, 3}:
[:, :, 1] =
 1  2

[:, :, 2] =
 3  3

[:, :, 3] =
 1  2
```
"""
generate(x::SqHist, dims::NTuple{N, Integer}) where {N} = generate(x.K, x.V, dims)
generate(x::SqHist, dims::Vararg{Integer, N}) where {N} = generate(x, dims)

"""
    generate!(A, u::AbstractArray{Float64}, x::SqHist)
    generate!(A, u::AbstractArray{Float64}, x::SqHistEquiprobable)
    generate!(A, x::SqHist)
    generate!(A, x::SqHistEquiprobable)

Populate the array `A` with random categories drawn from the discrete distribution which
corresponds to the squared histogram `x`. If provided, `u` is used as temporary storage
for uniform(0,1) draws; it must be an array which has equal size as `A`.
If `u` is not provided, it will be allocated; repeated callers are encouraged to
pre-allocate this temporary storage once in order to avoid unnecessary allocations.
"""
generate!(A, u, x::SqHist) = generate!(A, u, x.K, x.V)
generate!(A, x::SqHist) = generate!(A, x.K, x.V)

generate(x::SqHistEquiprobable) = generate(x.n)
generate(x::SqHistEquiprobable, dims::NTuple{N, Integer}) where {N} = generate(x.n, dims)
generate(x::SqHistEquiprobable, dims::Vararg{Integer, N}) where {N} = generate(x, dims)
generate!(A, u, x::SqHistEquiprobable) = generate!(A, u, x.n)
generate!(A, x::SqHistEquiprobable) = generate!(A, x.n)

"""
    vgenerate(x::SqHist)
    vgenerate(x::SqHistEquiprobable)

Generate a random category from the discrete distribution which corresponds to the
squared histogram `x`.
"""
vgenerate(x::SqHist) = vgenerate(x.K, x.V)

"""
    vgenerate(x::SqHist, dims::Tuple)
    vgenerate(x::SqHist, dims::Integer...)
    vgenerate(x::SqHistEquiprobable, dims::Tuple)
    vgenerate(x::SqHistEquiprobable, dims::Integer...)

Generate an array of random categories from the discrete distribution given the
squared histogram `x`.

See also: [`generate!`](@ref)

# Examples
```julia-repl
julia> p = [2/15, 7/15, 6/15]; x = SqHist(p);

julia> vgenerate(x, 1,2,3)
1×2×3 Array{Int64, 3}:
[:, :, 1] =
 2  3

[:, :, 2] =
 2  2

[:, :, 3] =
 1  2
```
"""
vgenerate(x::SqHist, dims::NTuple{N, Integer}) where {N} = vgenerate(x.K, x.V, dims)
vgenerate(x::SqHist, dims::Vararg{Integer, N}) where {N} = vgenerate(x, dims)

"""
    vgenerate!(A, u::AbstractArray{Float64}, x::SqHist)
    vgenerate!(A, u::AbstractArray{Float64}, x::SqHistEquiprobable)
    vgenerate!(A, x::SqHist)
    vgenerate!(A, x::SqHistEquiprobable)

Populate the array `A` with random categories drawn from the discrete distribution which
corresponds to the squared histogram `x`. If provided, `u` is used as temporary storage
for uniform(0,1) draws; it must be an array which has equal size as `A`.
If `u` is not provided, it will be allocated; repeated callers are encouraged to
pre-allocate this temporary storage once in order to avoid unnecessary allocations.
"""
vgenerate!(A, u, x::SqHist) = vgenerate!(A, u, x.K, x.V)
vgenerate!(A, x::SqHist) = vgenerate!(A, x.K, x.V)

vgenerate(x::SqHistEquiprobable) = vgenerate(x.n)
vgenerate(x::SqHistEquiprobable, dims::NTuple{N, Integer}) where {N} = vgenerate(x.n, dims)
vgenerate(x::SqHistEquiprobable, dims::Vararg{Integer, N}) where {N} = vgenerate(x, dims)
vgenerate!(A, u, x::SqHistEquiprobable) = vgenerate!(A, u, x.n)
vgenerate!(A, x::SqHistEquiprobable) = vgenerate!(A, x.n)
