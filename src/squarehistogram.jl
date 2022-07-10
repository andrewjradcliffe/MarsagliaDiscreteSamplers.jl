#
# Date created: 2022-07-05
# Author: aradclif
#
#
############################################################################################
# https://www.jstatsoft.org/article/view/v011i03
#### Marsaglia's Square Histogram (Method II in the above article)
# p ∈ ℝᴺ, ∑ᵢpᵢ = 1, a = 1/N
# K ∈ ℕᴺ, Kᵢ = i
# V ∈ ℝⁿ, Vᵢ = i * a
# Generate: j = ⌊N*U+1⌋; if U < V[j], return j, else return K[j]
# Theoretically, just one U ~ Uniform(0,1) is sufficient. In practice, it is also faster as
#     U=rand(); j=floor(Int, N * U + 1)
# is far fewer instructions than
#     j=rand(1:N); U=rand()
#### Robin Hood
## Motivation
# The frequency of the `else` statement being required is proportional to the "over-area",
# i.e. the part of the final squared histogram that lies above the division points.
# Or, to quantify it: ∑ᵢ (i * a) - V[i]
# The objective is to minimize the number of times the `else return K[j]` occurs
# by creating V such that each V[j] is as close to j/N as possible -- an NP-hard problem.
# Marsaglia's suggestion of the Robin Hood method is a good solution which is 𝒪(NlogN).
## Thoughts
# The reduction in `else` statements leads to faster sampling -- 𝒪(1) regardless --
# as it would certainly lead to a more predictable instruction pipeline due to minimizing
# the number of times the `else` branch is executed (or used, even if executed).
# Does it justify the 𝒪(NlogN) construction cost for the alias tables? -- given
# that we could pay 𝒪(N) instead to construct an inferior table (but no limit on how terrible).
#     - An analysis based on the objective function above would be informative;
#       can be verified with Monte Carlo: compute 𝒻(𝐱) = ∑ᵢ (i * a) - V[i]
#       using the V produced by each procedure.
#     - Number of samples to be drawn is an orthogonal decision variable;
#       one surmises that increasing number of samples favor better table.
## Algorithm
# repeat these two steps N - 1 times
#     1. find the smallest probability, pᵢ, and the largest probability, pⱼ
#     2. set K[i] = j; V[i] = (i - 1) * a + pᵢ; replace pⱼ with pⱼ - (a - pᵢ)
## Numerical stability
# Replacing pⱼ is the only point at which stability is a real concern. There are a few
# options for the order of operations and parenthesis. Unsurprisingly, Marsaglia gives
# the most stable form: pⱼ = pⱼ - (a - pᵢ)
# But it is worthwhile to show that this is the most stable form.
#     First, consider that it should be the case that pᵢ ≤ a, hence 0 ≤ (a - pᵢ) ≤ 1/n.
#     (a - pᵢ) may be a small number, but (a - pᵢ) is always well-defined since it occurs
#     at eps(a)≡ulp(a).
# It is worth noting that the (a - pᵢ) operation becomes unstable when pᵢ ≤ eps(a)/4, assuming
# the worst case, a=1. However, this has the opposite relationship to n: increasing n will
# result in a subtraction which takes place at smaller values, hence, the (floating)
# points are more densely packed (i.e. distance to nearest float is smaller).
# It is reassuring to note that eps(.5) = 2⁻⁵³, hence, even for a vector of length 2,
# the (a - pᵢ) is stable for pᵢ > 2⁻⁵⁵.
#     The subsequent subtraction, i.e. pⱼ - (a - pᵢ), will occur at eps(pⱼ)≡ulp(pⱼ). Thus,
#     the operation will be unstable when pⱼ - c * ulp(pⱼ), c ≤ 1/4 (for pⱼ = 1, the worst case).
# That is, unstable when: (a - pᵢ) ≤ eps(pⱼ)/4
# If pᵢ ≈ 0, (a - pᵢ) ≈ 1/n ⟹ 1/n ≤ eps(pⱼ)/4 is unstable
# As pⱼ is at most 1, the worst case will be eps(pⱼ)/4 = 2⁻⁵⁴, i.e. 1/n ≤ 2⁻⁵⁴.
# ∴ in the worse case, instability begins at n ≥ 2⁵⁴ if pᵢ ≈ 0.
# ∴ in general, expect instability if (a - pᵢ) ≤ 2⁻⁵⁴.
# The above assumed Float64 with 53 bits of precision; the general form in terms of precision
# replaces 2⁻⁵⁴ → 2⁻ᵖ⁻¹.
# These are very permissive bounds; one is likely to run into other issues well before
# the algorithm becomes numerically unstable.
## Numerical stability, revisit
# Oddly, Marsaglia's implementation in TplusSQ.c uses pⱼ + pᵢ - a, which has slightly
# worse numerical stability. (pⱼ + pᵢ) becomes unstable when pᵢ ≤ eps(pⱼ)/2, which
# corresponds to 2⁻ᵖ at pⱼ = 1.
# It is possible to find cases where both suffer roundoff, which is ≤ 2⁻ᵖ⁺¹ for pⱼ + pᵢ - a
# and ≤ 2⁻ᵖ for pⱼ - (a - pᵢ).
# Provided that one is working with Float64, it most likely does not matter.
# However, if p is provided as Float32, it may be preferable to forcibly promote to Float64
# just to ensure stability; naturally, Float16 is largely unsuitable and needs promotion.
# ##
# # Comparison of stability
# # f1 is unstable at n = 10, qⱼ = .999999999999
# # However, f2 and f3 are unstable at n = 10, qⱼ = 5
# f1(qⱼ, qᵢ, a) = qⱼ - (a - qᵢ)
# f2(qⱼ, qᵢ, a) = qⱼ + qᵢ - a
# ff2(qⱼ, qᵢ, a) = @fastmath qⱼ + qᵢ - a
# f3(qⱼ, qᵢ, a) = (qⱼ + qᵢ) - a
# n = 10
# a = 1 / n
# qⱼ = .999999999999
# qᵢ = (1 - qⱼ) / n
# f1(qⱼ, qᵢ, a)
# f2(qⱼ, qᵢ, a)
# ff2(qⱼ, qᵢ, a)
# f3(qⱼ, qᵢ, a)
# f1(big(qⱼ), big(qᵢ), big(a))
# f2(big(qⱼ), big(qᵢ), big(a))
# f3(big(qⱼ), big(qᵢ), big(a))
# # Another example
# f1(1.0, eps()/2, .0003)
# f2(1.0, eps()/2, .0003)
# f2(1.0, big(eps()/2), .0003)
# isfeq(qⱼ, qᵢ, a) = f1(qⱼ, qᵢ, a) == f2(qⱼ, qᵢ, a)
# isfeq_big(qⱼ, qᵢ, a) = f1(qⱼ, qᵢ, a) == Float64(f2(big(qⱼ), big(qᵢ), big(a)))
# count(i -> isfeq(qⱼ, i, a), qᵢ:qᵢ:n*qᵢ)
# isfeq.(qⱼ, qᵢ:qᵢ:n*qᵢ, a)
# isfeq_big.(qⱼ, qᵢ:qᵢ:n*qᵢ, a)

@inline _sqhist_init(T::Type{<:AbstractFloat}, n::Int) = Vector{Int}(undef, n), Vector{T}(undef, n), Vector{T}(undef, n)

"""
    sqhist_robinhood!(K::Vector{Int}, V::Vector{<:AbstractFloat}, q, p)

Construct into `K` and `V` a square histogram using the Robin Hood method,
given a vector of probabilities, `p`, which satisfies ∑ᵢpᵢ = 1 and 0 ≤ pᵢ ≤ 1 ∀i.
`q` is temporary storage, which will be overwritten, of the same type and size as `p`.

# Examples
```jldoctest
julia> p = [2/15, 7/15, 6/15];

julia> K = similar(p, Int); V = similar(p, Float64); q = zeros(3)
3-element Vector{Float64}:
 0.0
 0.0
 0.0

julia> sqhist_robinhood!(K, V, q, p)
([2, 3, 3], [0.13333333333333333, 0.6000000000000001, 1.0])

julia> q
3-element Vector{Float64}:
 0.3333333333333333
 0.3333333333333333
 0.3333333333333334
```
"""
function sqhist_robinhood!(K::Vector{Int}, V::Vector{T}, q::Vector{T}, p::Vector{T}) where {T<:AbstractFloat}
    n = length(p)
    a = inv(n)
    # initialize
    @inbounds for i ∈ eachindex(K, V, p, q)
        K[i] = i
        V[i] = i * a
        q[i] = p[i]
    end
    for _ = 1:n-1
        qᵢ, i = findmin(q)
        qⱼ, j = findmax(q)
        K[i] = j
        V[i] = (i - 1) * a + qᵢ
        q[j] = qⱼ - (a - qᵢ)
        q[i] = a
    end
    K, V
end

"""
    sqhist_robinhood!(K::Vector{Int}, V::Vector{<:AbstractFloat}, p)

Construct into `K` and `V` a square histogram using the Robin Hood method,
given a vector of probabilities, `p`, which satisfies ∑ᵢpᵢ = 1 and 0 ≤ pᵢ ≤ 1 ∀i.
Repeat callers are encouraged to pre-allocate temporary storage once, and call through
`sqhist_robinhood!(K, V, q, p)`

# Examples
```jldoctest
julia> p = [2/15, 7/15, 6/15];

julia> K = similar(p, Int); V = similar(p, Float64);

julia> sqhist_robinhood!(K, V, p)
([2, 3, 3], [0.13333333333333333, 0.6000000000000001, 1.0])
```
"""
sqhist_robinhood!(K, V, p) = sqhist_robinhood!(K, V, similar(p), p)

"""
    sqhist_robinhood(p::Vector{<:AbstractFloat})

Construct a square histogram using the Robin Hood method, given a vector of probabilities,
`p`, which satisfies ∑ᵢpᵢ = 1 and 0 ≤ pᵢ ≤ 1 ∀i. This algorithm requires 𝒪(N²) time, but
produces a square histogram which has fewer column pivots than `sqhist`, potentially
leading to gains in generation efficiency due to decreased probability of the `else`
branch occurring. As all generation algorithms are vectorized, this will have no effect
as both branches will be evaluated unconditionally. Consequently, there is little reason
to utilize this algorithm.

See also: [`sqhist_robinhood!`](@ref)

# Examples
```jldoctest
julia> p = [2/15, 7/15, 6/15];

julia> sqhist_robinhood(p)
([2, 3, 3], [0.13333333333333333, 0.6000000000000001, 1.0])
```
"""
function sqhist_robinhood(p::Vector{T}) where {T<:AbstractFloat}
    n = length(p)
    K, V, q = _sqhist_init(promote_type(T, Float64), n)
    sqhist_robinhood!(K, V, q, p)
end

function vsqhist_robinhood!(K::Vector{Int}, V::Vector{T}, q::Vector{T}, p::Vector{T}) where {T<:AbstractFloat}
    n = length(p)
    a = inv(n)
    @turbo for i ∈ eachindex(K, V, p, q)
        K[i] = i
        V[i] = i * a
        q[i] = p[i]
    end
    for _ = 1:n-1
        # qᵢ, i = vfindmin(q)
        # qⱼ, j = vfindmax(q)
        ((qᵢ, i), (qⱼ, j)) = vfindextrema(q)
        K[i] = j
        V[i] = (i - 1) * a + qᵢ
        q[j] = qⱼ - (a - qᵢ)
        q[i] = a
    end
    K, V
end

function vsqhist_robinhood(p::Vector{T}) where {T<:AbstractFloat}
    n = length(p)
    K, V, q = _sqhist_init(promote_type(T, Float64), n)
    vsqhist_robinhood!(K, V, q, p)
end

################

"""
    sqhist!(K::Vector{Int}, V::Vector{<:AbstractFloat}, L::Vector{Int}, S::Vector{Int}, q, p)
    sqhist!(K::Vector{Int}, V::Vector{<:AbstractFloat}, q, p)
    sqhist!(K::Vector{Int}, V::Vector{<:AbstractFloat}, p)

Construct into `K` and `V` a square histogram given a vector of probabilities, `p`,
which satisfies ∑ᵢpᵢ = 1 and 0 ≤ pᵢ ≤ 1 ∀i.
`L` and `S` are temporary storage, which will be overwritten, of the same size and type as `K`.
`q` is temporary storage, which will be overwritten, of the same type and size as `p`.
Repeat callers are encouraged to pre-allocate temporary storage once, and call through
`sqhist!(K, V, L, S, q, p)`; otherwise, `L`, `S` and `q` will be allocated as needed.

# Examples
```jldoctest
julia> p = [2/15, 7/15, 6/15];

julia> K = zeros(Int, 3); V = zeros(3); L = zeros(Int, 3); S = zeros(Int, 3); q = zeros(3);

julia> sqhist!(K, V, L, S, q, p)
([3, 2, 2], [0.13333333333333333, 0.6666666666666666, 0.8666666666666667])

julia> L, S, q
([2, 3, 0], [3, 0, 0], [0.3333333333333333, 0.33333333333333337, 0.3333333333333333])
```
"""
function sqhist!(K::Vector{Int}, V::Vector{T}, larges::Vector{Int}, smalls::Vector{Int}, q::Vector{T}, p::Vector{T}) where {T<:AbstractFloat}
    n = length(p)
    checkbounds(larges, n)
    checkbounds(smalls, n)
    a = inv(n)
    kl = 0
    ks = 0
    # initialize
    @inbounds for i ∈ eachindex(K, V, p, q)
        K[i] = i
        V[i] = i * a
        pᵢ = p[i]
        q[i] = pᵢ
        if pᵢ > a
            larges[kl+=1] = i
        else
            smalls[ks+=1] = i
        end
    end
    @inbounds while kl > 0 && ks > 0
        j = larges[kl]; kl -= 1
        i = smalls[ks]; ks -= 1
        qᵢ = q[i]
        qⱼ = q[j]
        K[i] = j
        V[i] = (i - 1) * a + qᵢ
        q[j] = qⱼ - (a - qᵢ)
        q[i] = a
        if q[j] > a
            larges[kl+=1] = j
        else
            smalls[ks+=1] = j
        end
    end
    K, V
end

"""
    sqhist!(K::Vector{Int}, V::Vector{<:AbstractFloat}, q, p)

Construct into `K` and `V` a square histogram given a vector of probabilities, `p`,
which satisfies ∑ᵢpᵢ = 1 and 0 ≤ pᵢ ≤ 1 ∀i.
`q` is temporary storage, which will be overwritten, of the same type and size as `p`.
Repeat callers are encouraged to pre-allocate temporary storage once, and call through
`sqhist!(K, V, L, S, q, p)`

# Examples
```jldoctest
julia> p = [2/15, 7/15, 6/15];

julia> K = similar(p, Int); V = similar(p, Float64); q = zeros(3)
3-element Vector{Float64}:
 0.0
 0.0
 0.0

julia> sqhist!(K, V, q, p)
([3, 2, 2], [0.13333333333333333, 0.6666666666666666, 0.8666666666666667])

julia> q
3-element Vector{Float64}:
 0.3333333333333333
 0.33333333333333337
 0.3333333333333333
```
"""
sqhist!(K, V, q, p) = (n = length(p); sqhist!(K, V, Vector{Int}(undef, n), Vector{Int}(undef, n), q, p))

"""
    sqhist!(K::Vector{Int}, V::Vector{<:AbstractFloat}, p)

Construct into `K` and `V` a square histogram given a vector of probabilities, `p`,
which satisfies ∑ᵢpᵢ = 1 and 0 ≤ pᵢ ≤ 1 ∀i.
Repeat callers are encouraged to pre-allocate temporary storage once, and call through
`sqhist!(K, V, L, S, q, p)`

# Examples
```jldoctest
julia> p = [2/15, 7/15, 6/15];

julia> K = similar(p, Int); V = similar(p, Float64);

julia> sqhist!(K, V, p)
([3, 2, 2], [0.13333333333333333, 0.6666666666666666, 0.8666666666666667])
```
"""
sqhist!(K, V, p) = sqhist!(K, V, similar(p), p)

"""
    sqhist(p::Vector{<:AbstractFloat})

Construct a square histogram, given a vector of probabilities,
`p`, which satisfies ∑ᵢpᵢ = 1 and 0 ≤ pᵢ ≤ 1 ∀i. This algorithm requires 𝒪(N) time, but
may produce a square histogram which has frequent column pivots, leading to a reduction
in generation efficiency due to increased probability of the `else` branch occurring.
As all generation algorithms are vectorized, this will have no effect, as both branches
will be evaluated unconditionally. Consequently, there is no penalty to using this algorithm.

See also: [`sqhist!`](@ref)

# Examples
```jldoctest
julia> p = [2/15, 7/15, 6/15];

julia> sqhist(p)
([3, 2, 2], [0.13333333333333333, 0.6666666666666666, 0.8666666666666667])
```
"""
function sqhist(p::Vector{T}) where {T<:AbstractFloat}
    n = length(p)
    K, V, q = _sqhist_init(promote_type(T, Float64), n)
    sqhist!(K, V, q, p)
end

# An optimized case, but why it would be used is questionable
function sqhist!(K::Vector{Int}, V::Vector{T}) where {T<:AbstractFloat}
    n = length(K)
    a = inv(n)
    @inbounds @simd for i ∈ eachindex(K, V)
        K[i] = i
        V[i] = i * a
    end
    K, V
end
sqhist!(K::Vector{Int}, V::Vector{<:AbstractFloat}, n::Int) = sqhist!(resize!(K, n), resize!(V, n))

################
# Interface -- for convenience of holding items together

"""
    SqHist(p)

Construct a square histogram, given a vector of probabilities, `p`, which satisfies
∑ᵢpᵢ = 1 and 0 ≤ pᵢ ≤ 1, using an 𝒪(N) algorithm. The resultant `K` and `V` are
wrapped in a struct for subsequent usage.

See also: [`sqhist`](@ref)

# Examples
julia> p = [2/15, 7/15, 6/15];

julia> x = SqHist(p)
SqHist{Int64, Float64}([3, 2, 2], [0.13333333333333333, 0.6666666666666666, 0.8666666666666667])
"""
struct SqHist{Ti<:Integer, Tv<:AbstractFloat}
    K::Vector{Ti}
    V::Vector{Tv}
end
SqHist(p::Vector{<:AbstractFloat}) = SqHist(sqhist(p)...)

Base.length(x::SqHist) = length(x.K)
Base.getindex(x::SqHist, i) = x.K[i], x.V[i]
Base.getindex(x::SqHist, ::Colon) = x.K, x.V
Base.:(==)(x::SqHist, y::SqHist) = x.K == y.K && x.V == y.V
Base.hash(x::SqHist, h::UInt) = hash(x.K, hash(x.V, h))
Base.resize!(x::SqHist, n::Int) = (resize!(x.K, n); resize!(x.V, n); x)

"""
    SqHistEquiprobable(n::Int)

Construct an efficient representation of the square histogram which corresponds to
`n` ≥ 1 categories of equal probability mass.

# Examples
julia> x = SqHistEquiprobable(3)
SqHistEquiprobable{Int64}(3)

julia> x[:]
(1:3, 0.3333333333333333:0.3333333333333333:1.0)

julia> K, V = map(collect, x[:])    # the literal version
([1, 2, 3], [0.3333333333333333, 0.6666666666666666, 1.0])
"""
struct SqHistEquiprobable{T<:Integer}
    n::T
    function SqHistEquiprobable{T}(x) where {T<:Integer}
        x > 0 || throw(ArgumentError("n must be > 0"))
        new(x)
    end
end
SqHistEquiprobable(x::T) where {T<:Integer} = SqHistEquiprobable{T}(x)

Base.convert(::Type{SqHistEquiprobable{T}}, x::SqHistEquiprobable) where {T<:Integer} =
    SqHistEquiprobable{T}(x.n)
Base.convert(::Type{SqHistEquiprobable{T}}, x::SqHistEquiprobable{T}) where {T<:Integer} = x

Base.length(x::SqHistEquiprobable) = x.n
Base.getindex(x::SqHistEquiprobable, i::Integer) = (checkbounds(1:x.n, i); (i, i / x.n))
Base.getindex(x::SqHistEquiprobable, i) = (checkbounds(1:x.n, i); a = inv(x.n); (i, i .* a))
Base.getindex(x::SqHistEquiprobable, ::Colon) = (a = inv(x.n); i = 1:x.n; (i, i .* a))
Base.:(==)(x::SqHistEquiprobable, y::SqHistEquiprobable) = x.n == y.n
Base.hash(x::SqHistEquiprobable, h::UInt) = hash(x.n, h)

# maybe?

"""
    sqhist!(x::SqHist, L::Vector{Int}, S::Vector{Int}, q, p)
    sqhist!(x::SqHist, q, p)
    sqhist!(x::SqHist, p)

Overwrite the `K`, `V` in `x` with a square histogram constructed from the vector of
probabilities, `p`, which satisfies ∑ᵢpᵢ = 1 and 0 ≤ pᵢ ≤ 1 ∀i.
`L` and `S` are temporary storage, which will be overwritten, of the same size and type as `K`.
`q` is temporary storage, which will be overwritten, of the same type and size as `p`.
Repeat callers are encouraged to pre-allocate temporary storage once, and call through
`sqhist!(x, L, S, q, p)`; otherwise, `L`, `S`, and `q` will be allocated as needed.

# Examples
```jldoctest
julia> x = SqHist([2/15, 7/15, 6/15])
SqHist{Int64, Float64}([3, 2, 2], [0.13333333333333333, 0.6666666666666666, 0.8666666666666667])

julia> sqhist!(x, [10/15, 1/15, 4/15])
SqHist{Int64, Float64}([1, 1, 1], [0.3333333333333333, 0.39999999999999997, 0.9333333333333333])
```
"""
sqhist!(x::SqHist, L, S, q, p) = (sqhist!(x.K, x.V, L, S, q, p); x)
sqhist!(x::SqHist, q, p) = (sqhist!(x.K, x.V, q, p); x)
sqhist!(x::SqHist, p) = (sqhist!(x.K, x.V, p); x)
