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


