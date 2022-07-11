# MarsagliaDiscreteSamplers

## Installation

```julia
using Pkg
Pkg.add("MarsagliaDiscreteSamplers")
```

## Usage

Given a discrete distribution defined by a probability table, `p ∈ ℝᴺ : ∑ᵢpᵢ = 1, 0 ≤ pᵢ ≤ 1`,
construct a square histogram in `𝒪(N)` time and space, then draw samples in `𝒪(1)` time.
Each draw is an index. In the simplest case, this index is identical to the category the
probability mass of which is `pᵢ`; this corresponds to `x ~ Categorical(p)`.

Any arbitrary discrete distribution may be sampled in this manner, as its probability mass function can be represented by a vector `p`. However, the draw(s) will need to be converted to whatever entity is indexed.
The categorical distribution is the trivial case, as the conversion function is the identity function.
If we take a slightly more complicated example, using the binomial distribution, e.g. `Binomial(3, 0.3)`, the probability mass function for which is

`p(i, n, p) = binomial(n, i) * p^i * (1-p)^(n-i)`

This yields the probability vector

`p = [1 * 0.7^3, 3 * 0.3 * 0.7^2, 3 * 0.3^2 * 0.7, 1 * 0.3^3]`

The draws returned by the sampler will be ∈ {1,2,3,4}, as there are 4 "categories" of outcome defined by the pmf of Bin(3, 0.3). It just happens to be the case that a 1-indexed numbering scheme for categories produces this. In other words, the first category corresponds to 0, the second to 1, etc. The conversion function is simply `f(x) = x - 1`.

At first glance, this may seem like a disadvantage, but it is in fact an advantage, as the sampler efficiency is not tied to any particular distribution. Consider a motivating example: `Binomial(100, 0.99)`,
the probability mass function for which is negligible for all but a small subset of its support.
One can sample from said distribution by providing the pmf evaluated on the subset of the support with non-negligible mass. Being generous, let us use the range `80:100`, such that `p ∈ ℝ²¹` which will return `x ∈ 1:21`. The conversion function is then `f(x) = x + 79`.
More generally, if the support of the binomial is restricted to `l:u` using lower bound, `l ≥ 0`, and upper bound, `u ≤ n`, this produces a `p ∈ ℝᵘ⁻ˡ⁺¹`, returns `x ∈ 1:(u-l+1)` and the conversion function to the true support is `f(x) = x + l - 1`. One could view this as the composition of functions, `h = f ∘ g`, `f(x) = x - 1` (unrestricted support), `g(x) = x + l` (adjustment for restricted support).
This applies to the Poisson, negative binomial, geometric, hypergeometric, etc. distributions.

Let us consider another example, applicable to any probability mass function which admits a sparse representation. The categorical distribution is selected for simplicity, but the approach is the same for any arbitrary pmf.
Consider `p ∈ ℝ¹⁰⁰⁰⁰`, with `p₁ = 0.3, p₂ = 0.2, p₁₀₀₀₀ = 0.5, and pᵢ = 0, i = 3,…,9999`.
It would be much more efficient to sample from `p′ = [p₁, p₂, p₁₀₀₀₀]`, and re-index the x's returned. To this end define `I = [1, 2, 10000]`. Now, `x ∈ 1:3` and one converts to the true category with `f(x) = I[x]`. Once one has the true category, `y = f(x)`, additional conversions (defined on the original support) can be applied.

## Examples
### `Binomial(3, 0.3)`

<details>
 <summaryClick me! ></summary>
<p>

The intent is to demonstrate the speed which can be achieved with Marsaglia's square histogram method.
```julia
julia> using MarsagliaDiscreteSamplers, BenchmarkTools

julia> pmf(i, n, p) = binomial(n, i) * p^i * (1-p)^(n-i);

julia> n = 3; p = 0.3;

julia> 𝓅 = pmf.(0:n, n, p)
4-element Vector{Float64}:
 0.3429999999999999
 0.4409999999999999
 0.189
 0.027

julia> z = SqHist(𝓅);

julia> x = generate(z, 10^6);

julia> f(x) = x - 1

julia> y = f.(x);

julia> @benchmark generate!($x, $z)
BenchmarkTools.Trial: 3117 samples with 1 evaluation.
 Range (min … max):  1.380 ms …   4.139 ms  ┊ GC (min … max): 0.00% … 34.68%
 Time  (median):     1.516 ms               ┊ GC (median):    0.00%
 Time  (mean ± σ):   1.597 ms ± 350.273 μs  ┊ GC (mean ± σ):  5.63% ± 11.95%

 Memory estimate: 7.63 MiB, allocs estimate: 2.

julia> u = similar(x, Float64); # pre-allocate temporary storage for repeated calls

julia> @benchmark generate!($x, $u, $z)
BenchmarkTools.Trial: 5988 samples with 1 evaluation.
 Range (min … max):  825.778 μs …  1.622 ms  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     827.810 μs              ┊ GC (median):    0.00%
 Time  (mean ± σ):   829.109 μs ± 16.532 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

 Memory estimate: 0 bytes, allocs estimate: 0.
 
julia> using Distributions

julia> d = Binomial(n, p)
Binomial{Float64}(n=3, p=0.3)

julia> w = rand(d, 10^6);

julia> @benchmark rand!($d, $w)
BenchmarkTools.Trial: 105 samples with 1 evaluation.
 Range (min … max):  47.509 ms …  49.507 ms  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     47.686 ms               ┊ GC (median):    0.00%
 Time  (mean ± σ):   47.743 ms ± 266.732 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%
 
 Memory estimate: 0 bytes, allocs estimate: 0.
 
julia> using Plots

julia> gr(size=(1200,800))

julia> p1 = histogram(y, label="Marsaglia square histogram method");

julia> p2 = histogram(w, label="Distributions");

julia> savefig(plot(p1, p2), joinpath(pwd(), "binomials_$(n)_$(p).pdf"))
```
</p>
</details>

### `Binomial(100, 0.99)`

<details>
 <summaryClick me! ></summary>
<p>

It is often beneficial, but not strictly necessary, to use only the points at which the pmf has non-negligible mass.
```julia
julia> using MarsagliaDiscreteSamplers, SpecialFunctions, BenchmarkTools

julia> pmf(i, n, p) = exp(loggamma(n+1) - loggamma(i+1) - loggamma(n-i+1)) * p^i * (1-p)^(n-i);

julia> n = 100; p = 0.99;

julia> lb, ub = 80, 100;

julia> 𝓅 = pmf.(lb:ub, n, p)
21-element Vector{Float64}:
 2.3986500044707484e-20
 5.86336667759521e-19
 1.3449991122629671e-17
 2.887696889220166e-16
 5.785706981616233e-15
 1.0781835128094341e-13
 1.8617471122347987e-12
 2.965955744318942e-11
 4.337710276066538e-10
 5.790112143783152e-9
 7.006035693977694e-8
 7.621950919821446e-7
 7.381693771261788e-6
 6.286345663268353e-5
 0.000463450802621816
 0.0028977871237615923
 0.014941714856895579
 0.060999165807532
 0.18486481882487396
 0.36972963764971967
 0.3660323412732292

julia> z = SqHist(𝓅);

julia> x = generate(z, 10^6);

julia> f(x, lb) = x - 1 + lb

julia> y = f.(x, lb);

julia> @benchmark generate!($x, $z)
BenchmarkTools.Trial: 2795 samples with 1 evaluation.
 Range (min … max):  1.372 ms …   3.968 ms  ┊ GC (min … max): 0.00% … 39.03%
 Time  (median):     1.524 ms               ┊ GC (median):    0.00%
 Time  (mean ± σ):   1.782 ms ± 588.526 μs  ┊ GC (mean ± σ):  6.31% ± 12.23%

 Memory estimate: 7.63 MiB, allocs estimate: 2.

julia> u = similar(x, Float64); # pre-allocate temporary storage for repeated calls

julia> @benchmark generate!($x, $u, $z)
BenchmarkTools.Trial: 5984 samples with 1 evaluation.
 Range (min … max):  826.627 μs …  1.621 ms  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     828.592 μs              ┊ GC (median):    0.00%
 Time  (mean ± σ):   829.513 μs ± 13.763 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

 Memory estimate: 0 bytes, allocs estimate: 0.

julia> using Distributions

julia> d = Binomial(n, p)
Binomial{Float64}(n=100, p=0.99)

julia> w = rand(d, 10^6);

julia> @benchmark rand!($d, $w)
BenchmarkTools.Trial: 116 samples with 1 evaluation.
 Range (min … max):  42.929 ms …  45.213 ms  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     43.037 ms               ┊ GC (median):    0.00%
 Time  (mean ± σ):   43.097 ms ± 321.675 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%
 
 Memory estimate: 0 bytes, allocs estimate: 0.
 
julia> using Plots

julia> gr(size=(1200,800))

julia> p1 = histogram(y, label="Marsaglia square histogram method");

julia> p2 = histogram(w, label="Distributions");

julia> savefig(plot(p1, p2), joinpath(pwd(), "binomials_$(n)_$(p).pdf"))
```
</p>
</details>

### `Categorical([0.3, 0.2, 0, …, 0, 0.5])`

<details>
 <summaryClick me! ></summary>
<p>

Far more efficient to sample from sparse pmf's using an approach such as below, though,
one does not begin to truly realize the gains until `n` becomes larger than can be held in the L2 cache.
```julia
julia> using MarsagliaDiscreteSamplers, BenchmarkTools

julia> n = 10^4

julia> p = [0.3; 0.2; fill(0.0, n - 3); 0.5];

julia> p′ = [0.3, 0.2, 0.5];

julia> I = [1, 2, n];

julia> z = SqHist(p′);

julia> x = generate(z, 10^6);

julia> f(x) = I[x];

julia> y = f.(x);

julia> @benchmark generate!($x, $z)
BenchmarkTools.Trial: 3162 samples with 1 evaluation.
 Range (min … max):  1.349 ms …   4.729 ms  ┊ GC (min … max): 0.00% … 69.94%
 Time  (median):     1.499 ms               ┊ GC (median):    0.00%
 Time  (mean ± σ):   1.574 ms ± 347.685 μs  ┊ GC (mean ± σ):  5.59% ± 12.02%

 Memory estimate: 7.63 MiB, allocs estimate: 2.

julia> u = similar(x, Float64); # pre-allocate temporary storage for repeated calls

julia> @benchmark generate!($x, $u, $z)
BenchmarkTools.Trial: 6117 samples with 1 evaluation.
 Range (min … max):  808.113 μs …  1.616 ms  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     809.980 μs              ┊ GC (median):    0.00%
 Time  (mean ± σ):   811.404 μs ± 16.157 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

 Memory estimate: 0 bytes, allocs estimate: 0.

julia> using Distributions

julia> d = Categorical(p);

julia> w = rand(d, 10^6);

julia> @benchmark rand!($d, $w)
BenchmarkTools.Trial: 664 samples with 1 evaluation.
 Range (min … max):  7.000 ms …  10.171 ms  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     7.467 ms               ┊ GC (median):    0.00%
 Time  (mean ± σ):   7.527 ms ± 303.820 μs  ┊ GC (mean ± σ):  0.03% ± 0.72%

 Memory estimate: 312.69 KiB, allocs estimate: 8.
 
julia> using Plots

julia> gr(size=(1200,800))

# Difficult to distinguish on plots

julia> p1 = histogram(y, label="Marsaglia square histogram method", bins=n);

julia> p2 = histogram(w, label="Distributions", bins=n);

julia> savefig(plot(p1, p2), joinpath(pwd(), "categoricals_$(n).pdf"))

julia> function unsafe_countcategory!(v::AbstractArray, A::AbstractArray)
           @inbounds for i ∈ eachindex(A)
               v[A[i]] += 1
           end
           v
       end;

julia> unsafe_countcategory(A, n::Int) = unsafe_countcategory!(zeros(Int, n), A);

julia> unsafe_countcategory(y, n)[I]
3-element Vector{Float64}:
 0.300233
 0.199834
 0.499933

julia> unsafe_countcategory(w, n)[I]
3-element Vector{Float64}:
 0.299828
 0.199998
 0.500174
```
</p>
</details>

### Comparison to results from article itself

<details>
 <summaryClick me! ></summary>
<p>

In Section 6 of the article, 3 tables of comparisons are given. The following replicates those
tables, using what is ultimately a simplification of "Method II" -- proceed straight to
the square histogram method, skipping the alias table. It turns out that for a vectorized sampler,
skipping the alias table portion of Method II yields faster code, as the square histogram step can be
made branchless. Further testing is needed to determine under what circumstances the branch (alias table or square histogram) may be worthwhile.
```julia
julia> using MarsagliaDiscreteSamplers, Distributions, BenchmarkTools

julia> n_sample = 10^8; A = Vector{Int}(undef, n_sample); U = similar(A, Float64);

# Time required to draw 10^8 samples from Binomial

julia> for n ∈ [20, 100, 1000, 10000, 100000]
           println("n = ", n)
           for ρ ∈ (.1, .4)
               println("\t p = ", ρ)
               d = Binomial(n, ρ)
               p = map(n -> pdf(d, n), 0:n)
               K, V = sqhist(p)
               @btime generate!($A, $U, $K, $V)
           end
       end
n = 20
         p = 0.1
  168.665 ms (0 allocations: 0 bytes)
         p = 0.4
  168.530 ms (0 allocations: 0 bytes)
n = 100
         p = 0.1
  169.132 ms (0 allocations: 0 bytes)
         p = 0.4
  168.981 ms (0 allocations: 0 bytes)
n = 1000
         p = 0.1
  169.180 ms (0 allocations: 0 bytes)
         p = 0.4
  169.197 ms (0 allocations: 0 bytes)
n = 10000
         p = 0.1
  190.279 ms (0 allocations: 0 bytes)
         p = 0.4
  190.923 ms (0 allocations: 0 bytes)
n = 100000
         p = 0.1
  318.831 ms (0 allocations: 0 bytes)
         p = 0.4
  315.848 ms (0 allocations: 0 bytes)
  
# Time required to draw 10^8 samples from Poisson

julia> for λ ∈ [1, 10, 25, 100, 250, 1000]
           println("λ = ", λ)
           d = Poisson(λ)
           p = map(n -> pdf(d, n), 0:max(1.5λ, 100))
           K, V = sqhist(p)
           @btime generate!($A, $U, $K, $V)
       end

λ = 1
  169.456 ms (0 allocations: 0 bytes)
λ = 10
  169.094 ms (0 allocations: 0 bytes)
λ = 25
  168.813 ms (0 allocations: 0 bytes)
λ = 100
  168.890 ms (0 allocations: 0 bytes)
λ = 250
  169.144 ms (0 allocations: 0 bytes)
λ = 1000
  169.457 ms (0 allocations: 0 bytes)

  
# Time required to draw 10^8 samples from Hypergeometric
# -- admittedly, not entirely clear that these are the parameterizations from the article
# This assumes that N1 is number of draws, N2 the population size and K the number of success states

julia> for (N1, N2, K) ∈ [(20, 20, 20), (100, 100, 20), (100, 100, 100), (100, 1000, 100),
                          (1000, 1000, 100), (1000, 1000, 1000), (1000, 10000, 100),
                          (1000, 10000, 1000), (10000, 10000, 1000), (10000, 10000, 10000)]
           println("N1 = ", N1, " N2 = ", N2, " K = ", K)
           𝑠, 𝑓, 𝑛 = K, N2 - K, N1
           d = Hypergeometric(𝑠, 𝑓, 𝑛)
           p = map(n -> pdf(d, n), support(d))
           K, V = sqhist(p)
           @btime generate!($A, $U, $K, $V)
       end
N1 = 20 N2 = 20 K = 20
  62.897 ms (0 allocations: 0 bytes)
N1 = 100 N2 = 100 K = 20
  62.913 ms (0 allocations: 0 bytes)
N1 = 100 N2 = 100 K = 100
  62.929 ms (0 allocations: 0 bytes)
N1 = 100 N2 = 1000 K = 100
  170.922 ms (0 allocations: 0 bytes)
N1 = 1000 N2 = 1000 K = 100
  62.888 ms (0 allocations: 0 bytes)
N1 = 1000 N2 = 1000 K = 1000
  62.919 ms (0 allocations: 0 bytes)
N1 = 1000 N2 = 10000 K = 100
  170.938 ms (0 allocations: 0 bytes)
N1 = 1000 N2 = 10000 K = 1000
  171.234 ms (0 allocations: 0 bytes)
N1 = 10000 N2 = 10000 K = 1000
  62.914 ms (0 allocations: 0 bytes)
N1 = 10000 N2 = 10000 K = 10000
  62.889 ms (0 allocations: 0 bytes)
  
# It could alternatively be N1 the number of draws, N2 the number of failures and K the
# number of success states -- this gives a more reasonable range of support.
# Alas, the files attached to the paper seem to be missing the hypergeometric distribution

julia> for (N1, N2, K) ∈ [(20, 20, 20), (100, 100, 20), (100, 100, 100), (100, 1000, 100),
                          (1000, 1000, 100), (1000, 1000, 1000), (1000, 10000, 100),
                          (1000, 10000, 1000), (10000, 10000, 1000), (10000, 10000, 10000)]
           println("N1 = ", N1, " N2 = ", N2, " K = ", K)
           𝑠, 𝑓, 𝑛 = K, N2, N1
           d = Hypergeometric(𝑠, 𝑓, 𝑛)
           p = map(n -> pdf(d, n), support(d))
           K, V = sqhist(p)
           @btime generate!($A, $U, $K, $V)
       end
N1 = 20 N2 = 20 K = 20
  170.708 ms (0 allocations: 0 bytes)
N1 = 100 N2 = 100 K = 20
  170.384 ms (0 allocations: 0 bytes)
N1 = 100 N2 = 100 K = 100
  170.881 ms (0 allocations: 0 bytes)
N1 = 100 N2 = 1000 K = 100
  170.931 ms (0 allocations: 0 bytes)
N1 = 1000 N2 = 1000 K = 100
  170.873 ms (0 allocations: 0 bytes)
N1 = 1000 N2 = 1000 K = 1000
  170.993 ms (0 allocations: 0 bytes)
N1 = 1000 N2 = 10000 K = 100
  170.966 ms (0 allocations: 0 bytes)
N1 = 1000 N2 = 10000 K = 1000
  171.030 ms (0 allocations: 0 bytes)
N1 = 10000 N2 = 10000 K = 1000
  171.167 ms (0 allocations: 0 bytes)
N1 = 10000 N2 = 10000 K = 10000
  192.990 ms (0 allocations: 0 bytes)
```
</p>
</details>
