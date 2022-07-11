#
# Date created: 2022-07-05
# Author: aradclif
#
#
############################################################################################

p = normalize1!(rand(2^6));
K, V = sqhist(p);
A = Vector{Int}(undef, 2^11);
u = rand(length(A));
K2, V2 = sqhist_robinhood(p);
Av = view(A, 1:100);
uv = view(u, 1:100);

M = ones(Int, 100,2^10);
Mv = view(M, 1, :);

@benchmark generate!($A, $u, $K, $V)
@benchmark generate_v2!($A, $u, $K, $V)
@benchmark vgenerate!($A, $u, $K, $V)
@benchmark vgenerate_v2!($A, $u, $K, $V)
@benchmark generate!($Av, $uv, $K, $V)
@benchmark generate_v2!($Av, $uv, $K, $V)
@benchmark vgenerate!($Av, $uv, $K, $V)
@benchmark vgenerate_v2!($Av, $uv, $K, $V)

@benchmark generate!($Mv, $u, $K, $V)
@benchmark generate_v2!($Mv, $u, $K, $V)
@benchmark vgenerate!($Mv, $u, $K, $V)

generate!(Mv, u, K, V);
Mv0 = deepcopy(Mv);
vgenerate!(Mv, u, K, V);
Mv0 == Mv

A1 = generate!(similar(A), u, K, V);
A2 = generate_v2!(similar(A), u, K, V);
A1 == A2

u2 = rand(2^10);
@benchmark generate_option2!($A, $u2, $K, $V)

for i = 1:20
    for j = -1:1
        n = (1 << i) + j
        u = rand(n)
        A1 = similar(u, Int)
        A2 = similar(u, Int)
        generate!(A1, u, K, V)
        generate2!(A2, u, K, V)
        @test A1 == A2
    end
end

@code_native generate!(A, u, K, V)
@code_native generate2!(A, u, K, V)

v = 1
@benchmark vfill!($A, $v)

n = 2^10
@benchmark generate!($A, $u, $n)
@benchmark vgenerate!($A, $u, $n)
ur = 1:n
@benchmark rand!($A, $ur)

@benchmark generate($n)
@benchmark generate_v2($n)


@benchmark generate!($Mv, $u, $n)
@benchmark vgenerate!($Mv, $u, $n)



@benchmark fill!($Mv, $v)
@benchmark vfill!($Mv, $v)

################
x = SquareHistogram(p)
x[2:3]
x[:]
@timev x[2:3]
@timev K[2:3]
@timev V[2:3]
@timev Ref(K)


################
# Comparison: non-chunked vs. chunked (smaller U(0,1) buffer)
# Ultimately, non-chunked is slightly faster, even at small array sizes. As the array
# size increases, non-chunked pulls ahead -- at 2^20, it is 11% faster.
# In essence, if you have the memory to spare, then non-chunked is preferable.
## Thoughts
# Chunked could be offered as an option. It does impose the limitation that u must be ::Array{Float64}
# Furthermore, it makes it more difficult to guarantee safety, as non-unit stride A's
# would be a mess.
# Realistically, for repeated usage, it seems quite reasonable to ask callers to provide u
# in an appropriate size. This makes it easy to guarantee safety. The only time I can imagine
# it being perceived as a problem is if one want to sample an absurd number of draws in one
# pass, but then, the memory required is not really so bad -- only 2x that required to store
# the output.
include("generation_variants.jl")

A = rand(2^20);
u1 = rand(length(A));
u2 = rand(2^10);

@benchmark option1!($A, $u1)
@benchmark voption1!($A, $u1)
@benchmark option2!($A, $u2)
@benchmark voption2!($A, $u2)
@benchmark option2_v2!($A, $u2)
@benchmark option3!($A, $u2)
@benchmark option4!($A, $u2)

A = Vector{Int}(undef, 2^15);
u = rand(length(A));
uc = rand(2^10);

for i = 1:10
    n = 1 << i
    p = normalize1!(rand(n))
    K, V = sqhist(p)
    println("generate!, n = ", n)
    @btime generate!($A, $u, $K, $V)
    println("generate_option2!, n = ", n)
    @btime generate_option2!($A, $uc, $K, $V)
    println("vgenerate!, n = ", n)
    @btime vgenerate!($A, $u, $K, $V)
end


# Example: Bin(3, 0.3)
pmf(i, n, p) = binomial(n, i) * p^i * (1-p)^(n-i)
n = 3
p = 0.3
𝓅 = pmf.(0:n, n, p)

z = SqHist(𝓅)
x = generate(z, 10^8);
f(x) = x - 1
y = f.(x);

u = similar(x, Float64);
@benchmark generate!($x, $z)
@benchmark generate!($x, $u, $z)

using Distributions
d = Binomial(n, p)
w = rand(d, 10^8);
@benchmark rand!($d, $w)

using Plots
gr(size=(1200,800))
p1 = histogram(y, label="Marsaglia sampler");
p2 = histogram(w, label="Distributions");

p3 = plot(p1, p2);
savefig(p3, joinpath(pwd(), "binomials.pdf"))

# Example: Bin(100, 0.99)
using SpecialFunctions
pmf(i, n, p) = exp(loggamma(n + 1) - loggamma(i + 1) - loggamma(n - i + 1)) * p^i * (1-p)^(n-i)
pmf2(i, n, p) = exp(loggamma(n + 1) - loggamma(i + 1) - loggamma(n - i + 1)) * exp(i * log(p)) * exp((n - i) * log1p(-p))
pmf3(i, n, p) = exp(loggamma(n + 1) - loggamma(i + 1) - loggamma(n - i + 1) + i * log(p) + (n - i) * log1p(-p))

n = 100
p = 0.99
lb, ub = 80, 100
𝓅 = pmf.(lb:ub, n, p)

z = SqHist(𝓅)
x = generate(z, 10^8);
f(x) = x - 1
g(x, lb) = x + lb
g(x) = g(x, lb)
h = f ∘ g
f(x, lb) = x - 1 + lb
y = f.(x, lb); #h.(x);

@benchmark generate!($x, $z)
u = similar(x, Float64);
@benchmark generate!($x, $u, $z)

using Distributions
d = Binomial(n, p)
w = rand(d, 10^8);
@benchmark rand!($d, $w)

using Plots
gr(size=(1200,800))
p1 = histogram(y, label="Marsaglia sampler");
p2 = histogram(w, label="Distributions");

p3 = plot(p1, p2);
savefig(p3, joinpath(pwd(), "binomials_$(n).pdf"))

# Example: Categorical([0.3, 0.2, 0, …, 0, 0.5])
n = 10^4
p = [0.3; 0.2; fill(0.0, n - 3); 0.5];
p′ = [0.3, 0.2, 0.5];
I = [1, 2, n];
f(x) = I[x]
z = SqHist(p′)

x = generate(z, 10^8);
y = f.(x); getindex.(Ref(I), x);

@benchmark generate!($x, $z)
u = similar(x, Float64);
@benchmark generate!($x, $u, $z)

using Distributions
d = Categorical(p)
w = rand(d, 10^8);
@benchmark rand!($d, $w)

using Plots
gr(size=(1200,800))
p1 = histogram(y, label="Marsaglia sampler", bins=n);
p2 = histogram(w, label="Distributions", bins=n);

p3 = plot(p1, p2);
savefig(p3, joinpath(pwd(), "categoricals_$(n).pdf"))


#### Comparison to results from paper itself
using Distributions

l2cache = 1280 * 10^3
l2cache ÷ 2^4 # 80000

#
n_sample = 10^8 # 10^8
A = Vector{Int}(undef, n_sample);
U = similar(A, Float64);

d = Poisson(100.)
p = map(n -> pdf(d, n), 0:200)
K, V = sqhist(p);
@benchmark generate!($A, $U, $K, $V)
# 2^16 / (57.096 * 1e-6)
t = SqHist(p)
@benchmark generate!($A, $U, $t)

d = Binomial(100, .345)
p = map(n -> pdf(d, n), 0:100)
K, V = sqhist(p);
@benchmark generate!($A, $U, $K, $V)

for λ ∈ [1, 10, 25, 100, 250, 1000]
    println("λ = ", λ)
    d = Poisson(λ)
    p = map(n -> pdf(d, n), 0:max(1.5λ, 100))
    K, V = sqhist(p)
    @btime generate!($A, $U, $K, $V)
end


println("Time required to draw ", n_sample, " samples")
for n ∈ [20, 100, 1000, 10000, 100000]
    println("n = ", n)
    for ρ ∈ (.1, .4)
        println("\t p = ", ρ)
        d = Binomial(n, ρ)
        p = map(n -> pdf(d, n), 0:n)
        K, V = sqhist(p)
        @btime generate!($A, $U, $K, $V)
    end
end

println("Time required to draw ", n_sample, " samples")
for (N1, N2, K) ∈ [(20, 20, 20), (100, 100, 20), (100, 100, 100), (100, 1000, 100),
                   (1000, 1000, 100), (1000, 1000, 1000), (1000, 10000, 100),
                   (1000, 10000, 1000), (10000, 10000, 1000), (10000, 10000, 10000)]
    println("N1 = ", N1, " N2 = ", N2, " K = ", K)
    𝑠, 𝑓, 𝑛 = K, N2 - K, N1
    d = Hypergeometric(𝑠, 𝑓, 𝑛)
    p = map(n -> pdf(d, n), support(d))
    K, V = sqhist(p)
    # @btime generate!($A, $U, $K, $V)
    println(𝑠, " ", 𝑓, " ", 𝑛, "\t support: ", support(d))
end

for (N1, N2, K) ∈ [(20, 20, 20), (100, 100, 20), (100, 100, 100), (100, 1000, 100),
                   (1000, 1000, 100), (1000, 1000, 1000), (1000, 10000, 100),
                   (1000, 10000, 1000), (10000, 10000, 1000), (10000, 10000, 10000)]
    println("N1 = ", N1, " N2 = ", N2, " K = ", K)
    𝑠, 𝑓, 𝑛 = K, N2, N1
    d = Hypergeometric(𝑠, 𝑓, 𝑛)
    p = map(n -> pdf(d, n), support(d))
    K, V = sqhist(p)
    @btime generate!($A, $U, $K, $V)
    # println(𝑠, " ", 𝑓, " ", 𝑛, "\t support: ", support(d))
end

N, K, n = 100, 20, 100
𝑠, 𝑓, 𝑛 = K, N - K, n

N1, N2, K = 20,20,20
N1, N2, K = 100, 100, 20
N1, N2, K = 100,100,100
N1, N2, K = 100, 1000, 100
𝑠, 𝑓, 𝑛 = K, N1 - K, N2
d = Hypergeometric(𝑠, 𝑓, 𝑛)

gg(N1, N2, K) = K, N2 - K, N1
gg(N1, N2, K)
d = Hypergeometric(ans...)
