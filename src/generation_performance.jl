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
# It could be offered as an option. It does impose the limitation that u must be ::Array{Float64}
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
