#
# Date created: 2022-07-05
# Author: aradclif
#
#
############################################################################################

p = normalize1!(rand(2^10));
K, V = sqhist(p);
A = Vector{Int}(undef, 2^10);
u = rand(length(A));
K2, V2 = sqhist_robinhood(p);
Av = view(A, 1:100);
uv = view(u, 1:100);


@benchmark generate!($A, $u, $K, $V)
@benchmark generate_v2!($A, $u, $K, $V)
@benchmark vgenerate!($A, $u, $K, $V)
@benchmark vgenerate_v2!($A, $u, $K, $V)
@benchmark generate!($Av, $uv, $K, $V)
@benchmark generate_v2!($Av, $uv, $K, $V)
@benchmark vgenerate!($Av, $uv, $K, $V)
@benchmark vgenerate_v2!($Av, $uv, $K, $V)

A1 = generate!(similar(A), u, K, V);
A2 = generate_v2!(similar(A), u, K, V);
A1 == A2

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
