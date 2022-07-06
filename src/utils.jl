#
# Date created: 2022-07-05
# Author: aradclif
#
#
############################################################################################

function vfindextrema(f::F, initmin::Iₘᵢₙ, initmax::Iₘₐₓ, A::AbstractArray{T, N}) where {F, Iₘᵢₙ, Iₘₐₓ, T, N}
    Tₒ = Base.promote_op(f, T)
    mn, mx = initmin(Tₒ), initmax(Tₒ)
    i_mn, i_mx = firstindex(A), firstindex(A)
    @turbo for i ∈ eachindex(A)
        v = f(A[i])
        newmin = v < mn
        newmax = v > mx
        i_mn = ifelse(newmin, i, i_mn)
        i_mx = ifelse(newmax, i, i_mx)
        mn = ifelse(newmin, v, mn)
        mx = ifelse(newmax, v, mx)
    end
    ((mn, i_mn), (mx, i_mx))
end
vfindextrema(f, A) = vfindextrema(f, typemax, typemin, A)
vfindextrema(A) = vfindextrema(identity, A)

# Sometimes more performant than fill!
function vfill!(A::AbstractArray, v::Real)
    @turbo for i ∈ eachindex(A)
        A[i] = v
    end
    A
end
