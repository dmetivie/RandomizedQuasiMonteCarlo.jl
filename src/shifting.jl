struct Shifter{F<:AbstractMatrix{<:Real}} <: SamplerShifting{Matrixvariate,Continuous}
    points::F
end

function prepare_shift(x)
    return Shifter(x)
end

function shift!(rng::AbstractRNG, sampler::Shifter, x::AbstractMatrix{<:Real}, s::AbstractVector{<:Real})
    rand!(rng, s)
    for i in 1:size(x, 1)
        x[i, :] = sampler.points[i, :] + s
    end
    frac!(x)
end

function shift!(rng::AbstractRNG, sampler::Shifter, x::AbstractMatrix{<:Real})
    d = size(x, 2)
    s = zeros(d)
    shift!(rng, sampler, x, s)
end

function shift!(rng::AbstractRNG, sampler::Shifter, x::AbstractMatrix{<:Real})
    d = size(x, 2)
    s = zeros(d)
    shift!(rng, sampler, x, s)
end

shift!(sampler::Shifter, x::AbstractMatrix{<:Real}, s::AbstractVector{<:Real}) = shift!(default_rng(), sampler, x, s)
shift!(sampler::Shifter, x::AbstractMatrix{<:Real}) = shift!(default_rng(), sampler, x)


# function shift!(x, s)
#     rand!(s)
#     for i in 1:size(x, 1)
#         x[i, :] += s
#     end
#     frac!(x)
# end

# function shift!(x)
#     d = size(x, 2)
#     s = zeros(d)
#     shift!(x, s)
# end

function frac!(x)
    for i in eachindex(x)
        x[i] -= floor(x[i])
    end
end