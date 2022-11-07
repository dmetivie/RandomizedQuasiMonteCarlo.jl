"""
    shift(x::AbstractArray) 
Cranley Patterson Rotation i.e. `y = (x .+ U) mod 1` where `U ∼ 𝕌([0,1]ᵈ)` and `x` is a `n×d` matrix
"""
function shift(x::AbstractArray)
    y = copy(x)
    shift!(y)
    return y
end

function shift!(x::AbstractMatrix{<:Real})
    d = size(x, 2)
    s = zeros(d)
    shift!(x, s)
end

function shift!(x::AbstractMatrix{<:Real}, s::AbstractVector{<:Real})
    rand!(s)
    for i in 1:size(x, 1)
        x[i, :] += s
    end
    frac!(x)
end


function frac!(x)
    for i in eachindex(x)
        x[i] -= floor(x[i])
    end
end

# ! Tentative to use some kind of genereic API. Not super satysfing (intuitive) right now
struct Shifter{F<:AbstractMatrix{<:Real}} <: SamplerShifting{Matrixvariate,Continuous}
    points::F
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

shift!(sampler::Shifter, x::AbstractMatrix{<:Real}, s::AbstractVector{<:Real}) = shift!(default_rng(), sampler, x, s)
shift!(sampler::Shifter, x::AbstractMatrix{<:Real}) = shift!(default_rng(), sampler, x)

function digital_shift(x::AbstractArray, b::Integer; M=32)
    n, d = size(x)
    bits = points2bits(x, b; M = M)
    y = copy(x)
    for s in 1:d
        digital_shift!(@view(bits[:, s, :]), b)
        for i in 1:n    
            y[i, s] = bits2unif(bits[i, s, :], b)
        end
    end
    return y
end

function digital_shift!(rng::AbstractRNG, random_bits::AbstractMatrix{<:Integer}, b::Integer)
    n, M = size(random_bits)
    DS = rand(rng, 0:b-1, M)
    for i in 1:n
        random_bits[i,:] = (random_bits[i, :] + DS) .% b
    end
end

digital_shift!(random_bits, b) = digital_shift!(default_rng(), random_bits, b)
