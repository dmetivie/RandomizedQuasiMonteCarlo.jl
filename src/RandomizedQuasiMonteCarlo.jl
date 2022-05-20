module RandomizedQuasiMonteCarlo

using DelimitedFiles
using LinearAlgebra: I, UnitLowerTriangular
using Random
using Distributions: Matrixvariate, Continuous
using Random: default_rng
import Distributions: Sampleable
import Random: _rand!, rand!

include("sampler.jl")
include("conversion.jl")
include("sobol.jl")

include("scrambling.jl")
include("shifting.jl")

# RandomizationMethod
# struct SamplerNestedUniform{F<:AbstractArray{Bool,3},T<:AbstractArray{I,3} where {I<:Integer}} <: RandomizationMethod
#     Bit::F
#     Index::T
# end
# Sampler(RNG::Type{<:AbstractRNG}, die::RandomizationMethod, r::Random.Repetition) =
#     SamplerSimple(die, Sampler(RNG, 1:die.nsides, r))

# Sampler(RNG::Type{<:AbstractRNG}, die::RandomizationMethod, ::Val{1}) = SamplerDie1(...)
# Sampler(RNG::Type{<:AbstractRNG}, die::RandomizationMethod, ::Val{Inf}) = SamplerDieMany(...)

# rand(rng::AbstractRNG, sp::SamplerSimple{Die}) = rand(rng, sp.data)

# abstract type SamplerScrambling <: Sampler{AbstractFloat} end
# struct NestedUniformScrambler{F<:AbstractArray{Bool,3},T<:AbstractArray{I,3} where {I<:Integer}} <: SamplerScrambling
#     Bit::F
#     Index::T
# end

# function rand(rng::AbstractRNG, sampler::SamplerScrambling)
#     # ... generate a single matrix sample to x
#     x = zeros(size(sampler))
#     random_bits = similar(sampler.Bit)
#     scramble!(rng, x, random_bits, sampler)
#     return x
# end

# function rand!(rng::AbstractRNG, x::DenseMatrix{T}, sampler::SamplerScrambling) where {T<:Real}
#     # ... generate a single matrix sample to x
#     random_bits = similar(sampler.Bit)
#     scramble!(rng, x, random_bits, sampler)
# end



# Sampler(RNG::Type{<:AbstractRNG}, scrambling::SamplerScrambling, r::Random.Repetition) =
#     SamplerDie(die, Sampler(RNG, 1:die.nsides, r))
# # the `r` parameter will be explained later on

# rand(rng::AbstractRNG, sp::SamplerDie) = rand(rng, sp.sp)


# function _rand!(rng::AbstractRNG, sampler::SamplerScrambling, x::DenseMatrix{T}) where {T<:Real}
#     # ... generate a single matrix sample to x
#     random_bits = similar(sampler.Bit)
#     scramble!(rng, x, random_bits, sampler)
# end

export scramble!, shift!
export NestedUniformScrambler, LinearMatrixScrambler, Shifter
export sobol_pts2bits
export which_permutation
export nested_uniform_scramble_bit!, linear_matrix_scramble_bit!
export bits2unif, unif2bits

end
