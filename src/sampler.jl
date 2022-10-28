#* Scrambling types
abstract type SamplerScrambling{F<:Matrixvariate,S<:Continuous} <: Sampleable{F,S} end

Base.size(sampler::SamplerScrambling) = (size(sampler.Bit, 1), size(sampler.Bit, 2)) # the size of each matrix sample

function _rand!(rng::AbstractRNG, sampler::SamplerScrambling, A::DenseMatrix{T}) where {T<:Real}
    # ... generate a single matrix sample to x
    random_bits = similar(sampler.Bit)
    scramble!(rng, A, random_bits, sampler)
end

function rand!(rng::AbstractRNG, s::Sampleable{Matrixvariate}, A::AbstractMatrix)
    size(A) == size(s) ||
        throw(DimensionMismatch("Output size inconsistent with sample length."))
    _rand!(rng, s, A)
end

function _rand!(rng::AbstractRNG, sampler::SamplerScrambling, A::AbstractArray{<:AbstractMatrix{T}}) where {T<:Real}
    # ... generate multiple matrix sample to x
    random_bits = similar(sampler.Bit)
    for i in eachindex(A)
        scramble!(rng, A[i], random_bits, sampler)
    end
end

function rand!(rng::AbstractRNG, s::Sampleable{Matrixvariate}, A::AbstractArray{<:AbstractMatrix{T}}) where {T<:Real}
    (size(A[firstindex(A)], 1), size(A[firstindex(A)], 2)) == size(s) ||
        throw(DimensionMismatch("Output size inconsistent with sample length."))
    _rand!(rng, s, A)
end

#* Shifting types
abstract type SamplerShifting{F<:Matrixvariate,S<:Continuous} <: Sampleable{F,S} end
Base.size(sampler::SamplerShifting) = size(sampler.points) # the size of each matrix sample

function _rand!(rng::AbstractRNG, sampler::SamplerShifting, A::DenseMatrix{T}) where {T<:Real}
    # ... generate a single matrix sample to A
    #... potentially some initialization
    shift!(rng, sampler, A)
    return A #? is it not bad to return the value here for performance ? 
end


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

