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

function rand!(rng::AbstractRNG, s::Sampleable{Matrixvariate}, A::AbstractMatrix)
    size(A) == size(s) ||
        throw(DimensionMismatch("Output size inconsistent with sample length."))
    _rand!(rng, s, A)
end