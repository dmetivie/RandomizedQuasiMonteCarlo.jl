struct NestedUniformScrambler{F<:AbstractArray{Bool,3},T<:AbstractArray{I,3} where {I<:Integer}} <: SamplingAlgorithm
    Bit::F
    Index::T
end
"""
Owen (1995) Nested Uniform Scrambling
"""
function scramble!(random_points, random_bits, sampler::NestedUniformScrambler)
    nested_uniform_scramble_bit!(random_bits, sampler.Bit, sampler.Index)
    for s in 1:size(random_points, 2), i in 1:size(random_points, 1)
        random_points[i, s] = bits2unif(random_bits[i, s, :])
    end
end

function nested_uniform_scramble_bit!(random_bits, thebits, indices)
    # in place Scramble Sobol' bits; nested uniform.
    #
    s, m, n = size(indices)
    M = size(random_bits, 3)
    @assert m ≥ 1 "We need m ≥ 1" # m=0 causes awkward corner case below.  Caller handles that case specially.

    for j in 1:s
        theperms = getpermset2(m)          # Permutations to apply to bits 1:m
        for k in 1:m                             # Here is where we want m > 0 so the loop works ok
            random_bits[:, j, k] = (thebits[:, j, k] + theperms[k][1 .+ indices[j, k, :]]) .% 2   # permutation by adding a bit modulo 2
        end
    end
    if M > m     # Paste in random entries for bits after m'th one
        random_bits[:, :, (m+1):M] = rand(n * s * (M - m)) .> 0.5
    end
end

function getpermset2(J)
    # Get 2^(j-1) random binary permutations for j=1 ... J
    # J will ordinarily be m when there are n=2^m points
    #
    # A nuisance is that m=0 gives J=0, requiring a list of length 0
    # that the for loop doesn't do as desired.
    # The caller will handle that corner case a different way.
    #
    # Caller has set the seed
    y = Vector{BitVector}(undef, J)
    for j in 1:J
        y[j] = rand(2^(j - 1)) .> 1 / 2 #? Bernoulli?
    end
    return y
end

struct LinearMatrixScrambler{F<:AbstractArray{Bool,3}} <: SamplingAlgorithm
    Bit::F
end
"""
Matousek (1998) Nested Uniform Scrambling
"""
#? Weird it should be faster than nested uniform Scrambling but here it is not at all.-> look for other implementation and paper
function scramble!(random_points, random_bits, sampler::LinearMatrixScrambler)
    linear_matrix_scramble_bit!(random_bits, sampler.Bit)
    for s in 1:size(random_points, 2), i in 1:size(random_points, 1)
        random_points[i, s] = bits2unif(random_bits[i, s, :])
    end
end

function linear_matrix_scramble_bit!(random_bits, thebits)
    # in place  bits; nested uniform.
    #
    n, s, M = size(thebits)
    m = Int(log2(n))
    @assert m ≥ 1 "We need m ≥ 1" # m=0 causes awkward corner case below.  Caller handles that case specially.

    for j in 1:s
        matousek_M, matousek_C = getmatousek2(m)     # Permutations to apply to bits 1:m
        for k in 2:m                     # not correct for m=1
            bitmat = @view thebits[:, j, 1:(k-1)]
            random_bits[:, j, k] = (thebits[:, j, k] + bitmat * matousek_M[k, 1:(k-1)]) .% 2
        end
        for k in 1:m
            random_bits[:, j, k] = (random_bits[:, j, k] .+ matousek_C[k]) .% 2
        end
    end
    if M > m     # Paste in random entries for bits after m'th one
        random_bits[:, :, (m+1):M] = rand(n * s * (M - m)) .> 0.5
    end
end

function getmatousek2(J)
    # Genereate the Matousek linear scramble in base 2 for one of the s components
    # We need a J x J bit matrix M and a length J bit vector C
    #
    matousek_M = Matrix{Bool}(I, J, J)  # Identity
    matousek_C = rand(J) .> 1 / 2
    for i in 2:J
        for j in 1:(i-1)
            matousek_M[i, j] = rand() > 1 / 2
        end
    end
    matousek_M, matousek_C
end
