struct NestedUniformScrambler{F<:AbstractArray{Bool,3},T<:AbstractArray{I,3} where {I<:Integer}} <: RandomizationMethod
    Bit::F
    Index::T
end
""" 
Owen (1995) Nested Uniform Scrambling
"""
function nested_uniform_scramble(points::AbstractArray; M=32)
    n, s = size(points)
    unrandomized_bits = BitArray(undef, n, s, M)
    random_bits = similar(unrandomized_bits)
    for i in eachindex(view(points, 1:n, 1:s))
        unrandomized_bits[i, :] = unif2bits(points[i]; M = M)
    end
    indices = which_permutation(unrandomized_bits)
    nested_uniform_scramble_bit!(random_bits, unrandomized_bits, indices)
    random_points = similar(points)
    for i in eachindex(view(random_points, 1:n, 1:s))
        random_points[i] = bits2unif(random_bits[i, :])
    end
    return random_points
end

function scramble!(random_points::AbstractArray, random_bits::AbstractArray{Bool,3}, sampler::NestedUniformScrambler)
    nested_uniform_scramble_bit!(random_bits, sampler.Bit, sampler.Index)
    for i in eachindex(view(random_points, 1:size(random_points, 1), 1:size(random_points, 2)))
        random_points[i] = bits2unif(random_bits[i, :])
    end
end
# Note that we use rand(Bool) instead of bitrand() which seems faster for small array (but longer for larger which for m ≃ 20 is not achieved)
function nested_uniform_scramble_bit!(random_bits::AbstractArray{Bool,3}, thebits::AbstractArray{Bool,3}, indices::AbstractArray{T,3} where {T<:Integer})
    # in place Scramble Sobol' bits; nested uniform.
    #
    n, m, s = size(indices)
    M = size(random_bits, 3)
    @assert m ≥ 1 "We need m ≥ 1" # m=0 causes awkward corner case below.  Caller handles that case specially.

    for j in 1:s
        theperms = getpermset2(m)          # Permutations to apply to bits 1:m
        for k in 1:m                             # Here is where we want m > 0 so the loop works ok
            random_bits[:, j, k] = thebits[:, j, k] .⊻ theperms[k, indices[:, k, j]]   # permutation by adding a bit modulo 2 here with xor operator (only for base 2)
            # random_bits[:, j, k] = (thebits[:, j, k] + theperms[k][indices[j, k, :]]) .% b   # permutation by adding a bit modulo b
        end
    end
    if M > m     # Paste in random entries for bits after m'th one
        random_bits[:, :, (m+1):M] = rand(Bool, n * s * (M - m))
    end
end

function getpermset2(J::Integer)
    # Get 2^(j-1) random binary permutations for j=1 ... J
    # J will ordinarily be m when there are n=2^m points
    #
    # A nuisance is that m=0 gives J=0, requiring a list of length 0
    # that the for loop doesn't do as desired.
    # The caller will handle that corner case a different way.
    #
    # Caller has set the seed
    y = BitMatrix(undef, J, 2^(J - 1))
    for j in 1:J
        nj = 2^(j - 1)
        y[j, 1:nj] = rand(Bool, nj)
    end
    return y
end
"""
Assign at each points (for every dim) a number that will tell which permutation to use. 
"""
function which_permutation(bits::AbstractArray{Bool,3})
    n, s = size(bits, 1), size(bits, 2)
    m = Int(log2(n)) #? Should it be always 32 instead of that?
    indices = zeros(Int, n, m, s)
    for j in 1:s
        indices[:, 1, j] = zeros(Int, n) # same permutation for all observations i
        for k in 2:m                     # Here is where we want m > 0 so the loop works ok
            bitmat = bits[:, j, 1:(k-1)] # slice on dim j to get a matrix
            indices[:, k, j] = bits2int(bitmat) # index of which perms to use at bit k for each i
        end
    end
    return indices .+ 1
end

struct LinearMatrixScrambler{F<:AbstractArray{Bool,3}} <: RandomizationMethod
    Bit::F
end
"""
Matousek (1998) Linear Matrix Scrambling Scrambling
"""
#? Weird it should be faster than nested uniform Scrambling but here it is not at all.-> look for other implementation and paper
function scramble!(random_points::AbstractArray, random_bits::AbstractArray{Bool,3}, sampler::LinearMatrixScrambler)
    linear_matrix_scramble_bit!(random_bits, sampler.Bit)
    for s in 1:size(random_points, 2), i in 1:size(random_points, 1)
        random_points[i, s] = bits2unif(random_bits[i, s, :])
    end
end

#TODO add sparsification
function linear_matrix_scramble_bit!(random_bits::AbstractArray{Bool,3}, thebits::AbstractArray{Bool,3})
    # in place  bits; nested uniform.
    #
    n, s, M = size(thebits)
    m = Int(log2(n))
    @assert m ≥ 1 "We need m ≥ 1" # m=0 causes awkward corner case below.  Caller handles that case specially.

    for j in 1:s
        # Permutations matrix and shift to apply to bits 1:m
        matousek_M, matousek_C = getmatousek2(m)
        for k in 2:m                     # not correct for m=1
            bitmat = @view thebits[:, j, 1:(k-1)]
            random_bits[:, j, k] = (thebits[:, j, k] + bitmat * matousek_M[k, 1:(k-1)]) .% 2
        end
        # Digital shift
        for k in 1:m
            random_bits[:, j, k] = random_bits[:, j, k] .⊻ matousek_C[k]
        end
    end
    if M > m     # Paste in random entries for bits after m'th one
        random_bits[:, :, (m+1):M] = rand(Bool, n * s * (M - m))
    end
end

function getmatousek2(J::Integer)
    # Genereate the Matousek linear scramble in base 2 for one of the s components
    # We need a J x J bit matrix M and a length J bit vector C
    #
    matousek_M = UnitLowerTriangular(BitArray(undef, J, J))  # Identity
    matousek_C = rand(Bool, J)
    for i in 2:J
        for j in 1:(i-1)
            matousek_M[i, j] = rand(Bool)
        end
    end
    matousek_M, matousek_C
end
