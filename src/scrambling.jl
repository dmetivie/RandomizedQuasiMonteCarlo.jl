struct NestedUniformScrambler{F<:AbstractArray{Bool,3},T<:AbstractArray{I,3} where {I<:Integer}} <: SamplerScrambling{Matrixvariate,Continuous}
    Bit::F
    Index::T
end
""" 
    nested_uniform_scramble(points::AbstractArray; M=32)
Return a scrambled version of the `points`. 
The scrambling method is Nested Uniform Scrambling which was introduced in Owen (1995).
M is the number of bits used for each points. One need M ≥ log(base, n). 
"""
function nested_uniform_scramble(rng::AbstractRNG, points::AbstractArray; M=32)
    n, d = size(points)
    m = log2i(n)
    unrandomized_bits = BitArray(undef, n, d, M)
    indices = zeros(Int, n, m, d)
    prepare_nested_uniform_scramble!(unrandomized_bits, indices, points)
    random_bits = similar(unrandomized_bits)
    nested_uniform_scramble_bit!(rng, random_bits, unrandomized_bits, indices)
    random_points = similar(points)
    for i in eachindex(view(random_points, 1:n, 1:d))
        random_points[i] = bits2unif(random_bits[i, :])
    end
    return random_points
end

nested_uniform_scramble(points::AbstractArray; M=32) = nested_uniform_scramble(Random.default_rng(), points::AbstractArray; M=M)

function prepare_nested_uniform_scramble(points::AbstractArray)
    n, d = size(points)
    m = log2i(n)
    unrandomized_bits = BitArray(undef, n, d, M)
    indices = zeros(Int, n, m, d)
    prepare_nested_uniform_scramble!(unrandomized_bits, indices, points)
    return NestedUniformScrambler(unrandomized_bits, indices)
end

function prepare_nested_uniform_scramble!(bits::AbstractArray{Bool,3}, indices::AbstractArray{<:Integer}, points::AbstractArray)
    points2bits!(bits, points)
    indices[:] = which_permutation(bits)
end

function scramble!(rng::AbstractRNG, random_points::AbstractArray, random_bits::AbstractArray{Bool,3}, sampler::NestedUniformScrambler)
    nested_uniform_scramble_bit!(rng, random_bits, sampler.Bit, sampler.Index)
    for i in eachindex(view(random_points, 1:size(random_points, 1), 1:size(random_points, 2)))
        random_points[i] = bits2unif(random_bits[i, :])
    end
end

scramble!(random_points::AbstractArray, random_bits::AbstractArray{Bool,3}, sampler::NestedUniformScrambler) = scramble!(default_rng(), random_points, random_bits, sampler)

# Note that we use rand(Bool) instead of bitrand() which seems faster for small array (but longer for larger which for m ≃ 20 is not achieved)
function nested_uniform_scramble_bit!(rng::AbstractRNG, random_bits::AbstractArray{Bool,3}, thebits::AbstractArray{Bool,3}, indices::AbstractArray{T,3} where {T<:Integer})
    # in place Scramble Sobol' bits; nested uniform.
    #
    n, m, d = size(indices)
    M = size(random_bits, 3)
    @assert m ≥ 1 "We need m ≥ 1" # m=0 causes awkward corner case below.  Caller handles that case specially.

    for s in 1:d
        theperms = getpermset(rng, m)          # Permutations to apply to bits 1:m
        for k in 1:m                             # Here is where we want m > 0 so the loop works ok
            random_bits[:, s, k] = thebits[:, s, k] .⊻ theperms[k, indices[:, k, s]]   # permutation by adding a bit modulo 2 here with xor operator (only for base 2)
        end
    end
    if M > m     # Paste in random entries for bits after m'th one
        random_bits[:, :, (m+1):M] = rand(rng, Bool, n * d * (M - m))
    end
end

function getpermset(rng::AbstractRNG, J::Integer)
    # Get 2^(j-1) random binary permutations for j=1 ... J
    # J will ordinarily be m when there are n=2^m points
    #
    # A nuisance is that m=0 gives J=0, requiring a list of length 0
    # that the for loop doesn't do as desired.
    # The caller will handle that corner case a different way.
    #
    y = BitMatrix(undef, J, 2^(J - 1))
    for j in 1:J
        nj = 2^(j - 1)
        y[j, 1:nj] = rand(rng, Bool, nj)
    end
    return y
end

getpermset(J::Integer) = getpermset(Random.GLOBAL_RNG, J::Integer)

which_permutation(bits::AbstractArray{Bool,3}) = which_permutation(bits, 2)

struct LinearMatrixScrambler{F<:AbstractArray{Bool,3}} <: SamplerScrambling{Matrixvariate,Continuous}
    Bit::F
end
""" 
    nested_uniform_scramble(points::AbstractArray; M=32)
Return a scrambled version of the `points`. 
The scrambling method is Linear Matrix Scrambling Scrambling which was introduced in Matousek (1998).
`M` is the number of bits used for each points. One need `M ≥ log(base, n)`. 
"""
function linear_matrix_scramble(points::AbstractArray; M=32)
    n, d = size(points)
    m = log2i(n)
    unrandomized_bits = BitArray(undef, n, d, M)
    points2bits!(unrandomized_bits, points)
    random_bits = similar(unrandomized_bits)
    linear_matrix_scramble_bit!(random_bits, unrandomized_bits)
    random_points = similar(points)
    for i in eachindex(view(random_points, 1:n, 1:d))
        random_points[i] = bits2unif(random_bits[i, :])
    end
    return random_points
end
"""

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
    n, d, M = size(thebits)
    m = log2i(n)
    @assert m ≥ 1 "We need m ≥ 1" # m=0 causes awkward corner case below.  Caller handles that case specially.

    for j in 1:d
        # Permutations matrix and shift to apply to bits 1:m
        matousek_M, matousek_C = getmatousek(m)
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
        random_bits[:, :, (m+1):M] = rand(Bool, n * d * (M - m))
    end
end

function getmatousek(J::Integer)
    # Genereate the Matousek linear scramble in base 2 for one of the d components
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

function points2bits!(bits, points)
    for i in eachindex(view(points, 1:size(bits, 1), 1:size(bits, 2)))
        bits[i, :] = unif2bits(points[i]; M=size(bits, 3))
    end
end

function points2bits!(bits, points, b::Integer)
    for i in eachindex(view(points, 1:size(bits, 1), 1:size(bits, 2)))
        bits[i, :] = unif2bits(points[i], b; M=size(bits, 3))
    end
end

function points2bits(x::AbstractArray, b::Integer; M=32)
    n, d = size(x)
    bits = zeros(Int, n, d, M)
    points2bits!(bits, x, b)
    return bits
end