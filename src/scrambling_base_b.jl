struct NestedUniformScrambler_b{F<:AbstractArray{<:Integer,3},T<:AbstractArray{<:Integer,3}} <: SamplerScrambling{Matrixvariate,Continuous}
    Bit::F
    Index::T
    base::Integer
end
""" 
Owen (1995) Nested Uniform Scrambling
"""
function nested_uniform_scramble(points::AbstractArray, b::Integer; M=32)
    n, s = size(points)
    @assert isinteger(log(b, n)) "n must be of the form n=bᵐ with m ≥ 0"
    unrandomized_bits = zeros(Int, n, s, M)
    indices = zeros(Int, n, Int(log(b, n)), s)
    prepare_nested_uniform_scramble!(unrandomized_bits, indices, points, b)
    random_bits = similar(unrandomized_bits)
    nested_uniform_scramble_bit!(random_bits, unrandomized_bits, indices, b)
    random_points = similar(points)
    for i in eachindex(view(random_points, 1:n, 1:s))
        random_points[i] = bits2unif(random_bits[i, :], b)
    end
    return random_points
end

function prepare_nested_uniform_scramble(points::AbstractArray, b::Integer)
    n, s = size(points)
    @assert isinteger(log(b, n)) "n must be of the form n=bᵐ with m ≥ 0"
    unrandomized_bits = zeros(Int, n, s, M)
    indices = zeros(Int, n, Int(log(b, n)), s)
    prepare_nested_uniform_scramble!(unrandomized_bits, indices, points, b)
    return NestedUniformScrambler_b(unrandomized_bits, indices, b)
end

function prepare_nested_uniform_scramble!(bits, indices, points::AbstractArray, b::Integer)
    for i in eachindex(view(points, 1:size(bits, 1), 1:size(bits, 2)))
        bits[i, :] = unif2bits(points[i], b; M=size(bits, 3))
    end
    indices[:] = which_permutation(unrandomized_bits, b)
end

function scramble!(rng::AbstractRNG, random_points::AbstractArray, random_bits::AbstractArray{<:Integer,3}, sampler::NestedUniformScrambler_b)
    nested_uniform_scramble_bit!(rng, random_bits, sampler.Bit, sampler.Index, sampler.base)
    for i in eachindex(view(random_points, 1:size(random_points, 1), 1:size(random_points, 2)))
        random_points[i] = bits2unif(random_bits[i, :], sampler.base)
    end
end

scramble!(random_points::AbstractArray, random_bits, sampler::NestedUniformScrambler_b) = scramble!(default_rng(), random_points, random_bits, sampler)

function nested_uniform_scramble_bit!(rng::AbstractRNG, random_bits::AbstractArray{<:Integer,3}, thebits::AbstractArray{<:Integer,3}, indices::AbstractArray{T,3} where {T<:Integer}, b::Integer)
    # in place nested uniform Scramble.
    #
    n, m, d = size(indices)
    M = size(random_bits, 3)
    @assert m ≥ 1 "We need m ≥ 1" # m=0 causes awkward corner case below.

    for s in 1:d
        theperms = getpermset(rng, m, b)          # Permutations to apply to bits 1:m
        for k in 1:m                             # Here is where we want m > 0 so the loop works ok
            random_bits[:, s, k] = (thebits[:, s, k] + theperms[k, indices[:, k, s]]) .% b   # permutation by adding a bit modulo b
        end
    end
    if M > m     # Paste in random entries for bits after m'th one
        random_bits[:, :, (m+1):M] = rand(rng, 0:b-1, n * d * (M - m))
    end
end

function getpermset(rng::AbstractRNG, J::Integer, b::Integer)
    # Get b^(j-1) random binary permutations for j=1 ... J
    # J will ordinarily be m when there are n=b^m points
    #
    y = zeros(Int, J, b^(J - 1))
    for j in 1:J
        nj = b^(j - 1)
        y[j, 1:nj] = rand(rng, 0:b-1, nj)
    end
    return y
end

getpermset(J::Integer, b::Integer) = getpermset(Random.GLOBAL_RNG, J, b)
"""
    which_permutation(bits::AbstractArray{<:Integer,3}, b)
Assign at each points (for every dim) a number that will tell which permutation to use. 
"""
function which_permutation(bits::AbstractArray{<:Integer,3}, b)
    n, d = size(bits, 1), size(bits, 2)
    m = Int(log(b, n)) #? Should it be always M = 32 instead of that?
    #! M=32 is probably useless since you check equidistribution up to 1/bᵐ (small m not big M)
    indices = zeros(Int, n, m, d)
    for j in 1:d
        which_permutation!(@view(indices[:,:,j]), bits[:,j,:], b)
    end
    return indices
end

function which_permutation!(indices::AbstractMatrix{<:Integer}, bits::AbstractMatrix{<:Integer}, b::Integer)
    indices[:, 1] = zeros(Int, size(indices, 1)) # same permutation for all observations i
    for k in 2:size(indices, 2)                     # Here is where we want m > 0 so the loop works ok
        bitmat = bits[:, 1:(k-1)] # slice on dim j to get a matrix
        indices[:, k] = bits2int(bitmat, b) # index of which perms to use at bit k for each i
    end
    indices .+= 1 # array indexing starts at 1
end

#TODO add choice of dims
function isequidistributed(x::AbstractMatrix, b::Integer; M=32, dims = :all)
    n, d = size(x)
    @assert isinteger(log(b, n)) "n must be of the form n=bᵐ with m ≥ 0"
    bits = zeros(Int, n, d, M)
    for i in 1:n, s in 1:d
        unif2bits!(@view(bits[i, s, :]), x[i, s], b)
    end
    indices = which_permutation(bits, b)
    bool = zeros(Bool, d)
    for s in 1:d
        bool[s] = size(unique(indices[:,:,s] , dims = 1), 1) == n÷b
    end
    prod(bool)
end


struct LinearMatrixScrambler_b{F<:AbstractArray{<:Integer,3}} <: SamplerScrambling{Matrixvariate,Continuous}
    Bit::F
    base::Integer
end
"""
Matousek (1998) Linear Matrix Scrambling Scrambling
"""
#? Weird it should be faster than nested uniform Scrambling but here it is not at all.-> look for other implementation and paper
function scramble!(rng::AbstractRNG, random_points::AbstractArray, random_bits::AbstractArray{<:Integer,3}, sampler::LinearMatrixScrambler_b)
    linear_matrix_scramble_bit!(rng, random_bits, sampler.Bit, sampler.base)
    for s in 1:size(random_points, 2), i in 1:size(random_points, 1)
        random_points[i, s] = bits2unif(random_bits[i, s, :], sampler.base)
    end
end

scramble!(random_points::AbstractArray, random_bits, sampler::LinearMatrixScrambler_b) = scramble!(default_rng(), random_points, random_bits, sampler)

function linear_matrix_scramble_bit!(rng::AbstractRNG, random_bits::AbstractArray{<:Integer,3}, thebits::AbstractArray{<:Integer,3}, b::Integer)
    # https://statweb.stanford.edu/~owen/mc/ Chapter 17.6 around equation (17.15).
    #
    n, d, M = size(thebits)
    m = Int(log(b, n))
    @assert m ≥ 1 "We need m ≥ 1" # m=0 causes awkward corner case below.  Caller handles that case specially.

    for j in 1:d
        # Permutations matrix and shift to apply to bits 1:m
        matousek_M, matousek_C = getmatousek(rng, m, b)

        # xₖ = (∑ₗ Mₖₗ aₗ + Cₖ) mod b where xₖ is the k element in base b
        # thebits (n×m) × matousek_M (m×m) .+ matousek_C' (1×m) 
        random_bits[:, j, 1:m] = (thebits[:, j, 1:m] * matousek_M' .+ matousek_C') .% b
    end

    # Paste in random entries for bits after m'th one
    if M > m
        random_bits[:, :, (m+1):M] = rand(rng, Bool, n * d * (M - m))
    end
end

function digital_shift!(rng::AbstractRNG, random_bits::AbstractVector{<:Integer}, b::Integer)
    random_bits[:] = (random_bits[:] + rand(rng, 0:b-1, length(random_bits))) .% b
end

digital_shift!(random_bits, b) = digital_shift!(Random.GLOBAL_RNG, random_bits, b)

function getmatousek(rng::AbstractRNG, J::Integer, b::Integer)
    # Genereate the Matousek linear scramble in base b for one of the s components
    # We need a J x J bit matrix M and a length J bit vector C
    #
    matousek_M = LowerTriangular(zeros(Int, J, J)) + Diagonal(rand(rng, 1:b-1, J)) # Mₖₖ ∼ U{1, ⋯, b-1}
    matousek_C = rand(rng, 0:b-1, J)
    for i in 2:J
        for j in 1:(i-1)
            matousek_M[i, j] = rand(rng, 0:b-1)
        end
    end
    matousek_M, matousek_C
end

getmatousek(J::Integer, b::Integer) = getmatousek(Random.GLOBAL_RNG, J, b)

# old version. Was it faster? Correct not sure.
# function linear_matrix_scramble_bit!(rng::AbstractRNG, random_bits::AbstractArray{<:Integer,3}, thebits::AbstractArray{<:Integer,3}, b::Integer)
#     # in place  bits
#     #
#     n, d, M = size(thebits)
#     m = Int(log(b, n))
#     @assert m ≥ 1 "We need m ≥ 1" # m=0 causes awkward corner case below.  Caller handles that case specially.

#     for j in 1:d
#         # Permutations matrix and shift to apply to bits 1:m
#         # matousek_M, matousek_C = getmatousek(rng, m, b)
#         for k in 2:m                     # not correct for m=1
#             bitmat = @view thebits[:, j, 1:(k-1)]
#             random_bits[:, j, k] = (thebits[:, j, k] + bitmat * matousek_M[k, 1:(k-1)]) .% b
#         end
#         # Digital shift
#         for k in 1:m
#             random_bits[:, j, k] = (random_bits[:, j, k] .+ matousek_C[k]) .% b
#         end

#     end
#     if M > m     # Paste in random entries for bits after m'th one
#         random_bits[:, :, (m+1):M] = rand(rng, Bool, n * d * (M - m))
#     end
# end