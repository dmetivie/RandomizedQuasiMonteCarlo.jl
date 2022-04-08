function sobol_pts2bits(m, s, M; file = joinpath(string(pkgdir(RandomizedQuasiMonteCarlo)),"src/data/fiftysobol.col"))
    # Get array of sobol' bits.
    # Calls int2bits and sobol_generating_matrix
    # M must be the number of columns of data in the file named fn

    n = 2^m
    y = BitArray(undef, (n, s, M))             # n obs x s dimensions x M bits
    a = sobol_generating_matrix(file, s, M)
    bitsi = BitArray(undef, M)
    bitsj = zeros(Int, M)
    for i in 1:n
        int2bits!(bitsi, i - 1, M)   # bits of integer i-1 used for observation i
        for j in 1:s
            bitsj[:] = a[j, :, :] * bitsi        # a[j,,] is the sobol' matrix for dimension j
            y[i, j, :] = bitsj .% 2
        end
    end
    return y
end

function sobol_generating_matrix(file, s, M)
    # main workhorse function for sobol_pts2bits
    # Return an array of sobol' matrices.
    # Not tuned for efficiency in time or space.
    # For instance, the matrices have many zeros.
    #
    col = CSV.read(file, DataFrame, header=false) # D.M. Consider the first line as first line and not a header
    @assert s ≤ nrow(col) "Not enough colunns in file fn = $fn. There should be at least s = $s."
    a = BitArray(undef, (s, M, M))
    for j in 1:s
        for k in 1:M
            a[j, :, k] = int2bits(col[j, k], M)
        end
    end
    return a
end

function sobol_indices(sobol_bits)
    n, s = size(sobol_bits, 1), size(sobol_bits, 2)
    m = Int(log2(n)) #? Should it be always 32 instead of that?
    indices = zeros(Int, s, m, n)
    for j in 1:s
        indices[j, 1, :] = zeros(Int, n) # same permutation for all observations i
        for k in 2:m                     # Here is where we want m > 0 so the loop works ok
            bitmat = sobol_bits[:, j, :] # slice on dim j to get a matrix
            bitmat = bitmat[:, 1:(k-1)]  #
            indices[j, k, :] = bits2int(bitmat) # index of which perms to use at bit k for each i
        end
    end
    return indices
end


"""
Utilities
"""
function bits2unif(bits::BitArray)
    # Turn sequence of bits into a point in [0,1)
    # First bits are highest order
    y = 0
    for j in length(bits):-1:1
        y = (y + bits[j]) / 2
    end
    y
end

function bits2int(b)
    # Convert bits b into integers.
    # Inverse of int2bits

    n, p = size(b)
    y = zeros(Int, n)

    for j in p:-1:1
        y = y * 2 + b[:, j]
    end
    return y
end

function int2bits(x::Int, M)
    # Convert an integer x into M bits
    # For large x (or small M) you get bits of x modulo 2^M
    # This does just one integer (not a vector of them).

    y = BitArray(undef, M)
    for j in 1:M
        y[j] = x % 2
        x = (x - y[j]) ÷ 2
    end
    return y
end

function int2bits!(y::BitVector, x::Int, M)
    # Convert an integer x into M bits
    # For large x (or small M) you get bits of x modulo 2^M
    # This does just one integer (not a vector of them).

    for j in 1:M
        y[j] = x % 2
        x = (x - y[j]) ÷ 2
    end
    return y
end