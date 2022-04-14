function sobol_pts2bits(m, s, M; file=joinpath(@__DIR__, "data", file_name(s)))
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

file_name(s) = s ≤ 50 ? "fiftysobol.col" : "sobol_Cs.col"

function sobol_generating_matrix(file, s, M)
    # main workhorse function for sobol_pts2bits
    # Return an array of sobol' matrices.
    # Not tuned for efficiency in time or space.
    # For instance, the matrices have many zeros.
    #
    col = readdlm(file, Int) #? Just export s line directly
    @assert s ≤ size(col, 1) "Not enough colunns in file fn = $fn. There should be at least s = $s."
    a = BitArray(undef, (s, M, M))
    for j in 1:s
        for k in 1:M
            a[j, :, k] = int2bits(col[j, k], M)
        end
    end
    return a
end