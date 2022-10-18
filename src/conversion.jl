
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

function unif2bits(y; M = 32, b = 2)
    # Turn a point in [0,1) into sequence of bits 
    # First bits are highest order
    bits = BitArray(zeros(M))
    for j in 1:M
        a = y - 1 / b^j
        bⱼ = a ≥ 0
        bits[j] = bⱼ
        if bⱼ
            y = a
        end
    end
    bits
end

#? TODO is it better to use this vectorized version or a bit(s)2int version?
function bits2int(b::BitMatrix)
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