"""
    bits2unif(bits::BitArray)
Convert a vector of M "bits" in base 2 into a number y∈[0,1[.
"""
function bits2unif(bits::BitArray)
    # Turn sequence of bits into a point in [0,1)
    # First bits are highest order
    y = 0//1
    for j in length(bits):-1:1
        y = (y + bits[j]) // 2
    end
    y
end

"""
    bits2unif(bits::AbstractArray{<:Integer}, b::Integer)
Convert a vector of M "bits" in base b into a number y∈[0,1[.
"""
function bits2unif(bits::AbstractArray{<:Integer}, b::Integer)
    # Turn sequence of bits into a point in [0,1)
    # First bits are highest order
    y = 0//1
    for j in length(bits):-1:1
        y = (y + bits[j]) // b
    end
    y
end

"""
    unif2bits(y, b::Integer; M=32)
Return the b-adic decomposition y = ∑ₖ yₖ/bᵏ a number y∈[0,1[ -> [y₁, ⋯, yₘ] 
"""
function unif2bits(y, b::Integer; M=32)
    bits = zeros(Int, M)
    unif2bits!(bits, y, b)
    return bits
end

"""
    unif2bits(y; M=32)
Return the binary decomposition y = ∑ₖ yₖ/2ᵏ a number y∈[0,1[ -> [y₁, ⋯, yₘ] 
"""
function unif2bits(y; M=32)
    bits = BitArray(undef, M)
    unif2bits!(bits, y, 2)
    return bits
end

function unif2bits!(bits::AbstractVector, y::Rational, b::Integer)
    for j in eachindex(bits), bb in b-1:-1:1 # see comment on nested loop to justify this writing
        a = y - bb // b^j
        bool = a > 0
        if a == 0
            bits[j] = bb
            break # it breaks from the nested loop (see here)[https://stackoverflow.com/questions/39796234/how-to-break-out-of-nested-for-loops-in-julia]
        elseif bool
            bits[j] = bb
            y = a
        end
    end
end

function unif2bits!(bits::AbstractVector, y, b::Integer)
    for j in eachindex(bits), bb in b-1:-1:1 # see comment on nested loop to justify this writing
        a = y - bb // b^j
        bool = a > 0
        if isapprox(a, 0, atol=3e-16) #! isequal(0.18518518518518515... -1/3^2-2/3^3 , 0) Should be true but is not! It fouls the binary expansion! Hence the isapprox with a hand tuned tolerence
            bits[j] = bb
            break # it breaks from the nested loop (see here)[https://stackoverflow.com/questions/39796234/how-to-break-out-of-nested-for-loops-in-julia]
        elseif bool
            bits[j] = bb
            y = a
        end
    end
end

#? TODO is it better to use this vectorized version or a bit(s)2int version?
"""
    bits2int(bits::BitMatrix)
Convert a matrix n×M of "bits" in base 2 into a vector of integers.
Inverse of int2bits
"""
function bits2int(bits::BitMatrix)
    # Convert bits into integers.
    # Inverse of int2bits

    n, p = size(bits)
    y = zeros(Int, n)

    for j in p:-1:1
        y = y * 2 + bits[:, j]
    end
    return y
end

"""
    bits2int(bit::AbstractMatrix{<:Integer}, b::Integer)
Convert a matrix n×M of "bits" in base b into a vector of integers.
Inverse of int2bits
"""
function bits2int(bit::AbstractMatrix{<:Integer}, b::Integer)

    #TODO @evalpoly faster?
    n, p = size(bit)
    y = zeros(Int, n)

    for j in p:-1:1
        y = y * b + bit[:, j]
    end
    return y
end
"""
    int2bits(x::Int, M)
 Convert an integer x into M bits
 For large x (or small M) you get bits of x modulo 2^M
"""
function int2bits(x::Int, M)
    y = BitArray(undef, M)
    int2bits!(y, x)
    return y
end

function int2bits!(y::BitVector, x::Int)
    for j in eachindex(y)
        y[j] = x % 2
        x = (x - y[j]) ÷ 2
    end
end