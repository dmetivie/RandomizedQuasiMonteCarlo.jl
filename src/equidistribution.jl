
#TODO add choice of dims
function isequidistributed(x::AbstractMatrix, b::Integer; M=32, dims=:all)
    n, d = size(x)
    m = logi(b, n)
    bits = points2bits(x, b, M=M)
    indices = which_permutation(bits, b)
    bool = zeros(Bool, d)
    for s in 1:d
        bool[s] = size(unique(indices[:, :, s], dims=1), 1) == n ÷ b
    end
    all(bool)
end

#TODO
# test if point set is a t,m,d net in base b
function istmdnet(x::AbstractArray, b::Integer)
    bits = points2bits(x, b, M=M)
    for i in 0:m
        j = m - i
        xᵢ = range(0, 1, step=1 / b^(i))
        xⱼ = range(0, 1, step=1 / b^(j))
    end
end