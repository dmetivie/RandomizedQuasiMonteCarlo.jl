function frac!(x)
    for i in eachindex(x)
        x[i] -= floor(x[i])
    end
end

function shift!(x, s)
    rand!(s)
    for i in 1:size(x, 1)
        x[i, :] += s
    end
    frac!(x)
end

function shift!(x)
    d = size(x, 2)
    s = rand(d)
    for i in 1:size(x, 1)
        x[i, :] += s
    end
    frac!(x)
end