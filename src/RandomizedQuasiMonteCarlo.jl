module RandomizedQuasiMonteCarlo

using DelimitedFiles
using LinearAlgebra: LowerTriangular, Diagonal, UnitLowerTriangular
using Random
using Distributions: Matrixvariate, Continuous
using Random: default_rng
import Distributions: Sampleable
import Random: _rand!, rand!

# See https://discourse.julialang.org/t/is-there-a-dedicated-function-computing-m-int-log-b-b-m/89776/10
function logi(b::Int, n::Int)
    m = round(Int, log(b, n))
    b^m == n || throw(ArgumentError("$n is not a power of $b"))
    return m
end

function log2i(n::Int)
    m = round(Int, log2(n))
    2^m == n || throw(ArgumentError("$n is not a power of 2"))
    return m
end

include("sampler.jl")
include("conversion.jl")
include("sobol.jl")

include("scrambling.jl")
include("scrambling_base_b.jl")
include("shifting.jl")

export scramble!, shift!
export NestedUniformScrambler, LinearMatrixScrambler, Shifter, NestedUniformScrambler_b, LinearMatrixScrambler_b
export sobol_pts2bits
export which_permutation, isequidistributed
export nested_uniform_scramble_bit!, linear_matrix_scramble_bit!
export bits2unif, unif2bits

@deprecate sobol_indices which_permutation

end
