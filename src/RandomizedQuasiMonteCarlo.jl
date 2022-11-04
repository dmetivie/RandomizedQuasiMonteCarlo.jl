module RandomizedQuasiMonteCarlo

using DelimitedFiles
using LinearAlgebra: LowerTriangular, Diagonal, UnitLowerTriangular
using Random
using Distributions: Matrixvariate, Continuous
using Random: default_rng
import Distributions: Sampleable
import Random: _rand!, rand!

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
