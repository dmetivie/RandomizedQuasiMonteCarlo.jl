module RandomizedQuasiMonteCarlo

using DelimitedFiles
using LinearAlgebra: I, UnitLowerTriangular
using Random
using Distributions: Matrixvariate, Continuous
using Random: default_rng
import Distributions: Sampleable
import Random: _rand!, rand!

include("sampler.jl")
include("conversion.jl")
include("sobol.jl")

include("scrambling.jl")
include("shifting.jl")

export scramble!, shift!
export NestedUniformScrambler, LinearMatrixScrambler, Shifter
export sobol_pts2bits
export which_permutation
export nested_uniform_scramble_bit!, linear_matrix_scramble_bit!
export bits2unif, unif2bits

@deprecate sobol_indices which_permutation

end
