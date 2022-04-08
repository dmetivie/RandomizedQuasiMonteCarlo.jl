module RandomizedQuasiMonteCarlo

using CSV, DataFrames
using LinearAlgebra: I

include("sobol.jl")
abstract type SamplingAlgorithm end
include("scrambling.jl")

export scramble!
export NestedUniformScrambler, LinearMatrixScrambler
export sobol_indices, sobol_pts2bits
export nested_uniform_scramble_bit!, linear_matrix_scramble_bit!
export bits2unif

end
