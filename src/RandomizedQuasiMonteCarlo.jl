module RandomizedQuasiMonteCarlo

using DelimitedFiles
using Random: rand!
using LinearAlgebra: I, UnitLowerTriangular

include("sobol.jl")
abstract type SamplingAlgorithm end
include("scrambling.jl")
include("shifting.jl")

export scramble!, shift!
export NestedUniformScrambler, LinearMatrixScrambler
export sobol_pts2bits
export which_permutation
export nested_uniform_scramble_bit!, linear_matrix_scramble_bit!
export bits2unif

end
