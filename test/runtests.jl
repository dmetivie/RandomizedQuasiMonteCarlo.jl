using RandomizedQuasiMonteCarlo
using Test

#? What to test execpt function? Discrepency?
@testset "RandomizedQuasiMonteCarlo.jl" begin
    m = 7
    N = 2^m # Number of points
    d = 2 # dimension

    u_uniform = rand(N, d) # i.i.d. uniform

    unrandomized_bits = sobol_pts2bits(m, d, 32)
    indices = which_permutation(unrandomized_bits) #32 bit version
    random_bits = similar(unrandomized_bits) # 32 bit version
    nus = NestedUniformScrambler(unrandomized_bits, indices)
    lms = LinearMatrixScrambler(unrandomized_bits)

    u_sob = dropdims(mapslices(bits2unif, unrandomized_bits, dims=3), dims=3)
    u_nus = copy(u_sob)
    u_lms = copy(u_sob)
    u_shift = copy(u_sob)

    scramble!(u_nus, random_bits, nus)
    scramble!(u_lms, random_bits, lms)
    shift!(u_shift)
end
