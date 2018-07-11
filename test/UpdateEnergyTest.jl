# Rose Teague 06/07/2018
# UpdateEnergyTest.jl
#
# Tests the function 'SiteEnergyFromDipoles' using eigenstate of the system.
push!(LOAD_PATH,"../src/")

using TheDancer
using Base.Test

N, Edisorder, Jdisorder, modelJ, B, dipolestrength, r = init!(3, 0.0, 0.0, 0.0, 298, 20)

"""
Test the on-site energy is correctly calculated as the potential energy felt by
a test charge on site i, in the presence of the electric fields established by
all other sites.

Testing in a symmetric 3-site environment, energy at the central position should
be 0 as effects from dipoles on either side cancel each other out
"""

function main()
    dipoles = [1.0,1.0,1.0]
    E = [0.0,0.0]
    S = [0.0,0.0,0.0]
    S,H = SiteEnergyFromDipoles(dipoles,S,E)

    ε = 1e-2
    @test 0.0 ≈ S[2] atol = ε
end

main()
