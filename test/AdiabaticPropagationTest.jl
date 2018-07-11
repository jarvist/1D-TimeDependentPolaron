# Rose Teague 06/07/2018
# AdiabaticPropagationTest.jl
#
# Runs the Adiabatic propagation of an eigenstate
push!(LOAD_PATH,"../src/")

using TheDancer
using Base.Test

N, Edisorder, Jdisorder, modelJ, B, dipolestrength, r = init!(3, 0.0, 0.0, 0.001, 298, 20)

"""
Test full adiabatic propagation for a simple 3-site system. With small transfer
integral (J=0.001eV) the electron density should end up all at the central
site.
"""
function main()
    S,E,H,psi,dipoles=prepare_model()
    density = abs.(psi.^2)

    for i in 1:200
        S,H,psi,density,dipoles = AdiabaticPropagation(S,dipoles,E,false)
    end

    ε = 1E-6
    @test  1 ≈ density[2] atol = ε
end

main()
