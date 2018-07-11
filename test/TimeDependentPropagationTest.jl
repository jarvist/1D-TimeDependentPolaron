# Rose Teague 06/07/2018
# AdiabaticPropagationTest.jl
#
# Runs the Adiabatic propagation of an eigenstate
push!(LOAD_PATH,"../src/")

using TheDancer
using Base.Test


N, Edisorder, Jdisorder, modelJ, B, dipolestrength, r = init!(10, 0.0, 0.0, 0.001, 298, 20)

"""
Test full time-dependent propagation for a simple 3-site system. With small transfer
integral (J=0.001eV).

Unitary propagation matrix should be unitary - test U'U = I

Prepare system in an eigenstate, unitary propagation should leave psi.^2 unchanged.
"""
function main()
    S,E,H,psi,dipoles=prepare_model()
    density = abs.(psi.^2)
    psi_1 = psi

    for i in 1:200
        psi,U=TimeDependentPropagation(psi,H,1,decompose=false,verbose=false,test=true)
        density = abs.(psi.^2)
        unit = sum(U'U)/N
        e = 1E-16
        @test 1 ≈ unit atol = e
    end

    diff = sum(density - abs.(psi_1.^2))

    ε = 1E-7
    @test  0 ≈ diff atol = ε
end

main()
