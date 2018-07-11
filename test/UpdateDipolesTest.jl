# Rose Teague 06/07/2018
# UpdateEnergyTest.jl
#
# Tests the function 'UpdateDipoles' using eigenstate of the system.
push!(LOAD_PATH,"../src/")

using TheDancer
using Base.Test

N, Edisorder, Jdisorder, modelJ, B, dipolestrength, r = init!(1, 0.0, 0.0, 0.01, 298, 20)

"""
Checks if dipoles are being updated correctly.

Tests against a single dipole with alpha polarisability being induced in a
unitary field. The dipole should then be given by alpha_test 
"""
function main()

    alpha_test = 2.63 #Polarisability volume of HCl
    Field = [1] # Unitary Field
    dipole = 0

    Dipoles = UpdateDipole(Field, dipole, 1, alpha_test)

    ε = 1e-2
    @test 2.63 ≈ Dipoles[1] atol = ε
end

main()
