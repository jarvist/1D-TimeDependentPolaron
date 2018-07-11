# Rose Teague 06/07/2018
# UpdateEnergyTest.jl
#
# Tests the function 'UpdateDipoles' using eigenstate of the system.
push!(LOAD_PATH,"../src/")

using TheDancer
using Base.Test

N, Edisorder, Jdisorder, modelJ, B, dipolestrength, r = init!(50, 0.0, 0.0, 0.01, 298, 20, 1)

function main()
    S,E,H,psi,dipoles=prepare_model()
    density = abs.(psi.^2)

    for i in 1:500
        S,H,psi,density,dipoles = AdiabaticPropagation(S,dipoles,E,false)
    end

    Field = FieldFromDensity(density)
    Dipoles = UpdateDipoles(Field, dipoles, 1)

    ε = 1e-6
    @test 0.001935443 ≈ maximum(Field) atol = ε
end

main()
