push!(LOAD_PATH,"../src/")

using TheDancer
using Base.Test
using Plots
gr()

N, Edisorder, Jdisorder, modelJ, B, dipolestrength, r = init!(3, 0.0, 0.0, 0.001, 298, 20)

"""
Test to compare old and new methods of calculating the dipole response.

Test will pass if the two methods are equal to within 1% of each other, else
test will fail.
"""
function New(dipoles, density, dampening)
    Field = FieldFromDensity(density)
    dipoles = UpdateDipole(Field, dipoles, dampening)
end

function main()
    S,E,H,psi,dipoles=prepare_model()
    density = abs.(psi.^2)

    dipoles_new = New(dipoles, density, 0.5)
    dipoles_old = dipoles_by_relaxation(dipoles, density, 0.5)

    diff = sum(abs.(dipoles_old - dipoles_new))
    percentage_diff = diff/sum(dipoles_new)

    ε = 1E-2
    @test  0 ≈ percentage_diff atol = ε


end

main()
