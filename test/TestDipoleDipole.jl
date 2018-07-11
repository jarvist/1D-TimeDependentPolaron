push!(LOAD_PATH,"../src/")

using TheDancer
using Base.Test

N, Edisorder, Jdisorder, modelJ, B, dipolestrength, r = init!(3, 0.0, 0.0, 0.001, 298, 20,0.2)


function main()

    densities = zeros(2)
    j=0
    for d_d in [false,true]
        j+=1
        S,E,H,psi,dipoles=prepare_model()
        density = abs.(psi.^2)
        for i in 1:200
            S,H,psi,density,dipoles  = AdiabaticPropagation(S,dipoles,E,false,d_d)
        end
        densities[j] = density[2]
    end


    diff = densities[1] - densities[2]
    println(diff)

    ε = 1E-5
    @test  0 ≈ diff atol = ε

end

main()
