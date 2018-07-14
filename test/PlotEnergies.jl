push!(LOAD_PATH,"../src/")

using TheDancer
using Base.Test
using Plots
gr()

N, Edisorder, Jdisorder, modelJ, B, dipolestrength, r = init!(50, 0.0, 0.0, 0.005, 298, 7.5)

function main()
    S,E,H,psi,dipoles=prepare_model()
    density = abs.(psi.^2)
    dampening = 0.07

    SCFcycles=200
    Unitarycycles=500
    TotCycles = SCFcycles+Unitarycycles

    Energies = SCFthenUnitary(dampening, SCFcycles, Unitarycycles,true)
    #
    # dt = 1000
    # for i in 1:100
    #     #S,H,psi,density,dipoles = AdiabaticPropagation(S,dipoles,density,E,dampening,false,false)
    #     S,H,psi,density,dipoles = UnitaryPropagation(dipoles,density,S,E,psi,dt,dampening)
    #     Energies[(i-1)*50+1:i*50] = eigvals(H)
    #
    # end

    plot(Energies'[1:TotCycles])
      # for i in 2:2
      #     plot!(Energies'[i*TotCycles+1:(i+1)*TotCycles])
      # end
     gui()
end


main()
