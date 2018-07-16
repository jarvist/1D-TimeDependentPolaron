push!(LOAD_PATH,"../src/")

using TheDancer
using Plots
gr()

N, Edisorder, Jdisorder, modelJ, B, dipolestrength, r = init!(50, 0.0, 0.0, 0.5, 298, 7.5)

function main()
    S,E,H,psi,dipoles=prepare_model()
    density = abs.(psi.^2)
    dampening = 0.05

    SCFcycles=200
    Unitarycycles=500
    TotCycles = SCFcycles+Unitarycycles

    Energies = SCFthenUnitary(dampening, SCFcycles, Unitarycycles,true)

    plot(Energies'[1:TotCycles])
      for i in 2:10
          plot!(Energies'[i*TotCycles+1:(i+1)*TotCycles])
      end
     gui()
end


main()
