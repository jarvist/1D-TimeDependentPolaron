push!(LOAD_PATH,"../src/") # load module from local directory

using TheDancer

N, Edisorder, Jdisorder, modelJ, B, dipolestrength, r = init!(50, 0.0, 0.0, 0.001, 298, 7.5)


function main()

    SCFcycles=200
    Unitarycycles=500

    for dampening in [0.07,0.025,0.05]
        SCFthenUnitary(dampening, SCFcycles, Unitarycycles, PNG=false) # PNG=true for .pngs for movie making
    end
end

main() # Party like it's C99!
