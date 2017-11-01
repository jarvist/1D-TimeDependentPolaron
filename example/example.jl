push!(LOAD_PATH,"../src/") # load module from local directory

using TD

function main()
    SCFcycles=50
    Unitarycycles=1000

    for dampening in [0.07,0.025,0.05]
        SCFthenUnitary(dampening, SCFcycles, Unitarycycles, PNG=false) # PNG=true for .pngs for movie making
    end
end

main() # Party like it's C99!
