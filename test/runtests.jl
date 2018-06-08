push!(LOAD_PATH,"../src/") # load module from local directory

using TheDancer 
using Base.Test

SCFthenUnitary(0.2, 10, 500) # This isn't really a test - it just runs the full GUI!

println("Tests finished succesfully.")

