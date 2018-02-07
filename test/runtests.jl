push!(LOAD_PATH,"../src/") # load module from local directory

using TD 
using Base.Test

SCFthenUnitary(0.2, 50, 50)

println("Tests finished succesfully.")

