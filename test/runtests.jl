push!(LOAD_PATH,"../src/") # load module from local directory

using TheDancer 
using Base.Test

#init!()
#SCFthenUnitary(0.2, 10, 500) # This isn't really a test - it just runs the full GUI!

model=TightBindingModel()
init!(model.N) # temporary kludge until every function takes TightBindingModel as argument 


UnitaryPropagation!(model, 0.1; slices=1)

println("Tests finished succesfully.")

