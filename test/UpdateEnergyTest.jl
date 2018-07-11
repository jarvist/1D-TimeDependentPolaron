# Rose Teague 06/07/2018
# UpdateEnergyTest.jl
#
# Tests the function 'SiteEnergyFromDipoles' using eigenstate of the system.
push!(LOAD_PATH,"../src/")

using TheDancer
using Base.Test

N, Edisorder, Jdisorder, modelJ, B, dipolestrength, r = init!(50, 0.0, 0.0, 0.001, 298, 20)
S,E,H,psi,dipoles=prepare_model()
