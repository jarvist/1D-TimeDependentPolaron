#= TD.jl
# Simple Julia codes to play with 1D Polaron propagation

# These codes simulate the formation of a Polaron in a 1D Tight-Binding model.
# The model is an |N> site model, where the sites are expected to be ~molecular~ units in real space.
# The response of the lattice is modelled as a dipole response, for the polarisation of the dielectric modes.
# Electronic structure is by tight-binding, with site energies perturbed by the response of the lattice,
# and a parameterised 'transfer integral' for the kinetic energy between nearest neighbour sites =#

module TheDancer
#=
Module to call all the functions required for polaron surface hopping simulation.
Set all system parameters in this module.
=#

println("\t\"He came riding fast, like a phoenix out of fire-flames.\" -- The Dancer, PJ Harvey" )

export randH, SiteEnergyFromDipoles, DipolesFromDensity, TimeDependentPropagation
export Plot_S_psi_density_dipoles, overlap, Decompose_H, Plot_H
export nondispersive_wavepacket, planewave, prepare_model, outputpng
export SCFthenUnitary, init

const UsePlots=true
using Plots
gr()

include("init.jl")
include("models.jl")
include("propagators.jl")
include("simulations.jl")

end
