# TheDancer.jl
# Simple Julia codes to play with 1D Polaron propagation
# By Jarvist Moore Frost and Rose Teague (2017--2018)

# These codes simulate the formation of a Polaron in a 1D Tight-Binding model.
# The model is an |N> site model, where the sites are expected to be ~molecular~ units in real space.
# The response of the lattice is modelled as a dipole response, for the polarisation of the dielectric modes.
# Electronic structure is by tight-binding, with site energies perturbed by the response of the lattice,
# and a parameterised 'transfer integral' for the kinetic energy between nearest neighbour sites =#

module TheDancer
println("\t\"He came riding fast, like a phoenix out of fire-flames.\" -- The Dancer, PJ Harvey" )

export nondispersive_wavepacket, planewave, randH, decompose_H
export SiteEnergyFromDipoles, FieldFromDensity, UpdateDipole, dipoles_by_relaxation
export AdiabaticPropagation, UnitaryPropagation, TimeDependentPropagation
export plot_S_psi_density_dipoles, plot_H, outputpng, overlap
export SCFthenUnitary, init!, prepare_model

const UsePlots=true
using Plots
gr()

include("initialise.jl")
include("models.jl")
include("propagators.jl")
include("simulations.jl")

end
