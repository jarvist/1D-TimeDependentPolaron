# TD.jl
# Simple Julia codes to play with 1D Polaron propagation

# These codes simulate the formation of a Polaron in a 1D Tight-Binding model. 
# The model is an |N> site model, where the sites are expected to be ~molecular~ units in real space. 
# The response of the lattice is modelled as a dipole response, for the polarisation of the dielectric modes.
# Electronic structure is by tight-binding, with site energies perturbed by the response of the lattice, 
# and a parameterised 'transfer integral' for the kinetic energy between nearest neighbour sites

module TheDancer 

println("\t\"He came riding fast, like a phoenix out of fire-flames.\" -- The Dancer, PJ Harvey" )

export randH, SiteEnergyFromDipoles, DipolesFromDensity, TimeDependentPropagation
export Plot_S_psi_density_dipoles, overlap, Decompose_H, Plot_H
export nondispersive_wavepacket, planewave, prepare_model, outputpng
export SCFthenUnitary

const UsePlots=true
if UsePlots
    using Plots # This is a meta-plotting package, wrapping around multiple backends
    gr() # High performance
    #unicodeplots() # Take it back to the 80s
else
    using UnicodePlots # Use UnicodePlots directly 
end

# CONSTANTS
# - genuine constants are fine
# - FIXME: put simulation parameters elsewhere; make things more functionally pure.

#Bolztmann's constant in units of eV - thereby all the potentials (of functional form or tabulated data) in units eV
const kB=8.6173324E-5
const hbar=1.0

const N=50 # Number of sites in model
             # --> Size of tridiagonal Hamiltonian constructed; larger value -> better statistics, slower runtime
const Edisorder=0.0 # Energetic disorder eV, Gaussian form
const Jdisorder=0.0 # Transfer integral disorder, eV. 

# Model setup
const J0=0.1
modelJ(θ) = J0*cos(θ*π/180.0).^2

const T=300
const B=1/(T*kB) #300K * k_B in eV

# This effectively reduces down to the 'alpha' parameter in the Frohlich polaron Hamiltonian
const dipolestrength=0.2 

include("models.jl")
include("propagators.jl")
include("simulations.jl")

end

