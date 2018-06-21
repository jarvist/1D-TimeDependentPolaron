push!(LOAD_PATH,"../")

include("testing_dipole_response.jl")
using Plots

#'''SetUp'''
const kB=8.6173324E-5
const hbar=1.0

const N=50 # Number of sites in model
             # --> Size of tridiagonal Hamiltonian constructed; larger value -> better statistics, slower runtime
const Edisorder=0.0 # Energetic disorder eV, Gaussian form
const Jdisorder=0.0 # Transfer integral disorder, eV.

# Model setup
const J0=0.001
modelJ(θ) = J0*cos(θ*π/180.0).^2

const T=300
const B=1/(T*kB) #300K * k_B in eV

# This effectively reduces down to the 'alpha' parameter in the Frohlich polaron Hamiltonian
const dipolestrength=0.2

#'''prepare'''
S,E,H,psi,density = prepare_model()
dipole = zeros(N)
Field = FieldFromDensity(density)
S,E,H = UpdateEnergy(dipole, Field, S, E)
dipole = UpdateDipole(Field,density)


psi=psi+nondispersive_wavepacket(20,4.0)
normalise = norm(psi)
psi = psi/normalise
density = conj(psi).*psi
#plot(real(density))
#gui()
Field = FieldFromDensity(density)
dipole = UpdateDipole(Field,density)
S,E,H = UpdateEnergy(dipole, Field, S, E)
norm_dipole = norm(dipole)
#plot(real(dipole/norm_dipole))
#gui()


#'''propagate'''
for i in 1:1000
    norm_dipole = norm(dipole)
    plot(real(dipole/norm_dipole))
    plot!(imag(psi), lims=(-1,1))
    plot!(real(psi))
    xlims!(1,N)
    ylims!(-1,1)
    plot!(real(density))
    gui()
    psi, density = Propagate(H,psi,1)
    Field = FieldFromDensity(real(density))
    dipole = UpdateDipole(Field, density)
    S,E,H = UpdateEnergy(dipole, Field, S, E)

end
