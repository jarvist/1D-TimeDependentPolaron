# WRH.jl
# Simple Julia codes to play with 1D Polaron propogation

# What's your name man? -- https://www.youtube.com/watch?v=SZXHoWwBcDc
println("Time is said to have only one dimension, and space to have three dimensions. - William Rowan Hamilton")

#Bolztmann's constant in units of eV - thereby all the potentials (of functional form or tabulated data) in units eV
const kB=8.6173324E-5

const N=10^1 # Length of tridiagonal Hamiltonian constructed; larger value -> better statistics, slower runtime

const Edisorder=0.0 # Energetic disorder eV, Gaussian form
const Jdisorder=0.0 # Transfer integral disorder, eV. 

# Model setup
const J0=0.8
modelJ(theta) = J0*cos(theta*pi/180.0).^2

const T=300
const B=1/(T*kB) #300K * k_B in eV

### Liberated from 'Sturm': https://github.com/jarvist/Teclo/blob/master/Sturm.jl ###
# Generate a random tridiagonal TightBinding Hamiltonian, in a form suitable for the Sturm sequence
# Given: 
#   SiteEnergy - scalar eV; reference for site energy of monomer
#   disorder - scalar eV ; amount of Gaussian / normal energetic disorder, for trace of Hamiltonian
#   modelJ(theta) - function, takes degrees, returns eV ; model for the transfer integral (e.g. E=(J0*cos(thetas*pi/180)).^2 )
#   B - scalar (units?); Thermodynamic (B)eta, used to populate Probability Density Function
#   Z - scalar (units?); Partition function, weighting for absolute Boltzmann populations
#   U - function(theta angle); Free energy function, used to generate Bolztmann populations
#   N - integar ; size of diagonal of Hamiltonian
function randH(SiteEnergy, Edisorder, Jdisorder, modelJ, N)
# Random Trace / diagonal elements
    S=SiteEnergy + Edisorder*randn(N)
# Random Off-diag elements
    E=modelJ(0) + Jdisorder*randn(N-1)
    return (S,E)
end

function main()
    # generates separate (S)ite (diagonal) and (E)-offdiagonal terms of Tight Binding Hamiltonian
    S,E=randH(5.0,Edisorder, Jdisorder, modelJ, N)

    #println("Full square matrix H");
    H=diagm(E,-1)+diagm(S)+diagm(E,1) #build full matrix from diagonal elements; for comparison
    println(H)

    println("Eigenvalues")
    println(eigvals(H))
    println("Min Eigenvalue")
    println(eigmin(H))
end

main() # Party like it's C99!

