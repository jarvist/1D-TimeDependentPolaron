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
#   Edisorder - scalar eV ; amount of Gaussian / normal energetic disorder, for trace of Hamiltonian
#   Jdisorder - scalar eV
#   modelJ(theta) - function, takes degrees, returns eV ; model for the transfer integral (e.g. E=(J0*cos(thetas*pi/180)).^2 )
#   N - integar ; size of diagonal of Hamiltonian
function randH(SiteEnergy, Edisorder, Jdisorder, modelJ, N)
# Random Trace / diagonal elements
    S=SiteEnergy + Edisorder*randn(N)
# Random Off-diag elements
    E=modelJ(0) + Jdisorder*randn(N-1)
    return (S,E)
end

function SiteEnergyFromDipoles(dipoles)
    S=zeros(N)
    for i in 1:N
        for j in 1:N
            if (j==i) 
                continue # avoid infinity self energies
            end
            S[i]+=dipoles[j]/(i-j)^3 # Contribution to site energy (1 e- at site) from dipoles
        end
#        @printf("Site: i %d SiteEnergy: S[i] %f\n",i,S[i])
    end
    S
end

const dampening=0.1

function DipolesFromDensity(dipoles,density)
    for i in 1:N
        for j in 1:N
            if (j==i)
                continue # avoid infinite self energies
            end
            dipoles[i]+=dampening*density[j]/(i-j)^1 # How much do the dipoles respond to the electron density?
        end
    end
    dipoles
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

    println("Eig Vecs")
    println(eigvecs(H))

    ## Testing
    dipoles=zeros(N)

    # Self consistent field loop
    for i in 1:100
        @printf("\n\tSCF loop: %d\n",i)
        S=SiteEnergyFromDipoles(dipoles)
        println("Site energies: ",S)
        H=diagm(E,-1)+diagm(S)+diagm(E,1) #build full Hamiltonian
        psi=eigvecs(H)[1,:] # gnd state
        println("Psi: ",psi)
        density=psi.^2
        println("Electron density: ",density)
        dipoles=DipolesFromDensity(dipoles,density)
        println("Dipoles: ",dipoles)
    end
end

main() # Party like it's C99!

