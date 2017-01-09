# WRH.jl
# Simple Julia codes to play with 1D Polaron propagation

# What's your name man? -- https://www.youtube.com/watch?v=SZXHoWwBcDc
println("Time is said to have only one dimension, and space to have three dimensions. - William Rowan Hamilton")

# These codes simulate the formation of a Polaron in a 1D Tight-Binding model. 
# The model is an |N> site model, where the sites are expected to be ~molecular~ units in real space. 
# The response of the lattice is modelled as a dipole response, for the polarisation of the dielectric modes.
# Electronic structure is by tight-binding, with site energies perturbed by the response of the lattice, 
# and a parameterised 'transfer integral' for the kinetic energy between nearest neighbour sites

#Bolztmann's constant in units of eV - thereby all the potentials (of functional form or tabulated data) in units eV
const kB=8.6173324E-5

const N=10^1 # Number of sites in model
             # --> Size of tridiagonal Hamiltonian constructed; larger value -> better statistics, slower runtime

const Edisorder=0.0 # Energetic disorder eV, Gaussian form
const Jdisorder=0.0 # Transfer integral disorder, eV. 

# Model setup
const J0=0.1
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

const dipolestrength=0.2 

function SiteEnergyFromDipoles(dipoles)
    S=zeros(N)
    for i in 1:N
        for j in 1:N
            if (j==i) 
                continue # avoid infinity self energies
            end
            S[i]+=dipolestrength*dipoles[j]/(i-j)^3 # Contribution to site energy (1 e- at site) from dipoles
        end
#        @printf("Site: i %d SiteEnergy: S[i] %f\n",i,S[i])
    end
    S
end

const dampening=0.5 # How much to update in each step

function DipolesFromDensity(dipoles,density)
    for i in 1:N
        relaxeddipole=0
        for j in 1:N
            if (j==i)
                continue # avoid infinite self energies
            end
            relaxeddipole+=density[j]/(i-j)^1 # How much do the dipoles respond to the electron density?
        end

        dipoles[i]+=dampening*(relaxeddipole-dipoles[i]) # Approaches the infinite time limit via Zeno's dichotomy
    end
    dipoles
end

function AdiabaticPropagation(dipoles,E)
    S=SiteEnergyFromDipoles(dipoles)

    H=diagm(E,-1)+diagm(S)+diagm(E,1) #build full Hamiltonian
    psi=eigvecs(H)[:,1] # gnd state
    
    density=psi.^2
    dipoles=DipolesFromDensity(dipoles,density)

    return S,psi,density,dipoles
end

const hbar=1.0
# Not currently unitary!
function TimeDependentPropagation(psi,H,dt) # propagate directly using full Hamiltonian=T+V
    psi=exp(-im*H*dt/hbar)*psi
    println("Pre normalised Norm of psi: ",norm(psi))
    psi/=norm(psi.^2) # Normalise propagated wavefunction
    return psi
end

function TimeDependentPropagation(psi,H,dt,E) # propagate using eigenvalue
    psi=exp(-im*E*dt/hbar)*psi
    println("Pre normalised Norm of psi: ",norm(psi))
#    psi/=norm(psi*psi') # Normalise propagated wavefunction
    return psi
end

using UnicodePlots # Take it back to the 80s

function Plot_S_psi_density_dipoles(S,psi,density,dipoles)
    println("Site energies: ",S)
    println("Psi: ",psi)
    println("Electron density: ",density)
    println("Dipoles: ",dipoles)

    myplot=lineplot(S,name="Site Energies",color=:red,width=80,ylim=[-1,1])
    lineplot!(myplot,psi,name="Psi",color=:green)
    lineplot!(myplot,density,name="Electon Density",color=:yellow)
    lineplot!(myplot,dipoles,name="Dipoles",color=:blue)
    print(myplot)
end

function main()
    # generates separate (S)ite (diagonal) and (E)-offdiagonal terms of Tight Binding Hamiltonian
    S,E=randH(5.0,Edisorder, Jdisorder, modelJ, N)

    #println("Full square matrix H");
    H=diagm(E,-1)+diagm(S)+diagm(E,1) #build full matrix from diagonal elements; for comparison
    psi=eigvecs(H)[:,1] # gnd state
 
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
    for i in 1:10
        @printf("\n\tSCF loop: %d\n",i)
        S,psi,density,dipoles = AdiabaticPropagation(dipoles,E)
        Plot_S_psi_density_dipoles(S,psi,density,dipoles)
    end

    dipoles=zeros(N)
    S,psi,density,dipoles = AdiabaticPropagation(dipoles,E)

    H=diagm(E,-1)+diagm(S)+diagm(E,1) #build full matrix from diagonal elements; for comparison
    psi=eigvecs(H)[:,1] # gnd state

    println("Hamiltonian: ")
    display(H) # Nb: display does the pretty-print which you see at the julia> command line.
    println("Eigvecs: ")
    display(eigvecs(H))

    dt=1 # Time step; not sure of units currently; hbar set to 1 above, energies in eV


    println("Psi: ",psi)
    #myplot=lineplot(psi,name="Psi",color=:red,width=80,ylim=[-1,1])
    for i in 1:20
        #psi=TimeDependentPropagation(psi,H,dt)
        
        psi=TimeDependentPropagation(psi,H,dt,eigvals(H)[1]) # Eigenvalue version
        
        println("TimeDependentPropagation Psi: ")
        display(psi)
        println("Psi.^2 : ")
        display(psi.^2)
        println()

        myplot=lineplot(real(psi),ylim=[-1,1],color=:red,width=80) # psi, wavefunction
        lineplot!(myplot,abs(psi.^2),color=:yellow) #psi^2, density
        print(myplot)
    end
end

main() # Party like it's C99!
