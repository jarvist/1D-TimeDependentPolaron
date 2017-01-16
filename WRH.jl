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

"""
    randH(SiteEnergy, Edisorder, Jdisorder, modelJ, N)

Liberated from 'Sturm': https://github.com/jarvist/Teclo/blob/master/Sturm.jl 

Generate a random tridiagonal TightBinding Hamiltonian, in a form suitable for the Sturm sequence

# Given: 
 *   SiteEnergy - scalar eV; reference for site energy of monomer
 *   Edisorder - scalar eV ; amount of Gaussian / normal energetic disorder, for trace of Hamiltonian
 *   Jdisorder - scalar eV
 *   modelJ(theta) - function, takes degrees, returns eV ; model for the transfer integral (e.g. E=(J0*cos(thetas*pi/180)).^2 )
 *   N - integar ; size of diagonal of Hamiltonian
"""
function randH(SiteEnergy, Edisorder, Jdisorder, modelJ, N)
# Random Trace / diagonal elements
    S=SiteEnergy + Edisorder*randn(N)
# Random Off-diag elements
    E=modelJ(0) + Jdisorder*randn(N-1)
    return (S,E)
end

const dipolestrength=0.2 

"Calculate site-energies for sites, from potential generated by dipoles at other sites."
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

"Step-forwards in time, and allow dipoles (dielectric response) of sites to respond to electron density."
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

"Build Hamiltonian from Dipoles (via SiteEnergyFromDipoles), diagonalise and update dipoles responding to ground state electron density. A wrapper function."
function AdiabaticPropagation(dipoles,E)
    S=SiteEnergyFromDipoles(dipoles)

    H=diagm(E,-1)+diagm(S)+diagm(E,1) #build full Hamiltonian
    psi=eigvecs(H)[:,1] # gnd state
    
    density=psi.^2
    dipoles=DipolesFromDensity(dipoles,density)

    return S,psi,density,dipoles
end

const hbar=1.0
"Warning - not currently unitary! Propagate Wavefunction directly with Hamiltonian and time dependent Schrodinger equation."
function TimeDependentPropagation(psi,H,dt;slices::Int=1,decompose::Bool=false) # propagate directly using full Hamiltonian=T+V
    # Decompose unitary evolution into this many slices
    dt=dt/slices

    if decompose
        S,J=Decompose_H(H) # split into diagonal and off-diag terms
        # Trotter decomposition
        U=expm(-im*J*dt/2*hbar)*expm(-im*S*dt/hbar)*expm(-im*J*dt/2*hbar)
    else
        U=expm(-im*H*dt/hbar)
    end

    for i=2:slices
        if decompose
            S,J=Decompose_H(H) # split into diagonal and off-diag terms
            # Trotter decomposition
            U*=expm(-im*J*dt/2*hbar)*expm(-im*S*dt/hbar)*expm(-im*J*dt/2*hbar)
        else
            U*=expm(-im*H*dt/hbar)
        end
    end

    @printf(" (matrix-squaring slices: %d) U: \n",slices)
    display(U)
    println("\n\tUU': (in any sane world this should be = Identity )\n")
    display(U*U') # Why you no unitary?
    psi=U*psi

    println("\nPre normalised Norm of psi: ",norm(psi))
    psi/=norm(abs(psi.^2)) # Normalise propagated wavefunction
    return psi
end

function TimeDependentPropagationDecompose(psi,H,dt) # propagate directly using full Hamiltonian=T+V
    S,J=Decompose_H(H) # split into diagonal and off-diag terms
    U=exp(-im*S*dt/hbar)*exp(-im*J*dt/hbar)

    @printf(" (matrix-squaring slices: %d) U: \n",slices)
    display(U)
    println("\n\tUU': (in any sane world this should be = Identity )\n")
    display(U*U') # Why you no unitary?
    psi=U*psi

    println("\nPre normalised Norm of psi: ",norm(psi))
    psi/=norm(psi.^2) # Normalise propagated wavefunction
    return psi
end

"Propagate wavefunction directly from eigenergy, and time dependent Schrodinger equation."
function TimeDependentPropagation(psi,H,dt,E) # propagate using eigenvalue
    psi=exp(-im*E*dt/hbar)*psi
    println("Pre normalised Norm of psi: ",norm(psi))
#    psi/=norm(psi*psi') # Normalise propagated wavefunction
    return psi
end

"Self-consistent response of the lattice with unitary (time dependent) evolution of the wavefunction. Nb: Doesn't work currently - need a working Hamiltonian based unitary operator!"
function UnitaryPropagation(dipoles,E,psi,dt;slices::Int=1)
    S=SiteEnergyFromDipoles(dipoles)

    H=diagm(E,-1)+diagm(S)+diagm(E,1) 
    
    En=eigvals(H)[1]
    println("Eigvals: ",En)
    #psi=TimeDependentPropagation(psi,H,dt,En)
    psi=TimeDependentPropagation(psi,H,dt,slices=slices,decompose=true)
 
    density=abs(psi.^2) # can be Complex!
    dipoles=DipolesFromDensity(dipoles,density)
    
    return S,psi,density,dipoles
end
 
using UnicodePlots # Take it back to the 80s

"Wrapper function to pretty-print and plot (UnicodePlots) relevant items of interest."
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

" Decompose Hamiltonian into Diagonal/S/PE and Off-diag/J/KE elements"
function Decompose_H(H)
    S=eye(N).*H # elementwise to select for just diagonal terms
    J=H-S
    return S,J
end

" Plot spectrum of (H)amiltonian, other useful info."
function Plot_H(H)
    # display -> use Julia's prettyprinting, a la the REPL
    display(H)
    println("Eigenvalues: ")
    display(eigvals(H))
    println("Min Eigenvalue: ",eigmin(H))
    
    myvecs=eigvecs(H)
    myvals=eigvals(H)

    println("Eig Vecs: ")
    display(myvecs)

    siteenergies=[H[i,i] for i in 1:size(H,2)]
    myplot=lineplot(siteenergies,name="Site Energies",color=:red,width=80)
    
    for i in 1:size(myvecs,2)
        vec=myvecs[i,:] + myvals
        println("Vec: ",vec)
        lineplot!(myplot,vec,name="Eigenvectors",color=:green)
    end
    print(myplot)
end

function main()
    # generates separate (S)ite (diagonal) and (E)-offdiagonal terms of Tight Binding Hamiltonian
    S,E=randH(5.0,Edisorder, Jdisorder, modelJ, N)

    #println("Full square matrix H");
    H=diagm(E,-1)+diagm(S)+diagm(E,1) #build full matrix from diagonal elements; for comparison
    psi=eigvecs(H)[:,1] # gnd state
 
    Plot_H(H)

    ## Testing
    dipoles=zeros(N)

    # Self consistent field loop
    for i in 1:1
        @printf("\n\tSCF loop: %d\n",i)
        S,psi,density,dipoles = AdiabaticPropagation(dipoles,E)
        Plot_S_psi_density_dipoles(S,psi,density,dipoles)
    end

    H=diagm(E,-1)+diagm(S)+diagm(E,1) #build full matrix from diagonal elements; for comparison
    psi=eigvecs(H)[:,1] # gnd state

    println("Hamiltonian: ")
    display(H) # Nb: display does the pretty-print which you see at the julia> command line.
    println("Eigvecs: ")
    display(eigvecs(H))

    dt=1 # Time step; not sure of units currently; hbar set to 1 above, energies in eV

    println("Psi: ",psi)
    #myplot=lineplot(psi,name="Psi",color=:red,width=80,ylim=[-1,1])
    for i in 1:25
        @printf("\n\tUnitary Propagation Loop: %d\n",i)
#        psi=eigvecs(H)[:,1] # gnd state
        S,psi,density,dipoles = UnitaryPropagation(dipoles,E,psi,dt,slices=i)
        Plot_S_psi_density_dipoles(S,real(psi),density,dipoles)
    end
end

main() # Party like it's C99!
