# TD.jl
# Simple Julia codes to play with 1D Polaron propagation
println("\t\"He came riding fast, like a phoenix out of fire-flames.\" -- The Dancer, PJ Harvey" )

# These codes simulate the formation of a Polaron in a 1D Tight-Binding model. 
# The model is an |N> site model, where the sites are expected to be ~molecular~ units in real space. 
# The response of the lattice is modelled as a dipole response, for the polarisation of the dielectric modes.
# Electronic structure is by tight-binding, with site energies perturbed by the response of the lattice, 
# and a parameterised 'transfer integral' for the kinetic energy between nearest neighbour sites

#Bolztmann's constant in units of eV - thereby all the potentials (of functional form or tabulated data) in units eV
const kB=8.6173324E-5

const N=50 # Number of sites in model
             # --> Size of tridiagonal Hamiltonian constructed; larger value -> better statistics, slower runtime

const Edisorder=0.0 # Energetic disorder eV, Gaussian form
const Jdisorder=0.0 # Transfer integral disorder, eV. 

# Model setup
const J0=0.1
modelJ(theta) = J0*cos(theta*pi/180.0).^2

const T=300
const B=1/(T*kB) #300K * k_B in eV

const hbar=1.0

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

# This effectively reduces down to the 'alpha' parameter in the Frohlich polaron Hamiltonian
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

# What fraction of the dipole response to update with each step of the electronic degree of freedom. 
# This is the rate of response of the dipole lattice c.f. an updated step of the electronic degree of freedom
# And can be imagined as a solution to a heavily damped Simple-Harmonic-Oscillator --> exponential (half life) solution
# dampening=0.025 

"Step-forwards in time, and allow dipoles (dielectric response) of sites to respond to electron density."
function DipolesFromDensity(dipoles,density,dampening)
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
function AdiabaticPropagation(dipoles,E,dampening)
    S=SiteEnergyFromDipoles(dipoles)

    H=diagm(E,-1)+diagm(S)+diagm(E,1) #build full Hamiltonian
    psi=eigvecs(H)[:,1] # gnd state
    
    density=psi.^2
    dipoles=DipolesFromDensity(dipoles,density,dampening)

    return S,psi,density,dipoles
end

"Self-consistent response of the lattice with unitary (time dependent) evolution of the wavefunction. "
function UnitaryPropagation(dipoles,E,psi,dt,dampening;slices::Int=1)
    S=SiteEnergyFromDipoles(dipoles)

    H=diagm(E,-1)+diagm(S)+diagm(E,1) 
    
    En=eigvals(H)[1]
    println("First Eigval (adibatic): ",En)

    #psi=eigvecs(H)[:,1] # gnd state; reproduces adiabtic state (eigval) above
    println("State energy: <psi|H|psi> = ",psi'*H*psi)
    
    psi=TimeDependentPropagation(psi,H,dt,slices=slices,decompose=false,verbose=false)
 
    density=abs(psi.^2) # can be Complex!

    dipoles=DipolesFromDensity(dipoles,density,dampening)
    
    return S,psi,density,dipoles
end
 
"""
    TimeDependentPropagation(psi,H,dt;slices::Int=1,decompose::Bool=false,verbose::Bool=false)

Propagate Wavefunction directly with Hamiltonian and time dependent Schrodinger equation.

Psi (N) is the wavefunction; H the Hamiltonian (NxN); dt the length of time to
propgate along; slices is how many slices to decompose the Unitary operator
into; verbose sets the display of info on the unitary nature of U, whether U*U'
~= I.  
"""
function TimeDependentPropagation(psi,H,dt;slices::Int=1,decompose::Bool=false,verbose::Bool=false) # propagate directly using full Hamiltonian=T+V
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

    psi=U*psi # OK; we've built out unitary time evolution operator U, now apply it

    if verbose # for debugging / introspection
        @printf(" (matrix-squaring slices: %d) U: \n",slices)
        display(U)
        println("\n\tUU': (in any sane world this should be = Identity )\n")
        display(U*U') # Why you no unitary?

        println("\nPre normalised Norm of psi: ",norm(psi))
    end

    normalisation=sqrt(sum(abs(psi.^2))) # Nb: for normalisation of WAVEFUNCTION; must SQUAREROOT the sum of density
    println("Normalising Psi by dividing by: ",normalisation)
    psi=psi/normalisation # Normalise propagated wavefunction
    return psi
end

"Propagate wavefunction directly from eigenergy, and time dependent Schrodinger equation."
function TimeDependentPropagation(psi,H,dt,E) # propagate using eigenvalue
    psi=exp(-im*E*dt/hbar)*psi
    println("Pre normalised Norm of psi: ",norm(psi))
#    psi/=norm(psi*psi') # Normalise propagated wavefunction
    return psi
end

const UsePlots=true
if UsePlots
    using Plots # This is a meta-plotting package, wrapping around multiple backends
    gr() # High performance
    #unicodeplots() # Take it back to the 80s
else
    using UnicodePlots # Use UnicodePlots directly 
end

"Wrapper function to pretty-print and plot (UnicodePlots) relevant items of interest."
function Plot_S_psi_density_dipoles(S,psi,density,dipoles;title="",verbose::Bool=false)
    if verbose
        println("Site energies: ",S)
        println("Psi: ",psi)
        println("Electron density: ",density," Sum: ",sum(density))
        println("Dipoles: ",dipoles," Sum: ",sum(dipoles))
    end

    # Using Plots interface
    if UsePlots
        plot(S,label="Site Energies",color=:red, width=2, size=(1024,768)) # makes it larger, but also slower!
        
        psioffset=-0.5 # Shift wavefunctions down
        plot!(real(psi)+psioffset,label="Re[Psi]",color=:green , width=2, fill=true, fill=(psioffset,0.3,:green))
        plot!(imag(psi)+psioffset,label="Im[Psi]",color=:pink, width=2, fill=true, fill=(psioffset,0.3,:pink)) 
        
        plot!(density,label="Electon Density",color=:orange, width=4, fill=true, fill=(0,0.3,:orange))

        plot!(dipoles,label="Dipoles",color=:blue, width=2)
        
        xaxis!("Site")
        yaxis!("Psi",[-1,1]) # Fix y-axis limits for animation
        
        title!(title)

        #gui() # Show figure as pop up window...
    else
        # Directly using UnicodePlots
        myplot=lineplot(S,name="Site Energies",color=:red,width=80,ylim=[-1,1])
        lineplot!(myplot,psi,name="Psi",color=:green)
        lineplot!(myplot,density,name="Electon Density",color=:yellow)
        lineplot!(myplot,dipoles,name="Dipoles",color=:blue)
        print(myplot)
    end
end

function overlap(psia,psib)
    overlaps=zeros(eltype(psib), size(psia,1))
    for i in 1:size(psia,1)
        overlap=dot(psia[:,i], psib)
        overlaps[i]=overlap
    end
    return overlaps
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

    # TODO: UPDATE THIS FOR PLOTS.jl
    siteenergies=[H[i,i] for i in 1:size(H,2)]
    myplot=lineplot(siteenergies,name="Site Energies",color=:red,width=80)
    
    for i in 1:size(myvecs,2)
        vec=myvecs[i,:] + myvals
        println("Vec: ",vec)
        lineplot!(myplot,vec,name="Eigenvectors",color=:green)
    end
    print(myplot)
end

" Shamelessly copied from Wikipedia:
https://en.wikipedia.org/wiki/Wave_packet "
function nondispersive_wavepacket(x0, lambda)
    psi=[ exp(-(x-x0)^2)*(cos(2*pi*(x-x0)/lambda) - im * sin(2*pi*(x-x0)/lambda) ) for x=1:N ]
    return psi
end

function planewave(lambda)
    k=2*pi/lambda
    psi=[ exp(-im* k*x) for x=1:N ]
    return psi
end

function prepare_model()
    # generates separate (S)ite (diagonal) and (E)-offdiagonal terms of Tight Binding Hamiltonian
    S,E=randH(5.0,Edisorder, Jdisorder, modelJ, N)
    H=diagm(E,-1)+diagm(S)+diagm(E,1) #build full matrix from diagonal elements; for comparison
    psi=eigvecs(H)[:,1] # 1=gnd state

    #    Plot_H(H) 
    #    TODO: UPDATE THIS FOR PLOTS.jl

    ## Construct initial dipoles structure 
    dipoles=zeros(N)
    #        dipoles=[ (x-N/2)^2*0.01 for x in 1:N ] # quadratic set of dipoles, to give energy slope across device 

    S,E,H,psi,dipoles
end

framecounter=0 # variable to keep track of which frame / plot for later movie we are in

function outputpng()
    global framecounter
    
    Plots.png(@sprintf("%05D.png",framecounter)) # Save plot to PNG file; with XXXXX.png filename
    framecounter=framecounter+1
    #println("Just output frame: $framecounter")
end

function SCFthenUnitary(dampening, SCFcycles, Unitarycycles)

    S,E,H,psi,dipoles=prepare_model()

    # Self consistent field loop; Adiabatic response of lattice + polaron
    # Sets up disorted lattice with polaron, before time-based propagation (should you want it)
    for i in 1:SCFcycles
        @printf("\n\tSCF loop: %d\n",i)
        S,psi,density,dipoles = AdiabaticPropagation(dipoles,E,dampening)
        Plot_S_psi_density_dipoles(S,psi,density,dipoles,title="Dampening: $dampening SCF (Adiabatic): Cycle $i / $SCFcycles")
        outputpng()
    end

    H=diagm(E,-1)+diagm(S)+diagm(E,1) #build full matrix from diagonal elements; for comparison
    psi=eigvecs(H)[:,2] # 1=gnd state, 2=1st excited state, etc.

    # Setup wavefunction for time-based propagation
    psi=psi+nondispersive_wavepacket(15,8.0)

    #        psi=nondispersive_wavepacket(20,8.0) # Centered on 20, with Width (speed?) 8.0
    #        psi=psi+nondispersive_wavepacket(40,-20.0) # Fight of the wavepackets!
    #        psi=planewave(8.0) # Plane wave, lambda=8.0 lattice units

    dt=1.0 # Time step; not sure of units currently; hbar set to 1 above, energies in eV

    println("Psi: ",psi)
    #myplot=lineplot(psi,name="Psi",color=:red,width=80,ylim=[-1,1])
    for i in 1:Unitarycycles
        @printf("\n\tUnitary Propagation Loop: %d\n",i)
        #        psi=eigvecs(H)[:,1] # gnd state
        S,psi,density,dipoles = UnitaryPropagation(dipoles,E,psi,dt,dampening,slices=1)
        Plot_S_psi_density_dipoles(S,psi,density,dipoles,title="Dampening: $dampening TDSE: Cycle $i / $Unitarycycles")
        outputpng()
 
        H=diagm(E,-1)+diagm(S)+diagm(E,1) #build full matrix from diagonal elements; for comparison
#        psi=eigvecs(H)[:,2] # 1=gnd state, 2=1st excited state, etc.
       
        ## Sort of terrible surface hopping implementation
        overlaps=abs(overlap(eigvecs(H),psi))
        closestAdiabaticState=indmax(overlaps) # index of maximum overlap vector
        println("Orbital overlaps; polaron c.f. complete set of states\n",overlaps) # calc and print overlaps of propagated function with full set of adiabatic states.
        println("Max overlap state $closestAdiabaticState with ",overlaps[closestAdiabaticState])

        # Should be minimum switching algorithm / partition function sampling
        # For now we just randomly jump
        if i%50==0
            println(" BORED! JUMPING STATE!")
            psi=eigvecs(H)[:,closestAdiabaticState+1] # TODO: +1 is a lie; but otherwise it just jumps to the ground state again and again 

            Plot_S_psi_density_dipoles(S,psi,density,dipoles,title="JUMPING JACK FLASH TO: $closestAdiabaticState")
            outputpng()
        end
    end
end

function main()
    SCFcycles=50
    Unitarycycles=1000

    for dampening in [0.07,0.025,0.05]
        SCFthenUnitary(dampening, SCFcycles, Unitarycycles)
    end
end

main() # Party like it's C99!
