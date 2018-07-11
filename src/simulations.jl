# simulations.jl
# Convenience functions which run and plot the simulation

function overlap(psia,psib)
    overlaps=zeros(eltype(psib), size(psia,1))
    for i in 1:size(psia,1)
        overlap=dot(psia[:,i], psib)
        overlaps[i]=overlap
    end
    return overlaps
end

" Decompose Hamiltonian into Diagonal/S/PE and Off-diag/J/KE elements"
function decompose_H(H)
    S=eye(N).*H # elementwise to select for just diagonal terms
    J=H-S
    return S,J
end

" Plot spectrum of (H)amiltonian, other useful info."
function plot_H(H)
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

framecounter=0 # variable to keep track of which frame / plot for later movie we are in

function outputpng()
    global framecounter

    Plots.png(@sprintf("%05D.png",framecounter)) # Save plot to PNG file; with XXXXX.png filename
    framecounter=framecounter+1
    #println("Just output frame: $framecounter")
end

function SCFthenUnitary(dampening, SCFcycles, Unitarycycles; PNG::Bool=false)

    S,E,H,psi,dipoles=prepare_model()

    # Self consistent field loop; Adiabatic response of lattice + polaron
    # Sets up distorted lattice with polaron, before time-based propagation (should you want it)
    for i in 1:SCFcycles
        @printf("\tSCF loop: %d\n",i)
        S,H,psi,density,dipoles = AdiabaticPropagation(S,dipoles,E)
        plot_model(S,psi,density,dipoles,title="Dampening: $dampening SCF (Adiabatic): Cycle $i / $SCFcycles")
        if PNG outputpng() end
    end

    H=diagm(E,-1)+diagm(S)+diagm(E,1) #build full matrix from diagonal elements; for comparison
    psi=eigvecs(H)[:,2] # 1=gnd state, 2=1st excited state, etc.

    # Setup wavefunction for time-based propagation
    psi=psi+nondispersive_wavepacket(15,4.0)

    #        psi=nondispersive_wavepacket(20,8.0) # Centered on 20, with Width (speed?) 8.0
    #        psi=psi+nondispersive_wavepacket(40,-20.0) # Fight of the wavepackets!
    #        psi=planewave(8.0) # Plane wave, lambda=8.0 lattice units

    dt=1000 # Time step
    # As energy is in eV; hbar=1; we believe unit is ħ/q = ~0.658 fs

    println("Psi: ",psi)
    #myplot=lineplot(psi,name="Psi",color=:red,width=80,ylim=[-1,1])
    for i in 1:Unitarycycles
        @printf("\n\tUnitary Propagation Loop: %d\n",i)
        #        psi=eigvecs(H)[:,1] # gnd state
        S,H,psi,density,dipoles = UnitaryPropagation(dipoles,S,E,psi,dt,dampening,slices=1)
        plot_model(S,psi,density,dipoles,title="Dampening: $dampening TDSE: Cycle $i / $Unitarycycles")
        if PNG outputpng() end

        H=diagm(E,-1)+diagm(S)+diagm(E,1) #build full matrix from diagonal elements; for comparison
#        psi=eigvecs(H)[:,2] # 1=gnd state, 2=1st excited state, etc.

        ## Prototype surface hopping implementation
        overlaps=abs(overlap(eigvecs(H),psi)) # Overlap of current time-dep ψ with the ADIABATIC set of Φ for this H.
        closestAdiabaticState=indmax(overlaps) # index of maximum overlap vector; i.e. this is the adiabatic state with max overlap

        println("Wavefunction overlap; polaron c.f. complete set of (adiabatic) states:")
        for j in overlaps
            @printf("%.2f ",j)
        end
        println("\nSum overlaps: ",sum(overlaps)," Mean: ",sum(overlaps)/N)

        println("Maximum overlap with state $closestAdiabaticState =",overlaps[closestAdiabaticState])

        # Should be minimum switching algorithm / partition function sampling
        # For now we just randomly jump
        if i%50==0
            println(" BORED! JUMPING STATE!")
            psi=eigvecs(H)[:,closestAdiabaticState+1] # TODO: +1 is a lie; but otherwise it just jumps to the ground state again and again

            plot_model(S,psi,density,dipoles,title="JUMPING JACK FLASH TO: $closestAdiabaticState")
            if PNG outputpng() end
        end
    end
end

"""
    UnitarySim(dampening,Unitarycycles; slices=5, dt=1.0, PNG::Bool=false)

Simple Unitary (time evolution) simulation.

slices - number of slices in Trotter decomposition of Hamiltonian for matrix exponentiation
dt - Time step. As energy is in eV; hbar=1; we believe unit is ħ/q = ~0.658 fs
"""
function UnitarySim(dampening,Unitarycycles; slices=5, dt=1.0, PNG::Bool=false)
    S,E,H,psi,dipoles=prepare_model()

#    dipoles=[ -(x-N/2)^2*0.04 for x in 1:N ] # quadratic set of dipoles, to give energy slope across simulation
    dipoles=[ (x-N/2)^3*0.0005 for x in 1:N ]

    # Setup wavefunction for time-based propagation
    psi=nondispersive_wavepacket(10,16.0)

    #        psi=nondispersive_wavepacket(20,8.0) # Centered on 20, with Width (speed?) 8.0
    #        psi=psi+nondispersive_wavepacket(40,-20.0) # Fight of the wavepackets!
    #        psi=planewave(8.0) # Plane wave, lambda=8.0 lattice units
    println("Psi: ",psi)

    for i in 1:Unitarycycles
        @printf("\n\tUnitary Propagation Loop: %d\n",i)
        #        psi=eigvecs(H)[:,1] # gnd state
        S,psi,density,dipoles = UnitaryPropagation(dipoles,E,psi,dt,dampening,slices=slices)
        plot_model(S,psi,density,dipoles,title="Dampening: $dampening TDSE: Cycle $i / $Unitarycycles")
        if PNG outputpng() end
    end
end
