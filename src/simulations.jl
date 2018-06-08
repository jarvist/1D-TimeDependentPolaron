# simulations.jl
# Convenience functions which run and plot the simulation 

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

        gui() # Show figure as pop up window...
    else
        # Directly using UnicodePlots
        myplot=lineplot(S,name="Site Energies",color=:red,width=80,ylim=[-1,1])
        lineplot!(myplot,real(psi),name="Psi",color=:blue)
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
        S,psi,density,dipoles = AdiabaticPropagation(dipoles,E,dampening)
        Plot_S_psi_density_dipoles(S,psi,density,dipoles,title="Dampening: $dampening SCF (Adiabatic): Cycle $i / $SCFcycles")
        if PNG outputpng() end
    end

    H=diagm(E,-1)+diagm(S)+diagm(E,1) #build full matrix from diagonal elements; for comparison
    psi=eigvecs(H)[:,2] # 1=gnd state, 2=1st excited state, etc.

    # Setup wavefunction for time-based propagation
    psi=psi+nondispersive_wavepacket(15,8.0)

    #        psi=nondispersive_wavepacket(20,8.0) # Centered on 20, with Width (speed?) 8.0
    #        psi=psi+nondispersive_wavepacket(40,-20.0) # Fight of the wavepackets!
    #        psi=planewave(8.0) # Plane wave, lambda=8.0 lattice units

    dt=1.0 # Time step
    # As energy is in eV; hbar=1; we believe unit is ħ/q = ~0.658 fs 

    println("Psi: ",psi)
    #myplot=lineplot(psi,name="Psi",color=:red,width=80,ylim=[-1,1])
    for i in 1:Unitarycycles
        @printf("\n\tUnitary Propagation Loop: %d\n",i)
        #        psi=eigvecs(H)[:,1] # gnd state
        S,psi,density,dipoles = UnitaryPropagation(dipoles,E,psi,dt,dampening,slices=1)
        Plot_S_psi_density_dipoles(S,psi,density,dipoles,title="Dampening: $dampening TDSE: Cycle $i / $Unitarycycles")
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

            Plot_S_psi_density_dipoles(S,psi,density,dipoles,title="JUMPING JACK FLASH TO: $closestAdiabaticState")
            if PNG outputpng() end
        end
    end
end

