# propagators.jl
# Let's get those ψ moving. 

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
    return dipoles
end

"Build Hamiltonian from Dipoles (via SiteEnergyFromDipoles), diagonalise and update dipoles responding to ground state electron density. A wrapper function."
function AdiabaticPropagation(dipoles,E,dampening)
    S=SiteEnergyFromDipoles(dipoles)

    H=diagm(E,-1)+diagm(S)+diagm(E,1) #build full Hamiltonian
    psi=eigvecs(H)[:,1] # gnd state
   
    println("Adiabtic State energy: <psi|H|psi> = ",psi'*H*psi)
    KE=diagm(E,-1)+diagm(E,+1)
    PE=diagm(S)
    println("KE: <psi|KE|psi> = ",psi'*KE*psi)
    println("PE: <psi|PE|psi> = ",psi'*PE*psi)
  
    density=psi.^2
    dipoles=DipolesFromDensity(dipoles,density,dampening)

    return S,psi,density,dipoles
end

"Self-consistent response of the lattice with unitary (time dependent) evolution of the wavefunction. "
function UnitaryPropagation(dipoles,E,psi,dt,dampening;slices::Int=1)
    S=SiteEnergyFromDipoles(dipoles)

    H=diagm(E,-1)+diagm(S)+diagm(E,+1) 
    
    En=eigvals(H)[1]
    println("First Eigval (adibatic): ",En)

#    psi=eigvecs(H)[:,1] # gnd state; reproduces adiabtic state (eigval) above
    println("State energy: <psi|H|psi> = ",psi'*H*psi)
    KE=diagm(E,-1)+diagm(E,+1)
    PE=diagm(S)
    println("KE: <psi|KE|psi> = ",psi'*KE*psi)
    println("PE: <psi|PE|psi> = ",psi'*PE*psi)
    
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
