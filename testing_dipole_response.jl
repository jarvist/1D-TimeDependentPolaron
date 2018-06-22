

println("\t\"He came riding fast, like a phoenix out of fire-flames.\" -- The Dancer, PJ Harvey" )


function randH(SiteEnergy, Edisorder, Jdisorder, modelJ, N)
# Random Trace / diagonal elements
    S=SiteEnergy + Edisorder*randn(N)
# Random Off-diag elements
    E=modelJ(0) + Jdisorder*randn(N-1)
    return (S,E)
end

function prepare_model()
    S,E=randH(5.0,Edisorder, Jdisorder, modelJ, N)
    H=diagm(E,-1)+diagm(S)+diagm(E,1) #build full matrix from diagonal elements; for comparison
    psi=eigvecs(H)[:,1] # 1=gnd state

    #    Plot_H(H)
    #    TODO: UPDATE THIS FOR PLOTS.jl

    ## Construct initial dipoles structure
    density = conj(psi).*psi

    return S,E,H,psi,density
end

function FieldFromDensity(density)
    Field=zeros(N)
    for i in 1:N
        Fsum=zero(eltype(density)) # Calculate total field at each
        for j in 1:N
            if (j==i)
                continue # avoid infinite self energies
            end
            Fsum+=density[j]/(j-i)^3 # How much do the dipoles respond to the electron density?
        end
        Field[i]=Fsum/r^3
    end
    return Field
end

" Decompose Hamiltonian into Diagonal/S/PE and Off-diag/J/KE elements"
function Decompose_H(H)
    S=eye(N).*H # elementwise to select for just diagonal terms
    J=H-S
    return S,J
end

function Propagate(H,psi,dt)
    U=eye(H) # identiy matrix same size + type as H

    S,J=Decompose_H(H) # split into diagonal and off-diag terms
    # Trotter decomposition
    U*=expm(-im*J*dt/2*hbar)*expm(-im*S*dt/hbar)*expm(-im*J*dt/2*hbar)
    psi = U*psi

    normalisation = norm(psi)
    psi = psi/normalisation

    density = conj(psi).*psi
    return psi, density
end

function UpdateDipole(Field,density)
    alpha = 10*(N*density)
    M = diagm(alpha)
    n=0
    n2=1
    for i in 1:N
        for j in 1:N
            if (j==i)
                continue
            end
            M[n+j] = 1/(r*(i-j))^3
        end
        n+=N
    end
    dipole = \(M,Field)
    return dipole

end

function UpdateEnergy(dipole,Field,S,E)
    V = dipole.*Field
    KE = S - V
    #E=?
    S=V

    H=diagm(E,-1)+diagm(S)+diagm(E,1) #build full matrix from diagonal elements; for comparison
    return S,E,H
end


function nondispersive_wavepacket(x0, λ)
    ψ=[ exp(-(x-x0)^2)*(cos(2π*(x-x0)/λ) - im * sin(2π*(x-x0)/λ) ) for x=1:N ]
    return ψ
end
