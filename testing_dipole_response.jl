

println("\t\"He came riding fast, like a phoenix out of fire-flames.\" -- The Dancer, PJ Harvey" )


function randH(SiteEnergy, Edisorder, Jdisorder, modelJ, N)
# Random Trace / diagonal elements
    S=SiteEnergy + Edisorder*randn(N)
# Random Off-diag elements
    E=modelJ(0) + Jdisorder*randn(N-1)
    return (S,E)
end

function prepare_model2()
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
            Fsum+=density[j]/(i-j)^3 # How much do the dipoles respond to the electron density?
        end
        Field[i]=Fsum/r^3
    end
    return Field
end

function FieldFromDipole(dipole)
    FieldDI=zeros(N)
    for i in 1:N
        FDIsum=zero(eltype(dipole)) # Calculate total field at each
        for j in 1:N
            if (j==i)
                continue # avoid infinite self energies
            end
            FDIsum+=dipole[j]/(i-j)^3 # How much do the dipoles respond to the electron density?
        end
        FieldDI[i]=FDIsum/r^3
    end
    return FieldDI
end

" Decompose Hamiltonian into Diagonal/S/PE and Off-diag/J/KE elements"
function Decompose_H(H)
    S=eye(N).*H # elementwise to select for just diagonal terms
    J=H-S
    return S,J
end

function Propagate(H,psi,dt)
    #U=eye(H) # identiy matrix same size + type as H

    #S,J=Decompose_H(H) # split into diagonal and off-diag terms
    #Trotter decomposition
    #U*=expm(-im*J*dt/2*hbar)*expm(-im*S*dt/hbar)*expm(-im*J*dt/2*hbar)
    #psi = U*psi
    psi=eigvecs(H)[:,1]
    En = eigvals(H)
    normalisation = norm(psi)
    psi = psi/normalisation

    density = conj(psi).*psi
    return psi, density
end

function UpdateDipole(Field,density,dipole)
    alpha = 1/221*ones(N)#(N*density)
    M = diagm(alpha)
    n=0
    n2=1
    for i in 1:N
        for j in 1:N
            if (j==i)
                #M[n+j] = alpha[i]
                continue
            end
            M[n+j] = 1/(r*(i-j))^3
        end
        n+=N
    end
    new_dipole = \(M,Field)
    dipole = 0.2*(new_dipole-dipole)
    return new_dipole

end

function UpdateEnergy(dipole,density,FieldDI, Field,S,E)
    Vpq = (-dipole.*(Field)+ density.*(FieldDI))
    norm_vpq = norm(Vpq)
    #plot(real(Vpq), label = "Vpq")
    Vpp = dipole.*(FieldDI)
    norm_vpp = norm(Vpp)
    #plot!(real(Vpp), label = "Vpp")
    Vqq = -density.*(Field)
    norm_vqq = norm(Vqq)
    #plot!(real(Vqq), label = "Vqq")
    #KE = S - V
    #E=?
    S=Vpq+Vpp+Vqq
    norm_S = norm(S)
    #plot!(real(S), label = "Total")
    println("norms: vpq", norm_vpq, "vpp", norm_vpp, "vqq", norm_vqq, "total", norm_S)

    H=diagm(E,-1)+diagm(S)+diagm(E,1) #build full matrix from diagonal elements; for comparison
    return S,E,H
end

function UpdateEnergy2(dipole,Field,S,E)
    S=zeros(N)
    for i in 1:N
        for j in 1:N
            if (j==i)
                continue # avoid infinity in self energies
            end
            S[i]+=dipole[j]/(r*(i-j))^3 # Contribution to site energy (1 e- at site) from dipoles
        end
#        @printf("Site: i %d SiteEnergy: S[i] %f\n",i,S[i])
    end
    H=diagm(E,-1)+diagm(S)+diagm(E,1)
    return S,E,H
end


function nondispersive_wavepacket(x0, λ)
    ψ=[ exp(-(x-x0)^2)*(cos(2π*(x-x0)/λ) - im * sin(2π*(x-x0)/λ) ) for x=1:N ]
    return ψ
end



"""
Think about which way the energy, dipoles and fields should be....
Think about magnitudes of dipoles and densities vs Fields (hence which have the most weight for energy)
"""
