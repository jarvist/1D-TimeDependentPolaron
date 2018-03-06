# models.jl
# Setup the 1D models 

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

" Shamelessly copied from Wikipedia:
https://en.wikipedia.org/wiki/Wave_packet "
function nondispersive_wavepacket(x0, λ)
    ψ=[ exp(-(x-x0)^2)*(cos(2π*(x-x0)/λ) - im * sin(2π*(x-x0)/λ) ) for x=1:N ]
    return ψ
end

function planewave(λ)
    k=2π/λ
    ψ=[ exp(-im* k*x) for x=1:N ]
    return ψ
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

    return S,E,H,psi,dipoles
end

