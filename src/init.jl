"""
Rose Teague & Jarvist Frost 04/07/2018
Script to initialise the system.
Define functions to prepare global constants and plotting styles

- init
- Plotting
--------------------------------------------------------------------------------

CONSTANTS

Atomic units:
-------------------------
kb = 1
hbar = 1
4*pi*e0 = 1
a_0 = 1
e = 1
-------------------------
distance = bohr (a_0)
energy = hartree (27.211385 eV)
time = hbar/E_h (2.418884326505(16)×10−17 s)
dipole moment = ea_0 (8.47835326(19)×10−30 C·m = 2.541746 Debye)
polarisability = e2 a_02/E_h (1.65E-41 C^2 m^2 /J)
temperature = E_h/kB (3.1577464E5 K)
"""

#-------------------------------------------------------------------------------
function init(num,Edis, Jdis, J, Temp, radius)
    #=
    Initialise global parameters to define the system
    ---------------------
    Input
    ---------------------
    N - (int)  Number of sites in model.
    Edisorder - (eV)  Energetic disorder, Gaussian form.
    Jdisorder - (eV)  Transfer integral disorder.
    J0 - (eV) Transfer Integral.
    T - (K) Temperature
    radius - (Bohr)
    ---------------------
    output
    ---------------------
    N - (int)
    Edisorder - (Hartree)
    Jdisorder - (Hartree)
    J0 - (Hartree)
    T - (E_h/kB)
    B - (1/E_h)
    r - (Bohr)
    =#

    global N=num
    global Edisorder=Edis
    global Jdisorder=Jdis
    global J0=J/27.211385
    global modelJ = J0*cos(π/180.0).^2
    global T=Temp/3.1577464E5
    global B=1/T
    global r=radius
    # This effectively reduces down to the 'alpha' parameter in the Frohlich polaron Hamiltonian
    # Contains the exponential decay of the dipole response (exp(-TimeStep*DecayRate))
    global dipolestrength=0.2

    return N, Edisorder, Jdisorder, modelJ, B, dipolestrength, r
end

#-------------------------------------------------------------------------------
#=
--------------------------------------------------------------------------------
Liberated from 'Sturm': https://github.com/jarvist/Teclo/blob/master/Sturm.jl

Generate a random tridiagonal TightBinding Hamiltonian, in a form suitable for
the Sturm sequence

# Given:
 *   SiteEnergy - scalar eV; reference for site energy of monomer
 *   Edisorder - scalar eV ; amount of Gaussian / normal energetic disorder,
        for trace of Hamiltonian
 *   Jdisorder - scalar eV
 *   modelJ(theta) - function, takes degrees, returns eV ; model for the
        transfer integral (e.g. E=(J0*cos(thetas*pi/180)).^2 )
 *   N - integar ; size of diagonal of Hamiltonian
 =#
function randH(SiteEnergy, Edisorder, Jdisorder, modelJ, N)
    #=
    Determine on-site and nearest neighbour energies for TightBinding
    Hamiltonian
    ----------------
    Inputs
    ----------------
    SiteEnergy - self-energy of sites (Hartree)
    Edisorder  - scale of disorder of on-site energies (Hartree)
    Jdisorder  - scale of disorder of transfer integral (Hartree)
    modelJ     - Transfer Integral for nearest neighbours (Hartree)
    N          - Number of sites in the system
    ----------------
    Outputs
    ----------------
    S - Diagonal elements of TightBinding Hamiltonian (Hartree)
    E - Off-diagonal elements of TightBinding Hamiltonian (Hartree)
    =#

    S=SiteEnergy + Edisorder*randn(N)
    E=modelJ + Jdisorder*randn(N-1)
    return (S,E)
end

#-------------------------------------------------------------------------------

function prepare_model()
    #=
    Convenience function to put together everything for initiating a 1D model.
    ------------
    Outputs
    ------------
    S       - Diagonal elements of TightBinding Hamiltonian (Hartree)
    E       - Off-Diagonal elements of TightBinding Hamiltonian (Hartree)
    H       - TightBinding Hamiltonian (Hartree)
    psi     - wavefunction
    dipoles - initial dipole matrix of zeros (e*a_0)
    =#

    S,E = randH(0.05, Edisorder, Jdisorder, modelJ, N)
    H = diagm(E,-1)+diagm(S)+diagm(E,1) #build full matrix from diagonal elements
    psi = eigvecs(H)[:,1] # 1=gnd state

    #    Plot_H(H)
    #    TODO: UPDATE THIS FOR PLOTS.jl

    ## Construct initial dipoles structure
    dipoles=zeros(N)

    return S,E,H,psi,dipoles
end

#-------------------------------------------------------------------------------
#Wrapper function to pretty-print and plot (UnicodePlots) relevant items of interest.
function Plotting(S,psi,density,dipoles;title="",verbose::Bool=false)
        #=
        Plot functions to display energy, wavefunction, electron density and dipole
        -----------------
        Inputs
        -----------------
        S - On-site energies (Hartrees)
        psi - Wavefunction
        density - Electron Density
        dipoles - Site Dipoles
        -----------------
        Output
        -----------------
        Plots of each cycle
        =#
    if verbose
        println("Site energies: ",S)
        println("Psi: ",psi)
        println("Electron density: ",density," Sum: ",sum(density))
        println("Dipoles: ",dipoles," Sum: ",sum(dipoles))
    end

    # Using Plots interface

    plot(S/norm(S),label="Site Energies",color=:red, width=2, yaxis = (" ", (-1,1)), size=(1024,768)) # makes it larger, but also slower!
    plot!(density,label="Electon Density",color=:orange, width=4, fill=true, fill=(0,0.3,:orange))
    plot!(dipoles/norm(dipoles),label="Dipoles",color=:blue, width=2)


    plt = twinx()
    psioffset=-0.0 # Shift wavefunctions down
    plot!(plt,real(psi)+psioffset,label="Re[Psi]", legend = :bottomright, color=:green , width=2, fill=true, fill=(0,0.3,:green),yaxis=(" ",(-0.5,1.5)))
    plot!(plt,imag(psi)+psioffset,label="Im[Psi]",color=:pink, width=2, fill=true, fill=(0,0.3,:pink))


    xaxis!("Site")
    title!(title)

    gui()

end
