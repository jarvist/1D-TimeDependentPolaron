push!(LOAD_PATH,"../")

include("testing_dipole_response.jl")
using Plots

"""
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

#'''SetUp'''
const kB=1 #8.6173324E-5 #eV
const hbar=1.0

const N=50 # Number of sites in model
             # --> Size of tridiagonal Hamiltonian constructed; larger value -> better statistics, slower runtime
const Edisorder=0.0 # Energetic disorder eV, Gaussian form
const Jdisorder=0.0 # Transfer integral disorder, eV.
const r = 20

# Model setup
const J0=0.0005 #(~0.1eV)
modelJ(θ) = J0*cos(θ*π/180.0).^2

const T=0.0032 #(~300K)
const B=1/(T*kB) #300K * k_B in eV

# This effectively reduces down to the 'alpha' parameter in the Frohlich polaron Hamiltonian
const dipolestrength=0.2

function TestNewDipole()
    #'''prepare'''
    S,E,H,psi,density = prepare_model2()
    dipole = zeros(N)
    Field = FieldFromDensity(density)
    FieldDI = FieldFromDipole(dipole)
    Field = Field+FieldDI
    S,E,H = UpdateEnergy2(dipole,Field,S,E)
    dipole = UpdateDipole(Field,density,dipole)


    for i in 1:100
        S,E,H = UpdateEnergy2(dipole,Field,S,E)
        psi = eigvecs(H)[:,1]
        density = psi.^2
        Field = FieldFromDensity(density)
        Field = Field+FieldDI
        dipole2 = UpdateDipole(Field,density,dipole)
        FieldDI = FieldFromDipole(dipole)
        dipole = dipole2
        norm_dipole = norm(dipole)
        norm_Field = norm(Field)
        norm_FieldDI = norm(FieldDI)
        norm_S = norm(S)
        norm_density=norm(density)
        xlims!(1,N)
        ylims!(-1,1)
        plot(real(dipole/norm_dipole), label = "dipole")
        plot!(Field, label="Field")
        plot!(real(density)/norm_density, label = "density")
        plot!(real(S/norm_S), label = "energy")
        gui()
    end


    psi=psi+nondispersive_wavepacket(20,4.0)
    normalise = norm(psi)
    psi = psi/normalise
    density = conj(psi).*psi
    #plot(real(density))
    #gui()
    Field = FieldFromDensity(density)
    Field = Field+FieldDI
    dipole = UpdateDipole(Field,density,dipole)
    FieldDI = FieldFromDipole(dipole)
    S,E,H = UpdateEnergy2(dipole,Field,S,E)


    #'''propagate'''
    for i in 1:1000
        norm_dipole = norm(dipole)
        plot(real(dipole/norm_dipole), label = "dipole")
        norm_Field = norm(Field)
        norm_FieldDI = norm(FieldDI)
        norm_S = norm(S)
        #plot!(imag(psi))
        #plot!(real(psi))
        xlims!(1,N)
        #ylims!(-1,1)
        norm_density=norm(density)
        plot!(real(density)/norm_density, label = "density")
        #plot!(FieldDI/norm_FieldDI, label = "FieldDI")
        #plot!(Field/norm_Field, label = "Field")
        #plot!((Field+FieldDI)/(norm_Field+norm_FieldDI), label = "Fields")
        plot!(real(S/norm_S), label = "energy")
        gui()
        psi, density = Propagate(H,psi,100)
        Field = FieldFromDensity(real(density))
        Field = Field+FieldDI
        dipole = UpdateDipole(Field, density,dipole)
        FieldDI = FieldFromDipole(dipole)
        S,E,H = UpdateEnergy2(dipole,Field,S,E)

    end
end

TestNewDipole()
