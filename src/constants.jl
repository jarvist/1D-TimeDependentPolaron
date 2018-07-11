# constants.jl  -  Physical constants

# 2016-07 - Updated accuracy from https://github.com/DataWookie/PhysicalConstants.jl
const h = 6.62606896e-34;                          # Js = kg m2 / s 
const hbar = const ħ = 1.05457162825e-34;          # Js = kg m2 / s 

const eV = const q = const ElectronVolt = 1.602176487e-19;    # kg m2 / s2 
const MassElectron = 9.10938188e-31;                          # kg

const ε0 = const VacuumPermittivity = 8.854187817e-12;        # A2 s4 / kg m3 
const VacuumPermeability = 1.25663706144e-6;                  # kg m / a2 s2 
const Debye = 3.33564095198e-30;                              # A s2 / m2 
const Gauss = 1e-4;                                           # kg / A s2 

const Rydberg = 2.17987196968e-18;                            # J = kg m2 / s2 
const Hartree = 2*Rydberg;                                    # J = kg m2 / s2
const Boltzmann = const kB =  1.3806504e-23;                  # J/K = kg m2 / K s2 

# convenience for conversions
const hbar_in_eV = const ħeV = hbar/eV #eV energy units
const kB_in_eV = kB/eV # in units of eV
const Hartree_in_eV = Hartree/eV # eV
const Rydberg_in_eV = Rydberg/eV 

# Dodgy unit conversions, done dirt cheap.
const Å=1E-10 # Angstrom

