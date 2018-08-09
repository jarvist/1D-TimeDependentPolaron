# PotentialEnergySurfaces.jl
# Define quadratic diabatic potential energy surfaces
# Determine adiabatic potential energy surfaces from diagonalising the
#     diabatic hamiltonian potential.

"""
Function to define quadratic potential energy surface
------------
R_0 = equilibrium bond length
e_0 = energy at equilibrium bond length
------------
dpl = quadratic function representing site to the left of the bond centre
dpr = quadratic function representing site to the right of the bond centre
"""
function diabatic_potential(K::Float64,R_0::Float64,e_0::Float64)
    dpl = function (R) return K*(R+R_0/2)^2 + e_0 end
    dpr = function (R) return K*(R-R_0/2)^2 + e_0 end
    return dpl, dpr
end


"""
Function to diagonalise the diabatic potential
------------
e_11 = local energy of site 1
e_22 = local energy of site 2
J = coupling between sites 1 and 2
------------
H_d = [e_11 J
        J   e_22]
"""
function H_diabatic(e_11::Float64,e_22::Float64,J::Float64)
    H_d = [e_11 J;J e_22]; return H_d
end


"""
Function to diagonalise the diabatic hamiltonian and obtain the
unitary transformation matrix
------------
H_d = diabatic Matrix
------------
U_nk = Unitary transformation matrix
"""
function diagonalise(H_d::Array, test::Bool=false)
    U_nk = eigvecs(H_d)
    if test
        U = U_nk'U_nk
        return U, U_nk
    else
        return U_nk
    end
end


"""
Function to generate adiabatic potential energy surfaces
--------------
dpl = diabatic potential surface from site left of centre
dpr = diabatic potential surface from site right of centre
--------------
ape = excited adiabatic potential energy surface
apg = ground adiabatic potential energy surface
"""
function adiabatic_potential(dpl::Function,dpr::Function,J_if::Float64)
    ape = function (R) return 0.5*(dpl(R)+dpr(R)) + 0.5*sqrt((dpl(R)-dpr(R))^2+4*J_if^2) end
    apg = function (R) return 0.5*(dpl(R)+dpr(R)) - 0.5*sqrt((dpl(R)-dpr(R))^2+4*J_if^2) end
    return apg, ape
end


"""
Function to generate adiabatic hamiltonian matrix
---------------
H_d = diabatic hamiltonian
---------------
H_ad = adiabatic hamiltonian
"""
function H_adiabatic(H_d::Array)
    lambdas = eigvals(H_d)
    H_ad = diagm(lambdas)
end
