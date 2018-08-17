# PotentialEnergySurfaces.jl
# Define quadratic diabatic potential energy surfaces
# Determine adiabatic potential energy surfaces from diagonalising the
#     diabatic hamiltonian potential.

# """
# Function to define quadratic potential energy surface
# ------------
# A_0 = real constant
# A_1 = real coefficient of R
# A_2 = real coefficient of R^2
# ------------
# dpl = quadratic function representing site to the left of the bond centre
# dpr = quadratic function representing site to the right of the bond centre
# """
# function diabatic_potential(A_0::Float64, A_1::Float64, A_2::Float64)
#     dpl = function (R) return A_0 + A_1*R + A_2*R^2 end
#     dpr = function (R) return dpl(-R) end
#     return dpl, dpr
# end

"""
Function to define linear potential energy surface
Simple avoided crossing
------------
alpha = on-site EPC constant
tau = neigbouring sites EPC constant
overlap = <phi_l(R)|phi_r(R)>
------------
dpl = linear function representing site to the left of the bond centre
dpr = linear function representing site to the right of the bond centre
"""
function diabatic_potential(alpha::Float64, tau::Float64, overlap::Function)
    dpl = function (R) return -alpha*R/2*(1-abs(overlap(R))^2)-tau*(overlap(R)+conj(overlap(R))) end
    dpr = function (R) return dpl(-R) end
    dpx = function (R) return -tau*(overlap(R)^2+1) end
    return dpl, dpr, dpx
end
function diabatic_potential(alpha::Float64, tau::Float64, overlap::Float64)
    dpl = function (R) return -alpha*R/2*(1-abs(overlap)^2)-tau*(overlap+conj(overlap)) end
    dpr = function (R) return dpl(-R) end
    dpx = function (R) return -tau*(overlap^2+1) end
    return dpl, dpr, dpx
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
function adiabatic_potential(dpl::Function,dpr::Function,J_lr::Float64)
    ape = function (R) return 0.5*(dpl(R)+dpr(R)) + 0.5*sqrt((dpl(R)-dpr(R))^2+4*abs(J_lr)^2) end
    apg = function (R) return 0.5*(dpl(R)+dpr(R)) - 0.5*sqrt((dpl(R)-dpr(R))^2+4*abs(J_lr)^2) end
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
