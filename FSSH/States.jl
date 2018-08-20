# States.jl
# Define electronic and nuclear states

"""
Function to define diabatic states as gaussian functions of width σ
--------------
R_n = bond length
σ = width of states
--------------
phi_l = gaussian centered on left site
phi_r = gaussian centered on right site
"""
function diabatic_state(R_n::Float64,σ::Float64=0.998)
    gaussian = function (R) return exp(-R^2/(2*σ^2)) end
    gaussian2(R) = gaussian(R)^2
    A = riemann(gaussian2, -10,10,1000)
    A = 1/sqrt(A)
    phi(R) = A*gaussian(R-R_n)
    return phi
end


"""
Overlap of diabatic states to determine coupling strength
--------------
phi_l = gaussian centered on left site
phi_r = gaussian centered on right site
--------------
J_lr = coupling strength to be used in diabatic Hamiltonian
"""
function overlap_phi(phi_l::Function, phi_r::Function)
    overlap_phi(R) = conj(phi_l(R))*phi_r(R)
    J_lr = riemann(overlap_phi, -10,10,1000)
    return J_lr
end



"""
Function for the derivative of a gaussian function
    (to be used in NACV when the derivative of phi w.r.t R is needed)
--------------
R_n = bond length
σ = width of states
--------------
d_phi = gaussian derivative
"""
function diabatic_derivative(σ::Float64=0.374)
    dphi = function (R) return (R)*exp(-(R)^2/(2*σ^2))/σ^2 end
end



"""
Function to define the adiabatic states from the diabatic states and
transformation Matrix
------------
U_nk = Unitary transformation matrix
phi_l = diabatic state on left site
phi_r = diabatic state on right site
------------
psi_g = adiabatic lower state
phi_e = adiabatic upper state
"""
function adiabatic_states(phi_l::Function, phi_r::Function, U_nk::Array)
    psi_g = function (R) return U_nk[1]*phi_l(R) + U_nk[2]*phi_r(R) end
    psi_e = function (R) return U_nk[3]*phi_l(R) + U_nk[4]*phi_r(R) end
    return psi_g, psi_e
end


"""
Function to define the non-adiabatic coupling vectors for all possible transisions
--------------
phi_l = diabatic state on left site
phi_r = diabatic state on right site
dphi = generic derivative of gaussian function
U_nk = Unitary transformation matrix
R_n = bond length
dR_l = change of nuclear coordinate on the left site
dR_r = change of nuclear coordinate on the right site
dt = time step this change has occured in
--------------
d_ge = Non-adiabatic coupling vector from lower to upper adiabatic states
d_eg = Non-adiabatic coupling vector from upper to lower adiabatic states
d_gg = Non-adiabatic coupling vector from lower to lower adiabatic states
d_ee = Non-adiabatic coupling vector from upper to upper adiabatic states
"""
function NACV(phi_l::Function, phi_r::Function, dphi::Function, U_nk::Array, R_n::Float64, dR_l, dR_r, dt)
    dphi_l = function (R) return dphi((R)) end
    dphi_r = function (R) return dphi((R-R_n)) end

    # perform integrations for each pair of diabatic states
    int_ll = function (R) return conj(phi_l(R))*dphi_l(R) end
    ll = riemann(int_ll,-10,10,1000)
    int_lr = function (R) return conj(phi_l(R))*dphi_r(R) end
    lr = riemann(int_lr,-10,10,1000)
    int_rl = function (R) return conj(phi_r(R))*dphi_l(R) end
    rl = riemann(int_rl,-10,10,1000)
    int_rr = function (R) return conj(phi_r(R))*dphi_r(R) end
    rr = riemann(int_rr,-10,10,1000)


    l_ge = dR_l*(conj(U_nk[1]*ll + U_nk[2]*rl)*U_nk[3])/dt
    r_ge = dR_r*(conj(U_nk[1]*lr + U_nk[2]*rr)*U_nk[4])/dt
    d_ge = l_ge+r_ge


    l_eg = dR_l*(conj(U_nk[3]*ll + U_nk[4]*rl)*U_nk[1])/dt
    r_eg = dR_r*(conj(U_nk[3]*lr + U_nk[4]*rr)*U_nk[2])/dt
    d_eg = l_eg+r_eg

    l_gg = dR_l*(conj(U_nk[1]*ll + U_nk[2]*rl)*U_nk[1])/dt
    r_gg = dR_r*(conj(U_nk[1]*lr + U_nk[2]*rr)*U_nk[2])/dt
    d_gg = l_gg+r_gg

    l_ee = dR_l*(conj(U_nk[3]*ll + U_nk[4]*rl)*U_nk[3])/dt
    r_ee = dR_r*(conj(U_nk[3]*lr + U_nk[4]*rr)*U_nk[4])/dt
    d_ee = l_ee+r_ee

    return d_ge,d_eg,d_gg,d_ee, [l_ge*dt/dR_l, r_ge*dt/dR_r, l_eg*dt/dR_l, r_eg*dt/dR_r]
end


"""
Function to determine adiabatic coefficients from diabatic coefficients
    and unitary transformation matrix
-------------
U_nk = Unitary Transformation Matrix
c1 = diabatic coefficient of left site
c2 = diabatic coefficient of right site
-------------
ag = adiabatic coefficient of ground state
ae = adiabatic coefficient of excited state
"""
function a_mn(U_nk::Array,cl::Complex{Float64},cr::Complex{Float64})
    ag = (cl*U_nk[4] - cr*U_nk[3])/(U_nk[1]*U_nk[4]-U_nk[2]*U_nk[3])
    ae = (cl*U_nk[2] - cr*U_nk[1])/(U_nk[2]*U_nk[3]-U_nk[1]*U_nk[4])
    return ag,ae
end


"""
Function to determine the probability of switching states
------------
a_m = adiabatic coefficient of state m
a_n = adiabatic coefficient of state n
d_mn = Non-adiabatic coupling vector for transition from state n to state m
dt = time step
------------
g_mn = probability of a transition from adiabatic surface m to adiabatic surface n
"""
function g_mn(a_m::Complex{Float64}, a_n::Complex{Float64}, d_mn::Float64, dt::Float64)
    g_mn =  -2*real(dt*(conj(a_n)*a_m*d_mn)/(conj(a_m)*a_m))
    return g_mn
end
