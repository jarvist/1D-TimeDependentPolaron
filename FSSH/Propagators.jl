# Propagators.jl
# Define functions to propagate the sites using the classical equation of motion
# and the wavefunctions using the time dependent Schrodinger equation

"""
Function to calculate the forces produced by the potential energy surface
-----------
pes = potential energy surface to be differentiated
-----------
f = force from potential energy surface as a function of site position
"""
function force(pes::Function)
    f = function (R) return -derivative(r->pes(r),R) end
    return f
end

"""
Define ordinary differential equation to be solved for classical equation of motion
-------------
du = [velocity, acceleration]
u = [position,velocity]
p = coefficient
t = time
F = force
-------------
ODE problem to be used in OrdinaryDiffEq Solver
"""
function acceleration(du::Array,u::Array,p,t,F::Function, M::Float64, K::Float64, R_0::Float64, ExtField::Float64)
    x = u[1]
    dx = u[2]
    du[1] = dx
    du[2] = 1.6e8/1.66*(-K*(x-R_0)+real(F(x))-ExtField)/M #1eV/(a.m.u*A) = 1.6e8/1.66 A/(0.1ns)^2
end

"""
Function to classically propagate nuclei
---------------
T = total time
x_0 = initial position
v_0 = initial velocity
dt = time step
F = function of force
---------------
x = array of positions
v = array of velocities
"""
function classical_propagation(T::Float64, x_0::Float64, R_0::Float64, v_0::Float64, dt::Float64, F::Function, M::Float64, K::Float64, ExtField::Float64)
    tspan = (0.0,T)
    myacceleration(du,u,p,t)=acceleration(du,u,p,t,F,M,K,R_0,ExtField)

    EoM = ODEProblem(myacceleration,[x_0,v_0],tspan)
    sol = solve(EoM, Euler(), dt = dt)
    x = sol[1,:]
    v = sol[2,:]
    return x,v
end



"""
Function to solve the time dependent Schrodinger equation to propagate the
adiabatic coefficients.
--------------------
a_ig = initial lower adiabatic state coeffient
a_ie = initial upper adiabatic state coeffient
d_ge = Non-adiabatic coupling vector from lower to upper adiabatic states
d_eg = Non-adiabatic coupling vector from upper to lower adiabatic states
d_gg = Non-adiabatic coupling vector from lower to lower adiabatic states
d_ee = Non-adiabatic coupling vector from upper to upper adiabatic states
R_e = Location of the electronic state on the nuclear coordinates
--------------------
a_fg = final lower adiabatic state coeffient
a_fe = final upper adiabatic state coeffient
"""
function new_a_mn(a_ig::Complex{Float64}, a_ie::Complex{Float64}, apg::Function, ape::Function, d_ge::Float64, d_eg::Float64, d_gg::Float64, d_ee::Float64, R_e::Float64, dt)
    Matrix = [-im*apg(R_e)/6.56e-26 -d_ge ; -d_eg -im*ape(R_e)/6.56e-26]

    lambdas = eigvals(Matrix)
    one = eigvecs(Matrix)[:,1]; two = eigvecs(Matrix)[:,2]
    C = (a_ig*two[2] - a_ie*two[1])/(one[1]*two[2] - one[2]*two[1])
    D = -(a_ig*one[2] - a_ie*one[1])/(one[1]*two[2] - one[2]*two[1])

    a_fg = C*one[1]*exp(lambdas[1]*dt) + D*two[1]*exp(lambdas[2]*dt)
    a_fe = C*one[2]*exp(lambdas[1]*dt) + D*two[2]*exp(lambdas[2]*dt)
    return a_fg, a_fe
end
