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
    f = function (R) return -derivative(r->pes(r/2),R) end
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
function acceleration(du::Array,u::Array,p,t,F::Function, M::Float64)
    x = u[1]
    dx = u[2]
    du[1] = dx
    du[2] = real(F(x)/M)
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
function classical_propagation(T::Float64, x_0::Float64, v_0::Float64, dt::Float64, F::Function, M::Float64)
    tspan = (0.0,T)
    myacceleration(du,u,p,t)=acceleration(du,u,p,t,F,M)

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
    Matrix = [-im*apg(R_e)-d_gg -d_ge ; -d_eg -im*ape(R_e)-d_ee]

    lambdas = eigvals(Matrix)
    one,two,three,four = eigvecs(Matrix)
    C = (a_ig*four - a_ie*three)/(one*four - two*three)
    D = -(a_ig*two - a_ie*one)/(one*four - two*three)

    a_fg = C*one*exp(lambdas[1]*dt) + D*three*exp(lambdas[2]*dt)
    a_fe = C*two*exp(lambdas[1]*dt) + D*four*exp(lambdas[2]*dt)
    return a_fg, a_fe
end
