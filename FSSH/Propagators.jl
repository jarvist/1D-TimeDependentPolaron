# Propagators.jl
# Define functions to propagate the sites using the classical equation of motion
# and the wavefunctions using the time dependent Schrodinger equation

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
function new_a_mn(a_ig::Float64, a_ie::Float64, apg::Function, ape::Function, d_ge::Float64, d_eg::Float64, d_gg::Float64, d_ee::Float64, R_e::Float64, dt)
    Matrix = [-im*apg(R_e)-d_gg -d_ge ; -d_eg -im*ape(R_e)-d_ee]
    
    lambdas = eigvals(Matrix)
    one,two,three,four = eigvecs(Matrix)
    C = (a_ig*four - a_ie*three)/(one*four - two*three)
    D = -(a_ig*two - a_ie*one)/(one*four - two*three)

    a_fg = C*one*exp(lambdas[1]*dt) + D*three*exp(lambdas[2]*dt)
    a_fe = C*two*exp(lambdas[1]*dt) + D*four*exp(lambdas[2]*dt)
    return a_fg, a_fe
end
