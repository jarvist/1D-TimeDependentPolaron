using Plots
gr()

module SurfaceHopping

export diabatic_potential,H_diabatic,diagonalise,adiabatic_potential,H_adiabatic
export diabatic_state, overlap_phi, diabatic_derivative, adiabatic_states, NACV, a_mn, g_mn
export force, acceleration, classical_propagation, new_a_mn
export riemann


using Distributions
using Calculus
using OrdinaryDiffEq



"""
Numeric integration, taken from http://mth229.github.io/integration.html
"""
function riemann(f::Function, a::Real, b::Real, n::Int; method="right")
    meth(f,l,r) = (1/6) * (f(l) + 4*(f((l+r)/2)) + f(r)) * (r-l)


    xs = a + (0:n) * (b-a)/n
    first = xs[1:end-1]; last = xs[2:end]; join = [first ; last];
    len = length(join)/2
    elements = reshape(join,(Int(len),2))
    pair = [elements[i,:] for i in 1:Int(len)]
    as = [meth(f,l,r) for (l,r) in pair]#meth(f, l, r)
    sum(as)

end

include("Propagators.jl")
include("States.jl")
include("PotentialEnergySurfaces.jl")


end
