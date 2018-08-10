# LetsGetTogether.jl
# Uses SurfaceHopping module to combine all functions and start the simulation


using SurfaceHopping

# Initialise system
PESg = true
cl = 1.0; cr = 0.0
R_0 = 1.442307; e_l = e_r = 0.0; K = 1.32041142; global M = 1.
R_n = 1.6; xr = R_n/2; xl = -xr; v_0 = 0.0
T = 500; dt = 0.01; ddt = 0.001
n = Int(dt/ddt +1); N = Int(T*n)

# Set up empty arrays
KE = zeros(N); PE = zeros(N); TE = zeros(N)

# Generate diabatic potential energy surfaces
dpl,dpr = diabatic_potential(K, R_0, e_l)
# Generate diabatic electonic states and form of their derivative
phi_l, phi_r = diabatic_states(R_n)
dphi = diabatic_derivative()
# Calculate coupling term
J_lr = overlap_phi(phi_l,phi_r)
#Generate adiabatic potential energy surfaces
apg,ape = adiabatic_potential(dpl,dpr,J_lr)

# plot(dpl)
# plot!(dpr)
# plot!(phi_l)
# plot!(phi_r)
# plot!(apg)
# plot!(ape)

# Determine which adiabatic potential energy surface to propagate nuclei on
if PESg
    f = force(apg)
else
    f = force(ape)
end
# Deteremine which site the electronic state is localised on
if cl^2 > cr^2
    R_e = xl
else
    R_e = xr
end
# Perform classical propagation of the nuclei from t -> t+dt
x,v = classical_propagation(dt, R_n, v_0, dt, f, M)
dx = x[end] - x[1]; dx_l = -dx/2; dx_r = -dx_l

# plot(x)
# plot!(v)

# Generate necessary parameterd for quantum propagation
H_d = H_diabatic(dpl(xl),dpr(xl),J_lr)
U_nk = diagonalise(H_d)
a_g, a_e = a_mn(U_nk, cl, cr); a_g = Complex(a_g); a_e = Complex(a_e)
d_ge, d_eg, d_gg, d_ee = NACV(phi_l, phi_r, dphi, U_nk, R_n, dx_l, dx_r, dt)
# Perform quantum propagation of the electronic state from t -> t+dt
for j in 1:n
    a_gn, a_en = new_a_mn(a_g, a_e, apg, ape, d_ge, d_eg, d_gg, d_ee, R_e, ddt)
    a_g = a_gn; a_e = a_en;
end

# Calculate the probability of hopping
if PESg
    g = g_mn(a_e, a_g, d_eg, dt) #g_eg
else
    g = g_mn(a_g, a_e, d_ge, dt) #g_ge
end
println(g)
