# load module from local directory

using SurfaceHopping
using Base.Test

# -----------------------------------------------------------------------------
H_d = H_diabatic(2.6,2.5,0.1)
U, U_nk = diagonalise(H_d,true)

ε = 1e-2
@test 1.0 ≈ U[1] atol = ε

# -----------------------------------------------------------------------------
H_ad = H_adiabatic(H_d)
H = *(U_nk',*(H_d,U_nk))

@test H[1] ≈ H_ad[1] atol = ε

# -----------------------------------------------------------------------------
phi_l,phi_r = diabatic_states(0.7,0.1)
phi_l2(R) = phi_l(R)^2; phi_r2(R) = phi_r(R)^2;

ε = 1e-4
integral = riemann(phi_l2,-10,10,1000)
@test 1.0 ≈ integral atol = ε

integral = riemann(phi_r2,-10,10,1000)
@test 1.0 ≈ integral atol = ε

# -----------------------------------------------------------------------------
psi_g, psi_e = adiabatic_states(phi_l, phi_r, U_nk)
psi_g2(R) = psi_g(R)^2; psi_e2(R) = psi_e(R)^2

integral = riemann(psi_g2,-10,10,1000)
@test 1.0 ≈ integral atol = ε

integral = riemann(psi_e2,-10,10,1000)
@test 1.0 ≈ integral atol = ε

# -----------------------------------------------------------------------------
dphi = diabatic_derivative()

d_ge, d_eg, d_gg, d_ee = NACV(phi_l, phi_r, dphi, U_nk, 2.6, -0.1, 0.1, 0.1)
#println(d_ge, "  ", d_eg, "  ", d_gg, "  ", d_ee)

# -----------------------------------------------------------------------------
cl = 1.0; cr = 0.0
ag, ae = a_mn(U_nk, cl, cr)
coeffs_diabatic = abs(cl)^2 + abs(cr)^2
@test 1.0 ≈ coeffs_diabatic atol = ε

coeffs_adiabatic = abs(ag)^2 + abs(ae)^2
@test 1.0 ≈ coeffs_adiabatic atol = ε

# -----------------------------------------------------------------------------
dpl, dpr = diabatic_potential(1.0,2.6,0.0)
apg, ape = adiabatic_potential(dpl, dpr, 2.0)
agn, aen = new_a_mn(ag, ae, apg, ape, d_ge, d_eg, d_gg, d_ee, 2.6, 0.1)

coeffs_adiabatic = abs(agn)^2 + abs(aen)^2
@test 1.0 ≈ coeffs_adiabatic atol = ε

# -----------------------------------------------------------------------------
println("Tests finished succesfully.")
