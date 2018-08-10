# LetsGetTogether.jl
# Uses SurfaceHopping module to combine all functions and start the simulation


using SurfaceHopping

# Initialise system
PESg = true; plotting = true
cl = Complex(1.0); cr = Complex(0.0)
R_0 = 1.442307; e_l = e_r = 0.0; K = 1.32041142; global M = 1.
R_n = 2.6; xr = R_n/2; xl = -xr; v_0 = 0.0; R_n0 = 2.6;
T = 5000.; dt = 0.001; ddt = 0.0001
n = Int(dt/ddt +1); N = Int(T)

# Set up empty arrays
KE = zeros(N); PE = zeros(N); TE = zeros(N); Ee = zeros(N)
V_change = zeros(N)

# Generate diabatic potential energy surfaces
dpl,dpr = diabatic_potential(K, R_0, e_l)
# Generate diabatic electonic states and form of their derivative
phi_l, phi_r = diabatic_states(R_n)
dphi = diabatic_derivative()
# Calculate coupling term
J_lr = overlap_phi(phi_l,phi_r)
#Generate adiabatic potential energy surfaces
apg,ape = adiabatic_potential(dpl,dpr,J_lr)


for i in 1:Int(T)

    phi_l_new(R) = phi_l(R-(R_n0-R_n)/2);phi_r_new(R) = phi_r(R+(R_n0-R_n)/2)
    J_lr_new = overlap_phi(phi_l_new,phi_r_new)
    apg,ape = adiabatic_potential(dpl,dpr,J_lr_new)

    # Determine which adiabatic potential energy surface to propagate nuclei on
    if PESg
        f = force(apg)
    else
        f = force(ape)
    end
    # Deteremine which site the electronic state is localised on
    if abs(cl)^2 > abs(cr)^2
        R_e = xl
    else
        R_e = xr
    end
    # Perform classical propagation of the nuclei from t -> t+dt

    x,v = classical_propagation(2*dt, R_n, v_0, dt, f, M)
    dx = x[end] - x[1]; dx_l = -dx/2; dx_r = -dx_l

    # Generate necessary parameterd for quantum propagation
    H_d = H_diabatic(dpl(xl),dpr(xl),J_lr_new)
    H_ad = H_adiabatic(H_d)
    U_nk = diagonalise(H_d)
    a_g, a_e = a_mn(U_nk, cl, cr); a_g = Complex(a_g); a_e = Complex(a_e)
    d_ge, d_eg, d_gg, d_ee = NACV(phi_l_new, phi_r_new, dphi, U_nk, R_n, dx_l, dx_r, 2*dt)
    # Perform quantum propagation of the electronic state from t -> t+dt
    for j in 1:2*n
        a_gn, a_en = new_a_mn(a_g, a_e, apg, ape, d_ge, d_eg, d_gg, d_ee, R_e, ddt)
        a_g = a_gn; a_e = a_en;
    end

    cl = (U_nk[1]*a_g + U_nk[3]*a_e); c2 = (U_nk[2]*a_g + U_nk[4]*a_e)

    # Calculate the probability of hopping
    if PESg
        g = g_mn(a_g, a_e, d_eg, dt);
        V_change[i] = ape(R_e) - apg(R_e)
         #FROM g TO e
    else
        g = g_mn(a_e, a_g, d_ge, dt);
        V_change[i] = apg(R_e) - ape(R_e)
         #FROM e TO g
    end

    Ee[i] = abs(a_g)^2*H_ad[1] + abs(a_e)^2*H_ad[4]
    KE[i] = 0.5*M*v[end]^2
    PE[i] = apg(R_e)
    if Ee[i] > V_change[i]
        println("swap possible ")
    end
    # Invoke surface hopping:
    chi = Distributions.Uniform()
    chi_rand = rand(chi)
    if chi_rand<g
        #println(chi_rand, "    ", g)
        if PE[i] + KE[i] > V_change[i]
            v_new = sqrt(complex(v[end]^2+2*(PE[i]-V_change[i])/M))
            v[end] = v_new
            if PESg
                println(PESg," SWAP ", g)
                PESg = false
            else
                println(PESg," SWAP ", g)
                PESg = true
            end
        end
    end
    # Update parameters for next iteration
    R_n = x[end]; v_0 = v[end]; xr = R_n/2; xl = -xr;

    if plotting
        psi_g, psi_e = adiabatic_states(phi_l_new, phi_r_new, U_nk)
        Psi_r(R) = real(a_g)*psi_g(R) + real(a_e)*psi_e(R)
        plot(apg, color=:blue, -1.5,1.5)
        plot!(ape, color=:blue)
        plot!(phi_l_new)
        plot!(phi_r_new)
        plot!(Psi_r)
        scatter!([xr],[apg(xr)], markershape=:circle)
        scatter!([xl],[apg(xl)], markershape=:circle)
        gui()
    end

end

TE = KE + PE
Ee;
plot(KE); plot!(PE); plot!(TE); plot!(V_change)
