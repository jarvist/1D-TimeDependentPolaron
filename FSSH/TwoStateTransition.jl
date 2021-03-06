
using Distributions
using Calculus
using OrdinaryDiffEq
using Plots
gr()

# Force constant from Wang = 14500 amu/ps^2 = 14500 * 1.66e-27(kg) / 10^(-24)(s^2) = 24.07
# equilibrium bond length for hydrogen = 7.5e-11m
# energy = 4.5 eV = 7.2e-19
# coupling constant = 3500 cm−1/Å = 3500 *100*10^10 = 3500 e-8 = 3.5 e-5


#
# ---------------
# constants used 04/08/18:
# J_if = 0.136056925
# en = 122.4512325 (but plotted with en = 0)
# R_0 = 1.4423076923076923
# K = 1.3204114285714283
# ----------------

const M = 1

"""
Function to define the non-adiabatic potential for a single site centered on R_0
Basic quadratic approximation is used with minimum energy e_0
"""
function non_adiabatic_potential(K, R_0, e_i0, e_f0)
    n_a_p_i = function (R) return (K*(R+R_0/2)^2 + e_i0) end
    n_a_p_f = function (R) return (K*(R-R_0/2)^2 + e_f0) end
    return n_a_p_i, n_a_p_f
end

"""
Function to define the adiabatic potentials from the non_adiabatic_potentials
"""
function adiabatic_potential(K, R_0, e_i0, e_f0, J_if)
    n_a_p_i, n_a_p_f = non_adiabatic_potential(K, R_0, e_i0, e_f0)
    a_p_m = function (R) return 0.5*(n_a_p_i(R)+n_a_p_f(R)) - 0.5*sqrt((n_a_p_i(R)-n_a_p_f(R))^2+4*J_if^2) end
    a_p_p = function (R) return 0.5*(n_a_p_i(R)+n_a_p_f(R)) + 0.5*sqrt((n_a_p_i(R)-n_a_p_f(R))^2+4*J_if^2) end
    return a_p_m, a_p_p
end

"""
Function to plot diabatic and adiabatic potentials for simple 2-atom system.
"""
function plot_potentials(K, R_0, e_i0, e_f0, J_if)
    plot(r -> non_adiabatic_potential(K, r, R_0, e_i0, e_f0)[1], -R_0, R_0, color = :red, label = "e_i")
    plot!(r -> non_adiabatic_potential(K, r, R_0, e_i0, e_f0)[2], -R_0, R_0, color = :red, label = "e_f")
    plot!(r -> adiabatic_potential(K, r,R_0, e_i0, e_f0, J_if)[1], -R_0, R_0, color = :blue, label = "e_+")
    plot!(r -> adiabatic_potential(K, r,R_0, e_i0, e_f0, J_if)[2], -R_0, R_0, color = :blue, label = "e_-")
    gui()
end

"""
Function to calculate the force field from the potential energy curves, V.
using the relation F = m*a = -∇V
"""
function force_from_V(a_p_m, a_p_p, plotting::Bool=true)
    F_p(x) = -derivative(R -> a_p_p(R),x)
    F_m(x) = -derivative(R -> a_p_m(R),x)

    if plotting
        plot(x->F_p(x), -R_0, R_0, label = "F_p")
        plot!(x->F_m(x), -R_0, R_0, label = "F_m")
        gui()
    end
    return F_m, F_p
end


"""
Function to calculate the second derivative of the potential (i.e the slope of the
force) at every point in space to be used to determine a classical propagation.
"""
function slope(F_m)
    slope_m(x) = derivative(R -> F_m(R),x)
end

"""
Function to define the equation of motion of adiabatic surface
"""
function Acceleration(du,u,p,t,F)
   x = u[1]
   dx = u[2]
   du[1] = dx
   du[2] = real(F(x/2)/M)
end

"""
choice of T and dt are important.
"""
function classical_propagation_2(x_0, v_0,T::Float64, dt, F)
    #F = function (R) return F_m(R/2) end
    tspan = (0.0,T)
    myAcceleration(du,u,p,t)=Acceleration(du,u,p,t,F)

    EoM = ODEProblem(myAcceleration,[x_0,v_0],tspan)
    sol = solve(EoM, Euler(), dt = dt)
    x = sol[1,:]
    v = sol[2,:]
    return x,v
end


"""
Function to trace motion of atom
x_0 = initial displacement from equilibrium bond length.
"""
function plot_trajectory(dx::Float64, v_0::Float64, K, R_0, e_i0, e_f0, m, dt,T::Float64)
    x = zeros(T/dt); v = zeros(T/dt); V1 = zeros(T/dt); V2 = zeros(T/dt); t = zeros(T/dt)
    x_f0 = (R_0+dx)/2; x_i0 = -x_f0; x_eq = R_0;
    x[1] = R_0 + dx; v[1] = v_0; t[1] = dt
    n_a_p_i = non_adiabatic_potential(K, R_0, e_i0, e_f0)[1]
    n_a_p_f = non_adiabatic_potential(K, R_0, e_i0, e_f0)[2]

    for i in 1:Int(T/dt-1)
        phis(R) = non_adiabatic_states(R, [-x[i]/2,x[i]/2])
        phi_1(R)= phis(R)[1]; phi_2(R) = phis(R)[2]
        J_if = overlap(phi_1,phi_2)
        a_p_m = adiabatic_potential(K, R_0, e_i0, e_f0, J_if)[1]
        a_p_p = adiabatic_potential(K, R_0, e_i0, e_f0, J_if)[2]

        F_m, F_p = force_from_V(a_p_m, a_p_p, false)

        x_new, v_new = classical_propagation_2(x[i], v[i], dt, dt, F_m)
        x[i+1],v[i+1] = x_new[end], v_new[end]
        V1 = a_p_m(-x[i+1]/2); V2 = a_p_p(-x[i+1]/2); t[i+1] = (i+1)*dt

        if i%10 == 0
            p1 = plot(n_a_p_i,-R_0,R_0, color = :blue)
            p1 = plot!(n_a_p_f,-R_0,R_0,color = :blue)
            p1 = plot!(a_p_m,-R_0,R_0,color = :red)
            p1 = plot!(a_p_p,-R_0,R_0,color = :red)
            p1 = scatter!([-x[i+1]/2], [a_p_m(-x[i+1]/2)],color = :black, markershape = :circle)
            p1 = scatter!([x[i+1]/2], [a_p_m(x[i+1]/2)],color = :red, markershape = :circle)
            gui()
        end

    end
    x1 = -x/2; x2 = x/2;
    KE = 0.5*v.^2
    PE = 2*V1
    TE = KE + PE
    # for i in 1:length(x)
    #     if i%100==0
    #         p1 = plot(r->a_p_m(r),-R_0,R_0)
    #         p1 = scatter!([x1[Int(i)]],[V1[Int(i)]],color = :black, markershape = :circle)
    #         gui()
    #     end
    # end

    p2 = plot(x, KE)
    p2 = plot!(x, PE)
    p2 = plot!(x, TE)
    p3 = plot(t, x)
    plot(p2,p3,layout = (1,2))
    gui()
    #println("min x = ",xmin, "F = ", f_xmin, "max x = ", xmax, "F = ", f_xmax )
    return x, v
end

"""
Function to calculate eigenstates
create Hamiltonian from the onsite energies (e_i and e_f) calculated and
coupling constants J.
find adiabatic states as eigenvectors of this Hamiltonian.
"""
# function adiabatic_states(x, a_p_m, a_p_p)
#     E_11 = a_p_m(x)
#     E_22 = a_p_p(x)
#     H = [E_11 0; 0 E_22]
#     psis = eigvecs(H)
#     return psis
# end

"""
Generate non adiabatic potentials as functions of R.
Find eigenvectors (adiabatic states) from the eigenstates.
Plot eigenstate as a function of nuclear coordinates.
    (Figures show excited state )
"""
function plot_states(K, R_0, dR, e_i0, e_f0, J_if)
    plot()
    x = R_0+dR
    R_n1 = R_e = x/2
    R_n2 = -x/2
    n_a_p_i = non_adiabatic_potential(K, R_0, e_i0, e_f0)[1]
    n_a_p_f = non_adiabatic_potential(K, R_0, e_i0, e_f0)[2]
    a_p_p = adiabatic_potential(K, R_0, e_i0, e_f0, J_if)[1]
    a_p_m = adiabatic_potential(K, R_0, e_i0, e_f0, J_if)[2]
    H_diabatic = [n_a_p_i(R_e) J_if; J_if n_a_p_f(R_e)]
    H2_diabatic = [n_a_p_i(-R_e) J_if; J_if n_a_p_f(-R_e)]
    #H_adiabatic = [a_p_m(R_e) 0; 0 a_p_p(R_e)p3 = ]
    cr_11, cr_12, cr_21, cr_22 = eigvecs(H_diabatic)
    cl_11, cl_12, cl_21, cl_22 = eigvecs(H2_diabatic)
    phis(R) = non_adiabatic_states(R, [R_n1, R_n2])
    psir_1(R) = cr_11*phis(R)[1] + cr_12*phis(R)[2]
    psir_2(R) = cr_21*phis(R)[1] + cr_22*phis(R)[2]
    psil_1(R) = cl_11*phis(R)[1] + cl_12*phis(R)[2]
    psil_2(R) = cl_21*phis(R)[1] + cl_22*phis(R)[2]
    # Psis(R) = adiabatic_states(R, a_p_m, a_p_p)
    # Phis(R) = non_adiabatic_states(R, [-x/2,x/2])
    # psi_1(R) = Psis(R)[:,1][1]*Phis(R)[1] + Psis(R)[:,1][2]*Phis(R)[2]
    # psi_2(R) = Psis(R)[:,2][1]*Phis(R)[1] + Psis(R)[:,2][2]*Phis(R)[2]
    # g_s_1 = zeros(201); g_s_2 = zeros(201); e_s_1 = zeros(201); e_s_2 = zeros(201); x = zeros(201);
    # for i in 0:0.1:20
    #     g_s_1[Int(i*10+1)] = psi_1(i)[1]
    #     g_s_2[Int(i*10+1)] = psi_1(i)[2]
    #     e_s_1[Int(i*10+1)] = psi_2(i)[1]
    #     e_s_2[Int(i*10+1)] = psi_2(i)[2]
    #     x[Int(i*10+1)] = i
    # end
    # scatter(x, g_s_1, markershape = :square, color = :blue, label = "ground_site_1")
    # scatter!(x, g_s_2, markershape = :square, color = :red, label = "ground_site_2")
    # scatter!(x, e_s_1, markershape = :diamond, color = :green, label = "excited_site_1")
    # scatter!(x, e_s_2, markershape = :diamond, color = :yellow, label = "excited_site_2")
    # gui()
    plot(r->psir_1(r), -R_0, R_0)
    plot!(r->psir_2(r)+1, -R_0, R_0)
    plot!(r->psil_1(r), -R_0, R_0)
    plot!(r->psil_2(r)+1, -R_0, R_0)
    gui()

end

"""
Non_adiabatic_states
N = number of sites
----------
Find a good estimate of the width of each state.
using exp(-(R-R_0)^2/2), A = 1.331335363800389
"""
function non_adiabatic_states(R,R_0s)
    phis = zeros(2)
    A = 1/0.354021770137867
    phi(R, R_0) = A*exp(-(R-R_0)^2/(0.01))
    for i in 1:length(R_0s)
        phis[i] = phi(R, R_0s[i])
    end
    return phis
end

"""
Overlap of diabatic states
"""
function overlap_phi(phi_1, phi_2)
    overlap_phi(R) = conj(phi_1(R))*phi_2(R)
    J_12 = riemann(overlap_phi, -10,10,1000)
    return J_12
end

"""
Function to calculate non-adiabatic coupling vector
using ADIABATIC states.
"""
function NACV(psi_i, psi_j, site)
    #grad_psi_1(R) = derivative(r->psi_1(r),R)
    grad_psi_j(R) = derivative(r->psi_j(site*r/2),R)
    overlap(R) = conj(psi_i(R))*grad_psi_j(R) #fix the time derivative of R.
    d_ij = riemann(overlap,-10,10,1000)
    # test(R) = conj(psi_i(R))*psi_i(R)
    # println(riemann(test,-10,10,1000))

    return d_ij

end

"""
Function to generate the full electronic wavefunction as a linear combination of
adiabatic_states by propagating the initial coeffients forwards in time (solving
the TDSE)
"""
function wavefunction_coefficients(d_12, d_21,a_i1, a_i2, a_p_m, a_p_p, R_e, dR, dt)
    # a_f1 = a_i1 - dt*im*a_p_m(R_e) - a_i2*NACV(psi_1, psi_2)*dR
    # a_f2 = a_i2 - dt*im*a_p_p(R_e) - a_i1*NACV(psi_2, psi_1)*dR
    Matrix = [-im*a_p_m(R_e) -d_12*dR/dt ; -d_21*dR/dt -im*a_p_p(R_e)]
    # println(M)
    lambdas = eigvals(Matrix)
    one,two,three,four = eigvecs(Matrix)
    C = (a_i1*four - a_i2*three)/(one*four - two*three)
    D = -(a_i1*two - a_i2*one)/(one*four - two*three)
    # println(C,"   ",D)
    # println(one,"    ", four)
    a_f1 = C*one*exp(lambdas[1]*dt) + D*three*exp(lambdas[2]*dt)
    a_f2 = C*two*exp(lambdas[1]*dt) + D*four*exp(lambdas[2]*dt)
    return a_f1, a_f2
end

"""
Funtion to generate the coefficients to link the diabatic and adiabatic states.
Psi = a1*psi_1 + a_2*psi_2
psi_1 = c_11*phi_1 + c_12*phi_2
psi_2 = c_21*phi_1 + c_22*phi_2

where psi_1 and psi_2 are the adiabatic states. for a 2 electron system these
are (0,1) and (1,0)

electron located at R_e; between 2 nuclei at R_n1 and R_n2
"""
function electronic_wavefunction(site, R_n, n_a_p_i, n_a_p_f, a_1, a_2, J_if, phis, plotting::Bool=true)
    R_e = site*R_n; R_n1 = -R_n; R_n2 = R_n
    H_diabatic = [n_a_p_i(R_e) J_if; J_if n_a_p_f(R_e)]
    phi_1(R) = phis(R)[1]; phi_2(R) = phis(R)[2]
    #H_adiabatic = [a_p_m(R_e) 0; 0 a_p_p(R_e)p3 = ]
    c_11, c_12, c_21, c_22 = eigvecs(H_diabatic)
    psi_1(R) = c_11*phi_1(R)[1] + c_12*phi_2(R)
    psi_2(R) = c_21*phi_1(R)[1] + c_22*phi_2(R)
    Psi(R) = a_1*psi_1(R) + a_2*psi_2(R)
    if plotting
        plot(r-> Psi(r), -10, 10)
    else
        return Psi, psi_1, psi_2
    end
end

"""
Function to calculate switching probabilities
Just considering state on V_11 to start.
"""
function g_12(d_12,dR, dt, a1, a2)
    g_12 = 2*real(conj(a2)*a1*dR*d_12/(dt*a1*conj(a1)))
    return g_12
end


"""
Function to calculate switching probabilities
Just considering state on V_22 to start.
"""
function g_21(d_21,dR, dt, a1, a2)
    g_21 = 2*real(conj(a1)*a2*dR*d_21/(dt*a2*conj(a2)))
    return g_21
end

"""
Function to define a change in parameters
"""
function r(R,Rn)
    r = R - Rn
    return r
end

"""
set up system parameters
plot initial state
propagate nuclei classically by dt
propagate wavefunction coefficients.

"""
function surface_hopping_2(R_0,x_0,v_0, n_a_p_i, n_a_p_f, a_1, a_2, e_i0, e_f0, phis)
# Initialise System
    PESm = true
    site = -1
    x1 = -x_0/2; x2 = x_0/2; x_n = x_0; v_n = v_0; PE = zeros(105000); KE = zeros(105000);TE = zeros(105000)
    # phis(R) = non_adiabatic_states(R, [x1, x2])
    phi_1(R) = phis(R)[1]
    phi_2(R) = phis(R)[2]
    J_if = overlap_phi(phi_1,phi_2)
    a_p_m = adiabatic_potential(K, R_0, e_i0, e_f0, J_if)[1]
    a_p_p = adiabatic_potential(K, R_0, e_i0, e_f0, J_if)[2]
    F_m, F_p = force_from_V(a_p_m, a_p_p, false)

    p1 = plot(r -> n_a_p_i(r),-R_0,R_0, color = :red)
    p1 = plot!(r -> n_a_p_f(r),-R_0,R_0, color = :red)
    p1 = plot!(r -> a_p_m(r),-R_0,R_0, color = :blue)
    p1 = plot!(r -> a_p_p(r),-R_0,R_0, color = :blue)
    p1 = scatter!([x_0/2],[a_p_m(x_0/2)], markershape = :circle, color = :black)
    p1 = scatter!([-x_0/2],[a_p_m(-x_0/2)], markershape = :circle, color = :black)
    gui()

    for j in 1:5000
        if PESm == true
            x,v = classical_propagation_2(real(x_n),real(v_n),0.01,0.0005,F_m)
        else
            x,v = classical_propagation_2(real(x_n),real(v_n),0.01,0.0005,F_p)
        end

        x1 = -x/2; x2 = x/2; Vg = [a_p_m(pos) for pos in x1];Ve = [a_p_p(pos) for pos in x1];

        KE[21*(j-1)+1:21*(j-1)+21] = 0.5*v.^2
        if PESm == true
            PE[21*(j-1)+1:21*(j-1)+21] = 2*Vg
            V_change = 2*(Ve-Vg)
        else
            PE[21*(j-1)+1:21*(j-1)+21] = 2*Ve
            V_change = 2*(Vg-Ve)
        end
        println("KE = ", KE[2*(j-1)+2], "PE = ", PE[2*(j-1)+2])
        TE[21*(j-1)+1:21*(j-1)+21] = KE[21*(j-1)+1:21*(j-1)+21]+PE[21*(j-1)+1:21*(j-1)+21]



        a1 = zeros(Complex,21); a2 = zeros(Complex,21); a1[1] = a_1; a2[1] = a_2
        x_n = x[end]; v_n = v[end]
        for i in 1:20
            phi_1n(R) = phi_1(R-(x_0/2-x2[i]))
            phi_2n(R) = phi_2(R+(x_0/2-x2[i]))
            phisn(R) = [phi_1n(R),phi_2n(R)]
            J_if = overlap_phi(phi_1n,phi_2n)
            a_p_m = adiabatic_potential(K, R_0, e_i0, e_f0, J_if)[1]
            a_p_p = adiabatic_potential(K, R_0, e_i0, e_f0, J_if)[2]
            F_m, F_p = force_from_V(a_p_m, a_p_p, false)
            Psi, psi_1, psi_2 = electronic_wavefunction(site, x2[i], n_a_p_i, n_a_p_f, a1[i], a2[i], J_if, phisn, false)
            d_12 = NACV(psi_1,psi_2, site); d_21 = NACV(psi_2, psi_1, site)
            a1[i+1], a2[i+1] = wavefunction_coefficients(d_12,d_21,a1[i],a2[i],a_p_m, a_p_p,site*x[i]/2,x[i+1]-x[i],0.0005)
        end
        Psi, psi_1, psi_2 = electronic_wavefunction(site, x2[end], n_a_p_i, n_a_p_f, a1[end], a2[end], J_if, phisn, false)

        if abs(Psi(x1[end])) > abs(Psi(x2[end]))
            site = -1
        else
            site = +1
        end

        if PESm == true
            d_12 = NACV(psi_1,psi_2,site)
            gm_12 = g_12(d_12,x[end]-x[end-1], 0.0005, a1[end], a2[end]); gp_21 = 0
        else
            d_21 = NACV(psi_2, psi_1,site)
            gp_21 = g_21(d_21,x[end]-x[end-1], 0.0005, a1[end], a2[end]); gm_12 =0
            println(gp_21)
        end

        chi = Distributions.Uniform()
        chi_rand = rand(chi)
        if chi_rand < 0.01
            println("small chi: ", chi_rand, "    ", max(real(gm_12), real(gp_21)))
        end
        if chi_rand < max(real(gm_12), real(gp_21))
            println(chi_rand, "    ", max(gm_12, gp_21))
            if PE[end] + KE[end] > V_change[end]
                v_0 = sqrt(complex(v[end]^2+2*(PE[end]-V_change[end])/M))

                if PESm
                    println(PESm," SWAP ", gm_12)
                    PESm = false
                    gm_12 = 0
                else
                    println(PESm," SWAP ", gp_21)
                    PESm = true
                    gp_21 = 0
                end
            end
        end
        #println(Psi(0.5))
        p1 = plot(r -> n_a_p_i(r),-R_0,R_0, color = :red)
        p1 = plot!(r -> n_a_p_f(r),-R_0,R_0, color = :red)
        p1 = plot!(r -> a_p_m(r),-R_0,R_0, color = :blue)
        p1 = plot!(r -> a_p_p(r),-R_0,R_0, color = :blue)
        p1 = plot!(r -> real(Psi(r))-1, -R_0, R_0, color = :cyan, linewidth = 4)
        p2 = scatter([x1[end]],[0], markershape = :circle, color = :black, xlims = (-R_0,R_0))
        p2 = scatter!([x2[end]],[0], markershape = :circle, color = :black)
        plot(p1,p2,layout = (2,1), size = (1200,800))
        gui()
    end
    TE = KE + PE
    plot(PE)

    gui()
end
# function test_wf_propagation()
#     plot_trajectory()
# end
"""
site = +/- 1 for nuclei on right or left

"""
function surface_hopping(T, dt, F_i, F_f, m, dx, R_i0, R_f0, n_a_p_i, n_a_p_f, J_if)

    PESm = true
    site = -1

    Total_E = zeros(T+2)
    xs = zeros(T+1)
    xs[1] = dx; x1_new = R_i0-dx/2; x2_new = R_f0 +dx/2; xe = x1_new; x_eq = R_f0-R_i0
    p1 = plot(r -> n_a_p_i(r),0,20)
    p1 = plot!(r -> n_a_p_f(r),0,20)
    t = 0; t_total = 0; gf_21 = 0; g1_12 = 0;
    x, v = classical_propagation_2(dx, 0., 100., 0.001)
    electronic_wavefunction = a_12

    while t+t_total<(T-2)
        #F(R) = F_m(R/2); M=1;
        x,v = classical_propagation(F_m, m, real(x), R_0, real(v), dt, slope_m)
        xs[Int(2+(t+t_total)/dt)] = real(x)
        xn_new = real(x)/2

        v_new = real(v)
        V1 = a_p_m(site*xn_new)
        V2 = a_p_p(site*xn_new)
        dR = (xs[Int(2+(t+t_total)/dt)] - xs[Int(1+(t+t_total)/dt)])
        if PESm
            Ve = V1; V_change = V2 - V1
            gi_12 = g_12(x2_new, xe, dR, dt, n_a_p_i, n_a_p_f, J_if)
            #println(x_new)
        else
            Ve = V2; V_change = V1 - V2
            gf_21 = g_21(x1_new, xe, dR, dt, n_a_p_i, n_a_p_f, J_if)
            #println(x_new)
        end
        t+=dt

        chi = Distributions.Uniform()
        chi_rand = rand(chi)
        if chi_rand < 0.01
            println(chi_rand, "    ", max(gi_12, gf_21))
        end
        K_E = 0.5*m*v_new^2
        P_E = Ve
        Total_E[Int(2+(t+t_total)/dt)] = K_E + P_E
        #println("chi = ", chi_rand, " g_12 = ", gi_12, " g_21 = ", gf_21)

        if chi_rand < max(gi_12, gf_21)

            if Total_E[Int(2+(t+t_total)/dt)] > V_change
                v_new2 = sqrt(v_new^2+2*(Ve-V_change)/m)

                if PESm
                    println(PESm," SWAP ", gi_12)
                    PESm = false
                    dx  = x2_new - R_f0
                    x, v = classical_propagation(F_f, m, dx, R_i0, R_f0, v_new2*sign(gi_12))
                    t_total += t
                    t = dt
                    gi_12 = 0
                else
                    println(PESm," SWAP ", gf_21)
                    PESm = true
                    dx  = x1_new - R_i0
                    x, v = classical_propagation(F_i, m, dx, R_i0, R_f0, v_new2*sign(gf_21))
                    t_total += t
                    t = dt
                    gf_21 = 0
                end
            end
        end

        p1 = plot(r -> n_a_p_i(r),0,20)
        p1 = plot!(r -> n_a_p_f(r),0,20)
        p1 = scatter!([x1_new],[V1], markershape = :circle, color = :black)
        p1 = scatter!([x2_new],[V2], markershape = :circle, color = :black)
        p1 = scatter!([xe],[Ve], markershape = :circle, color = :cyan, markersize = :3)
        gui()
    end

end


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
