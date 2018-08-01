
using QuadGK
using Distributions
using Calculus
using Plots
gr()

"""
Function to define the non-adiabatic potential for a single site centered on R_0
Basic quadratic approximation is used with minimum energy e_0
"""
function non_adiabatic_potential(R,R_0,e_0)
    return n_a_p = (R-R_0)^2 + e_0
end

"""
Function to define the adiabatic potentials from the non_adiabatic_potentials
"""
function adiabatic_potential(R, R_i0, e_i0, R_f0, e_f0, J_if)
    n_a_p_i = non_adiabatic_potential(R, R_i0, e_i0)
    n_a_p_f = non_adiabatic_potential(R, R_f0, e_f0)
    a_p_p = 0.5*(n_a_p_i+n_a_p_f) + 0.5*sqrt((n_a_p_i-n_a_p_f)^2+4*J_if^2)
    a_p_m = 0.5*(n_a_p_i+n_a_p_f) - 0.5*sqrt((n_a_p_i-n_a_p_f)^2+4*J_if^2)
    return a_p_p, a_p_m
end

"""
Function to plot diabatic and adiabatic potentials for simple 2-atom system.
"""
function plot_potentials(R_i0, e_i0, R_f0, e_f0, J_if)
    plot(r -> non_adiabatic_potential(r, R_i0, e_i0), 0, 20, color = :red, label = "e_i")
    plot!(r -> non_adiabatic_potential(r, R_f0, e_f0), 0, 20, color = :red, label = "e_f")
    plot!(r -> adiabatic_potential(r,R_i0, e_i0, R_f0, e_f0, J_if)[1], 0, 20, color = :blue, label = "e_+")
    plot!(r -> adiabatic_potential(r,R_i0, e_i0, R_f0, e_f0, J_if)[2], 0, 20, color = :blue, label = "e_-")
    gui()
end

"""
Function to calculate the force field from the potential energy curves, V.
using the relation F = m*a = -∇V
"""
function force_from_V(a_p_p,a_p_m)
    F_p(x) = -derivative(R -> a_p_p(R),x)
    F_m(x) = -derivative(R -> a_p_m(R),x)

    plot(x->F_p(x), 0, 20, label = "F_p")
    plot!(x->F_m(x), 0, 20, label = "F_m")
    gui()
    return F_p, F_m
end

"""
Function to classically propagate nuclei along lower potential energy surface
initial_position = x_0
initial_velocity = v_0
acceleration = F_m(initial_position)/mass
time_step = dt
"""
function classical_propagation(F_m, m, x_0, R_i0, R_f0, v_0)
    x_eq = R_f0 - R_i0
    slope = (F_m(x_0/2) - F_m(x_0/2+2))/2.0
    omega = sqrt(slope/m)
    phase = atan(v_0/(omega*(x_0-x_eq)))/omega
    #println(phase, "     ", (x_0-x_min)/(2*cos(omega*phase)), "    ", (v_0)/(2*omega*sin(omega*phase)))
    A = (x_0-x_eq)/(2*cos(omega*phase))
    x(t) = x_eq + A*(exp(-im*omega*(t-phase)) + exp(im*omega*(t-phase)))
    v(t) = omega*im*A*(-exp(-im*omega*(t-phase)) + exp(im*omega*(t-phase)))
    return x, v
end

"""
Function to trace motion of atom
x_0 = initial displacement from equilibrium bond length.
"""
function plot_trajectory(dx::Float64, v_0::Float64, R_i0, R_f0, F_i, m, dt,T, n_a_p_i, n_a_p_f)
    trace = zeros(T)
    f = zeros(T)
    kinetic_energy = zeros(T)
    potential_energy = zeros(T)
    x_i0 = R_i0 - dx; x_f0 = R_f0 + dx; x_eq = R_f0 - R_i0
    x, v = classical_propagation(F_i, m, x_i0, R_i0, R_f0, v_0)
    xs = zeros(T)
    #p1 = scatter([x_0], markershape = :circle, color = :black)
    f_xmin = f_xmax = F_i(x_i0)
    p1 = plot(r -> n_a_p_i(r),0,20)
    p1 = scatter!([R_i0-x_i0],[F_i(x_i0/2)], markershape = :circle, color = :black)
    for i in 1:T
        t = i*dt
        xs[i] = real(x(t))
        x1_new = R_i0-real(x(t)-x_eq)/2
        x2_new = R_f0 + real(x(t)-x_eq)/2
        v_new = real(v(t))
        V1 = n_a_p_i(x1_new)
        V2 = n_a_p_f(x2_new)
        kinetic_energy[i] = 0.5*m*v_new^2
        potential_energy[i] = V2+V1
        p1 = plot(r -> n_a_p_i(r),0,20)
        p1 = plot!(r -> n_a_p_f(r),0,20)
        p1 = scatter!([x1_new],[V1], markershape = :circle, color = :black)
        p1 = scatter!([x2_new],[V2], markershape = :circle, color = :black)
        gui()
    end

    p2 = plot(xs, [f, kinetic_energy, potential_energy])
    plot(p1,p2,layout = (2,1))
    gui()
    #println("min x = ",xmin, "F = ", f_xmin, "max x = ", xmax, "F = ", f_xmax )
    return f, kinetic_energy, potential_energy
end

"""
Function to calculate eigenstates
create Hamiltonian from the onsite energies (e_i and e_f) calculated and
coupling constants J.
find adiabatic states as eigenvectors of this Hamiltonian.
"""
function adiabatic_states(x, n_a_p_i, n_a_p_f, J_if)
    E_i = n_a_p_i(x)
    E_f = n_a_p_f(x)
    V = J_if
    H = [E_i V; V E_f]
    psis = eigvecs(H)
    return psis
end

"""
Generate non adiabatic potentials as functions of R.
Find eigenvectors (adiabatic states) from the eigenstates.
Plot eigenstate as a function of nuclear coordinates.
    (Figures show excited state )
"""
function plot_states(R_i0, e_i0, R_f0, e_f0, J_if)
    plot()
    n_a_p_i(R) = non_adiabatic_potential(R, R_i0, e_i0)
    n_a_p_f(R) = non_adiabatic_potential(R, R_f0, e_f0)
    Psis(R) = adiabatic_states(R, n_a_p_i, n_a_p_f, J_if)
    Psi_1(R) = Psis(R)[:,1]
    Psi_2(R) = Psis(R)[:,2]
    g_s_1 = zeros(201); g_s_2 = zeros(201); e_s_1 = zeros(201); e_s_2 = zeros(201); x = zeros(201);
    for i in 0:0.1:20
        g_s_1[Int(i*10+1)] = Psi_1(i)[1]
        g_s_2[Int(i*10+1)] = Psi_1(i)[2]
        e_s_1[Int(i*10+1)] = Psi_2(i)[1]
        e_s_2[Int(i*10+1)] = Psi_2(i)[2]
        x[Int(i*10+1)] = i
    end
    scatter(x, g_s_1, markershape = :square, color = :blue, label = "ground_site_1")
    scatter!(x, g_s_2, markershape = :square, color = :red, label = "ground_site_2")
    scatter!(x, e_s_1, markershape = :diamond, color = :green, label = "excited_site_1")
    scatter!(x, e_s_2, markershape = :diamond, color = :yellow, label = "excited_site_2")
    gui()
end

"""
Non_adiabatic_states
N = number of sites
----------
Find a good estimate of the width of each state.
"""
function non_adiabatic_states(R,R_0s)
    phis = zeros(2)
    phi(R, R_0) = exp(-(R-R_0)^2/(2))
    for i in 1:2
        phis[i] = phi(R, R_0s[i])
    end
    return phis
end

"""
Function to calculate non-adiabatic coupling vector
"""
function NACV(R_1, R_2, R_i2, dt)
    Phis(R) = non_adiabatic_states(R,[R_1,R_2])
    Phi_1(R) = Phis(R)[1]
    Phi_2(R) = Phis(R)[2]
    #grad_phi_1(R) = derivative(r->Phi_1(r),R)
    grad_phi_2(R) = derivative(r->Phis(r)[2],R)
    overlap(R) = dot(conj(Phi_1(R)),grad_phi_2(R))*(R_2 - R_i2)/dt #fix the time derivative of R.
    d_jk = (quadgk(overlap, 0, 20))[1]
    return d_jk
end

"""
Function to calculate switching probabilities
Just considering state on V_11 to start.
"""
function g_12(R_1, R_2, R_i2, dt, n_a_p_i, n_a_p_f, J_if)
    d_12 = NACV(R_1, R_2, R_i2, dt)
    psi_1 = adiabatic_states(R_1, n_a_p_i, n_a_p_f, J_if)[:,1]
    g_12 = 2*real(conj(psi_1[2])*psi_1[1]*R_2*d_12)*dt/(psi_1[1]*conj(psi_1[1]))
    return g_12
end


"""
Function to calculate switching probabilities
Just considering state on V_22 to start.
"""
function g_21(R_2, R_1, R_i1, dt, n_a_p_i, n_a_p_f, J_if)
    d_21 = NACV(R_2, R_1, R_i1, dt)
    psi_1 = adiabatic_states(R_1, n_a_p_i, n_a_p_f, J_if)[:,2]
    g_21 = 2*real(conj(psi_1[1])*psi_1[2]*R_2*d_21)*dt/(psi_1[2]*conj(psi_1[2]))
    return g_21
end

"""
"""
function surface_hopping(T, dt, F_i, F_f, m, x_0, n_a_p_i, n_a_p_f, J_if)

    PESi = true

    T_E = zeros(T)
    xs = zeros(T+1)
    xs[1] = x_0; x_new = x_0
    p1 = plot(r -> n_a_p_i(r),0,20)
    p1 = plot!(r -> n_a_p_f(r),0,20)
    p1 = scatter!([x_0],[F_i(x_0)], markershape = :circle, color = :black)
    t = 0; t_total = 0; gf_21 = 0; g1_12 = 0;
    xi, vi = classical_propagation(F_i, m, x_new, 5, 0.)

    while t+t_total<(T-2)
        if PESi
            x_new = xs[Int(2+(t+t_total)/dt)] = real(xi(t))
            v_new = real(vi(t))
            V = n_a_p_i(x_new)
            V_change = n_a_p_f(x_new)
            gi_12 = g_12(10, xs[Int(1+(t+t_total)/dt)], xs[Int(2+(t+t_total)/dt)], dt, n_a_p_i, n_a_p_f, J_if)
            t+=dt
            #println(x_new)
        else
            x_new = xs[Int(2+(t+t_total)/dt)] = real(xf(t))
            v_new = real(vf(t))
            V = n_a_p_f(x_new)
            V_change = n_a_p_i(x_new)
            #println( xs[Int(1+(t+t_total)/dt)], " ", xs[Int(2+(t+t_total)/dt)])
            gf_21 = g_21(5, xs[Int(1+(t+t_total)/dt)], xs[Int(2+(t+t_total)/dt)], dt, n_a_p_i, n_a_p_f, J_if)
            t+=dt
            #println(x_new)
        end
        chi = Distributions.Uniform()
        chi_rand = rand(chi)
        K_E = 0.5*m*v_new^2
        P_E = V
        T_E[Int(2+(t+t_total)/dt)] = K_E + P_E
        #println("chi = ", chi_rand, " g_12 = ", gi_12, " g_21 = ", gf_21)
        if chi_rand < max(gi_12, gf_21)

            if T_E[Int(2+(t+t_total)/dt)] > V_change
                v_new2 = sign(v_new)*sqrt(v_new^2+2*(V-V_change)/m)

                if PESi
                    println(PESi," SWAP ", gi_12)
                    PESi = false
                    xf, vf = classical_propagation(F_f, m, x_new, 10.,v_new2)
                    t_total += t
                    t = dt
                    gi_12 = 0
                else
                    println(PESi," SWAP ", gf_21)
                    PESi = true
                    xi, vi = classical_propagation(F_i, m, x_new, 5,v_new2)
                    t_total += t
                    t = dt
                    gf_21 = 0
                end
            end
        end

        p1 = plot(r -> n_a_p_i(r),0,20)
        p1 = plot!(r -> n_a_p_f(r),0,20)
        p1 = scatter!([x_new],[V], markershape = :circle, color = :black)
        gui()
    end

end
