
using QuadGK
using Distributions
using Calculus
using Plots
using OrdinaryDiffEq
gr()


"""
Function to define the non-adiabatic potential for a single site centered on R_0
Basic quadratic approximation is used with minimum energy e_0
"""
function non_adiabatic_potential(R, R_0, e_i0, e_f0)
    n_a_p_i = (R+R_0/2)^2 + e_i0
    n_a_p_f = (R-R_0/2)^2 + e_f0
    return n_a_p_i, n_a_p_f
end

"""
Function to define the adiabatic potentials from the non_adiabatic_potentials
"""
function adiabatic_potential(R, R_0, e_i0, e_f0, J_if)
    n_a_p_i, n_a_p_f = non_adiabatic_potential(R, R_0, e_i0, e_f0)
    a_p_m = 0.5*(n_a_p_i+n_a_p_f) - 0.5*sqrt((n_a_p_i-n_a_p_f)^2+4*J_if^2)
    a_p_p = 0.5*(n_a_p_i+n_a_p_f) + 0.5*sqrt((n_a_p_i-n_a_p_f)^2+4*J_if^2)
    return a_p_m, a_p_p
end

"""
Function to plot diabatic and adiabatic potentials for simple 2-atom system.
"""
function plot_potentials(R_0, e_i0, e_f0, J_if)
    plot(r -> non_adiabatic_potential(r, R_0, e_i0, e_f0)[1], 0, 20, color = :red, label = "e_i")
    plot!(r -> non_adiabatic_potential(r, R_0, e_i0, e_f0)[2], 0, 20, color = :red, label = "e_f")
    plot!(r -> adiabatic_potential(r,R_0, e_i0, e_f0, J_if)[1], 0, 20, color = :blue, label = "e_+")
    plot!(r -> adiabatic_potential(r,R_0, e_i0, e_f0, J_if)[2], 0, 20, color = :blue, label = "e_-")
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
        plot(x->F_p(x), -10, 10, label = "F_p")
        plot!(x->F_m(x), -10, 10, label = "F_m")
        gui()
    end
    return F_m, F_p
end

"""
Analytic solution to adiabatic force
"""
function adiabatic_force(R)
    f = -2*R - (R0*(R0*R+de/2))/(sqrt((R*R0 + de/2)^2 + 4*jif^2))
    return f
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
function Acceleration(du,u,p,t)
   x = u[1]
   dx = u[2]
   du[1] = dx
   du[2] = F(x)/M
end

function classical_propagation_2(x_0, v_0,T::Float64, dt)
    tspan = (0.0,T)
    EoM = ODEProblem(Acceleration,[x_0,v_0],tspan)
    sol = solve(EoM, Euler(), dt = dt)
    x = sol[1,:]
    v = sol[2,:]
    return x,v
end


"""
Function to classically propagate nuclei along lower potential energy surface
initial_position = x_0
initial_velocity = v_0
acceleration = F_m(initial_position)/mass
time_step = dt
"""
function classical_propagation(F_m, m, x_0, R_0, v_0, dt, slope_m, plotting::Bool=false)
    x_eq = R_0
    # s = (F_m(x_0/2) - F_m(x_0/2+1))/1.0
    omega = sign(slope_m(x_0/2))sqrt(abs(slope_m(x_0/2)/m))
    phase = atan(v_0/(omega*(x_0-x_eq)))/omega
    #println(phase, "     ", (x_0-x_eq)/(2*cos(omega*phase)), "    ", (v_0)/(2*omega*sin(omega*phase)))
    A = (x_0-x_eq)/(2*cos(omega*phase))
    x = x_eq + A*(exp(-im*omega*(dt-phase)) + exp(im*omega*(dt-phase)))
    v = omega*im*A*(-exp(-im*omega*(dt-phase)) + exp(im*omega*(dt-phase)))

    if plotting
        xs = zeros(700)
        vs = zeros(700)
        ts = zeros(700)
        V = zeros(700)
        x = x_0
        v = v_0

        for i in 1:700

            dt = 0.5
            t = i*0.5
            xs[i] = x
            vs[i] = v
            V[i] = apm(x/2)
            omega = sign(slope_m(x/2))sqrt(abs(slope_m(x/2)/m))
            phase = atan(v/(omega*(x-x_eq)))/omega
            A = (x-x_eq)/(2*cos(omega*phase))
            if x<0.01
                x=0
                v=0
            else
                x = real(x_eq + A*(exp(-im*omega*(dt-phase)) + exp(im*omega*(dt-phase))))
                v = real(omega*im*A*(-exp(-im*omega*(dt-phase)) + exp(im*omega*(dt-phase))))
            end
            ts[i] = t
        end
        plot!(ts, [xs,vs,V])
    else
        return x,v
    end

end

"""
Function to trace motion of atom
x_0 = initial displacement from equilibrium bond length.
"""
function plot_trajectory(dx::Float64, v_0::Float64, R_0, e_i0, e_f0, J_if, m, dt,T::Float64)
    x_f0 = (R_0+dx)/2; x_i0 = -x_f0; x_eq = R_0; x = R_0 + dx;
    n_a_p_i(R) = non_adiabatic_potential(R, R_0, e_i0, e_f0)[1]
    n_a_p_f(R) = non_adiabatic_potential(R, R_0, e_i0, e_f0)[2]
    a_p_m(R) = adiabatic_potential(R, R_0, e_i0, e_f0, J_if)[1]
    a_p_p(R) = adiabatic_potential(R, R_0, e_i0, e_f0, J_if)[2]
    F_i, F_f = force_from_V(n_a_p_i, n_a_p_f, false)
    F_m, F_p = force_from_V(a_p_m, a_p_p, false)
    # slope_m = slope(F_m)
    # x, v = classical_propagation(F_i, m, R_0+dx, R_0, v_0, )
    F(R) = -F_m(R/2); M = m
    x,v = classical_propagation_2(x,v_0,T,dt)
    #p1 = scatter([x_0], markershape = :circle, color = :black)
    p1 = plot(r -> n_a_p_i(r),-10,10)
    # p1 = plot!(r -> a_p_p(r),0,20, color = :red)
    # p1 = plot!(r -> a_p_m(r),0,20, color = :red)
    x1 = -x/2; x2 = x/2;
    force = [F(pos) for pos in x]
    V1 = [a_p_m(pos) for pos in x1]
    V2 = [a_p_p(pos) for pos in x1]
    KE = 0.5*m*v.^2
    PE = 2*V1
    TE = KE + PE
    t = range(0.,dt,length(x))
    for i in 1:length(x)
        if i%1000==0
            p1 = plot(r->a_p_m(r),-10,10)
            p1 = scatter!([x1[Int(i)]],[V1[Int(i)]],color = :black, markershape = :circle)
            gui()
        end
    end

    p2 = plot(x, [KE, PE,TE])
    p3 = plot(t, x)
    plot(p1,p2,p3,layout = (1,3))
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
function adiabatic_states(x, a_p_m, a_p_p)
    E_11 = a_p_m(x)
    E_22 = a_p_p(x)
    H = [E_11 0; 0 E_22]
    psis = eigvecs(H)
    return psis
end

"""
Generate non adiabatic potentials as functions of R.
Find eigenvectors (adiabatic states) from the eigenstates.
Plot eigenstate as a function of nuclear coordinates.
    (Figures show excited state )
"""
function plot_states(R_0, dR, e_i0, e_f0, J_if)
    plot()
    x = R_0+dR
    R_n1 = R_e = x/2
    R_n2 = -x/2
    n_a_p_i(R) = non_adiabatic_potential(R, R_0, e_i0, e_f0)[1]
    n_a_p_f(R) = non_adiabatic_potential(R, R_0, e_i0, e_f0)[2]
    a_p_p(R) = adiabatic_potential(R, R_0, e_i0, e_f0, J_if)[1]
    a_p_m(R) = adiabatic_potential(R, R_0, e_i0, e_f0, J_if)[2]
    H_diabatic = [n_a_p_i(R_e) J_if; J_if n_a_p_f(R_e)]
    #H_adiabatic = [a_p_m(R_e) 0; 0 a_p_p(R_e)p3 = ]
    c_11, c_12, c_21, c_22 = eigvecs(H_diabatic)
    phis(R) = non_adiabatic_states(R, [R_n1, R_n2])
    psi_1(R) = c_11*phis(R)[1] + c_12*phis(R)[2]
    psi_2(R) = c_21*phis(R)[1] + c_22*phis(R)[2]
    # Psis(R) = adiabatic_states(R, a_p_m, a_p_p)
    # Phis(R) = non_adiabatic_states(R, [-x/2,x/2])
    # psi_1(R) = Psis(R)[:,1][1]*Phis(R)[1] + Psis(R)[:,1][2]*Phis(R)[2]
    # psi_2(R) = Psis(R)[:,2][1]*Phis(R)[1] + Psis(R)[:,2][2]*Phis(R)[2]
    g_s_1 = zeros(201); g_s_2 = zeros(201); e_s_1 = zeros(201); e_s_2 = zeros(201); x = zeros(201);
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
    plot(r->psi_1(r), -10, 10)
    plot!(r->psi_2(r), -10, 10)
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
using ADIABATIC states.
"""
function NACV(psi_1, psi_2, dR, dt)
    #grad_psi_1(R) = derivative(r->psi_1(r),R)
    grad_psi_2(R) = derivative(r->psi_2(r),R)
    overlap(R) = dot(conj(psi_1(R)),grad_psi_2(R))*(dR/dt) #fix the time derivative of R.
    d_jk = riemann(psi_1, real(-10), real(10),1000,"simpsons")

    return d_jk

end

"""
Function to generate the full electronic wavefunction as a linear combination of
adiabatic_states by propagating the initial coeffients forwards in time (solving
the TDSE)
"""
function wavefunction_coefficients(a_i1,a_i2, a_p, R_e,n_a_p_i, n_a_p_f, J_if, R_n1, R_n2, dR, dt)
    a_f1 = a_i1 - dt*(im*a_p(R_e) - a_i2*NACV(R_e, n_a_p_i, n_a_p_f, J_if, R_n1, R_n2, dR, dt))
    a_f2 = a_i2 - dt*(im*a_p(R_e) - a_i1*NACV(R_e, n_a_p_f, n_a_p_i, J_if, R_n2, R_n1, dR, dt))
    return real(a_f1), real(a_f2)
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
function electronic_wavefunction(site, R_n, n_a_p_i, n_a_p_f, a_1, a_2, J_if, plotting::Bool=true)
    R_e = site*R_n; R_n1 = -R_n; R_n2 = R_n
    H_diabatic = [n_a_p_i(R_e) J_if; J_if n_a_p_f(R_e)]
    #H_adiabatic = [a_p_m(R_e) 0; 0 a_p_p(R_e)p3 = ]
    c_11, c_12, c_21, c_22 = eigvecs(H_diabatic)
    phis(R) = non_adiabatic_states(R, [R_n1, R_n2])
    psi_1(R) = c_11*phis(R)[1] + c_12*phis(R)[2]
    psi_2(R) = c_21*phis(R)[1] + c_22*phis(R)[2]
    Psi(R) = a_1*psi_1(R) + a_2*psi_2(R)
    if plotting
        plot(r-> Psi(r), -10, 10)
    else
        return Psi, psi_1, psi_2, c_11, c_12, c_21, c_22
    end
end

"""
Function to calculate switching probabilities
Just considering state on V_11 to start.
"""
function g_12(n_a_p_i, n_a_p_f, R_n, J_if, dR, dt, a1, a2)
    d_12 = NACV(-1, n_a_p_i, n_a_p_f, R_n, J_if, dR, dt)
    g_12 = 2*real(conj(a2)*a1*dR*d_12)*dt/(a1*conj(a1))
    return g_12
end


"""
Function to calculate switching probabilities
Just considering state on V_22 to start.
"""
function g_21(n_a_p_i, n_a_p_f, R_n, J_if, dR, dt, a1, a2)
    d_21 =  NACV(1, n_a_p_i, n_a_p_f, R_n, J_if, dR, dt)
    g_21 = 2*real(conj(a1)*a2*dR*d_21)*dt/(a2*conj(a2))
    return g_21
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
    x, v = classical_propagation(F_i, m, x1_new, R_i0, R_f0, 0.)
    electronic_wavefunction = a_12

    while t+t_total<(T-2)
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
  if method == "right"
     meth(f,l,r) = f(r) * (r-l)
  elseif method == "left"
     meth(f,l,r) = f(l) * (r-l)
  elseif method == "trapezoid"
     meth(f,l,r) = (1/2) * (f(l) + f(r)) * (r-l)
  elseif method == "simpsons"
     meth(f,l,r) = (1/6) * (f(l) + 4*(f((l+r)/2)) + f(r)) * (r-l)
  end

  xs = a + (0:n) * (b-a)/n
  as = [meth(f, l, r) for (l,r) in zip(xs[1:end-1], xs[2:end])]
  sum(as)
end
