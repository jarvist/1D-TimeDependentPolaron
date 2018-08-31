# LetsGetTogether.jl
# Uses SurfaceHopping module to combine all functions and start the simulation


using SurfaceHopping
using Base
using Distributions

# Initialise system
plotting = false;
# A
# global hbar_2ma2 = (1.05e-34)^2/(2*9.11e-31*(5.29e-11)^2)/100;
# global e2_4pie0 = (1.6e-19)^2*(8.987551787368e9)/100;
# A_0 = (-hbar_2ma2 +3*e2_4pie0/(1.0e-10))/1.6e-19; # [eV]
# A_1 = -(4*e2_4pie0/(R_0*1e-10))/1.6e-19; # [eV/A]
# A_2 = (2*e2_4pie0/(R_0*1e-10))/1.6e-19; # [eV/A^2]
alpha = 0.55; tau = 0.1;
K = 1.0; global M = 28.07; # 1amu/ps^2 = 0.00010375eV/A^2
#ExtField = 5.0
temp =9000
kb = 8.625000000000001e-5 # eV/Kelvin
R_0  = 4.0
global mean_last
#R_n = R_0 + sqrt(kb*temp/K); xr = abs(R_n); xl = 0; v_0 = sqrt(1.6e8/1.66*kb*temp/M);  # [A]
cycles = 1; T = 2e-2; dt = 1e-6; ddt = 1e-6 #
# dt = 1 = 0.1ns (use dt = 1e-10 to balance Angstrom)
n = Int(round(dt/ddt+1)); N = Int(round(T/(2*dt)))

# Set up empty arrays
KE = zeros(N); PE = zeros(N); TE = zeros(N); Ee = zeros(N);  En = zeros(N)
V_change = zeros(N); Empty = zeros(N); Empty2 = zeros(N); Empty3 = zeros(N); Amps = zeros(N);
Empty2 = complex(Empty2); runfreq = zeros(N); Amps2 = zeros(N); temporary = 0; FieldCount = 0
freqvsR = Dict{Float64,Int64}();RateField = Dict{Float64,Float64}();FieldFreqs = Dict{Float64,Int64}();

if Base.Filesystem.isdir("$temp, from0.1_0.01 field divisions_31_08")
    Base.Filesystem.rm("$temp, from0.1_0.01 field divisions_31_08", recursive = true)
end
Base.Filesystem.mkdir("$temp, from0.1_0.01 field divisions_31_08")
Base.Filesystem.cd("$temp, from0.1_0.01 field divisions_31_08")


# p1 = plot(apG, color = :blue, label = "Adiabatic PES for ground state", linewidth = 5)
# plot!(p1,apE, color = :blue, label = "Adiabatic PES for excited state", linewidth = 5)
# plot!(p1,dpL, color = :green, label = "Diabatic PES for H--H+", linewidth = 5)
# plot!(p1,dpR, color = :green, label = "Diabatic PES for H+--H", linewidth = 5)
# plot!(p1,title = "PES for avoided crossing", xlabel = "Bond Distance, Angstrom", ylabel = "Energy, eV", leg= ((0.37,0.98)))
#savefig(p1, "Potential Energy Surfaces")
phi_l = diabatic_state(0.);
phi_r(R) = diabatic_state(R);
dphi = diabatic_derivative()
# Calculate coupling term
J_lr(R) = overlap_phi(phi_l,phi_r(R))
# Generate diabatic potential energy surfaces
alpha_k(R) = alpha+0.0*R_0
dpL,dpR,dpX = diabatic_potential(alpha,alpha_k,tau,J_lr)
H_D(R) = H_diabatic(dpL(R),dpR(R),dpX(R))
#Generate adiabatic potential energy surfaces
#apg,ape = adiabatic_potential(dpl,dpr,J_lr)
H_aD(R) = H_adiabatic(H_D(R))
apG(R) = H_aD(R)[1]; apE(R) = H_aD(R)[4]


p2 = plot()
p3 = scatter()
Mean_pos = zeros(100)
for ExtField in -0.1:0.01:0.1
    FieldCount+=1
    # Generate diabatic electonic states and form of their derivative
    #alpha_k(R) = alpha+ExtField*R

    println("External Field = ", ExtField*1000," mV/Angstrom")
    for l in 1:cycles
        #temp = temp*10
        State_Amplitudes = Array{Float64,}(5000)
        cl = Complex(1.0); cr = Complex(0.0)
        PESg = true
        R_0 = 4.0;
        #dR = rand(Normal(0,2*sqrt(kb*temp/K))); dv = rand(Normal(0,2*sqrt(1.6e8/1.66*kb*temp/M)))
        dR = 2*sqrt(kb*temp/K); dv = 2*sqrt(1.6e8/1.66*kb*temp/M)
        R_n = R_0 + dR; xr = R_n; xl = 0; v_0 = 0+dv;
        global count; count = 0
        for i in 1:N
            phi_r_new = phi_r(R_n)
            dphi = diabatic_derivative()
            J_lr_new = J_lr(R_n)
            dpl,dpr,dpx = diabatic_potential(alpha,alpha+ExtField*abs(R_n),tau,J_lr_new)
            H_d(R) = H_diabatic(dpl(R),dpr(R),dpx(R))
            H_ad(R) = H_adiabatic(H_d(R))
            apg(R) = H_ad(R)[1]; ape(R) = H_ad(R)[4]
            #apg,ape = adiabatic_potential(dpl,dpr,J_lr_new)
            # Determine which adiabatic potential energy surface is the active layer
            if PESg
                f = force(apg)
            else
                f = force(ape)
            end
            # Deteremine which site the electronic state is localised on
            # if PESg & R_n>0
            #     R_e = R_n
            # elseif PESg & R_n>0
            #     R_e = R_n
            # end
            # Perform classical propagation of the nuclei from t -> t+dt
            x,v = classical_propagation(2*dt, R_n, sign(R_n)*R_0, v_0, dt, f, M, K, 0.)
            dx = x[end] - x[1]; dx_l = -dx/2; dx_r = -dx_l

            # Generate necessary parameterd for quantum propagation
            #H_d = H_diabatic(dpl(R_n),dpr(R_n),J_lr_new)
            #H_ad = H_adiabatic(H_d)
            U_nk = diagonalise(H_d(R_n))
            #a_g, a_e = a_mn(U_nk, cl, cr);
            as = *(U_nk',[cl cr]')
            a_g = Complex(as[1]); a_e = Complex(as[2]);
            Empty2[1] = a_g; Empty2[2] = a_e
            #a_g = [a_gn; a_e = [a_en]
            d_ge, d_eg, d_gg, d_ee,ds = NACV(phi_l, phi_r_new, dphi, U_nk, R_n, dx_l, dx_r, 2*dt)
            # Perform quantum propagation of the electronic state from t -> t+dt
            for j in 2:2*n+1
                a_gn, a_en = new_a_mn(Empty2[2*j-3], Empty2[2*j-2], apg, ape, d_ge, d_eg, d_gg, d_ee, R_n, ddt)
                Empty2[j*2-1] = a_gn; Empty2[2*j] = a_en
            end
            a_g = Empty2[4*n+1]; a_e = Empty2[4*n+2]
            as[1] = a_g; as[2] = a_e
            #Empty2[i] = a_g

            #cl = (U_nk[1]*a_g + U_nk[3]*a_e); cr = (U_nk[2]*a_g + U_nk[4]*a_e)
            cs = *(as',U_nk'); cl = cs[1]; cr = cs[2]
            Amps[i] = abs(cl)^2; Amps2[i] = abs(a_g)^2
            #Ee[i] = abs(a_g)^2*apg(R_n) + abs(a_e)^2*ape(R_n)
            # Calculate the probability of hopping
            if PESg
                Empty3[i] = g_mn(a_g, a_e, d_eg, 2*dt);
                En[i] = apg(R_n)
                Ee[i] = apg(R_n)
                V_change[i] = ape(R_n) - Ee[i]
                #FROM g TO e
            else
                Empty3[i] = g_mn(a_e, a_g, d_ge, 2*dt);
                En[i] = ape(R_n)
                Ee[i] = ape(R_n)
                V_change[i] = apg(R_n) - Ee[i]
                 #FROM e TO g
            end
            #Empty3[i] = g
            if Empty3[i]>=0
                g = Empty3[i]
            else
                g = 0
                Empty3[i] = 0
            end

            KE[i] = 1.66/1.6e8*0.5*M*v[end]^2
            PE[i] = 0.5*K*(abs(x[end])-R_0)^2

            # Invoke surface hopping:
            chi = Distributions.Uniform()
            chi_rand = rand(chi)/1e3
            if chi_rand<g
                #println(chi_rand, "    ", g)
                if PESg
                    A = M*d_eg
                    B = M*(ds[3]^2+ds[4]^2)
                    d = (ds[3]+ds[4])
                else
                    A = M*d_ge
                    B = M*(ds[1]^2+ds[2]^2)
                    d = (ds[1]+ds[2])
                end
                #println(V_change[i] ,"      ", KE[i], "      ", 1.66/1.6e8*A^2/(2*B))
                if  V_change[i] < KE[i]
                    runfreq[i]+=1
                    if get(freqvsR, round(abs(R_n),1)+10*(ExtField-1), 0)==0
                        freqvsR[round(abs(R_n),1)+10*(ExtField-1)]=1
                    else
                        freqvsR[round(abs(R_n),1)+10*(ExtField-1)]+=1
                    end
                    v_new = sign(real(v[end]))*sqrt(complex(v[end]^2-2*(1.6e8/1.66*V_change[i])/M))#v[end]+d*A/B*(-1+sqrt(complex(1-2*(1.6e8/1.66*V_change[i]*B/A^2))))
                    #println(v_new-v[end])
                    v[end] = -real(v_new)

                    x[end] = -x[end]
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
            Empty[i] = R_n
            R_n = x[end]; v_0 = v[end]; xr = R_n; xl = 0;

            # if plotting
            #     psi_g, psi_e = adiabatic_states(phi_l, phi_r_new, U_nk)
            #     Psi_r(R) = real(a_g)*psi_g(R) + real(a_e)*psi_e(R)
            #     # plot(apg, color=:blue, -1.5,1.5)
            #     # plot!(ape, color=:blue)
            #     plot(phi_l,-2,7)
            #     plot!(phi_r_new)
            #     plot!(Psi_r)
            #     scatter!([0],[0], markershape=:circle)
            #     scatter!([R_n],[0], markershape=:circle)
            #     ylims!((-1.5,1.5))
            #     gui()
            # end

        end
        println(l)
        if l > 1
            mean_last = mean([mean(Amps),mean_last])
        else
            mean_last = mean(Amps)
        end

        if l == cycles
            APG = [apG(pos) for pos in Empty]
            APE = [apE(pos) for pos in Empty]
            writedlm("APE, External Field = $ExtField", APE)
            writedlm("APG, External Field = $ExtField", APG)
            writedlm("Active Surface, External Field = $ExtField", En)
            writedlm("Amplitude of state, External Field = $ExtField", Amps)
        end
        temporary+=sum(Empty3)
    end

    Mean_pos[FieldCount] = mean_last
    RateField[ExtField] = temporary/cycles
    num = 0

    # Max = sort(collect(freqvsR), by = tuple -> last(tuple))


    # if length(freqvsR)>0
    #     while Max[end][1]<0
    #         Max = Max[1:end-1]
    #     end
    #     Empty[Int(ExtField*10+51)] = Max[end][1]
    #     println(Max[end][1])
    # else
    #     Empty[Int(ExtField*10+51)] = 0
    # end

    for tuple in freqvsR
          num+=tuple[2]
    end
   FieldFreqs[ExtField] = num
end
TE = KE + PE + Ee;
Ee; V_change; Empty; Empty2; Empty3; runfreq;freqvsR;Amps;Amps2;RateField;
Mean_pos;

writedlm("Num_Hops_per_Field.txt", FieldFreqs)
writedlm("Rate_Hops_per_Field.txt", RateField)
writedlm("Freq_Hops_Histogram.txt", freqvsR)

#plot(KE); plot!(PE); plot!(TE); plot!(Ee)
