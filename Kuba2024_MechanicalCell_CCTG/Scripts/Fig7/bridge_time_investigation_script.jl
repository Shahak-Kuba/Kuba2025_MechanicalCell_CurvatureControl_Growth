using Makie
using CairoMakie

# Time to bridge based on side length
Tb_Square = (L, kf, q₀) -> L./(4*kf*q₀)
Tb_Hex = (L, kf, q₀) -> (√3 .*L)./(4*kf*q₀)
Tb_Buenzli_2020 = (L,Tb₀,μ) -> Tb₀.*(L.^μ) 
# Analytical solution
S = LinRange(150,650,500) # range of side lengths
q₀ = 1/20; # initial cell density
kf = 87.842 # tissue formation rate
Tb_Square_2024 = Tb_Square(S,kf,q₀) 
Tb_Hex_2024 = Tb_Hex(S,kf,q₀)

# Model parameter values from Buenzli et al. 2020
Tb₀ = 0.0506460614771658
μ = 0.999093953616113 
Tb_2020 = Tb_Buenzli_2020(S,Tb₀, μ)

# Data extracted from Buenzli et al. 2020 - Figure 3
L_tested_2020 = [200,300,400,500,600]
Tb_approx_2020 = [11.1259,14.2681,17.9732,28.465,29.0009]
Tb_error_2020 = [0.131983,0.10861,0.874037,2.0047,0.629024]

function plotAnalytic_vs_Regression(S, Tb_2020, Tb_Square_2024, Tb_Hex_2024, L_tested_2020, Tb_approx_2020, Tb_error_2020)
    txtSize = 35;
    tickSize = 30;
    f = Figure(backgroundcolor=RGBf(1.0, 1.0, 1.0),
        size=(850, 850))
    ga = f[1, 1] = GridLayout()
    HexON = false
    if HexON
        gaxmain = Axis(ga[1, 1], height = 650, width=650, limits=(150,650,0,75),
                        xlabel="L [μm]", xlabelsize = txtSize, xticklabelsize = tickSize,
                        ylabel="Tb [days]", ylabelsize = txtSize, yticklabelsize = tickSize)
        
        
        CairoMakie.lines!(gaxmain,S,Tb_Square_2024,label=L"\text{Square:}\;T_{b}(L)", linewidth=4, linestyle=:solid, color=:blue)
        CairoMakie.lines!(gaxmain,S,Tb_Hex_2024,label=L"\text{Hex:}\;T_{b}(L)", linewidth=4, linestyle=:solid, color=:red)
        CairoMakie.lines!(gaxmain,S,Tb_2020,label=L"\text{Buenzli et al. 2020}", linewidth=4, linestyle=:dash, color=:black)
    else
        gaxmain = Axis(ga[1, 1], height = 650, width=650, limits=(150,650,5,40),
                        xlabel="L [μm]", xlabelsize = txtSize, xticklabelsize = tickSize,
                        ylabel="Tb [days]", ylabelsize = txtSize, yticklabelsize = tickSize)
        
        CairoMakie.lines!(gaxmain,S,Tb_Square_2024,label=L"\text{Square:}\;T_{\text{b}}(L)", linewidth=4, linestyle=:solid, color=:blue)
        CairoMakie.lines!(gaxmain,S,Tb_2020,label=L"\text{Buenzli et al. 2020}", linewidth=4, linestyle=:dash, color=:black)
        CairoMakie.errorbars!(gaxmain,L_tested_2020,Tb_approx_2020,Tb_error_2020, color=:black, linewidth=2, whiskerwidth = 12)
        CairoMakie.scatter!(gaxmain,L_tested_2020,Tb_approx_2020, color=:red, markersize=25, marker=:hline, label=L"\text{Experimental data}")
    end


    #Legend(f[1,1],[Analytic_Sol,Square_Sol,Hex_Sol], ["Analytic Circle", "Discrete Square","Discrete Hex"])
    axislegend(gaxmain, merge = true, unique = true, labelsize=txtSize, position=:lt)
    return f
end

f = plotAnalytic_vs_Regression(S, Tb_2020, Tb_Square_2024, Tb_Hex_2024,L_tested_2020, Tb_approx_2020, Tb_error_2020)
save("Scripts/Paper_Figures/Fig7_Side_Length_Compare.pdf", f)