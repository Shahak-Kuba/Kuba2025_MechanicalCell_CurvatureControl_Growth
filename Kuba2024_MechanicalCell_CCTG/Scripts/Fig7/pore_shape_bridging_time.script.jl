
using Kuba2024_MechanicalCell_CCTG
using Makie

# See parameter approximation document
# Calculating kf
KF = 8784.2;
Tb = 28.46
l = 500;
Ω₀ = l^2
P = l*4
q₀ = 1/20;
N1 = 120 #Int(P*q₀) # number of cells
kf = KF/N1
l_min = 5
l_max = 20

# set random seed number for reproducability 
seed = 99

####### 500μm side length square pores ############
# setting up simulation parameters
m = 2 # number of springs per cell
R₀ = 282.095  # shape radius μm
D = 0.00
kₛ = 7.5
Kₛ = 150
l₀ = 10.0
L₀ = ((l_max - l_min)/((kₛ/Kₛ)*((l_max^2 - l_min^2)/2 + l₀*(l_min - l_max)) - log(l_min/l_max)))
η = 1.0 
growth_dir = "inward" # Options: "inward", "outward"
domain_type = "2D"
Tmax = 28 # days
δt = 0.01
btypes = ["square", "hex"]  #Options: ["circle", "triangle", "square", "hex", "star","cross"]
dist_type = "Linear" #Options: ["Linear", "sigmoid", "2sigmoid", "exp",  "sine", "cosine", "quad", "cubic"]
q_lim = 0.2
ρ_lim = q_lim * m

## Cell Behaviours
prolif = false; death = false; embed = false;
α = 0.0;        β = 0.0;      Ot = 0.0;
event_δt = δt

# 2D simulations

sol_nonlinear_500 = Kuba2024_MechanicalCell_CCTG.GrowthSimulation(N1,m,R₀,D,Kₛ,L₀,kf,η,growth_dir,domain_type,Tmax,δt,btypes,"nonlinear",dist_type,
                    prolif, death, embed, α, β, Ot, event_δt, seed, 31);

Density_cmap =  :cool #:rainbow1
Stress_cmap = :winter 
Density_Range = (0.05,0.2)
axisTicks = [-600, -400, -200, 0, 200, 400, 600]

f1_nonlinear = Kuba2024_MechanicalCell_CCTG.plotResults2D(sol_nonlinear_500[1].u, sol_nonlinear_500[1].Density, Density_cmap, Density_Range,  L"\text{Density} \; q \; \text{[μm^{-1}]}", (350,350), axisTicks, N1, m, 10)
f2_nonlinear = Kuba2024_MechanicalCell_CCTG.plotResults2D(sol_nonlinear_500[2].u, sol_nonlinear_500[2].Density, Density_cmap, Density_Range,  L"\text{Density} \; q \; \text{[μm^{-1}]}", (350,350), axisTicks, N1, m, 10)

save("Scripts/Paper_Figures/Fig8_square_500.png", f1_nonlinear)
save("Scripts/Paper_Figures/Fig8_hex_500.png", f2_nonlinear)


####### 750μm side length square pores ############
# setting up simulation parameters
N2 = 150
R₀ = 423.1421876608172  # shape radius μm
Tmax = 48

# 2D simulations
sol_nonlinear_750 = Kuba2024_MechanicalCell_CCTG.GrowthSimulation(N2,m,R₀,D,Kₛ,L₀,kf,η,growth_dir,domain_type,Tmax,δt,btypes,"nonlinear",dist_type,
                    prolif, death, embed, α, β, Ot, event_δt, seed, 31);

f1_nonlinear = Kuba2024_MechanicalCell_CCTG.plotResults2D(sol_nonlinear_750[1].u, sol_nonlinear_750[1].Density, Density_cmap, Density_Range,  L"\text{Density} \; q \; \text{[μm^{-1}]}", (500,500), axisTicks, N2, m, 10)
f2_nonlinear = Kuba2024_MechanicalCell_CCTG.plotResults2D(sol_nonlinear_750[2].u, sol_nonlinear_750[2].Density, Density_cmap, Density_Range,  L"\text{Density} \; q \; \text{[μm^{-1}]}", (500,500), axisTicks, N2, m, 10)

save("Scripts/Paper_Figures/Fig8_square_750.png", f1_nonlinear)
save("Scripts/Paper_Figures/Fig8_hex_750.png", f2_nonlinear)                    

####### 1000μm side length square pores ############
# setting up simulation parameters
N3 = 204
R₀ = 564.1895835477563  # shape radius μm
Tmax = 48

# 2D simulations
sol_nonlinear_1000 = Kuba2024_MechanicalCell_CCTG.GrowthSimulation(N3,m,R₀,D,Kₛ,L₀,kf,η,growth_dir,domain_type,Tmax,δt,btypes,"nonlinear",dist_type,
                    prolif, death, embed, α, β, Ot, event_δt, seed, 31);


f1_nonlinear = Kuba2024_MechanicalCell_CCTG.plotResults2D(sol_nonlinear_1000[1].u, sol_nonlinear_1000[1].Density, Density_cmap, Density_Range,  L"\text{Density} \; q \; \text{[μm^{-1}]}", (650,650), axisTicks, N3, m, 10)
f2_nonlinear = Kuba2024_MechanicalCell_CCTG.plotResults2D(sol_nonlinear_1000[2].u, sol_nonlinear_1000[2].Density, Density_cmap, Density_Range,  L"\text{Density} \; q \; \text{[μm^{-1}]}", (650,650), axisTicks, N3, m, 10)

save("Scripts/Paper_Figures/Fig8_square_1000.png", f1_nonlinear)
save("Scripts/Paper_Figures/Fig8_hex_1000.png", f2_nonlinear)       


f = Kuba2024_MechanicalCell_CCTG.plotMultiAreaVsTime(sol_nonlinear_1000[1].t,sol_nonlinear_750[1].t,sol_nonlinear_500[1].t,sol_nonlinear_1000[1].Ω,sol_nonlinear_1000[2].Ω,sol_nonlinear_750[1].Ω,sol_nonlinear_750[2].Ω,sol_nonlinear_500[1].Ω,sol_nonlinear_500[2].Ω,N3,N2,N1,kf)
save("Scripts/Paper_Figures/Fig7_area_compare.pdf", f)


