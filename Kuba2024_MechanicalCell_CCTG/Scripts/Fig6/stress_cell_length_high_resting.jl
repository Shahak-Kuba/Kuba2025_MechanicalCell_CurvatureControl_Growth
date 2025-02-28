
using Kuba2024_MechanicalCell_CCTG
using Makie

# Approximating experimental kf: see parameter approximation section of the paper for details
KF = 8784.2;
Tb = 28.46
l = 500;
Ω₀ = l^2
P = l*4
q₀ = 1/20; 
N = Int(P*q₀) # number of cells
kf = KF/N
l_min = 5
l_max = 20


# simulation parameters
m = 1 # number of springs per cell
R₀ = 282.095  # radius of a circle with the same area as the selected boundary btype.
D = 0.00 # diffusivity in the continuum model
kₛ = 7.5 # spring constant (hookean)
Kₛ = 100  # spring constant (nonlinear)
l₀ = 15.0 # resting spring length (hookean)
L₀ = l₀ # resting spring length (nonlinear)
η = 1.0 # viscosity
growth_dir = "inward" # growth direction
domain_type = "2D" # domain type
Tmax = 24 # days
δt = 0.01 # numerical time step
btypes = ["square"]  #Options: ["circle", "triangle", "square", "hex", "star","cross"]
dist_type = "Linear" #Options: ["Linear", "sigmoid", "2sigmoid", "exp",  "sine", "cosine", "quad", "cubic"]

## Cell Behaviours (ignore this for this paper)
prolif = false; death = false; embed = false;
α = 0.0;        β = 0.0;      Ot = 0.0;
event_δt = δt

# simulations
sol_hookean = Kuba2024_MechanicalCell_CCTG.GrowthSimulation(N,m,R₀,D,kₛ,l₀,kf,η,growth_dir,domain_type,Tmax,δt,btypes,"hookean",dist_type,
                    prolif, death, embed, α, β, Ot, event_δt, 1, 5);
sol_nonlinear  = Kuba2024_MechanicalCell_CCTG.GrowthSimulation(N,m,R₀,D,Kₛ,L₀,kf,η,growth_dir,domain_type,Tmax,δt,btypes,"nonlinear",dist_type,
                    prolif, death, embed, α, β, Ot, event_δt, 1, 5);
sol_hookean2  = Kuba2024_MechanicalCell_CCTG.GrowthSimulation(N,m2,R₀,D,kₛ,l₀,kf,η,growth_dir,domain_type,Tmax,δt,btypes,"hookean",dist_type,
                    prolif, death, embed, α, β, Ot, event_δt, 1, 31);
sol_nonlinear2  = Kuba2024_MechanicalCell_CCTG.GrowthSimulation(N,m2,R₀,D,Kₛ,L₀,kf,η,growth_dir,domain_type,Tmax,δt,btypes,"nonlinear",dist_type,
                    prolif, death, embed, α, β, Ot, event_δt, 1, 31);

# Plotting                  
Density_cmap =  :cool #:rainbow1
Stress_cmap = :winter 
geo = 1
Stress_Range_Hookean = (-0.1, 0.1)
Stress_Range_Nonlinear = (-0.01, 0.01)
f_length_Hookean = Kuba2024_MechanicalCell_CCTG.plotForceLawCompareStairs(sol_hookean[geo].Density, sol_hookean[geo].t,1)
f2_hookean = Kuba2024_MechanicalCell_CCTG.plotStress2D_Quadrant(sol_hookean2[geo].u, sol_hookean2[geo].ψ, Stress_cmap, Stress_Range_Hookean, L"σ/E \; \text{[-]}", (300,300))
f_length_Nonlinear = Kuba2024_MechanicalCell_CCTG.plotForceLawCompareStairs(sol_nonlinear[geo].Density, sol_nonlinear[geo].t,1)
f2_nonlinear= Kuba2024_MechanicalCell_CCTG.plotStress2D_Quadrant(sol_nonlinear2[geo].u, sol_nonlinear2[geo].ψ, Stress_cmap, Stress_Range_Nonlinear, L"σ/E \; \text{[-]}", (300,300))
save("Scripts/Paper_Figures/Fig6_High_Hookean_Cell_Length.pdf", f_length_Hookean)
save("Scripts/Paper_Figures/Fig6_High_Hookean_Stress.pdf", f2_hookean)
save("Scripts/Paper_Figures/Fig6_High_Nonlinear_Cell_Length.pdf", f_length_Nonlinear)
save("Scripts/Paper_Figures/Fig6_High_Nonlinear_Stress.pdf", f2_nonlinear)