
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
m = 2 # number of springs per cell
R₀ = 282.095  # radius of a circle with the same area as the selected boundary btype.
D = 0.00 # diffusivity in the continuum model
kₛ = 7.5 # spring constant (hookean)
Kₛ = kₛ / 0.2^2 # spring constant (nonlinear)
l₀ = 10.0 # resting spring length (hookean)
L₀ = l₀  # resting spring length (nonlinear)
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
                    prolif, death, embed, α, β, Ot, event_δt, 1, 31);

sol_nonlinear = Kuba2024_MechanicalCell_CCTG.GrowthSimulation(N,m,R₀,D,Kₛ,L₀,kf,η,growth_dir,domain_type,Tmax,δt,btypes,"nonlinear",dist_type,
                    prolif, death, embed, α, β, Ot, event_δt, 1, 31);

# Plotting
Density_cmap =  :cool 
Stress_cmap = :winter 
geo = 1
Density_Range = (0.05,0.1)
Stress_Range_Hookean = (-2, 2)
Stress_Range_Nonlinear = (20, 50)
f1_hookean = Kuba2024_MechanicalCell_CCTG.plotResults2D_Quadrant(sol_hookean[geo].u, sol_hookean[geo].Density, Density_cmap, Density_Range,  L"q \; \text{[μm^{-1}]}", (280,280), N, m, 10)
f2_nonlinear = Kuba2024_MechanicalCell_CCTG.plotResults2D_Quadrant(sol_nonlinear[geo].u, sol_nonlinear[geo].Density, Density_cmap, Density_Range,  L"q \; \text{[μm^{-1}]}", (280,280), N, m, 10)
save("Scripts/Paper_Figures/Fig5_square_hookean_cell_traj.pdf", f1_hookean)
save("Scripts/Paper_Figures/Fig5_square_nonlinear_cell_traj.pdf", f2_nonlinear)




