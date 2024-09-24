
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
N = Int(P*q₀) # number of cells
kf = KF/N
l_min = 5
l_max = 20

# set random seed number for reproducability 
seed = 2

# See parameter approximation document
# Calculating kf
KF = 8784.2;
Tb = 28.46
l = 500;
Ω₀ = l^2
P = l*4
q₀ = 1/20; 
N = 50 #Int(P*q₀) # number of cells
kf = KF/N
l_min = 10
l_max = 20


# setting up simulation parameters
m = 2 # number of springs per cell
R₀ = 80 #282.095  # shape radius μm
D = 0.00
kₛ = 1
Kₛ = 15
l₀ = 15
L₀ = ((l_max - l_min)/((kₛ/Kₛ)*((l_max^2 - l_min^2)/2 + l₀*(l_min - l_max)) - log(l_min/l_max)))
η = 1.0 
growth_dir = "inward" # Options: "inward", "outward"
domain_type = "2D"
Tmax = 8 # days
δt = 0.01
btypes = ["PerturbedCircle"]  #Options: ["circle", "triangle", "square", "hex", "star","cross"]
dist_type = "Linear" #Options: ["Linear", "sigmoid", "2sigmoid", "exp",  "sine", "cosine", "quad", "cubic"]
q_lim = 0.2
ρ_lim = q_lim * m
restoring_force = "nonlinear"

## Cell Behaviours
prolif = false; death = false; embed = false;
α = 0.0;        β = 0.0;      Ot = 0.003;
event_δt = δt

# 2D simulations
sol_hookean, 🥔, 🌻 = Kuba2024_MechanicalCell_CCTG.GrowthSimulation(N,m,R₀,D,kₛ,l₀,kf,η,growth_dir,domain_type,Tmax,δt,btypes,"hookean",dist_type,
                    prolif, death, embed, α, β, Ot, event_δt, seed, 31);

sol_nonlinear, 🥔, 🌻 = Kuba2024_MechanicalCell_CCTG.GrowthSimulation(N,m,R₀,D,Kₛ,L₀,kf,η,growth_dir,domain_type,Tmax,δt,btypes,"nonlinear",dist_type,
                    prolif, death, embed, α, β, Ot, event_δt, seed, 31);

Density_cmap =  :cool #:rainbow1
Stress_cmap = :winter 

geo = 1

Density_Range = (0.05,0.2)
Stress_Range_Hookean = (-2, 2)
Stress_Range_Nonlinear = (-2, 2)

f1_hookean = Kuba2024_MechanicalCell_CCTG.plotResults2D_Quadrant(sol_hookean[geo].u, sol_hookean[geo].Density, Density_cmap, Density_Range,  L"\text{Density} \; q \; \text{[1/μm]}", (180,160), N, m, 10)

f1_nonlinear = Kuba2024_MechanicalCell_CCTG.plotResults2D_Quadrant(sol_nonlinear[geo].u, sol_nonlinear[geo].Density, Density_cmap, Density_Range,  L"\text{Density} \; q \; \text{[1/μm]}", (180,160), N, m, 10)

save("haversianPore_hookean_cell_traj.pdf", f1_hookean)

save("haversianPore_nonlinear_cell_traj.pdf", f1_nonlinear)




