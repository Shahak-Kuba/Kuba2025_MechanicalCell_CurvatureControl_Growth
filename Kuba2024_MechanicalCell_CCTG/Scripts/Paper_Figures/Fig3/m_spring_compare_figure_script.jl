using Kuba2024_MechanicalCell_CCTG
using Makie

# Shared variables
R₀ = 56.41895835477563;
D = 25;
kf = 20;
growth_dir = "inward";
Tmax = 22.0; # days
btype = "square"; #Options: ["circle", "triangle", "square", "hex", "star","cross"]

# Discrete Simulation Variables
# set random seed number for reproducability 
seed = 88;

# setting up simulation parameters
N = 20; # number of cells
m1 = 1; # number of springs per cell
m2 = 4;
l₀ = 10.0;
η = 1.0 ;
δt = 0.01;
dist_type = "Linear"; #Options: ["Linear", "sigmoid", "2sigmoid", "exp",  "sine", "cosine", "quad", "cubic"]

## Cell Behaviours
prolif = false; death = false; embed = false;
α = 0.0001;        βv = 0.001;      γv = 0.01;
event_δt = δt;

# Continuum simulation variavbles
Av = 0.0;

# Generating results
Discrete_Solution_m1, Discrete_Solution_m2, Continuum_Solution = Kuba2024_MechanicalCell_CCTG.ComparisonSim(N,m1,m2,R₀,D,l₀,kf,η,growth_dir,Tmax,δt,btype,dist_type, 
                                                                            prolif, death, embed, α, βv, γv, event_δt, seed, Av);

indicies = [1,1,5,5,11,11]
num_cols = 2
f1 = Kuba2024_MechanicalCell_CCTG.DiscVSContDensity_plot_all(Discrete_Solution_m1, m1, Discrete_Solution_m2, m2, Continuum_Solution, indicies, num_cols)
save("Scripts/Paper_Figures/Fig3_m_springs_compare_density_profiles.pdf",f1)

Density_cmap =  :cool #:rainbow1
Density_Range = (0.05,0.15)

f2 = Kuba2024_MechanicalCell_CCTG.plotResults2D_Fig3(Discrete_Solution_m1.u, Discrete_Solution_m1.Density, Density_cmap, Density_Range,  L"\text{Density} \; q \; [\text{μm^{-1}}]", (60,60), N, m1, 10)
save("Scripts/Paper_Figures/Fig3_square_infill_m1_springs.pdf",f2)
f3 = Kuba2024_MechanicalCell_CCTG.plotResults2D_Fig3(Discrete_Solution_m2.u, Discrete_Solution_m2.Density, Density_cmap, Density_Range,  L"\text{Density} \; q \; [\text{μm^{-1}}]", (60,60), N, m2, 10)
save("Scripts/Paper_Figures/Fig3_square_infill_m4_springs.pdf",f3)
