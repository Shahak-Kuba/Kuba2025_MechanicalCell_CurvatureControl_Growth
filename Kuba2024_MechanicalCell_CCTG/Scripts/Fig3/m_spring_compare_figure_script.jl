using Kuba2024_MechanicalCell_CCTG
using Makie

# Shared variables
btype = "square"; #Options: ["circle", "triangle", "square", "hex", "star","cross"]
R₀ = 56.41895835477563; # radius of a circle with the same area as the selected boundary btype.
D = 25; # diffusivity in the continuum model
kf = 20; # tissue formation rate
growth_dir = "inward"; # growth direction
Tmax = 22.0; # simulated days

# setting up simulation parameters
N = 20; # number of cells
m1 = 1; # number of springs per cell
m2 = 4; # number of springs per cell
l₀ = 10.0; # resting length of cells
η = 1.0 ; # viscosity
δt = 0.01; # numerical time step
dist_type = "Linear"; # distribution of cells along the interface

## Cell Behaviours (ignore this for this paper)
prolif = false; death = false; embed = false; # no cell behaviours
α = 0.00;        βv = 0.00;      γv = 0.00; # no cell behaviours
event_δt = δt;


# Generating results
Discrete_Solution_m1, Discrete_Solution_m2, Continuum_Solution = Kuba2024_MechanicalCell_CCTG.ComparisonSim(N,m1,m2,R₀,D,l₀,kf,η,growth_dir,Tmax,δt,btype,dist_type, 
                                                                            prolif, death, embed, α, βv, γv, event_δt, 1, 0.0);


indicies = [1,1,5,5,11,11]; # indecies of the recorded time points for the density profiles
num_cols = 2; # number of columns in the plot
f1 = Kuba2024_MechanicalCell_CCTG.DiscVSContDensity_plot_all(Discrete_Solution_m1, m1, Discrete_Solution_m2, m2, Continuum_Solution, indicies, num_cols)
save("Scripts/Paper_Figures/Fig3_m_springs_compare_density_profiles.pdf",f1)

Density_cmap =  :cool 
Density_Range = (0.05,0.15)

f2 = Kuba2024_MechanicalCell_CCTG.plotResults2D_Fig3(Discrete_Solution_m1.u, Discrete_Solution_m1.Density, Density_cmap, Density_Range,  L"\text{Density} \; q \; [\text{μm^{-1}}]", (60,60), N, m1, 10)
save("Scripts/Paper_Figures/Fig3_square_infill_m1_springs.pdf",f2)
f3 = Kuba2024_MechanicalCell_CCTG.plotResults2D_Fig3(Discrete_Solution_m2.u, Discrete_Solution_m2.Density, Density_cmap, Density_Range,  L"\text{Density} \; q \; [\text{μm^{-1}}]", (60,60), N, m2, 10)
save("Scripts/Paper_Figures/Fig3_square_infill_m4_springs.pdf",f3)
