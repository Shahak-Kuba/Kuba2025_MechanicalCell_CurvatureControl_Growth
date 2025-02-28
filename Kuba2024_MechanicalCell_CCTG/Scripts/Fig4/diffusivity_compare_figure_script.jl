using Kuba2024_MechanicalCell_CCTG
using Makie

# Shared variables
btype = "square"; #Options: ["circle", "triangle", "square", "hex", "star","cross"]
R₀ = 56.41895835477563; # radius of a circle with the same area as the selected boundary btype.
D_array = [1, 100, 10000]; # diffusivity array in the continuum model
kf = 87.84; # tissue formation rate
growth_dir = "inward"; # growth direction
Tmax = 5; # simulated days

# setting up simulation parameters
N = 20; # number of cells
m = 10; # number of springs per cell
l₀ = 10.0; # resting length of cells
η = 1.0 ; # viscosity
δt = 0.01; # numerical time step
dist_type = "Linear"; # distribution of cells along the interface

## Cell Behaviours (ignore this for this paper)
prolif = false; death = false; embed = false;
α = 0.0001;        βv = 0.001;      γv = 0.01;
event_δt = δt;

# Generating results
Discrete_Solution, Continuum_Solution = Kuba2024_MechanicalCell_CCTG.ComparisonSim_Density(N,m,R₀,D_array,l₀,kf,η,growth_dir,Tmax,δt,btype,dist_type, 
                                                                            prolif, death, embed, α, βv, γv, event_δt, 1, 0.0);

# Plotting
cmap = :cool
xbound = 60
ybound = 60
Cbar_min = 0.05
Cbar_max = 0.15

f4 = Kuba2024_MechanicalCell_CCTG.DiscVSContShape_plot_all(Discrete_Solution,m,Continuum_Solution, xbound, ybound, cmap, Cbar_min, Cbar_max)
save("Scripts/Paper_Figures/Fig4_Disc_VS_Cont_D.pdf", f4)