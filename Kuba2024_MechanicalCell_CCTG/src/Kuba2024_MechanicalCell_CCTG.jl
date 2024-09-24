module Kuba2024_MechanicalCell_CCTG

    # PACKAGES USED for solving equations
    using Base
    using DifferentialEquations
    using LinearAlgebra
    using Random
    using ElasticArrays
    using QuadGK
    using Roots
    # PACKAGES FOR DATA SMOOTHING
    using Loess
    # PACKAGES USED for benchmarking
    using BenchmarkTools
    # PACKAGES USED for plotting
    using Plots
    using Makie
    using CairoMakie
    using ColorSchemes
    using Colors
    # PACKAGES USED for misc
    using Printf
    using JLD2
    import FilePaths
    using CSV
    #using Tables
    using DataFrames

    # DEVELOPED SIMULATION CODE

    # discrete simulation code
    include("Discrete/GeneralEquations.jl")
    include("Discrete/DataStructs.jl")
    include("Discrete/ModifierFncs.jl")
    include("Discrete/Misc.jl")
    include("Discrete/PoreBoundaries.jl")
    include("Discrete/Model/CellMechanics.jl")
    include("Discrete/Model/TissueSecretion.jl")
    include("Discrete/Model/AnalyticSolution.jl")
    include("Discrete/ProblemSetup.jl")
    include("Discrete/TissueGrowthODEproblem.jl")
    include("Discrete/PostSimulation.jl")
    include("Discrete/GrowthSimulation.jl")

    # continuum limit simulation code
    include("Continuum/Semi-Implicit_FD/FD_ContinuumSolvers.jl")
    include("Continuum/Semi-Implicit_FD/FD_SolverFncs.jl")
    include("Continuum/FVM_K-T/FVM_ContinuumSolver.jl")
    include("Continuum/FVM_K-T/FVM_SolverFncs.jl")
    include("Continuum/ComparisonSimulation.jl") 

    # Plotting
    include("Plotting/GeneralPlotting.jl")
    include("Plotting/EmbeddedPlotting.jl")
    include("Plotting/HistogramPlotting.jl")
    include("Plotting/InterfaceAnimation.jl")
    include("Plotting/PlottingFncs1D.jl")
    include("Plotting/PlottingFncs2D.jl")
    include("Plotting/ComparisonPlottingFncs.jl") 
    include("Plotting/PlottingFncsPDE.jl")
end # module Kuba2024_MechanicalCell_CCTG
