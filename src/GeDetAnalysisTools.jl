__precompile__()

module GeDetAnalysisTools

    using Dates
    using HDF5
    using StatsBase
    using RecipesBase
    using Statistics
    using StaticArrays
    using ProgressMeter
    using LinearAlgebra
    using Distributions
    using ParallelProcessingTools
    using GeDetPulseShapeAnalysisToolsTmp
    using GeDetSpectrumAnalyserTmp
    using SIS3316
    using CompressedStreams
    using DataFrames
    using LaTeXStrings

    import Plots: @layout
    import Plots: plot, plot!
    import Plots: xlims, xlims!, ylims, ylims!
    import Plots: annotate!
    import Plots: png, pdf, eps


    const GAT = GeDetAnalysisTools
    export GAT

    const USER_DATA_PATH   = ENV["GEDET_USER_DATA_PATH"]   # e.g.: export GEDET_USER_DATA_PATH="/remote/ceph/group/gedet/data/lab"
    const USER_OUTPUT_PATH = ENV["GEDET_USER_OUTPUT_PATH"] # e.g.: export GEDET_USER_OUTPUT_PATH="/remote/ceph/user/l/lhauert/analysis_directory"

    include("types/DAQ.jl")
    include("types/Detector.jl")
    include("types/Measurement.jl")
    include("types/Dataset.jl")
    include("types/MFunction.jl")

    include("data_conversion/0_data_conversion.jl")

    include("hdf5_data/hdf5_data.jl")

    include("analysis_tools/analysis_tools.jl")
    
    include("plot_recipes/plot_recipes.jl")
    include("plot_recipes/plotting.jl")

end # module
