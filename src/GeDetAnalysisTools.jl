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
    using ParallelProcessingTools
    using GeDetPulseShapeAnalysisToolsTmp
    using GeDetSpectrumAnalyserTmp

    using Plots: @layout

    const USER_DATA_PATH   = ENV["GEDET_USER_DATA_PATH"]   # e.g.: export GEDET_USER_DATA_PATH="/remote/ceph/group/gedet/data/lab"
    const USER_OUTPUT_PATH = ENV["GEDET_USER_OUTPUT_PATH"] # e.g.: export GEDET_USER_OUTPUT_PATH="/remote/ceph/user/l/lhauert/analysis_directory"

    include("types/DAQ.jl")
    include("types/Detector.jl")
    include("types/Measurement.jl")

    include("hdf5_data/hdf5_data_functions.jl")

    include("plot_recipes.jl")

    include("calibration.jl")

    include("analysis_tools/analysis_tools.jl")

end # module
