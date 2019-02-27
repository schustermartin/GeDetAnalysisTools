__precompile__()

module GeDetAnalysisTools

    using CompressedStreams
    using DataFrames
    using Dates
    using Distributed
    using Distributions
    using DSP
    using GeDetPulseShapeAnalysisToolsTmp
    using HDF5
    using LaTeXStrings
    using LegendDataTypes
    using LegendHDF5IO
    using LegendTextIO
<<<<<<< HEAD
    using LinearAlgebra
    using ParallelProcessingTools
    using ProgressMeter
    using RadiationSpectra 
    using RecipesBase
    using SIS3316
    using SpecialFunctions
    using StaticArrays
    using Statistics
    using StatsBase
    using Unitful
=======
    using DSP
    using JSON
>>>>>>> ee75b43707463d95af43866ad13690db0b958f05

    import Plots: @layout, grid
    import Plots: plot, plot!
    import Plots: xlims, xlims!, ylims, ylims!
    import Plots: annotate!
    import Plots: png, pdf, eps
    import Plots: Measures

    import LegendHDF5IO: readdata, writedata

    import Base: show, display, print, println
    import Base: sort, sort!
    import Base: getindex, size, length

    const GAT = GeDetAnalysisTools
    export GAT
    export exists

    export r, phi_side, phi_top, z

    const USER_DATA_PATH   = ENV["GEDET_USER_DATA_PATH"]   # e.g.: export GEDET_USER_DATA_PATH="/remote/ceph/group/gedet/data/lab"
    const USER_OUTPUT_PATH = ENV["GEDET_USER_OUTPUT_PATH"] # e.g.: export GEDET_USER_OUTPUT_PATH="/remote/ceph/user/l/lhauert/analysis_directory"

    include("types/DAQ.jl")
    include("types/Detector.jl")
    include("types/Measurement.jl")
    include("types/Dataset.jl")
    include("types/MFunction.jl")
    
    include("isotopes/isotopes.jl")
    
    include("data_conversion/0_data_conversion.jl")

    include("hdf5_data/hdf5_data.jl")

    include("analysis_tools/analysis_tools.jl")
    
    include("plot_recipes/plot_recipes.jl")
    include("plot_recipes/plotting.jl")
    
    include("filters/electronics_filter.jl")
    include("filters/read_parameters.jl")

end # module
