module GeDetAnalysisTools

    using CompressedStreams
    using DataFrames
    using Dates
    using Distributed
    using Distributions
    using DSP
    using HDF5
    import HDF5
    using JSON
    using LaTeXStrings
    using LegendDataTypes
    using LegendHDF5IO
    using LegendTextIO
    using LinearAlgebra
    using Measures
    using ParallelProcessingTools
    using ProgressMeter
    using RadiationDetectorSignals
    using RadiationSpectra
    using RadiationSpectra_beforeBAT 
    using RecipesBase
    using SIS3316Digitizers
    using SpecialFunctions
    using StaticArrays
    using Statistics
    using StatsBase
    using Unitful
    using IntervalSets

    import Plots: @layout, grid
    import Plots: plot, plot!
    import Plots: xlims, xlims!, ylims, ylims!
    import Plots: annotate!
    import Plots: png, pdf, eps
    import Plots: vline!, hline!

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

    ########## Compatibility Hacks #############
    function HDF5.create_dataset(group::HDF5.Group, name::String, type::Type, tuple, str::AbstractString, chsize) 
        HDF5.create_dataset(group, name, type, tuple, chunk = chsize)
    end
    d_open = HDF5.open_dataset
    g_open = HDF5.open_group

    d_create = HDF5.create_dataset
    g_create = HDF5.create_group

    d_write = HDF5.write_dataset
    o_delete = HDF5.delete_object
    ############################################
    include("types/DAQ.jl")
    include("types/Detector.jl")
    include("types/Measurement.jl")
    include("types/Dataset.jl")
    include("types/MFunction.jl")

    include("isotopes/isotopes.jl")

    include("data_conversion/0_data_conversion.jl")

    include("hdf5_data/hdf5_data.jl")

    include("analysis_tools/event_flagging/0_event_flagging.jl")

    include("analysis_tools/0_0_analysis_tools.jl")

    include("plot_recipes/plot_recipes.jl")
    include("plot_recipes/plotting.jl")

    include("filters/electronics_filter.jl")
    include("filters/read_parameters.jl")

    include("utils/utils.jl")

end # module
