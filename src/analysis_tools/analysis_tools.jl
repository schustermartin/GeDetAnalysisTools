include("math_functions.jl")

# general
include("single_segment_selection.jl")

# standard calibration & cross-talk correction
include("standard_calibration/0_standard_calibration.jl")

#pile-up-rejection
include("baseline_slope.jl")

# Source analysis
include("determine_source_facing_segments.jl")
include("superpulses.jl")
include("crystal_axes.jl")
include("segment_boundaries.jl")
