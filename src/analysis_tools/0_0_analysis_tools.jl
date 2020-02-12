include("0_1_math_functions.jl")

# standard calibration & cross-talk correction
include("standard_calibration/0_0_standard_calibration.jl")


# general
include("single_segment_selection.jl")
include("1_0_determine_rise_times.jl")
include("2_0_determine_maxmin_pulse_amplitudes.jl")
include("3_0_determine_AoverE.jl")
include("4_0_determine_tail_slopes.jl")



# Source analysis
include("determine_source_facing_segments.jl")
include("superpulses.jl")
include("crystal_axes.jl")
include("segment_boundaries.jl")
