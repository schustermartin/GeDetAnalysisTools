include("0_1_baseline_analysis.jl")
include("0_2_exponential_tail_decay.jl")

include("1_0_daq_quick_cal_with_daq_energies.jl")
include("1_1_daq_determine_tau_decay_constants.jl")

include("2_0_determine_measured_pulse_amplitudes.jl")
include("2_1_determine_core_precalibration_factor.jl")
include("2_2_determine_core_calibration_factor.jl")
include("2_3_determine_calibration_matrix.jl")
include("2_4_recalculate_energies.jl")
include("2_5_determine_single_channel_indicies.jl")
include("2_6_determine_tau_decay_constants.jl")

include("3_0_quality_check_background_lines.jl")


include("full_chain.jl")
