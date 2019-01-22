function full_chain_standard_calibration(m::Measurement; overwrite=false, overwrite_init_tdcs::Bool=false)::Nothing

	if overwrite_init_tdcs || !exists(m, "Results/init_tau_decay_constants")
		write_analysis_result_dataset(m, "init_tau_decay_constants", Float32[ 50 for ichn in eachindex(1:m.detector.n_channels)])
	    write_analysis_result_dataset(m, "init_tau_decay_constants_err", Float32[ -1 for ichn in eachindex(1:m.detector.n_channels)])
	end

	if overwrite || !exists(m, "Processed_data/measured_pulse_amplitudes") 
		tdcs = read_analysis_result_dataset(m, "init_tau_decay_constants")
		println("Determing measured pulse amplitudes: $(m.name)")
		determine_measured_pulse_amplitudes(m, tdcs)
	end

	if overwrite || !exists(m, "Results/core_precalibration_factor") 
		c0_pre, h_core, h_peaks, h_pcf = determine_core_precalibration_factor_with_mpas(m);
		write_analysis_result_dataset(m, "core_precalibration_factor", c0_pre);
	else 
		c0_pre = read_analysis_result_dataset(m, "core_precalibration_factor");
	end

	if overwrite || !exists(m, "Results/core_calibration_factor") 
		c0, h0, peak_fits, c_fit = determine_core_calibration_factor_with_mpas(m, c0_pre);
		write_analysis_result_dataset(m, "core_calibration_factor", c0)
	else
		c0 = read_analysis_result_dataset(m, "core_calibration_factor")
	end

	if overwrite || !exists(m, "Results/calibration_matrix") 
		c, ratio_hists, cal_ratio_hists, cal_peak_fits, ct_ratio_hists, ct_peak_fits = determine_calibration_matrix_with_mpas(m, c0);
		write_analysis_result_dataset(m, "calibration_matrix", c)
	else
		c = read_analysis_result_dataset(m, "calibration_matrix")
	end

	if overwrite || !exists(m, "Processed_data/energies") 
		recalculate_energies(m, c);
	end

	if overwrite || !exists(m, "Results/background_photon_lines_fit_parameters_core")
		quality_check(m)
		println("Quality check done: $(m.name)")
	end
	
	if overwrite || !exists(m, "Processed_data/single_segment_indices") 
		determine_single_channel_indices(m, c);
	end
	
	if overwrite || !exists(m, "Processed_data/tau_decay_constants") 
		println("Determing individual decay time constants: $(m.name)")
		determine_individual_decay_time_constants(m, c);
	end
	
	if overwrite || !exists(m, "Results/tau_decay_constants") 
		tdcs, tdcs_err, hists, fit_results = determine_decay_time_constants(m; energy_range=200:3000)
		write_analysis_result_dataset(m, "tau_decay_constants", tdcs);
		write_analysis_result_dataset(m, "tau_decay_constants_err", tdcs_err);
	else
		tdcs = read_analysis_result_dataset(m, "tau_decay_constants")
	end

	return nothing
end