function full_chain_standard_calibration(	m::Measurement; overwrite=false, overwrite_init_tdcs::Bool=false,
											skip_dmpa = false,
											precal_nbins::Int = 6000, precal_photon_lines = [175.5, 609.312, 668, 785., 911.204, 1120.287, 1460.830, 1764.494, 2614.533],
											peak_threshold = 30.,
											cal_photon_lines = [609.312, 1460.830], 
											α = 0.005, rtol=1e-3 , min_n_peaks = length(precal_photon_lines), 
											max_n_peaks = 4 * length(precal_photon_lines), peak_sigma = 3.0,
											ssidcs_ΔE = 3,
											quality_check_photon_lines = cal_photon_lines,
											fit_individual_decay_time_constants = false,
											c0_pre = missing)::Nothing

	println("now $(m.name)")
	if (!exists(m, "Processed_data/tau_decay_constants")) && fit_individual_decay_time_constants
		println("Determing individual decay time constants: $(m.name)")
		determine_individual_decay_time_constants(m);
	end

	if exists(m, "Results/tau_decay_constants") && overwrite_init_tdcs == true
		tdcs = read_analysis_result_dataset(m, "tau_decay_constants")
		tdcs_err = read_analysis_result_dataset(m, "tau_decay_constants")
		write_analysis_result_dataset(m, "init_tau_decay_constants", tdcs)
	    write_analysis_result_dataset(m, "init_tau_decay_constants_err", tdcs_err)
		println("using fitted tau decay times.")
	end

	if  !exists(m, "Results/init_tau_decay_constants") && overwrite_init_tdcs == false
		if exists(m, "Processed_data/tau_decay_constants")
			println("using fitted tau decay times from daq data")
			# tdcs_daq, tdcs_daq_err, hists, fit_results = daq_determine_decay_time_constants(m)#, photon_lines = precal_photon_lines)
			tdcs_daq, hists, fit_results = daq_determine_decay_time_constants(m)#, photon_lines = precal_photon_lines)
			write_analysis_result_dataset(m, "init_tau_decay_constants", tdcs_daq)
			# write_analysis_result_dataset(m, "init_tau_decay_constants_err", tdcs_daq_err)
			write_analysis_result_dataset(m, "daq_tau_decay_constants", tdcs_daq)
			# write_analysis_result_dataset(m, "daq_tau_decay_constants_err", tdcs_daq_err)
		    write_analysis_result_dataset(m, "init_tau_decay_constants_err", Float32[ -1 for ichn in eachindex(1:m.detector.n_channels)])
			println(tdcs_daq)
		else
			write_analysis_result_dataset(m, "init_tau_decay_constants", Float32[ 50 for ichn in eachindex(1:m.detector.n_channels)])
		    write_analysis_result_dataset(m, "init_tau_decay_constants_err", Float32[ -1 for ichn in eachindex(1:m.detector.n_channels)])
			println("using hard initial tau decay times: [50, 50, 50, 50, 50]")
		end
	end

	if (overwrite || !exists(m, "Processed_data/measured_pulse_amplitudes")) && !skip_dmpa
		tdcs = read_analysis_result_dataset(m, "init_tau_decay_constants")
		println("Determing measured pulse amplitudes: $(m.name)")
		determine_measured_pulse_amplitudes(m, tdcs)
	end

	if !ismissing(c0_pre) write_analysis_result_dataset(m, "core_precalibration_factor", c0_pre); end
	if (overwrite || !exists(m, "Results/core_precalibration_factor")) && ismissing(c0_pre)
		c0_pre, h_core, h_deconv = determine_core_precalibration_factor_with_mpas(m, nbins = precal_nbins, photon_lines = precal_photon_lines, α = α, rtol=rtol, min_n_peaks = min_n_peaks, max_n_peaks = max_n_peaks, peak_sigma = peak_sigma, peak_threshold = peak_threshold );
		write_analysis_result_dataset(m, "core_precalibration_factor", c0_pre);
	else
		c0_pre = read_analysis_result_dataset(m, "core_precalibration_factor");
	end

	if overwrite || !exists(m, "Results/core_calibration_factor")
		c0, h0, peak_fits, c_fit = determine_core_calibration_factor_with_mpas(m, c0_pre, photon_lines = cal_photon_lines);
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
		quality_check(m, quality_check_photon_lines = quality_check_photon_lines)
		println("Quality check done: $(m.name)")
	end

	if overwrite || !exists(m, "Processed_data/single_segment_indices")
		determine_single_channel_indices(m, c, ΔE=ssidcs_ΔE);
	end


	if (overwrite || !exists(m, "Results/tau_decay_constants"))
		# tdcs, tdcs_err, hists, fit_results = determine_decay_time_constants(m; energy_range=200:3000)
		tdcs, hists, fit_results = determine_decay_time_constants(m; energy_range=200:3000)
		write_analysis_result_dataset(m, "tau_decay_constants", tdcs);
		write_analysis_result_dataset(m, "tau_decay_constants_err", Float32[ -1 for ichn in eachindex(1:m.detector.n_channels)]);
	elseif fit_individual_decay_time_constants
		tdcs = read_analysis_result_dataset(m, "tau_decay_constants")
	end

	return nothing
end
