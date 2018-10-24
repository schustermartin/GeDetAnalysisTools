function event_range_iterator(n_events::Int, chunk_n_events::Int)::Array{UnitRange{Int64},1}
	evt_ranges = [ idx:idx+chunk_n_events-1 for idx in 1:chunk_n_events:n_events ]
	if last(evt_ranges[end]) > n_events evt_ranges[end] = first(evt_ranges[end]):n_events end
	return evt_ranges
end



function determine_superpulse(m::Measurement, channel::Int, energyrange::AbstractRange; T::Type=Float32)::AbstractArray{T, 2} 
	n_total_events::Int  = get_number_of_events(m)

	input_filenames = gather_absolute_paths_to_hdf5_input_files(m)
	h5f = h5open(input_filenames[1], "r")
	g_daq = g_open(h5f, "DAQ_Data")
	g_pd = g_open(h5f, "Processed_data")
	g_results = g_open(h5f, "Processed_data")
	d_energies = d_open(g_pd, "energies")
	d_daq_pulses = d_open(g_daq, "daq_pulses")
	n_samples, n_channel, n_events = size(d_daq_pulses)

	energies = T.(read(d_energies))

	elimit_min::T = first(energyrange)
	elimit_max::T = last( energyrange)

	superpulse = zeros(T, n_samples, n_channel)

	core::Int = 1

	c::Array{T, 2} = read_calibration_matrix(input_filenames[1])

	bl::Int = m.daq.baseline_length
	bl_inv::T = 1 / bl
	sr::T = m.daq.sampling_rate
	tau_decay_constants = Float32.([50f0 for i in 1:n_channel])
    decay_factors = zeros(T, n_channel)
    for chn in 1:n_channel
        decay_factors[chn] = exp( - (1 / sr) / (tau_decay_constants[chn] * (10^(-6.0)) ) )
    end

    n_individual_pulses::Int = 0
	chunk_n_events::Int = 1000
	evt_ranges = event_range_iterator(n_events, chunk_n_events)
	@showprogress for evt_range in evt_ranges
		chunk_pulses::Array{Float32, 3} = d_daq_pulses[: , :, evt_range]
		for i in 1:length(evt_range)
			ievt::Int = evt_range[1] + i - 1
			if ss_event(energies[1, ievt], energies[channel, ievt], elimit_min, elimit_max)
				n_individual_pulses += 1
				tdc_pulse = PulseShapeAnalysis.baseline_substraction_and_decay_correction(Float32.(chunk_pulses[:, :, i]), bl, bl_inv, decay_factors); 
				superpulse += PulseShapeAnalysis.calibrate_pulses(tdc_pulse, c)
			end
		end
	end
	close(h5f)
	return superpulse ./ n_individual_pulses
end