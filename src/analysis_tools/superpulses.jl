function event_range_iterator(n_events::Int, chunk_n_events::Int)::Array{UnitRange{Int64},1}
	evt_ranges = [ idx:idx+chunk_n_events-1 for idx in 1:chunk_n_events:n_events ]
	if last(evt_ranges[end]) > n_events evt_ranges[end] = first(evt_ranges[end]):n_events end
	return evt_ranges
end


function determine_superpulse(m::Measurement, channels::Vector{Int}, energyrange::AbstractRange; T::Type=Float32)
	n_total_events::Int  = get_number_of_events(m)

	inputfiles = gather_absolute_paths_to_hdf5_input_files(m)
	tau_decay_constants = read_analysis_result_dataset(inputfiles[1], "tau_decay_constants")
	c::Array{T, 2} = read_analysis_result_dataset(m, "calibration_matrix")
	core::Int = 1
	elimit_min::T = first(energyrange)
	elimit_max::T = last( energyrange)
	superpulses = [zeros(T, 1, 1) for chn in channels]
    n_individual_pulses::Vector{Int} =[ 0 for channel in channels]
	
	h5f = h5open(inputfiles[1], "r")
	try 
		g_daq = g_open(h5f, "DAQ_Data")
		g_pd = g_open(h5f, "Processed_data")
		d_energies = d_open(g_pd, "energies")
		d_daq_pulses = d_open(g_daq, "daq_pulses")
		n_samples, n_channel, n_events = size(d_daq_pulses)
		superpulses = [zeros(T, n_samples, n_channel) for chn in channels]

		energies = T.(read(d_energies))

		bl::Int = m.daq.baseline_length
		bl_inv::T = 1 / bl
		sr::T = m.daq.sampling_rate
	    decay_factors = zeros(T, n_channel)
	    for chn in 1:n_channel
	        decay_factors[chn] = exp( - (1 / sr) / (tau_decay_constants[chn] * (10^(-6.0)) ) )
	    end

		chunk_n_events::Int = 1000
		# n_events = 10000 # debugging
		evt_ranges = event_range_iterator(n_events, chunk_n_events)

		@showprogress for evt_range in evt_ranges
			chunk_pulses::Array{Float32, 3} = d_daq_pulses[: , :, evt_range]
			for i in 1:length(evt_range)
				ievt::Int = evt_range[1] + i - 1
				if elimit_max >= energies[core, ievt] >= elimit_min
					for ichn in eachindex(channels)
						if elimit_max >= energies[channels[ichn], ievt] >= elimit_min
							n_individual_pulses[ichn] += 1
							tdc_pulses = GeDetPulseShapeAnalysisToolsTmp.baseline_substraction_and_decay_correction(Float32.(chunk_pulses[:, :, i]), bl, bl_inv, decay_factors); 
							superpulses[ichn] += GeDetPulseShapeAnalysisToolsTmp.calibrate_pulses(tdc_pulses, c)
						end
					end
				end
			end
		end
		close(h5f)
	catch err
		close(h5f)
		error(err)
	end

	for ichn in eachindex(channels)
		superpulses[ichn] = superpulses[ichn] ./ n_individual_pulses[ichn]
	end
	return superpulses, n_individual_pulses
end

