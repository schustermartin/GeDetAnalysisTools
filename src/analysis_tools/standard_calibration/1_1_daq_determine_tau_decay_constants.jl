function determine_decay_time_constants_distributions_with_daq_ss_selection(  m, cm_daq;
                                                        energy_limits=[500, 1500],
                                                        take_every_n_sample_for_fit = 1 )

    T = Float32
    energy_limits = T.(energy_limits)
    n_channel = get_number_of_channel(m)
    n_segments = n_channel - 1
    new_pulse_format = is_new_pulse_format(m)
    n_total_events = get_number_of_events(m)
    single_segment_abs_diff::T = 3

    decay_time_constants_distributions = [Histogram(30:0.5:70, :left) for ichn in 1:n_channel] # in microseconds
    decay_time_constants_vs_energy = [Histogram( ((energy_limits[1]:1:energy_limits[2]), 30:0.5:70), :left) for ichn in 1:n_channel] # in microseconds
    
    inputfiles = gather_absolute_paths_to_hdf5_input_files(m)

    time_array_of_pulses = get_time_array_of_pulses(m.daq)
    decay_window_start_index = m.daq.n_samples-m.daq.decay_window_length+1
    xarr = time_array_of_pulses[decay_window_start_index:take_every_n_sample_for_fit:end]
    dataarr::Array{Float64, 1} = zeros(Float64, length(xarr)) # must be Float64 due LsqFit...
    init_params::Array{Float64, 1} = zeros(Float64, 2) # must be Float64 due LsqFit...

    init_decay_rates = Float64.(50.0 / 1e6)

    tdc::T = 0
    core::UInt8 = 1

    fitf = RadiationSpectra.FitFunction(  exponential_decay )
     
    @info "Tau Decay Fitting: $n_total_events events to process"

    for (fi, f) in enumerate(inputfiles)
        h5f = h5open(f, "r")
        try
            g = g_open(h5f, "DAQ_Data")
            d_daq_pulses = d_open(g, "daq_pulses")
            daq_energies = calibrate_energies(T.(read(d_open(g, "daq_energies"))), cm_daq)
            n_events = size(d_daq_pulses, 3)
            chunksize = n_channel, chunksize_daq[3]
            # n_events = 2000 # debugging

            evt_ranges = event_range_iterator(n_events, chunk_n_events)
            if new_pulse_format
                @fastmath @inbounds begin
                    @showprogress for evt_range in evt_ranges
                        chunk_pulses::Array{T, 3} = d_daq_pulses[: , :, evt_range]
                        for event in 1:length(evt_range)
                            ssidx = get_single_segment_channel_index_abs(daq_energies[:, first(evt_range) + event - 1], single_segment_abs_diff)
                            if ssidx > 1
                                daq_pulses_evt = T.(chunk_pulses[:,:,event])  
                                core_energy::T = daq_energies[core, event]
                                if energy_limits[1] <= core_energy <= energy_limits[2] # ~keV
                                    for ichn in UInt8[core, ssidx]
                                        baseline::T = mean(T, daq_pulses_evt[1:m.daq.baseline_length, ichn])
                                        dataarr[:] = Float64.(daq_pulses_evt[decay_window_start_index:take_every_n_sample_for_fit:end, ichn] .- baseline) # must be Float64 due LsqFit...
                                        init_params[:] = [dataarr[1], init_decay_rates[ichn]] # must be Float64 due LsqFit...
                                        fitf.initial_parameters = init_params
                                        RadiationSpectra.lsqfit!(fitf, xarr, dataarr)
                                        tdc = 1f6 * fitf.parameters[2]
                                        # tdc = 1f6 * T( GeDetSpectrumAnalyserTmp.LSQFIT( xarr, dataarr, exponential_decay, init_params, estimate_uncertainties=false).parameters[2] ) 
                                        # tdc = 1e6 * LsqFit.curve_fit(exponential_decay, xarr, dataarr, init_params).param[2] 
                                        push!(decay_time_constants_distributions[ichn], tdc )
                                        push!(decay_time_constants_vs_energy[ichn], (core_energy, tdc))
                                    end
                                end
                            end
                        end
                    end
                end
            else # old pulse format
                @warn "old pulse format not yet tested"
                @fastmath @inbounds begin
                    @showprogress for evt_range in evt_ranges
                        chunk_pulses::Array{T, 3} = d_daq_pulses[: , :, evt_range]
                        for event in 1:length(evt_range)
                            ssidx = get_single_segment_channel_index_abs(daq_energies[:, first(evt_range) + event - 1], single_segment_abs_diff)
                            if ssidx > 1
                                daq_pulses_evt = T.(transpose(chunk_pulses[:,:,event]))  
                                core_energy::T = daq_energies[core, event]
                                if energy_limits[1] <= core_energy <= energy_limits[2] # ~keV
                                    for ichn in UInt8[core, ssidx]
                                        baseline::T = mean(T, daq_pulses_evt[1:m.daq.baseline_length, ichn])
                                        dataarr[:] = Float64.(daq_pulses_evt[decay_window_start_index:take_every_n_sample_for_fit:end, ichn] .- baseline) # must be Float64 due LsqFit...
                                        init_params[:] = [dataarr[1], init_decay_rates[ichn]] # must be Float64 due LsqFit...
                                        fitf.initial_parameters = init_params
                                        RadiationSpectra.lsqfit!(fitf, xarr, dataarr)
                                        tdc = 1f6 * fitf.parameters[2]
                                        # tdc = 1f6 * T( GeDetSpectrumAnalyserTmp.LSQFIT( xarr, dataarr, exponential_decay, init_params, estimate_uncertainties=false).parameters[2] ) 
                                        # tdc = 1e6 * LsqFit.curve_fit(exponential_decay, xarr, dataarr, init_params).param[2] 
                                        push!(decay_time_constants_distributions[ichn], tdc )
                                        push!(decay_time_constants_vs_energy[ichn], (core_energy, tdc))
                                    end
                                end
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
    end
    return decay_time_constants_distributions, decay_time_constants_vs_energy
end