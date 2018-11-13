function determine_individual_decay_time_constants( m, c; take_every_n_sample_for_fit = 25 )
    T = Float32
    n_channel = get_number_of_channel(m)
    n_segments = n_channel - 1
    new_pulse_format = is_new_pulse_format(m)
    n_total_events = get_number_of_events(m)

    inputfiles = gather_absolute_paths_to_hdf5_input_files(m)

    time_array_of_pulses = get_time_array_of_pulses(m.daq)
    decay_window_start_index = m.daq.n_samples-m.daq.decay_window_length+1
    xarr = time_array_of_pulses[decay_window_start_index:take_every_n_sample_for_fit:end]
    dataarr::Array{Float64, 1} = zeros(Float64, length(xarr)) # must be Float64 due LsqFit...
    init_params::Array{Float64, 1} = zeros(Float64, 2) # must be Float64 due LsqFit...

    init_decay_rate = Float64(50.0 / 1e6)

    tdc::T = 0
    core::UInt8 = 1

    @info "Tau Decay Fitting: $n_total_events events to process"

    for (fi, f) in enumerate(inputfiles)
        h5f = h5open(f, "r+")
        try
            g_daq = g_open(h5f, "DAQ_Data")
            g_pd  = g_open(h5f, "Processed_data")
            d_daq_pulses = d_open(g_daq, "daq_pulses")
            chunksize_daq = get_chunk(d_daq_pulses)
            n_samples, n_channel, n_events = new_pulse_format ? size(d_daq_pulses) : (size(d_daq_pulses, 2), size(d_daq_pulses, 1), size(d_daq_pulses, 3))
            chunk_n_events = get_chunk(d_daq_pulses)[end]
            # n_events = 3000 # debugging

            chunksize = n_channel, chunksize_daq[3]

            if exists(g_pd, "tau_decay_constants")
                o_delete(g_pd, "tau_decay_constants") 
            end
            d_tau_decay_constants = d_create(g_pd, "tau_decay_constants", Float32, ((n_channel, n_events),(n_channel, n_events)), "chunk", chunksize)

            evt_ranges = event_range_iterator(n_events, chunksize[2])
            if new_pulse_format
                @fastmath @inbounds begin
                    @showprogress for evt_range in evt_ranges
                    # @onthreads 1:2 for evt_range in workpartition(evt_ranges, 2, Base.Threads.threadid())
                        chunk_pulses::Array{T, 3} = d_daq_pulses[: , :, evt_range]
                        tdcs::Array{T, 2} = zeros(T, n_channel, length(evt_range))
                        for event in 1:length(evt_range)
                        # @onthreads 1:2 for event in workpartition(1:length(evt_range), 2, Base.Threads.threadid())
                            daq_pulses_evt = T.(chunk_pulses[:,:,event])  
                            for ichn in eachindex(1:n_channel)
                                baseline::T = mean(T, daq_pulses_evt[1:m.daq.baseline_length, ichn])
                                dataarr[:] = Float64.(daq_pulses_evt[decay_window_start_index:take_every_n_sample_for_fit:end, ichn] .- baseline) # must be Float64 due LsqFit...
                                init_params[:] = [dataarr[1], init_decay_rate] # must be Float64 due LsqFit...
                                tdc = 1f6 * T( GeDetSpectrumAnalyserTmp.LSQFIT( xarr, dataarr, exponential_decay, init_params, estimate_uncertainties=false).parameters[2] ) 
                                tdcs[ichn, event] = tdc
                            end
                        end
                        for i in eachindex(tdcs)
                            if isnan(tdcs[i])
                                tdcs[i]::T = 0
                            end
                            if isinf(tdcs[i])
                                tdcs[i]::T = 0
                            end
                        end
                        d_tau_decay_constants[:, evt_range] = tdcs
                    end
                end
            else # old pulse format
                @warn "old pulse format not yet tested"
                @fastmath @inbounds begin
                    @showprogress for evt_range in evt_ranges
                        chunk_pulses::Array{T, 3} = d_daq_pulses[: , :, evt_range]
                        tdcs::Array{T, 2} = zeros(T, n_channel, length(evt_range))
                        for event in 1:length(evt_range)
                            daq_pulses_evt = T.(transpose(chunk_pulses[:,:,event]))  
                            for ichn in eachindex(1:n_channel)
                                baseline::T = mean(T, daq_pulses_evt[1:m.daq.baseline_length, ichn])
                                dataarr[:] = Float64.(daq_pulses_evt[decay_window_start_index:take_every_n_sample_for_fit:end, ichn] .- baseline) # must be Float64 due LsqFit...
                                init_params[:] = [dataarr[1], init_decay_rate] # must be Float64 due LsqFit...
                                tdc = 1f6 * T( GeDetSpectrumAnalyserTmp.LSQFIT( xarr, dataarr, exponential_decay, init_params, estimate_uncertainties=false).parameters[2] ) 
                                tdcs[ichn, event] = tdc
                            end
                        end
                        for i in eachindex(tdcs)
                            if isnan(tdcs[i])
                                tdcs[i]::T = 0
                            end
                            if isinf(tdcs[i])
                                tdcs[i]::T = 0
                            end
                        end
                        d_tau_decay_constants[:, evt_range] = tdcs
                    end
                end
            end
          
            close(h5f)
        catch err
            close(h5f)
            error(err)
        end
    end
    return nothing
end

function determine_individual_decay_time_constants( m; take_every_n_sample_for_fit = 25 )
    c = read_analysis_result_dataset(m, "calibration_matrix")
    determine_individual_decay_time_constants(m, c, take_every_n_sample_for_fit=take_every_n_sample_for_fit)
end

function determine_decay_time_constants(m; energy_range=200:3000, create_plots=true)
    n_total_events = get_number_of_events(m)
    n_channel = get_number_of_channel(m)
    inputfiles = gather_absolute_paths_to_hdf5_input_files(m)
    core::UInt8 = 1
    fitted_tau_decay_constants::Array{Float64, 1} = Float64[ 0 for chn in 1:n_channel ]
    hists = Histogram[ Histogram(1:0.5:150, :left) for chn in 1:n_channel ]
    for (fi, f) in enumerate(inputfiles)
        h5f = h5open(f, "r+")
        try
            g_pd  = g_open(h5f, "Processed_data")
            d_energies = d_open(g_pd, "energies")
            d_tau_decay_constants = d_open(g_pd, "tau_decay_constants")
            d_ssidcs = d_open(g_pd, "single_segment_indices")
            T::Type = eltype(d_energies)
            energy_range = T.(energy_range)
            n_channel, n_events = size(d_energies)
            chunk_n_events = get_chunk(d_energies)[end]

            evt_ranges = event_range_iterator(n_events, chunk_n_events)
            @fastmath @inbounds begin
                @showprogress for evt_range in evt_ranges
                    chunk_energies::Array{T, 2} = d_energies[:, evt_range]
                    chunk_ssi::Array{UInt8, 1} = d_ssidcs[evt_range]
                    chunk_tdcs::Array{T, 2} = d_tau_decay_constants[:, evt_range]
                    for event in 1:length(evt_range)
                        ssi::UInt8 = chunk_ssi[event]
                        if ssi > 0
                            core_energy::T = chunk_energies[core, event]
                            if first(energy_range) <= core_energy <= last(energy_range)
                                push!(hists[core], chunk_tdcs[core, event])
                                push!(hists[ssi], chunk_tdcs[ssi, event])
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

    fit_results = GeDetSpectrumAnalyserTmp.Fit[]
    for ichn in eachindex(1:n_channel)
        h = hists[ichn]
        mp = midpoints(h.edges[1])
        μ::Float64 = mean(h)
        idx_μ = StatsBase.binindex(h, μ)
        fitrange = (μ - 10 * step(h)):(μ + 10 * step(h))
        σ::Float64 = 0.5 * std(h, mean=μ)
        A::Float64 = sum(h.weights)
        p0::Array{Float64, 1} = [ A, σ, μ ]
        fr = GeDetSpectrumAnalyserTmp.fit(h, fitrange, scaled_cauchy,p0)
        push!(fit_results, fr)
    end

    if create_plots
        det = m.detector
        p = histogramdisplay(hists, det)
        for (i, fr) in enumerate(fit_results)
            plotrange = fr.parameters[3] - 4 * fr.parameters[2], fr.parameters[3] + 4 * fr.parameters[2]
            plot!(p, fr, subplot=det.channel_display_order[i], xlims=plotrange)
        end
        savefig(m, p, "2_5_tau_decay_constants", "tau_decay_constants", fmt=:png); p = 0;
    end
    tdcs::Array{Float64, 1} = [ fr.parameters[3] for fr in fit_results ]

    return tdcs, hists, fit_results
end