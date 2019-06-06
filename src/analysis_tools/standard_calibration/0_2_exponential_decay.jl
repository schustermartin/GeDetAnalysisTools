function determine_individual_decay_time_constants( m; take_every_n_sample_for_fit = 25 )
    T::DataType = Float32
    n_channel::Int = get_number_of_channel(m)
    n_segments::Int = n_channel - 1
    new_pulse_format::Bool = is_new_pulse_format(m)
    n_total_events::Int = get_number_of_events(m)

    inputfiles = gather_absolute_paths_to_hdf5_input_files(m)

    time_array_of_pulses = get_time_array_of_pulses(m.daq)
    decay_window_start_index = m.daq.n_samples-m.daq.decay_window_length+1
    xarr::Vector{T} = time_array_of_pulses[decay_window_start_index:take_every_n_sample_for_fit:end]
    dataarr::Vector{T} = zeros(T, length(xarr)) # must be Float64 due LsqFit...
    init_params::Vector{T} = zeros(T, 2) # must be Float64 due LsqFit...

    init_decay_rate = Float64(50.0 / 1e6)

    core::UInt8 = 1
    tdc_factor::T = 1e6

    fitf = RadiationSpectra.FitFunction{T}(  exponential_decay, 1, 2 )
    set_fitranges!(fitf, ( (xarr[1], xarr[end]), ))

    @info "Tau Decay Fitting: $n_total_events events to process"

    for (fi, f) in enumerate(inputfiles)
        h5f = h5open(f, "r+")
        try
            g_daq = g_open(h5f, "DAQ_Data")
            
            g_pd = exists(h5f, "Processed_data") ? g_open(h5f, "Processed_data") : g_create(h5f, "Processed_data")
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
                                dataarr[:] = T.(daq_pulses_evt[decay_window_start_index:take_every_n_sample_for_fit:end, ichn] .- baseline) # must be Float64 due LsqFit...
                                # init_params[:] = [dataarr[1], init_decay_rate] # must be Float64 due LsqFit...
                                set_initial_parameters!(fitf, [dataarr[1], init_decay_rate])
                                RadiationSpectra.lsqfit!(fitf, xarr, dataarr)
                                tdcs[ichn, event] = tdc_factor * fitf.fitted_parameters[2]
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
                @error("Old pulse format: Not supported anymore. Reconvert data")
            end

            close(h5f)
        catch err
            close(h5f)
            error(err)
        end
    end
    return nothing
end
