function determine_measured_pulse_amplitudes(m::Measurement, tau_decay_constants::Array{<:Real, 1})
    inputfiles = gather_absolute_paths_to_hdf5_input_files(m)
    T = Float32
    sampling_rate = T(m.daq.sampling_rate)
    bl = Int(m.daq.baseline_length )
    bl_inv = T(1 / bl)
    decay_window_start_index = m.daq.n_samples - m.daq.decay_window_length + 1
    decay_window_length_inv = T(1 / m.daq.decay_window_length) 
    tau_decay_constants = T.(tau_decay_constants)
    decay_factors = T[exp( - (1 / sampling_rate) / (tdc * (10^(-6.0)) ) ) for tdc in tau_decay_constants]

    for f in inputfiles
        new_pulse_format = is_new_pulse_format(f)
        if !new_pulse_format error("Measurement $(m.name) does not have new pulse format. Old one ist not implemented yet.") end
        h5f = h5open(f, "r+")
        try
            g_daq = g_open(h5f, "DAQ_Data")
            d_daq_pulses = d_open(g_daq, "daq_pulses")
            chunksize_daq = get_chunk(d_daq_pulses)
            n_samples, n_channel, n_events = new_pulse_format ? size(d_daq_pulses) : (size(d_daq_pulses, 2), size(d_daq_pulses, 1), size(d_daq_pulses, 3))
            chunksize = n_channel, chunksize_daq[3]
            g_pd = exists(h5f, "Processed_data") ? g_open(h5f, "Processed_data") : g_create(h5f, "Processed_data")
            if exists(g_pd, "measured_pulse_amplitudes")
                o_delete(g_pd, "measured_pulse_amplitudes") 
            end
            d_measured_pulse_amplitudes = d_create(g_pd, "measured_pulse_amplitudes", Float32, ((n_channel, n_events),(n_channel, n_events)), "chunk", chunksize)
            
            # n_events = 32000 # debugging

            @fastmath function process_chunk(chunk_pulses::Array{<:Real, 3}, evt_range::AbstractRange, bl::Int, bl_inv::T, decay_factors::Array{T, 1}, decay_window_start_index::Int)::Array{T, 2} where {T<:AbstractFloat}
                mpas = Array{T, 2}(undef, n_channel, length(evt_range))
                @inbounds for i in eachindex(1:length(evt_range))
                    tdc_pulses::Array{T, 2} = GeDetPulseShapeAnalysisToolsTmp.baseline_substraction_and_decay_correction(Float32.(chunk_pulses[:, :, i]), bl, bl_inv, decay_factors); 
                    mpas[:, i] = [GeDetPulseShapeAnalysisToolsTmp.fastmean(tdc_pulses[decay_window_start_index:end, ichn], decay_window_length_inv) for ichn in eachindex(1:n_channel) ]
                end
                return mpas
            end

            @fastmath @inbounds begin
                @showprogress for evt_range in event_range_iterator(n_events, chunksize[2])
                    d_measured_pulse_amplitudes[:, evt_range] = process_chunk(d_daq_pulses[:, :, evt_range], evt_range, bl, bl_inv, decay_factors, decay_window_start_index)
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

function determine_baseline_information(m::Measurement)
    inputfiles = gather_absolute_paths_to_hdf5_input_files(m)
    T = Float32
    sampling_rate = T(m.daq.sampling_rate)
    bl = Int(m.daq.baseline_length )
    bl_inv = T(1 / bl)
    slope_length::Int = 100
    slope_w1::UnitRange = 1:slope_length
    slope_w2::UnitRange = bl - slope_length:bl
    slope_length_inv::T = 1 / slope_length

    for f in inputfiles[1:1]
        new_pulse_format = is_new_pulse_format(f)
        if !new_pulse_format error("Measurement $(m.name) does not have new pulse format. Old one ist not implemented yet.") end
        h5f = h5open(f, "r+")
        try
            g_daq = g_open(h5f, "DAQ_Data")
            d_daq_pulses = d_open(g_daq, "daq_pulses")
            daq_pulse_T = eltype(d_daq_pulses)
            chunksize_daq = get_chunk(d_daq_pulses)
            n_samples, n_channel, n_events = new_pulse_format ? size(d_daq_pulses) : (size(d_daq_pulses, 2), size(d_daq_pulses, 1), size(d_daq_pulses, 3))
            chunksize = n_channel, chunksize_daq[3]
            g_pd = exists(h5f, "Processed_data") ? g_open(h5f, "Processed_data") : g_create(h5f, "Processed_data")
            if exists(g_pd, "baseline_offset") o_delete(g_pd, "baseline_offset") end
            if exists(g_pd, "baseline_slope") o_delete(g_pd, "baseline_slope") end
            if exists(g_pd, "baseline_rms") o_delete(g_pd, "baseline_rms") end
            d_baseline_offset = d_create(g_pd, "baseline_offset", Float32, ((n_channel, n_events),(n_channel, n_events)), "chunk", chunksize)
            d_baseline_slope  = d_create(g_pd, "baseline_slope",  Float32, ((n_channel, n_events),(n_channel, n_events)), "chunk", chunksize)
            d_baseline_rms = d_create(g_pd, "baseline_rms",    Float32, ((n_channel, n_events),(n_channel, n_events)), "chunk", chunksize)
            
            # n_events = 2000 # debugging

            # @fastmath @inbounds begin
            begin
                chunk_offsets = Array{T, 2}(undef, n_channel, chunksize[2])
                chunk_slopes = Array{T, 2}(undef, n_channel, chunksize[2])
                chunk_rms = Array{T, 2}(undef, n_channel, chunksize[2])
                chunk_pulses = Array{T, 3}(undef, n_samples, n_channel, chunksize[2])
                # @showprogress for evt_range in event_range_iterator(n_events, chunksize[2])
                for evt_range in event_range_iterator(n_events, chunksize[2])
                    @info evt_range
                    levtr::Int = length(evt_range)
                    chunk_pulses[:, :, 1:levtr]  = d_daq_pulses[:, :, evt_range]
                    for i in 1:levtr
                        for ichn in 1:n_channel
                            baseline::Vector{T} = T.( chunk_pulses[1:bl, ichn, i] )
                            slope::T = slope_length_inv * (mean(baseline[slope_w2]) - mean(baseline[slope_w1]))
                            μ::T = mean(baseline)
                            σ::T = std(baseline, mean = μ)
                            chunk_slopes[ichn, i] = slope
                            chunk_offsets[ichn, i] = μ
                            chunk_rms[ichn, i] = σ
                        end
                    end
                    d_baseline_offset[:, evt_range] = chunk_offsets[:, 1:levtr]
                    d_baseline_slope[:, evt_range]  = chunk_slopes[:, 1:levtr]
                    d_baseline_rms[:, evt_range]    = chunk_rms[:, 1:levtr]
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

