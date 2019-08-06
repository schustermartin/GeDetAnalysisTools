function determine_rise_times(m::Measurement, rise_time_windows::Vector{<:Tuple{<:Real, <:Real}}; 
                                waveform_filter::Function = wv -> wv,
                                overwrite::Bool = false, debug::Bool = false)
    T::DataType = Float32
    n_risetimes::Int = length(rise_time_windows)
    rise_time_windows = Vector{Tuple{T, T}}(rise_time_windows)

	n_total_events::Int = get_number_of_events(m)
	inputfiles::Vector{AbstractString} = gather_absolute_paths_to_hdf5_input_files(m)
    tau_decay_constants::Vector{T} = T.(read_analysis_result_dataset(inputfiles[1], "init_tau_decay_constants"))
    n_channel::Int = length(tau_decay_constants)
    c::Array{T, 2} = T.(read_analysis_result_dataset(m, "calibration_matrix"))
    n_rs_values::Int = 2
    
    bl::Int = m.daq.baseline_length
    bl_inv::T = T(1 / bl)
    sampling_rate::T = T(m.daq.sampling_rate)
    Δt::T = T(inv(sampling_rate))
    decay_factors::Vector{T} = zeros(T, n_channel)
    for chn in 1:n_channel
        decay_factors[chn] = exp( - (1 / sampling_rate) / (tau_decay_constants[chn] * T(1e-6) ) )
    end

    @inbounds for fn in inputfiles
        h5f::HDF5File = h5open(fn, "r+")
        try
            g_daq::HDF5Group = g_open(h5f, "DAQ_Data")
            g_pd::HDF5Group = g_open(h5f, "Processed_data")
            d_energies::HDF5Dataset = d_open(g_pd, "energies")
            d_daq_pulses::HDF5Dataset = d_open(g_daq, "daq_pulses")
            
            chunk_size::NTuple{3, Int} = get_chunk(d_daq_pulses)
            n_samples::Int, _n_channel::Int, n_events::Int = size(d_daq_pulses)
            if debug n_events = chunk_size[3] * 2 end
            @info "$n_events events"
            evt_ranges::Vector{UnitRange{Int}} = event_range_iterator(n_events, chunk_size[3])
            
            ds_risetimes::Vector{HDF5Dataset} = HDF5Dataset[]
            for irs in 1:n_risetimes
                ds_name = "rise_time_$(rise_time_windows[irs][1])-$(rise_time_windows[irs][2])"
                if exists(g_pd, ds_name) o_delete(g_pd, ds_name) end
                ds::HDF5Dataset = d_create(g_pd, ds_name, Float32, ((n_rs_values, n_channel, n_events), (n_rs_values, n_channel, n_events)), "chunk", (n_rs_values, n_channel, chunk_size[3]))
                push!(ds_risetimes, ds)
            end

            @showprogress for evt_range in evt_ranges
                chunk_pulses::Array{T, 3} = d_daq_pulses[: , :, evt_range]
                chunk_energies::Array{T, 2} = d_energies[:, evt_range]
                chunk_rise_times::Array{T, 4} = Array{T, 4}(undef, n_rs_values, n_channel, length(evt_range), n_risetimes )
                for i in 1:length(evt_range)
                    evt_energies::Vector{T} = chunk_energies[:, i]
                    waveforms::Array{T, 2} = GeDetPulseShapeAnalysisToolsTmp.baseline_substraction_and_decay_correction(chunk_pulses[:, :, i], bl, bl_inv, decay_factors); 
                    GeDetPulseShapeAnalysisToolsTmp.calibrate_pulses!(waveforms, c)
                    
                    for ichn in 1:n_channel
                        processed_waveform::Vector{T} = waveform_filter(waveforms[:, ichn])
                        
                        for irs in 1:n_risetimes
                            rs::T, rs_start::T, rs_stop::T = GeDetPulseShapeAnalysisToolsTmp.determine_risetime(processed_waveform, 
                                                                                                                rise_time_windows[irs][1], rise_time_windows[irs][2], 
                                                                                                                evt_energies[ichn], Δt)
                            chunk_rise_times[:, ichn, i, irs] = [rs, rs_start]
                        end
                        
                    end
                end
                for irs in 1:n_risetimes
                    ds_risetimes[irs][:, :, evt_range] = chunk_rise_times[:, :, :, irs]
                end
            end
        catch err
            close(h5f)
            error(err)
        end
        close(h5f)
    end
    nothing
end

function get_risetime_dataset(m::Measurement, rise_time_window::Tuple{<:Real, <:Real})::Array{<:AbstractFloat, 3} 
    T = Float32
    rs_min::T = rise_time_window[1]
    rs_max::T = rise_time_window[2]
    ds_name = "rise_time_$(rs_min)-$(rs_max)"
    inputfiles::Vector{AbstractString} = gather_absolute_paths_to_hdf5_input_files(m)
    n_total_events::Int = get_number_of_events(m)
    n_total_channel::Int = get_number_of_channel(m)
    n_rs_values::Int = 2
    rsds::Array{T, 3} = Array{T, 3}(undef, n_rs_values, n_total_channel, n_total_events)
    last_evt_idx::Int = 0
    for fn in inputfiles
        h5f::HDF5File = h5open(fn, "r")
        try
            g_pd::HDF5Group = g_open(h5f, "Processed_data")
            d_rs::HDF5Dataset = d_open(g_pd, ds_name)
            n_events::Int = size(d_rs, 3)
            range = last_evt_idx+1:last_evt_idx+n_events
            rsds[:,:,range] = read(d_rs)
            last_evt_idx += n_events
        catch err
            close(h5f)
            error(err)
        end
        close(h5f)
    end
    return rsds
end