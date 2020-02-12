@fastmath function determine_risetime(data::Array{T, 1}, threshold_min::T, threshold_max::T, pulse_energy::T, Δt::T)::Tuple{T, T, T} where {T <: Real}
	@inbounds begin
		tmin = threshold_min * pulse_energy
		tmax = threshold_max * pulse_energy
		i_middle_point = 0
		i = 1
		while data[i] < 0.5 * (threshold_max + threshold_min) * pulse_energy
			i += 1
		end

		i_middle_point = i

		i_rt_start_low = 0
		i_rt_start_up  = 0
		i_rt_stop_low  = 0
		i_rt_stop_up   = 0
		i = i_middle_point
		while data[i] > tmin
			i -= 1
		end
		i_rt_start_up  = i + 1
		i_rt_start_low = i

		#y=mx+t
		m_start = data[i_rt_start_up] - data[i_rt_start_low]
		t_start = data[i_rt_start_low]
		x_start = (tmin - t_start) / m_start
		rt_start = (i_rt_start_low - 1 + x_start) * Δt # Δt = inv(sampling_rate)

		i = i_middle_point
		while data[i] < tmax
			i += 1
		end
		i_rt_stop_up  = i
		i_rt_stop_low = i - 1
		m_stop = data[i_rt_stop_up] - data[i_rt_stop_low]
		t_stop = data[i_rt_stop_low]
		x_stop = (tmax - t_stop) / m_stop
		rt_stop = (i_rt_stop_low - 1 + x_stop) * Δt
		
		return rt_stop - rt_start, rt_start, rt_stop
	end
end


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
                if (exists(g_pd, ds_name)&&overwrite) o_delete(g_pd, ds_name) end
            end
            filter!(x-> !exists(g_pd, "rise_time_$(x[1])-$(x[2])"), rise_time_windows)
            n_risetimes = length(rise_time_windows)
            @show(n_risetimes, rise_time_windows)
            if n_risetimes > 0
                for irs in 1:n_risetimes
                    ds_name = "rise_time_$(rise_time_windows[irs][1])-$(rise_time_windows[irs][2])"
                    ds::HDF5Dataset = d_create(g_pd, ds_name, Float32, ((n_rs_values, n_channel, n_events), (n_rs_values, n_channel, n_events)), "chunk", (n_rs_values, n_channel, chunk_size[3]))
                    push!(ds_risetimes, ds)
                end

                @showprogress for evt_range in evt_ranges
                    chunk_pulses::Array{T, 3} = d_daq_pulses[: , :, evt_range]
                    chunk_energies::Array{T, 2} = d_energies[:, evt_range]
                    chunk_rise_times::Array{T, 4} = Array{T, 4}(undef, n_rs_values, n_channel, length(evt_range), n_risetimes )
                    for i in 1:length(evt_range)
                        evt_energies::Vector{T} = chunk_energies[:, i]
                        waveforms::Array{T, 2} = baseline_substraction_and_decay_correction(chunk_pulses[:, :, i], bl, bl_inv, decay_factors);
                        calibrate_pulses!(waveforms, c)

                        for ichn in 1:n_channel
                            processed_waveform::Vector{T} = waveform_filter(waveforms[:, ichn])

                            for irs in 1:n_risetimes
                                rs::T, rs_start::T, rs_stop::T = determine_risetime(processed_waveform,
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
            else
                @info("nothing to do")
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

function determine_pulse_level_time_stamps(m::Measurement, pulse_levels::Vector{<:Real}; 
                                waveform_filters::Union{Missing, Vector{<:Tuple{<:Interval, Function}}} = missing,
                                overwrite::Bool = false, debug::Bool = false)
    T::DataType = Float32
    n_pulse_levels::Int = length(pulse_levels)
    pulse_levels = T.(pulse_levels)

    n_total_events::Int = get_number_of_events(m)
	inputfiles::Vector{AbstractString} = gather_absolute_paths_to_hdf5_input_files(m)
    tau_decay_constants::Vector{T} = T.(read_analysis_result_dataset(inputfiles[1], "init_tau_decay_constants"))
    n_channel::Int = length(tau_decay_constants)
    c::Array{T, 2} = T.(read_analysis_result_dataset(m, "calibration_matrix"))
    
    bl::Int = m.daq.baseline_length
    bl_inv::T = T(1 / bl)
    sampling_rate::T = T(m.daq.sampling_rate)
    Δt::T = T(inv(sampling_rate))
    decay_factors::Vector{T} = zeros(T, n_channel)
    for chn in 1:n_channel
        decay_factors[chn] = exp( - (1 / sampling_rate) / (tau_decay_constants[chn] * T(1e-6) ) )
    end

    some_filters::Bool = !ismissing(waveform_filters)

    n_filtered_events::Int = 0
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
            for ipl in 1:n_pulse_levels
                ds_name = "pulse_level_$(pulse_levels[ipl])_time_stamps"
                if exists(g_pd, ds_name) o_delete(g_pd, ds_name) end
                ds::HDF5Dataset = d_create(g_pd, ds_name, Float32, ((n_channel, n_events), (n_channel, n_events)), "chunk", (n_channel, chunk_size[3]))
                push!(ds_risetimes, ds)
            end


            @showprogress for evt_range in evt_ranges
                chunk_pulses::Array{T, 3} = d_daq_pulses[: , :, evt_range]
                chunk_energies::Array{T, 2} = d_energies[:, evt_range]
                chunk_rise_times::Array{T, 3} = Array{T, 3}(undef, n_channel, length(evt_range), n_pulse_levels )
                for i in 1:length(evt_range)
                    evt_energies::Vector{T} = chunk_energies[:, i]
                    waveforms::Array{T, 2} = baseline_substraction_and_decay_correction(chunk_pulses[:, :, i], bl, bl_inv, decay_factors); 
                    calibrate_pulses!(waveforms, c)
                    
                    filter_function = begin
                        f = identity
                        for ifilter in eachindex(waveform_filters)
                            if evt_energies[1] in waveform_filters[ifilter][1]
                                f = waveform_filters[ifilter][2]
                                n_filtered_events += 1
                            end
                        end
                        f 
                    end
                    for ichn in 1:n_channel
                        processed_waveform::Vector{T} = filter_function(waveforms[:, ichn])
                        
                        for ipl in 1:n_pulse_levels
                            rs::T, rs_start::T, rs_stop::T = determine_risetime(processed_waveform, 
                                                                                                                pulse_levels[ipl], T(1), 
                                                                                                                evt_energies[ichn], Δt)
                            if isnan(rs_start)
                                rs_start = T(-1)# error("Got NaN in ristetime determination for Energie $(evt_energies[ichn]) and pulse level $(pulse_levels[ipl])")
                            end
                            chunk_rise_times[ichn, i, ipl] = rs_start
                        end
                        
                    end
                end
                for ipl in 1:n_pulse_levels
                    ds_risetimes[ipl][:, evt_range] = chunk_rise_times[:, :, ipl]
                end
            end
        catch err
            close(h5f)
            error(err)
        end
        close(h5f)
        # @info n_filtered_events
    end
    nothing
end

function get_pulse_level_time_stamps(m::Measurement)
    T = Float32
    ds_name_regex = r"pulse_level_[0-9]*.?[0-9]*_time_stamps"
    inputfiles::Vector{AbstractString} = gather_absolute_paths_to_hdf5_input_files(m)
    n_total_events::Int = get_number_of_events(m)
    n_total_channel::Int = get_number_of_channel(m)
    ds_names = begin
        h5f::HDF5File = h5open(inputfiles[1], "r")
        _ds_names = try
            g_pd::HDF5Group = g_open(h5f, "Processed_data")
            filter(dsn -> occursin(ds_name_regex, dsn), names(g_pd))
        catch err
            close(h5f)
            error(err)
        end
        close(h5f)
        _ds_names
    end
    pulse_levels = T[parse(T, n[(length("pulse_level_")+1):end-(length("_time_stamps"))]) for n in ds_names]
    
    time_stamps::Vector{Array{T, 2}} = [Array{T, 2}(undef, n_total_channel, n_total_events) for ts in eachindex(ds_names)]
    last_evt_idx::Int = 0
    for fn in inputfiles
        h5f = h5open(fn, "r")
        try
            g_pd = g_open(h5f, "Processed_data")
            datasets = [d_open(g_pd, dsn) for dsn in ds_names]
            n_events = size(datasets[1], 3)
            range = last_evt_idx+1:last_evt_idx+n_events
            for (ids, ds) in enumerate(datasets)
                time_stamps[ids][:] = read(ds)
            end
            last_evt_idx += n_events
        catch err
            close(h5f)
            error(err)
        end
        close(h5f)
    end
    return time_stamps, pulse_levels
end
