function determine_AoverE( m::Measurement;
                            waveform_filters::Union{Missing, Vector{<:Tuple{<:Interval, Function}}} = missing,
                            overwrite::Bool = false, debug::Bool = false)
    T::DataType = Float32
    some_filters::Bool = !ismissing(waveform_filters)

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

    n_filtered_events::Int = 0
    @inbounds for fn in inputfiles
        h5f::HDF5File = h5open(fn, "r+")
        try
            g_daq::HDF5Group = g_open(h5f, "DAQ_Data")
            g_pd::HDF5Group = g_open(h5f, "Processed_data")
            d_daq_pulses::HDF5Dataset = d_open(g_daq, "daq_pulses")
            d_energies::HDF5Dataset = d_open(g_pd, "energies")

            chunk_size::NTuple{3, Int} = get_chunk(d_daq_pulses)
            n_samples::Int, _n_channel::Int, n_events::Int = size(d_daq_pulses)
            if debug n_events = chunk_size[3] * 2 end
            @info "$n_events events"
            evt_ranges::Vector{UnitRange{Int}} = event_range_iterator(n_events, chunk_size[3])

            if exists(g_pd, "AoverE") o_delete(g_pd, "AoverE") end
            ds_AoverE::HDF5Dataset = d_create(g_pd, "AoverE", Float32, ((n_channel, n_events), (n_channel, n_events)), "chunk", (n_channel, chunk_size[3]))

            @showprogress for evt_range in evt_ranges
                chunk_energies::Array{T, 2} = d_energies[:, evt_range]
                chunk_pulses::Array{T, 3} = d_daq_pulses[: , :, evt_range]
                chunk_AoverE::Array{T, 2} = Array{T, 2}(undef, n_channel, length(evt_range) )
                for i in 1:length(evt_range)
                    evt_energies::Vector{T} = chunk_energies[:, i]
                    waveforms::Array{T, 2} = GeDetPulseShapeAnalysisToolsTmp.baseline_substraction_and_decay_correction(chunk_pulses[:, :, i], bl, bl_inv, decay_factors);
                    GeDetPulseShapeAnalysisToolsTmp.calibrate_pulses!(waveforms, c)
                    filter_function = begin
                        f = identity
                        if some_filters
                            for ifilter in eachindex(waveform_filters)
                                if evt_energies[1] in waveform_filters[ifilter][1]
                                    f = waveform_filters[ifilter][2]
                                    n_filtered_events += 1
                                end
                            end
                        end
                        f 
                    end
                    for ichn in 1:n_channel
                        processed_waveform::Vector{T} = filter_function(waveforms[:, ichn])
                        I::Vector{T} = diff(processed_waveform)
                        A::T = sum(diff(processed_waveform))
                        E::T = maximum(I)
                        chunk_AoverE[ichn, i] = A / E
                    end
                end
                ds_AoverE[:, evt_range] = chunk_AoverE
            end
        catch err
            close(h5f)
            error(err)
        end
        close(h5f)
    end
    @show n_filtered_events
    nothing
end