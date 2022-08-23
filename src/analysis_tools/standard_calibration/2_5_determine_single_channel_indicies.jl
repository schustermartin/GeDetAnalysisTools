function determine_single_channel_indices( m, c; ΔE::Real = 3.0, ΔE_high_energy::Real  = 23.0 )
    n_total_events = get_number_of_events(m)
    inputfiles = gather_absolute_paths_to_hdf5_input_files(m)
    core::UInt8 = 1
    for (fi, f) in enumerate(inputfiles)
        h5f = h5open(f, "r+")
        try
            g_pd  = g_open(h5f, "Processed_data")
            d_energies = d_open(g_pd, "energies")
            T::Type = eltype(d_energies)
            d::T = T(ΔE)
            d_high_energy::T = T(ΔE_high_energy)
            n_channel, n_events = size(d_energies)
            chunk_n_events = HDF5.get_chunk(d_energies)[end]
            # n_events = 3000 # debugging
            if exists(g_pd, "single_segment_indices")
                o_delete(g_pd, "single_segment_indices")
            end
            d_single_segment_indices  = d_create(g_pd, "single_segment_indices", UInt8, ((n_events,),(n_events,)), "chunk", (chunk_n_events,) )

            if n_channel > 1
                evt_ranges = event_range_iterator(n_events, chunk_n_events)
                @fastmath @inbounds begin
                    @showprogress for evt_range in evt_ranges
                        chunk_energies::Array{T, 2} = d_energies[:, evt_range]
                        chunk_ssi::Array{UInt8, 1} = zeros(UInt8, length(evt_range))
                        for event in 1:length(evt_range)
                            chunk_ssi[event] = get_single_segment_channel_index_abs(chunk_energies[:, event], d, d_high_energy)
                        end
                        d_single_segment_indices[evt_range] = chunk_ssi
                    end
                end
            else
                d_single_segment_indices[:] = 0x01
            end

            close(h5f)
        catch err
            close(h5f)
            error(err)
        end
    end
    return nothing
end

function determine_single_channel_indices( m; ΔE::Real = 3.0 )
    c = read_analysis_result_dataset(m, "calibration_matrix")
    determine_single_channel_indices(m, c, ΔE=ΔE)
end
