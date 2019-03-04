function recalculate_energies(m::Measurement, c::Array{<:Real, 2}; create_plots=true)::Nothing
    inputfiles = gather_absolute_paths_to_hdf5_input_files(m)
    T = get_eltype_of_dataset(inputfiles[1], "Processed_data", "measured_pulse_amplitudes")
    n_channel = size(c, 1)

    c_transpose::Array{T, 2} = c'

    if create_plots
        hists = Histogram[ Histogram(0:1:6000, :left) for i in eachindex(1:n_channel)]
        h_seg_sum = Histogram(0:1:6000, :left)
    end

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
            d_measured_pulse_amplitudes = d_open(g_pd, "measured_pulse_amplitudes")
            if exists(g_pd, "energies") o_delete(g_pd, "energies") end
            d_energies = d_create(g_pd, "energies", Float32, ((n_channel, n_events),(n_channel, n_events)), "chunk", chunksize )

            @fastmath @inbounds begin
                @showprogress for evt_range in event_range_iterator(n_events, chunksize[2])
                    chunk_mpas::Array{T, 2} = d_measured_pulse_amplitudes[:, evt_range]
                    chunk_energies::Array{T, 2} = Array{T, 2}(undef, n_channel, length(evt_range))
                    for i in eachindex(1:length(evt_range))
                        chunk_energies[:, i] = c_transpose * chunk_mpas[:, i]
                        if create_plots
                            push!(h_seg_sum, sum(chunk_energies[2:end, i]))
                        end
                    end
                    if create_plots
                        for ichn in eachindex(1:n_channel)
                            append!(hists[ichn], chunk_energies[ichn, :])
                        end
                    end
                    d_energies[:, evt_range] = chunk_energies
                end
            end


            close(h5f)
        catch err
            close(h5f)
            error(err)
        end
    end

    if create_plots
        p = energyhistogramdisplay(hists, m.detector, yscale=:log10, size=(2200,1100))
        savefig(m, p, "2_4_recalculated_energies", "energy_spectra", fmt=:png); p = 0;

        p = plot(hists[1], st=:step, label="Core", xlabel="E / keV", size=(1600,900), yscale=:log10)
        plot!( h_seg_sum, st=:step, yscale=:log10, ylabel="Summed Segments"  )
        plts = []
        photon_lines::Vector{Float64} = [609.312, 911.204, 1120.287, 1460.830, 1764.494, 2614.533]
        for pl in photon_lines
            p_tmp = deepcopy(p)
            xl = Float64[pl - 15, pl + 15]
            idx1 = StatsBase.binindex(hists[1], xl[1])
            idx2 = StatsBase.binindex(hists[1], xl[2])
            ymax1 = maximum(hists[1].weights[idx1:idx2])
            ymax2 = maximum(h_seg_sum.weights[idx1:idx2])
            ymax = max(ymax1, ymax2)
            p_tmp = plot(p_tmp, xlims=xl, xticks=3, legend=false, yscale=:identity, ylims=Float64[0, ymax])
            push!(plts, p_tmp)
        end
        pfull = plot(p, plts..., layout= (@layout [a{0.5h}; b c d e f g]), size=(1920,1080))
        savefig(m, pfull, "2_4_recalculated_energies", "core_and_sum_segment_energy_spectrum", fmt=:png); p = 0;
    end

    return nothing
end

function recalculate_energies(m::Measurement; create_plots=true)::Nothing
    c = read_analysis_result_dataset(m, "calibration_matrix")
    recalculate_energies(m, c; create_plots=create_plots)
end
