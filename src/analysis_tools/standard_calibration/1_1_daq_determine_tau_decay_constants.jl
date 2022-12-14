
function daq_determine_decay_time_constants(m; photon_lines = [609.312, 911.204, 1120.287, 1460.830, 1764.494, 2614.533], energy_range=200:3000, create_plots=true, ssidcs_ΔE = 10.0, nbins = 10000)
    n_total_events = get_number_of_events(m)
    n_channel = get_number_of_channel(m)
    multi_channel_det::Bool = n_channel > 1
    inputfiles = gather_absolute_paths_to_hdf5_input_files(m)
    core::UInt8 = 1
    fitted_tau_decay_constants::Array{Float64, 1} = Float64[ 0 for chn in 1:n_channel ]
    hists = Histogram[ Histogram(20:0.1:80, :left) for chn in 1:n_channel ]
    d_energies = get_quick_calibrated_daq_energies(m, photon_lines = photon_lines, nbins = nbins)
    total_evt_start_idx::Int = 0

    for (fi, f) in enumerate(inputfiles)
        h5f = h5open(f, "r+")
        try
            g_pd  = g_open(h5f, "Processed_data")
            d_tau_decay_constants = d_open(g_pd, "tau_decay_constants")

            d_energies = get_quick_calibrated_daq_energies(m, photon_lines = photon_lines, nbins = nbins)

            T::Type = eltype(d_energies)
            energy_range = T.(energy_range)
            n_events_in_file = size(d_tau_decay_constants)[end]
            chunk_n_events = HDF5.get_chunk(d_tau_decay_constants)[end]
            # n_events = 2000
            evt_ranges = event_range_iterator(n_events_in_file, chunk_n_events)
            healthy_event_indices::Vector{Int64} = get_event_indices(GAT.get_event_flags(m), :healthy)
            @fastmath @inbounds begin
                @showprogress for evt_range in evt_ranges
                    chunk_energies::Array{T, 2} = d_energies[:, total_evt_start_idx .+ evt_range]
                    chunk_tdcs::Array{T, 2} = d_tau_decay_constants[:, evt_range]
                    for event in 1:length(evt_range)
                        abs_evt_idx = evt_range[1] - 1 + event
                        if abs_evt_idx in healthy_event_indices
                            if multi_channel_det
                                ssi::UInt8 = get_single_segment_channel_index_abs(d_energies[:,event], T(ssidcs_ΔE))
                                if ssi > 0
                                    core_energy::T = chunk_energies[core, event]
                                    if first(energy_range) <= core_energy <= last(energy_range)
                                        push!(hists[core], chunk_tdcs[core, event])
                                        push!(hists[ssi], chunk_tdcs[ssi, event])
                                    end
                                end
                            else
                                core_energy = chunk_energies[core, event]
                                if first(energy_range) <= core_energy <= last(energy_range)
                                    push!(hists[core], chunk_tdcs[core, event])
                                end
                            end
                        end
                    end
                end
            end
            total_evt_start_idx += n_events_in_file
            close(h5f)
        catch err
            close(h5f)
            error(err)
        end
    end

    TF::DataType = Float64
    # fit_results = GeDetSpectrumAnalyserTmp.Fit[]
    fit_results = RadiationSpectra_beforeBAT.FitFunction[]
    for ichn in eachindex(1:n_channel)
        fitf = RadiationSpectra_beforeBAT.FitFunction{TF}( scaled_cauchy, 1, 3 )
        try
            h = hists[ichn]
            mp = midpoints(h.edges[1])
            μ::Float64 = mean(h)
            idx_μ = StatsBase.binindex(h, μ)
            # fitrange = (μ - 10 * step(h)):(μ + 10 * step(h))
            fitrange = ((μ - 10 * step(h)),(μ + 10 * step(h)))
            σ::Float64 = 0.5 * std(h, mean=μ)
            A::Float64 = sum(h.weights)
            p0::Array{Float64, 1} = [ A, σ, μ ]
            set_fitranges!(fitf, (((μ - 10 * step(h)),(μ + 10 * step(h))), ))
            set_initial_parameters!(fitf, p0)
            RadiationSpectra_beforeBAT.lsqfit!(fitf, h)
            push!(fit_results, fitf)
            # push!(fit_results, fr)
        catch err
            h = hists[ichn]
            set_fitranges!(fitf, ((35, 65), ))

            set_initial_parameters!(fitf, TF[findmax(h.weights)[1],1.5, collect(h.edges[1])[findmax(h.weights)[2]]])
            # println(T[findmax(h.weights)[1],1.5, collect(h.edges[1])[findmax(h.weights)[2]]])

            RadiationSpectra_beforeBAT.lsqfit!(fitf, h)
            push!(fit_results, fitf)
        end
    end
    if create_plots
        det = m.detector
        p = histogramdisplay(hists, det)
        for (i, fr) in enumerate(fit_results)
            plotrange = fr.fitted_parameters[3] - 4 * fr.fitted_parameters[2], fr.fitted_parameters[3] + 4 * fr.fitted_parameters[2]
            plot!(p, fr, subplot=det.channel_display_order[i])#, bin_width = StatsBase.binvolume(hists[1], 1))#, xlims=plotrange)
        end
        # savefig(m, p, "2_5_tau_decay_constants", "tau_decay_constants", fmt=:png); p = 0;
    end
    tdcs::Array{Float64, 1} = Float64.([ fr.fitted_parameters[3] for fr in fit_results ])

    return tdcs, hists, fit_results
end
