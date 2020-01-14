function determine_decay_time_constants(m; energy_range=200:3000, create_plots=true)
    n_total_events = get_number_of_events(m)
    n_channel = get_number_of_channel(m)
    inputfiles = gather_absolute_paths_to_hdf5_input_files(m)
    core::UInt8 = 1
    fitted_tau_decay_constants::Array{Float64, 1} = Float64[ 0 for chn in 1:n_channel ]
    hists = Histogram[ Histogram(10:0.5:150, :left) for chn in 1:n_channel ]

    for (fi, f) in enumerate(inputfiles)
        h5f = h5open(f, "r+")
        try
            g_pd  = g_open(h5f, "Processed_data")
            d_energies = d_open(g_pd, "energies")
            d_tau_decay_constants = d_open(g_pd, "tau_decay_constants")
            d_ssidcs = d_open(g_pd, "single_segment_indices")
            d_pile_up_flag = d_open(g_pd, "event_flags")
            T::Type = eltype(d_energies)
            energy_range = T.(energy_range)
            n_channel, n_events = size(d_energies)
            chunk_n_events = get_chunk(d_energies)[end]

            # n_events = 1000
            evt_ranges = event_range_iterator(n_events, chunk_n_events)
            @fastmath @inbounds begin
                # @showprogress for evt_range in evt_ranges
                for evt_range in evt_ranges
                    chunk_energies::Array{T, 2} = d_energies[:, evt_range]
                    chunk_ssi::Array{UInt8, 1} = d_ssidcs[evt_range]
                    chunk_tdcs::Array{T, 2} = d_tau_decay_constants[:, evt_range]
                    for event in get_event_indices(d_pile_up_flag[evt_range], :healthy)#1:length(evt_range)
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

    # for i in eachindex(hists)
    #     hists[i] = normalize!(float(hists[i])) 
    # end

    T = Float64
    # fit_results = GeDetSpectrumAnalyserTmp.Fit[]
    fit_results = RadiationSpectra.FitFunction[]
    for ichn in eachindex(1:n_channel)
        fitf = RadiationSpectra.FitFunction{T}( scaled_cauchy, 1, 3 )
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
            set_fitranges!(fitf, (fitrange, ))
            set_initial_parameters!(fitf, p0)
            RadiationSpectra.lsqfit!(fitf, h)
            # fr = GeDetSpectrumAnalyserTmp.fit(h, fitrange, scaled_cauchy, p0)
            push!(fit_results, fitf)
            # push!(fit_results, fr)
        catch err
            h = hists[ichn]
            set_fitranges!(fitf, ((35, 65), ))
            set_initial_parameters!(fitf, [maximum(h.weights), 1, 50])
            RadiationSpectra.lsqfit!(fitf, h)
            push!(fit_results, fitf)
        end
    end
    if create_plots
        det = m.detector
        p = histogramdisplay(hists, det)
        for (i, fr) in enumerate(fit_results)
            plotrange = fr.fitted_parameters[3] - 4 * fr.fitted_parameters[2], fr.fitted_parameters[3] + 4 * fr.fitted_parameters[2]
            bw = StatsBase.binvolume(hists[i], 1)
            bw = 1.0 # in case of lsqfit! take always 1
            plot!(p, fr, subplot=det.channel_display_order[i], bin_width = bw )#, xlims=plotrange)
        end
        savefig(m, p, "2_6_tau_decay_constants", "tau_decay_constants", fmt=:png); p = 0;
    end
    tdcs::Array{Float64, 1} = [ fr.fitted_parameters[3] for fr in fit_results ]

    if create_plots
        p = plot(tdcs, st=:scatter, title="Tau Decay Constants", ylabel="τ / μs", xlabel="Channel Index", size=(1200,600), label="τ's")
        tdcs_init = read_analysis_result_dataset(m, "init_tau_decay_constants")
        tdcs_init_err = read_analysis_result_dataset(m, "init_tau_decay_constants_err")
        for i in eachindex(tdcs_init_err) if tdcs_init_err[i] < 0 tdcs_init_err[i] = 0 end end
        plot!(tdcs_init, label="Init τ's", st=:scatter)
        savefig(m, p, "2_6_tau_decay_constants", "scatter_tau_decay_constants", fmt=:png); p = 0;
    end


    return tdcs, hists, fit_results
end
