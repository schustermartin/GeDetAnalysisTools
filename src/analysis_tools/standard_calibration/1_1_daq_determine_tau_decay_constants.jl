
function daq_determine_decay_time_constants(m; photon_lines = [609.312, 911.204, 1120.287, 1460.830, 1764.494, 2614.533], energy_range=200:3000, create_plots=true, ssidcs_ΔE = 20.0)
    n_total_events = get_number_of_events(m)
    n_channel = get_number_of_channel(m)
    inputfiles = gather_absolute_paths_to_hdf5_input_files(m)
    core::UInt8 = 1
    fitted_tau_decay_constants::Array{Float64, 1} = Float64[ 0 for chn in 1:n_channel ]
    hists = Histogram[ Histogram(20:0.5:80, :left) for chn in 1:n_channel ]

    for (fi, f) in enumerate(inputfiles)
        h5f = h5open(f, "r+")
        try
            g_pd  = g_open(h5f, "Processed_data")
            d_tau_decay_constants = d_open(g_pd, "tau_decay_constants")
            d_energies = get_quick_calibrated_daq_energies(m, photon_lines = photon_lines)
            # d_ssidcs = d_open(g_pd, "single_segment_indices")
            T::Type = eltype(d_energies)
            energy_range = T.(energy_range)
            n_channel, n_events = size(d_energies)
            chunk_n_events = 1000#get_chunk(d_energies)[end]

            evt_ranges = event_range_iterator(n_events, chunk_n_events)
            @fastmath @inbounds begin
                # @showprogress for evt_range in evt_ranges
                for evt_range in evt_ranges
                    chunk_energies::Array{T, 2} = d_energies[:, evt_range]
                    # chunk_ssi::Array{UInt8, 1} = d_ssidcs[evt_range]
                    chunk_tdcs::Array{T, 2} = d_tau_decay_constants[:, evt_range]
                    for event in 1:length(evt_range)
                        ssi::UInt8 = get_single_segment_channel_index_abs(d_energies[:,event],T(ssidcs_ΔE))
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

    # fit_results = GeDetSpectrumAnalyserTmp.Fit[]
    fit_results = RadiationSpectra.FitFunction[]
    for ichn in eachindex(1:n_channel)
        fitf = RadiationSpectra.FitFunction( scaled_cauchy )
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
            fitf.fitrange = fitrange
            fitf.initial_parameters = p0
            RadiationSpectra.lsqfit!(fitf, h)
            # fr = GeDetSpectrumAnalyserTmp.fit(h, fitrange, scaled_cauchy, p0)
            push!(fit_results, fitf)
            # push!(fit_results, fr)
        catch err
            h = hists[ichn]
            fitf.fitrange = (35, 65)
            @show maximum(h.weights)
            fitf.initial_parameters = Float64[findmax(h.weights)[1],1.5, collect(h.edges[1])[findmax(h.weights)[2]]]
            try
                RadiationSpectra.lsqfit!(fitf, h, estimate_uncertainties = true)
            catch err
                RadiationSpectra.lsqfit!(fitf, h, estimate_uncertainties = false)
                fitf.uncertainties = eltype(fitf.uncertainties)[-1, -1, -1, -1]
            end
            # fr = GeDetSpectrumAnalyserTmp.Fit(scaled_cauchy, 0.0:100.0, Float64[1,1,0], Float64[-1,-1,-1], Float64[1,1,0], missing)
            push!(fit_results, fitf)
        end
    end
    if create_plots
        det = m.detector
        p = histogramdisplay(hists, det)
        for (i, fr) in enumerate(fit_results)
            plotrange = fr.parameters[3] - 4 * fr.parameters[2], fr.parameters[3] + 4 * fr.parameters[2]
            plot!(p, fr, subplot=det.channel_display_order[i])#, xlims=plotrange)
        end
        # savefig(m, p, "2_5_tau_decay_constants", "tau_decay_constants", fmt=:png); p = 0;
    end
    tdcs::Array{Float64, 1} = [ fr.parameters[3] for fr in fit_results ]
    # for fr in fit_results GeDetSpectrumAnalyserTmp.estimate_uncertainties!(fr) end
    tdcs_err::Array{Float64, 1} = [ fr.uncertainties[3] for fr in fit_results ]

    # if create_plots
    #     p = plot(tdcs, yerr=tdcs_err, st=:scatter, title="Tau Decay Constants", ylabel="τ / μs", xlabel="Channel Index", size=(1200,600), label="τ's")
    #     tdcs_init = read_analysis_result_dataset(m, "init_tau_decay_constants")
    #     tdcs_init_err = read_analysis_result_dataset(m, "init_tau_decay_constants_err")
    #     for i in eachindex(tdcs_init_err) if tdcs_init_err[i] < 0 tdcs_init_err[i] = 0 end end
    #     plot!(tdcs_init, yerr=tdcs_init_err, label="Init τ's", st=:scatter)
    #     # savefig(m, p, "2_5_tau_decay_constants", "scatter_tau_decay_constants", fmt=:png); p = 0;
    # end


    return tdcs, tdcs_err, hists, fit_results
end



# function determine_decay_time_constants_distributions_with_daq_ss_selection(  m, cm_daq;
#                                                         energy_limits=[500, 1500],
#                                                         take_every_n_sample_for_fit = 1 )
#
#     T = Float32
#     energy_limits = T.(energy_limits)
#     n_channel = get_number_of_channel(m)
#     n_segments = n_channel - 1
#     new_pulse_format = is_new_pulse_format(m)
#     n_total_events = get_number_of_events(m)
#     single_segment_abs_diff::T = 3
#
#     decay_time_constants_distributions = [Histogram(30:0.5:70, :left) for ichn in 1:n_channel] # in microseconds
#     decay_time_constants_vs_energy = [Histogram( ((energy_limits[1]:1:energy_limits[2]), 30:0.5:70), :left) for ichn in 1:n_channel] # in microseconds
#
#     inputfiles = gather_absolute_paths_to_hdf5_input_files(m)
#
#     time_array_of_pulses = get_time_array_of_pulses(m.daq)
#     decay_window_start_index = m.daq.n_samples-m.daq.decay_window_length+1
#     xarr = time_array_of_pulses[decay_window_start_index:take_every_n_sample_for_fit:end]
#     dataarr::Array{Float64, 1} = zeros(Float64, length(xarr)) # must be Float64 due LsqFit...
#     init_params::Array{Float64, 1} = zeros(Float64, 2) # must be Float64 due LsqFit...
#
#     init_decay_rates = Float64.(50.0 / 1e6)
#
#     tdc::T = 0
#     core::UInt8 = 1
#
#     fitf = RadiationSpectra.FitFunction(  exponential_decay )
#
#     @info "Tau Decay Fitting: $n_total_events events to process"
#
#     for (fi, f) in enumerate(inputfiles)
#         h5f = h5open(f, "r")
#         try
#             g = g_open(h5f, "DAQ_Data")
#             d_daq_pulses = d_open(g, "daq_pulses")
#             daq_energies = calibrate_energies(T.(read(d_open(g, "daq_energies"))), cm_daq)
#             n_events = size(d_daq_pulses, 3)
#             chunksize = n_channel, chunksize_daq[3]
#             # n_events = 2000 # debugging
#
#             evt_ranges = event_range_iterator(n_events, chunk_n_events)
#             if new_pulse_format
#                 @fastmath @inbounds begin
#                     @showprogress for evt_range in evt_ranges
#                         chunk_pulses::Array{T, 3} = d_daq_pulses[: , :, evt_range]
#                         for event in 1:length(evt_range)
#                             ssidx = get_single_segment_channel_index_abs(daq_energies[:, first(evt_range) + event - 1], single_segment_abs_diff)
#                             if ssidx > 1
#                                 daq_pulses_evt = T.(chunk_pulses[:,:,event])
#                                 core_energy::T = daq_energies[core, event]
#                                 if energy_limits[1] <= core_energy <= energy_limits[2] # ~keV
#                                     for ichn in UInt8[core, ssidx]
#                                         baseline::T = mean(T, daq_pulses_evt[1:m.daq.baseline_length, ichn])
#                                         dataarr[:] = Float64.(daq_pulses_evt[decay_window_start_index:take_every_n_sample_for_fit:end, ichn] .- baseline) # must be Float64 due LsqFit...
#                                         init_params[:] = [dataarr[1], init_decay_rates[ichn]] # must be Float64 due LsqFit...
#                                         fitf.initial_parameters = init_params
#                                         RadiationSpectra.lsqfit!(fitf, xarr, dataarr)
#                                         tdc = 1f6 * fitf.parameters[2]
#                                         # tdc = 1f6 * T( GeDetSpectrumAnalyserTmp.LSQFIT( xarr, dataarr, exponential_decay, init_params, estimate_uncertainties=false).parameters[2] )
#                                         # tdc = 1e6 * LsqFit.curve_fit(exponential_decay, xarr, dataarr, init_params).param[2]
#                                         push!(decay_time_constants_distributions[ichn], tdc )
#                                         push!(decay_time_constants_vs_energy[ichn], (core_energy, tdc))
#                                     end
#                                 end
#                             end
#                         end
#                     end
#                 end
#             else # old pulse format
#                 @warn "old pulse format not yet tested"
#                 @fastmath @inbounds begin
#                     @showprogress for evt_range in evt_ranges
#                         chunk_pulses::Array{T, 3} = d_daq_pulses[: , :, evt_range]
#                         for event in 1:length(evt_range)
#                             ssidx = get_single_segment_channel_index_abs(daq_energies[:, first(evt_range) + event - 1], single_segment_abs_diff)
#                             if ssidx > 1
#                                 daq_pulses_evt = T.(transpose(chunk_pulses[:,:,event]))
#                                 core_energy::T = daq_energies[core, event]
#                                 if energy_limits[1] <= core_energy <= energy_limits[2] # ~keV
#                                     for ichn in UInt8[core, ssidx]
#                                         baseline::T = mean(T, daq_pulses_evt[1:m.daq.baseline_length, ichn])
#                                         dataarr[:] = Float64.(daq_pulses_evt[decay_window_start_index:take_every_n_sample_for_fit:end, ichn] .- baseline) # must be Float64 due LsqFit...
#                                         init_params[:] = [dataarr[1], init_decay_rates[ichn]] # must be Float64 due LsqFit...
#                                         fitf.initial_parameters = init_params
#                                         RadiationSpectra.lsqfit!(fitf, xarr, dataarr)
#                                         tdc = 1f6 * fitf.parameters[2]
#                                         # tdc = 1f6 * T( GeDetSpectrumAnalyserTmp.LSQFIT( xarr, dataarr, exponential_decay, init_params, estimate_uncertainties=false).parameters[2] )
#                                         # tdc = 1e6 * LsqFit.curve_fit(exponential_decay, xarr, dataarr, init_params).param[2]
#                                         push!(decay_time_constants_distributions[ichn], tdc )
#                                         push!(decay_time_constants_vs_energy[ichn], (core_energy, tdc))
#                                     end
#                                 end
#                             end
#                         end
#                     end
#                 end
#             end
#             close(h5f)
#         catch err
#             close(h5f)
#             error(err)
#         end
#     end
#     return decay_time_constants_distributions, decay_time_constants_vs_energy
# end
