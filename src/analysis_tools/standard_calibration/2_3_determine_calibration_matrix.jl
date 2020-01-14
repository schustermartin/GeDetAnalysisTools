function determine_calibration_matrix_with_mpas(m::Measurement, core_calibration_constant::Real; create_plots = true, ct_peak_min_shift = 0.25)
    inputfiles = gather_absolute_paths_to_hdf5_input_files(m)
    n_channel = get_number_of_channel(m)
    n_segments = n_channel - 1
    core = 1

    T = Float32

    c0 = T(core_calibration_constant)

    c = zeros(Float64, n_channel, n_channel)
    for i in 1:n_channel c[i,i] = 1 end
    c[1, 1] = inv(c0)

    mpas::Array{T, 2} = transpose(get_measured_pulse_amplitudes(m))
    n_events = size(mpas, 1)

    ratios = Array{T, 2}(undef, n_events, n_segments )

    mpas *= c0

    for i in eachindex(1:n_events)
        if mpas[i, core] == 0 mpas[i, core] = 1 end
    end

    core_mpa_inv::Array{T, 1} = inv.(mpas[:, core])

    for iseg in 1:n_segments
        ichn = iseg + 1
        ratios[:, iseg] = mpas[:, ichn] .* core_mpa_inv[:]
    end

    ratio_hists = [ Histogram(-1.0:0.01:2, :left ) for iseg in eachindex(1:n_segments)]
    for iseg in eachindex(1:n_segments)
        append!(ratio_hists[iseg], ratios[:, iseg])
    end

    cal_peak_rough_pos = zeros(Float64, n_segments)
    plts = []
    for iseg in 1:n_channel-1
        h = ratio_hists[iseg]
        mp = (h.edges[1][1:length(h.edges[1])-1] .+ 0.5 * step(h.edges[1]))
        ct_peak_idx = findmax(h.weights[1:StatsBase.binindex(h, 0.2)])[2]
        ct_peak_pos = mp[ct_peak_idx]
        start_idx = StatsBase.binindex(h, ct_peak_pos + ct_peak_min_shift)
        cal_peak_idx = findmax(h.weights[start_idx:end])[2] - 1
        # if iseg == 4 ## <- modify for distorted pulses
        #     start_idx = StatsBase.binindex(h, ct_peak_pos + 0.6)
        #     cal_peak_idx = findmax(h.weights[start_idx:end])[2] - 1
        # end
        cal_peak_rough_pos[iseg] = mp[start_idx + cal_peak_idx]
        if create_plots
            p = plot(h, yscale=:log10, st=:step, label="", title="mpa seg$(iseg) / mpa core", xlims=[-0.1, cal_peak_rough_pos[iseg] > 1.1 ? cal_peak_rough_pos[iseg] * 1.1 : 1.1], legend=false)
            plot!(p, [cal_peak_rough_pos[iseg]],st=:vline, label="cal peak pos")
            push!(plts, p)
        end
    end
    if create_plots
        p = plot(plts..., size=(1920, 1080))
        savefig(m, p, "2_3_calibration_matrix", "1_mpa_seg_over_core_ratio_plots", fmt=:png ); p=0;
    end

    cal_ratio_hists = [ fit(Histogram, ratios[:, iseg], cal_peak_rough_pos[iseg] - 0.05:0.001:cal_peak_rough_pos[iseg] + 0.05, closed=:left ) for iseg in 1:size(ratios, 2)]
    cal_peak_pos = zeros(Float64, n_segments)
    # cal_peak_fits = GeDetSpectrumAnalyserTmp.Fit[]
    cal_peak_fits = RadiationSpectra.FitFunction[]

    plts = []
    for iseg in eachindex(1:n_channel-1)
        h = cal_ratio_hists[iseg]
        mp = (h.edges[1][1:length(h.edges[1])-1] .+ 0.5 * step(h.edges[1]))
        cal_peak_idx = findmax(h.weights)[2]
        init_fit_params = [ h.weights[cal_peak_idx] * 2π * step(h.edges[1]), 2 * step(h.edges[1]), mp[cal_peak_idx] ]
        fitrange = ((mp[cal_peak_idx] - 6 * step(h.edges[1])), (mp[cal_peak_idx] + 7.5 * step(h.edges[1])))
        fitf = RadiationSpectra.FitFunction{T}( scaled_cauchy, 1, 3 )
        set_initial_parameters!(fitf, init_fit_params)
        set_fitranges!(fitf, (fitrange,))
        RadiationSpectra.lsqfit!(fitf, h)
        push!(cal_peak_fits, fitf)
        # push!(cal_peak_fits, fr)
        cal_peak_pos[iseg] = fitf.fitted_parameters[3]
        # cal_peak_pos[iseg] = fr.parameters[3]
        if create_plots
            p = plot(h, st=:step, label="", title="seg$iseg")
            bw = StatsBase.binvolume(h, 1)
            plot!(p, fitf, label="", bin_width = bw)
            push!(plts, p)
        end
    end
    if create_plots
        p = plot(plts..., size=(1920, 1080))
        savefig(m, p, "2_3_calibration_matrix", "2_mpa_radio_cal_peaks_fits", fmt=:png ); p=0;
    end

    for iseg in eachindex(1:n_segments)
        ichn = iseg + 1
        c[ichn, ichn] = cal_peak_pos[iseg] * c[1, 1]
    end

    # Crosstralk Values
    ratios = transpose(ratios)
    ct_ratio_hists = Array{Array{<:Histogram, 1}}([ [Histogram(-0.3:0.001:0.3, :left) for i in eachindex(1:n_segments)] for j in eachindex(1:n_segments) ]) ## <- modify for distorted pulses
    for ichn_src in 2:n_channel
        iseg_src = ichn_src - 1
        ss_ratio_limits::Tuple{T, T} = cal_peak_pos[iseg_src] - 0.03, cal_peak_pos[iseg_src] + 0.03 # ss: single segment
        for i in eachindex(1:n_events) # n_events
            if ss_ratio_limits[1] <= ratios[iseg_src, i] <= ss_ratio_limits[2] # if single segment event
                for ichn_tar in 2:n_channel
                    iseg_tar = ichn_tar - 1
                    if ichn_tar == ichn_src continue end
                    push!(ct_ratio_hists[iseg_src][iseg_tar], ratios[iseg_tar, i])
                end
            end
        end
    end
    # ct_peak_fits = Array{Array{GeDetSpectrumAnalyserTmp.Fit, 1}}([])
    ct_peak_fits = Array{Array{RadiationSpectra.FitFunction, 1}}([])
    for ichn_src in 2:n_channel
        iseg_src = ichn_src - 1
        tmpfits = Array{RadiationSpectra.FitFunction, 1}([])
        plts = []
        for ichn_tar in 2:n_channel
            iseg_tar = ichn_tar - 1
            h = ct_ratio_hists[iseg_src][iseg_tar]
            mp = (h.edges[1][1:length(h.edges[1])-1] .+ 0.5 * step(h.edges[1]))
            ct_peak_idx = findmax(h.weights)[2]
            ct_peak_pos = mp[ct_peak_idx]

            init_fit_params = [ h.weights[ct_peak_idx] * 2π * step(h.edges[1]), 3 * step(h.edges[1]), mp[ct_peak_idx] ]
            # println("init fit params: $ichn_src, $ichn_tar", init_fit_params)
            # init_fit_params = [ h.weights[ct_peak_idx] *  step(h.edges[1]), 2 * step(h.edges[1]), mp[ct_peak_idx] ]
            # fitrange = (mp[ct_peak_idx] - 3 * step(h.edges[1])):step(h.edges[1]):(mp[ct_peak_idx] + 3 * step(h.edges[1]))
            fitrange = ((mp[ct_peak_idx] - 7 * step(h.edges[1])), (mp[ct_peak_idx] + 7 * step(h.edges[1])))## <- modify for distorted pulses
            fitf = RadiationSpectra.FitFunction{T}( scaled_cauchy, 1, 3 )
            set_initial_parameters!(fitf, init_fit_params)
            set_fitranges!(fitf, (fitrange,))
            RadiationSpectra.lsqfit!(fitf, h)
            push!(tmpfits, fitf)
            ct_peak_pos = fitf.fitted_parameters[3]
            # push!(tmpfits, fr)
            # ct_peak_pos = fr.parameters[3]

            if ichn_tar != ichn_src
                c[ichn_src, ichn_tar] = ct_peak_pos * c[1, 1]
                if create_plots
                    p = plot(h, st=:step, legend=false, title="mpa ratio: seg$(iseg_tar) / core", )
                    bw = StatsBase.binvolume(h, 1)
                    plot!(p, fitf, bin_width = bw)
                    # plot!(p, fr)
                    push!(plts, p)
                end
            elseif create_plots
                p = plot( [0,1], [0,1], legend=false, color=:red, title="mpa ratio: seg$(iseg_tar) / core")
                plot!(p,  [0,1], [1,0], legend=false, color=:red )
                push!(plts, p)
            end
        end
        if create_plots
            p = plot(plts..., size=(1920 ,1080))
            savefig(m, p, "2_3_calibration_matrix/ct_mpa_ratio_peak_fits_$(m.name)", "mca_ratio_ct_peak_fit_sseg$iseg_src", fmt=:png); p = 0;
        end
        push!(ct_peak_fits, tmpfits)
    end

    c = inv(c)
    c[1, 1] = c0

    if create_plots
        p = plot(c, st=:heatmap, aspect_ratio=1, size=(900,900));
        savefig(m, p, "2_3_calibration_matrix", "3_calibration_matrix", fmt=:png ); p=0;
    end

    return c, ratio_hists, cal_ratio_hists, cal_peak_fits, ct_ratio_hists, ct_peak_fits
end
