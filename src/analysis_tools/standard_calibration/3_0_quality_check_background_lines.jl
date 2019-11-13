function quality_check(	m::Measurement;	overwrite=false, quality_check_photon_lines=[])
	energies = get_energies(m)
	n_channel::Int, n_events::Int = size(energies)
	multi_channel_det::Bool = n_channel > 1
	T = Float64
	daqT = exists(m, "Results/DAQ_Lifetime") ? GAT.read_analysis_result_dataset(m, "DAQ_Lifetime") : missing

	df = DataFrame( Isotope = String[],
					Energy = Float64[],
					Type = Symbol[]	)

	# push!(df, ["Pb214",  351.932, :BackgroundPhotonLine ])
	push!(df, ["Bi214",  609.312, :BackgroundPhotonLine ])
	push!(df, ["K40",   1460.830, :BackgroundPhotonLine ])
	# push!(df, ["Bi214", 1764.494, :BackgroundPhotonLine ])
	push!(df, ["Tl208", 2614.533, :BackgroundPhotonLine ])
	photon_lines::Vector{Float64} = sort!(vcat([e for e in df[!, :Energy]], quality_check_photon_lines))

	photon_lines_fit_parameters_core::Array{Float64, 2} = zeros(Float64, 6, length(photon_lines))
	if multi_channel_det photon_lines_fit_parameters_sumseg::Array{Float64, 2} = zeros(Float64, 6, length(photon_lines)) end

	h_core    = float(Histogram(0:1.0:3000, :left))
	if multi_channel_det 
		h_sumsegs = float(Histogram(0:1.0:3000, :left))
		sum_seg_energies = zeros(eltype(energies), n_events )
		for i in eachindex(sum_seg_energies)
			sum_seg_energies[i] = sum(energies[2:end, i])
		end
		append!(h_sumsegs, sum_seg_energies) 
	end
	append!(h_core, energies[1, :])
	if !ismissing(daqT)
		h_core.weights = h_core.weights ./ daqT
		if multi_channel_det h_sumsegs.weights = h_sumsegs.weights ./ daqT end
	end
	p_spectra = plot(h_core, st=:step, label="core", xlims=[0, 3000], xticks=collect(0:500:3000), ylabel = "cts / s")
	if multi_channel_det plot!(h_sumsegs, st=:step, label="Summed Segments") end
	p_peaks = []
	p_fr_results = []
	for (ipl, pl) in enumerate(photon_lines)
		photon_lines_fit_parameters_core[1, ipl] = pl
		if multi_channel_det photon_lines_fit_parameters_sumseg[1, ipl] = pl end
		idx_pl = searchsortedfirst(h_core.edges[1], pl)
		maximum_counts_core = maximum(h_core.weights[(idx_pl-11):(idx_pl+11)])
		maximum_counts_segs = multi_channel_det ? maximum(h_sumsegs.weights[(idx_pl-11):(idx_pl+11)]) : maximum_counts_core
		p = plot(h_core, st=:step, lw=2., legend=false, xlims=[pl - 10, pl + 10], ylims=[0, max(maximum_counts_core, maximum_counts_segs)])
		if multi_channel_det plot!(h_sumsegs, st=:step, lw=2.) end
		scale::Float64 = maximum_counts_core
		σ::Float64 = 1
		μ::Float64 = pl
		coeff_0::Float64 = mean(h_sumsegs.weights[(idx_pl-20):(idx_pl-10)])
		coeff_1::Float64 = 0
		# fitrange = pl - 10:pl + 10
		fitrange = (pl - 10, pl + 10)
		fitf_core = RadiationSpectra.FitFunction{T}( RadiationSpectra.Gauss_plus_linear_background, 1, 5 )
		set_fitranges!(fitf_core, (fitrange,))
		set_initial_parameters!(fitf_core, [scale, σ, μ, coeff_0, coeff_1])
		RadiationSpectra.lsqfit!(fitf_core, h_core)
		# fr_core = GeDetSpectrumAnalyserTmp.fit(h_core, fitrange, RadiationSpectra.Gauss_plus_linear_background, Float64[scale, σ, μ, coeff_0, coeff_1])
		# println(RadiationSpectra.get_fitted_parameters(fitf_core))
		if multi_channel_det
			scale = maximum_counts_segs
			fitf_segs = RadiationSpectra.FitFunction{T}( RadiationSpectra.Gauss_plus_linear_background, 1, 5 )
			set_fitranges!(fitf_segs, (fitrange,))
			set_initial_parameters!(fitf_segs, [scale, σ, μ, coeff_0, coeff_1])
			RadiationSpectra.lsqfit!(fitf_segs, h_sumsegs)
			# println(RadiationSpectra.get_fitted_parameters(fitf_segs))
		# fr_segs = GeDetSpectrumAnalyserTmp.fit(h_sumsegs, fitrange, RadiationSpectra.Gauss_plus_linear_background, Float64[scale, σ, μ, coeff_0, coeff_1])
		end

		# plot!(fr_segs)
		# plot!(fr_core)
		if multi_channel_det plot!(fitf_segs) end
		plot!(fitf_core)
		push!(p_peaks, p)
		p_fr_result = plot(legend=false,grid=false,foreground_color_subplot=:white)
		# pars = fr_core.parameters
		pars = collect(get_fitted_parameters(fitf_core))
		photon_lines_fit_parameters_core[2:end, ipl] = pars

        # df[:FitParameters][1][1] = pars
	# 	df[ipl, 4][:, 1] = pars

		annotate!(0.4, 0.8, L"Core")
		if multi_channel_det annotate!(0.7, 0.8, L"SumSegs") end
		annotate!(0.2, 0.65, L"\mu:")
		annotate!(0.2, 0.50, L"\sigma:")
		annotate!(0.2, 0.35, L"A:")
		annotate!(0.4, 0.65, "$(round(pars[3], digits=1))")
		annotate!(0.4, 0.50, "$(round(pars[2], digits=1))")
		annotate!(0.4, 0.35, "$(round(pars[1], digits=1))")
		if multi_channel_det
			pars = collect(get_fitted_parameters(fitf_segs) )
			photon_lines_fit_parameters_sumseg[2:end, ipl] = pars

		# 	df[ipl, 4][:, end]  = pars
			annotate!(0.7, 0.65, "$(round(pars[3], digits=1))")
			annotate!(0.7, 0.50, "$(round(pars[2], digits=1))")
			annotate!(0.7, 0.35, "$(round(pars[1], digits=1))")
		end
		push!(p_fr_results, p_fr_result)
	end
	p_core_summed_segments = plot(p_peaks..., layout=(1, length(photon_lines)))
	p_fr_results = plot(p_fr_results..., layout=(1, length(photon_lines)))

	write_analysis_result_dataset(m, "background_photon_lines_fit_parameters_core", photon_lines_fit_parameters_core)
	if multi_channel_det write_analysis_result_dataset(m, "background_photon_lines_fit_parameters_sumseg", photon_lines_fit_parameters_sumseg) end

	p_summary = plot(p_spectra, p_core_summed_segments, p_fr_results, layout=(@layout [a{0.4h}; b{0.4h}; c]), size=(1600,1000))
	savefig(m, p_summary, "3_0_data_quality_check_background", "background_lines", fmt=:png)

	return nothing
end
