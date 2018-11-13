function quality_check(	m::Measurement;	overwrite=false)
	energies = get_energies(m)
	n_channel::Int, n_events::Int = size(energies)

	df = DataFrame( Isotope = String[],
					Energy = Float64[],
					Type = Symbol[]	)

	push!(df, ["Pb214",  351.932, :BackgroundPhotonLine ])
	push!(df, ["Bi214",  609.312, :BackgroundPhotonLine ])
	push!(df, ["K40",   1460.830, :BackgroundPhotonLine ])
	push!(df, ["Bi214", 1764.494, :BackgroundPhotonLine ])
	push!(df, ["Tl208", 2614.533, :BackgroundPhotonLine ])
	photon_lines::Vector{Float64} = [e for e in df[:Energy]]

	photon_lines_fit_parameters_core::Array{Float64, 2} = zeros(Float64, 6, length(photon_lines))
	photon_lines_fit_parameters_sumseg::Array{Float64, 2} = zeros(Float64, 6, length(photon_lines))

	h_core    = Histogram(0:1:3000, :left)
	h_sumsegs = Histogram(0:1:3000, :left)
	sum_seg_energies = zeros(eltype(energies), n_events )
	for i in eachindex(sum_seg_energies)
		sum_seg_energies[i] = sum(energies[2:end, i])
	end
	append!(h_core, energies[1, :])
	append!(h_sumsegs, sum_seg_energies)
	p_spectra = plot(h_core, st=:step, label="core", xlims=[0, 3000], xticks=collect(0:500:3000))
	plot!(h_sumsegs, st=:step, label="Summed Segments")
	p_peaks = []
	p_fr_results = []
	for (ipl, pl) in enumerate(photon_lines)
		photon_lines_fit_parameters_core[1, ipl] = pl
		photon_lines_fit_parameters_sumseg[1, ipl] = pl
		idx_pl = searchsortedfirst(h_core.edges[1], pl)
		maximum_counts_core = maximum(h_core.weights[(idx_pl-11):(idx_pl+11)])
		maximum_counts_segs = maximum(h_sumsegs.weights[(idx_pl-11):(idx_pl+11)])
		p = plot(h_core, st=:step, lw=2., legend=false, xlims=[pl - 10, pl + 10], ylims=[0, max(maximum_counts_core, maximum_counts_segs)])
		plot!(h_sumsegs, st=:step, lw=2.)
		scale::Float64 = maximum_counts_core
		σ::Float64 = 1
		μ::Float64 = pl
		coeff_0::Float64 = 0
		coeff_1::Float64 = 0
		fitrange = pl - 10:pl + 10
		fr_core = GeDetSpectrumAnalyserTmp.fit(h_core, fitrange, gauss_plus_first_order_polynom, Float64[scale, σ, μ, coeff_0, coeff_1])
		scale = maximum_counts_segs
		fr_segs = GeDetSpectrumAnalyserTmp.fit(h_sumsegs, fitrange, gauss_plus_first_order_polynom, Float64[scale, σ, μ, coeff_0, coeff_1])
		plot!(fr_segs)
		plot!(fr_core)
		push!(p_peaks, p)
		p_fr_result = plot(axis=false, grid=false)
		pars = fr_core.parameters
		photon_lines_fit_parameters_core[2:end, ipl] = pars

        # df[:FitParameters][1][1] = pars
	# 	df[ipl, 4][:, 1] = pars

		annotate!(0.4, 0.8, L"Core")
		annotate!(0.7, 0.8, L"SumSegs")
		annotate!(0.2, 0.65, L"\mu:")
		annotate!(0.2, 0.50, L"\sigma:")
		annotate!(0.2, 0.35, L"A:")
		annotate!(0.4, 0.65, "$(round(pars[3], digits=1))")
		annotate!(0.4, 0.50, "$(round(pars[2], digits=1))")
		annotate!(0.4, 0.35, "$(round(pars[1], digits=1))")
		pars = fr_segs.parameters
		photon_lines_fit_parameters_sumseg[2:end, ipl] = pars

	# 	df[ipl, 4][:, end]  = pars
		annotate!(0.7, 0.65, "$(round(pars[3], digits=1))")
		annotate!(0.7, 0.50, "$(round(pars[2], digits=1))")
		annotate!(0.7, 0.35, "$(round(pars[1], digits=1))")
		push!(p_fr_results, p_fr_result)
	end
	p_core_summed_segments = plot(p_peaks..., layout=(1, length(photon_lines)))
	p_fr_results = plot(p_fr_results..., layout=(1, length(photon_lines)))

	write_analysis_result_dataset(m, "background_photon_lines_fit_parameters_core", photon_lines_fit_parameters_core)
	write_analysis_result_dataset(m, "background_photon_lines_fit_parameters_sumseg", photon_lines_fit_parameters_sumseg)

	p_summary = plot(p_spectra, p_core_summed_segments, p_fr_results, layout=(@layout [a{0.4h}; b{0.4h}; c]), size=(1600,1000))
	savefig(m, p_summary, "3_0_data_quality_check_background", "background_lines", fmt=:png)

	return nothing
end
