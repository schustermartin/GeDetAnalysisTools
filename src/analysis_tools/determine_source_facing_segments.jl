function fit_source_peak_and_determine_sfc(m::Measurement, peakenergy::T, sno_limit::T=T(0.05); min_amplitude::T=T(100.0), overwrite = false) where {T<:Real} ##::Array{Int, 1}
	if exists(m, "Results/source_facing_channels_$(peakenergy)keV") && !overwrite
		println("`Results/source_facing_channels_$(peakenergy)keV` already determined for $(m.name). Skipping...")
		return nothing
	end
	
	source_facing_channels = Int[]
	energies::Array{T, 2} = get_energies(m)
	n_channel = size(energies, 1)
	ss_selection_diff::T = 20
	core::Int = 1

	ehists = [Histogram(-10:1.0:10, :left) for chn in 1:n_channel]

	for ievt in 1:size(energies, 2)
		ssidx = get_single_segment_channel_index_abs(energies[:, ievt], ss_selection_diff)
		if ssidx > 0
			push!(ehists[core],  energies[core, ievt]  - peakenergy )
			push!(ehists[ssidx], energies[ssidx, ievt] - peakenergy )
		end
	end

	# par_names = ["Scale", "σ", "Offset"] #, "Slope"]
	# init_fit_parameters = T[ min_amplitude, 1.0, 0.0] #, 0.]
	@fastmath function f(x::T, par::Array{T, 1})::T
		if par[1] < 0 || par[2] < 0.6 || par[2] > 5.0 || par[3] < 0 return NaN end
		return par[1] / (sqrt(2 * π * par[2]^2)) * exp(-0.5 * (x^2) / (par[2]^2)) + par[3]# + par[4] * x
	end
	f(x::Array{T, 1}, par::Array{T, 1})::Array{T, 1} = T[f(v, par) for v in x]

	fit_functions = RadiationSpectra.FitFunction[]
	for ichn in 1:n_channel
		fitf = RadiationSpectra.FitFunction{T}( f, 1, 3  )
		set_fitranges!(fitfunc, ((-5, 5),) )
		p0 = (
			scale = min_amplitude, 
			σ = 1.0, 
			offset = mean(ehists[ichn].weights[1:4])
		)
		set_initial_parameters!(fitf, p0)
		push!(fit_functions, fitf)
	end
	
	σ_core::T = 1.
	offset_core::T = 0.
	for ichn in 1:n_channel
		RadiationSpectra.lsqfit!(fit_functions[ichn], ehists[ichn])
		fitted_pars = collect(get_fitted_parameters(fit_functions[ichn]))
		# funcs[ichn].par = fit_result.parameters
		if ichn == 1
			σ_core = abs(fitted_pars[2])
			offset_core = fitted_pars[3]
		end
		σ = abs(fitted_pars[2])
		scale = fitted_pars[1] # scale
		offset = fitted_pars[3]
		bg_area = 2 * 2σ_core * offset_core
		sno = scale / (2 * 2σ * offset)
		sno_to_core = scale / bg_area
		# println(ichn, "\t", sno_to_core)
		if σ < 3.0 && σ > 0.4 && sno_to_core > sno_limit && scale >= min_amplitude
			# info("Channel $(chn):\tσ=$(signif(σ, 2))\t-\tSNO (2σ) = $(sno)")
			if ichn >= 1 push!(source_facing_channels, ichn) end
		end
	end

	write_analysis_result_dataset(m, "source_facing_channels_$(peakenergy)keV", source_facing_channels)
	println("created result dataset: `source_facing_channels_$(peakenergy)keV`")
	plt = histogramdisplay(ehists, m.detector, size=(2000,1200))
	for i in eachindex(fit_functions)
		c = in(i, source_facing_channels) ? "red" : "black"
		plot!(plt.subplots[m.detector.channel_display_order[i]], fit_functions[i], lc=c, lw=2)
	end
	savefig(m, plt, "4_0_source_facing_segment_determination", "source_facing_segment_determination_$(peakenergy)keV", fmt=:png)

	# return source_facing_channels, ehists, funcs, frs
	return source_facing_channels, ehists, fit_functions
end
