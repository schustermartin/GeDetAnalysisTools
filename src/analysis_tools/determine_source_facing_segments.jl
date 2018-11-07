function fit_source_peak_and_determine_sfc(m::Measurement, peakenergy::T, sno_limit::T=T(0.05); min_amplitude::T=T(100.0)) where {T<:Real} ##::Array{Int, 1}
	source_facing_channels = Int[]
	energies::Array{T, 2} = get_energies(m)
	n_channel = size(energies, 1)
	ss_selection_diff::T = 3
	core::Int = 1

	ehists = [Histogram(-10:0.5:10, :left) for chn in 1:n_channel]

	for ievt in 1:size(energies, 2)
		ssidx = get_single_segment_channel_index_abs(energies[:, ievt], ss_selection_diff)
		if ssidx > 0 
			push!(ehists[core],  energies[core, ievt]  - peakenergy )
			push!(ehists[ssidx], energies[ssidx, ievt] - peakenergy )
		end
	end

	par_names = ["Scale", "σ", "Offset", "Slope"]
	init_fit_parameters = T[ 0.0, 50.0, 0., 0.]
	f(x::T, par::Array{T, 1})::T = @fastmath par[1] / (sqrt(2 * π * par[2]^2)) * exp(-0.5 * (x^2) / (par[2]^2))  + par[3] + par[4] * x
	f(x::Array{T, 1}, par::Array{T, 1})::Array{T, 1} = T[f(v, par) for v in x]
	funcs = [MFunction("Gauss fixed at zero plus first order polynomial", init_fit_parameters, par_names, f) for ichn in 1:n_channel]

	frs = []

 	fit_results = GeDetSpectrumAnalyserTmp.Fit[]
	σ_core::T = 1.
	offset_core::T = 0.
	for ichn in 1:n_channel
		fit_result = GeDetSpectrumAnalyserTmp.fit(ehists[ichn], -10:10, funcs[ichn].f, funcs[ichn].par, estimate_uncertainties=false)
		funcs[ichn].par = fit_result.parameters
		if ichn == 1
			σ_core = abs(fit_result.parameters[2])
			offset_core = fit_result.parameters[3]
		end
		σ = abs(fit_result.parameters[2])
		scale = fit_result.parameters[1] # scale
		offset = fit_result.parameters[3]
		bg_area = 2 * 2σ_core * offset_core 
		sno = scale / (2 * 2σ * offset)
		sno_to_core = scale / bg_area
		# println(ichn, "\t", sno_to_core)
		if σ < 2.5 && σ > 0.4 && sno_to_core > sno_limit && scale >= min_amplitude
			# info("Channel $(chn):\tσ=$(signif(σ, 2))\t-\tSNO (2σ) = $(sno)")
			if ichn != 1 push!(source_facing_channels, ichn) end
		end
		push!(frs, fit_result)
	end

	write_analysis_result_dataset(m, "source_facing_channels", source_facing_channels)
	return source_facing_channels, ehists, funcs, frs
end