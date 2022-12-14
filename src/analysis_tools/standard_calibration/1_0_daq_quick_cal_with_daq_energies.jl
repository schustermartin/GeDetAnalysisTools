# function determine_daq_core_calibration_constant(m::Measurement; photon_lines=[609.312, 911.204, 1120.287, 1460.830, 1764.494, 2614.533])
# 	daq_core_energies = get_daq_energies(m)[1, :]
# 	h = fit(Histogram, daq_core_energies, nbins=10000, closed=:left)
# 	c, pcg_hist = RadiationSpectra.determine_calibration_constant_through_peak_ratios(h, photon_lines )
# 	# c, pcg_hist = GeDetSpectrumAnalyserTmp.determine_calibration_constant_through_peak_ratios(h)
# 	return c, pcg_hist
# end

function determine_calibration_matrix_with_daq_energies(m::Measurement; photon_lines=[609.312, 911.204, 1120.287, 1460.830, 1764.494, 2614.533], nbins = 10000)
	@fastmath @inbounds begin
		daq_core_energies = transpose(get_daq_energies(m))
		T = Float32
		n_events::Int, n_channel::Int = size(daq_core_energies)
		n_segments::Int = n_channel - 1
		h_core = fit(Histogram, daq_core_energies[:, 1], nbins=nbins, closed=:left)
		@info photon_lines
		try_cnt = 1
		sigma_modifiers=[0.0,0.25,0.5,1,1.5,2.0,-0.25,-0.5,-1,-1.5,-2.0]
		while try_cnt<=11
		c0_daq, pcg_hist = try
				try_cnt >1 ?  @info("try $try_cnt: now using σ = $(3.5+sigma_modifiers[try_cnt])") : nothing
				RadiationSpectra_beforeBAT.determine_calibration_constant_through_peak_ratios(h_core, photon_lines, min_n_peaks = 20, threshold = 40.0, α = 0.5, σ = 3.5+sigma_modifiers[try_cnt], rtol = 7e-3 )
			catch err
				try_cnt +=1
				println("daq data quick calibration failed with σ = $(3.5+sigma_modifiers[try_cnt-1]). Retrying with modified σ..." )
				nothing, nothing
			end
			c0_daq != nothing ? try_cnt = 12 : nothing
		end
		if c0_daq == nothing
			@info "daq precal factor determination failed, processing with 1e-3"
			c0_daq= 1e-3
		end
		@info c0_daq
		c0_daq, core_peak_fits, core_c0_fit = RadiationSpectra_beforeBAT.determine_calibration_constant_through_peak_fitting(h_core, photon_lines, c0_daq)

		ratios = Array{T, 2}(undef, size(daq_core_energies, 1), size(daq_core_energies, 2) - 1 )

		daq_core_energies *= c0_daq

		for i in eachindex(daq_core_energies)
			if daq_core_energies[i] == 0
				daq_core_energies[i] = 1
			end
		end

		daq_core_energies_inv::Array{T, 1} = inv.(daq_core_energies[:, 1])

		for iseg in 1:size(ratios, 2)
			ichn = iseg + 1
			ratios[:, iseg] = daq_core_energies[:, ichn] .* daq_core_energies_inv[:]
		end

		ratio_hists = [ Histogram(-1.0:0.01:2, :left ) for iseg in 1:size(ratios, 2) ]
		for iseg in 1:size(ratios, 2)
			append!(ratio_hists[iseg], ratios[:, iseg])
		end

		cal_peak_pos = zeros(Float64, n_segments)
		ct_peak_pos  = zeros(Float64, n_segments)
		for iseg in 1:n_channel-1
			h = ratio_hists[iseg]
			mp = (h.edges[1][1:length(h.edges[1])-1] .+ 0.5 * step(h.edges[1]))
			ct_peak_idx = findmax(h.weights)[2]
	        ct_peak_pos[iseg] = mp[ct_peak_idx]
	        start_idx = StatsBase.binindex(h, ct_peak_pos[iseg] + 0.25)
	        cal_peak_idx = findmax(h.weights[start_idx:end])[2]
	        cal_peak_pos[iseg] = mp[start_idx + cal_peak_idx]
		end

		# Segment Calibration
		cal_ratio_hists = [ fit(Histogram, ratios[:, iseg], cal_peak_pos[iseg] - 0.05:0.001:cal_peak_pos[iseg] + 0.05, closed=:left ) for iseg in 1:size(ratios, 2)]
		cal_peak_pos = zeros(Float64, n_segments)
		for iseg in 1:n_channel-1
			h = cal_ratio_hists[iseg]
			mp = (h.edges[1][1:length(h.edges[1])-1] .+ 0.5 * step(h.edges[1]))
			cal_peak_idx = findmax(h.weights)[2]
			# init_fit_params = [ h.weights[cal_peak_idx] * 2π * step(h.edges[1]), 2 * step(h.edges[1]), mp[cal_peak_idx] ]
			fitrange = (mp[cal_peak_idx] - 8 * step(h.edges[1])), (mp[cal_peak_idx] + 8 * step(h.edges[1]))
			# println("test")
			fitf = RadiationSpectra_beforeBAT.FitFunction{Float64}( scaled_cauchy, 1, 3 )
			p0 = (
				scale =  h.weights[cal_peak_idx] * 2π * step(h.edges[1]),
				σ = 2.5 * step(h.edges[1]),
				μ = mp[cal_peak_idx]
			)
			set_initial_parameters!(fitf, p0)
			set_fitranges!(fitf, ( ((mp[cal_peak_idx] - 7 * step(h.edges[1])), (mp[cal_peak_idx] + 7 * step(h.edges[1]))), ))
			# fitf.initial_parameters = init_fit_params
			# fitf.fitrange = fitrange

			RadiationSpectra_beforeBAT.lsqfit!(fitf, h)
	        # fr = GeDetSpectrumAnalyserTmp.fit(h, fitrange, scaled_cauchy, init_fit_params)
	        # cal_peak_pos[iseg] = fr.parameters[3]
	        cal_peak_pos[iseg] = fitf.fitted_parameters[3]
		end

		# Calibration Matrix
		c = zeros(T, n_channel, n_channel)
		for i in 1:n_channel c[i,i] = 1 end
		c[1, 1] = inv(c0_daq)
		for ichn in 2:n_channel
			iseg = ichn - 1
			c[ichn, ichn] = cal_peak_pos[iseg] * c[1, 1]
		end

		# Crosstralk Values
	 	ratios = transpose(ratios)
        ct_ratios = [ [Histogram(-0.3:0.001:0.3, :left) for i in 2:n_channel] for j in 2:n_channel ]
		for ichn_src in 2:n_channel
	        iseg_src = ichn_src - 1
	        ss_ratio_limits::Tuple{T, T} = cal_peak_pos[iseg_src] - 0.03, cal_peak_pos[iseg_src] + 0.03 # ss: single segment
        	for i in axes(ratios, 2) # n_events
        		if ss_ratio_limits[1] <= ratios[iseg_src, i] <= ss_ratio_limits[2] # if single segment event
			        for ichn_tar in 2:n_channel
			        	iseg_tar = ichn_tar - 1
			        	if ichn_tar == ichn_src continue end
			        	push!(ct_ratios[iseg_src][iseg_tar], ratios[iseg_tar, i])
			        end
        		end
                # if ss_ratio_limits[1] <= ratios[i, iseg_src] <= ss_ratio_limits[2] push!(ct_ratios, ratios[i, iseg_tar]) end
	        end
	    end
	    for ichn_src in 2:n_channel
	    	iseg_src = ichn_src - 1
    		for ichn_tar in 2:n_channel
    			iseg_tar = ichn_tar - 1
    			if ichn_tar == ichn_src continue end
    			h = ct_ratios[iseg_src][iseg_tar]
    			mp = (h.edges[1][1:length(h.edges[1])-1] .+ 0.5 * step(h.edges[1]))
    			ct_peak_idx = findmax(h.weights)[2]
	        	ct_peak_pos = mp[ct_peak_idx]
	        	c[ichn_src, ichn_tar] = ct_peak_pos * c[1, 1]
    		end
	    end

		c = inv(c)
		return c
	end
end


"""
	calibrate_energies(energies::Array{<:Real, 2}, c::Array{<:Real, 2})::Array{Real, 2}

The dimensions of `energies` should be (channel, events).

`c` should be the already converted cross-talk and calibration matrix.
"""
function calibrate_energies(energies::Array{<:AbstractFloat, 2}, c::Array{<:Real, 2})::Array{Real, 2}
	T = eltype(energies)
	c = T.(c)
	n_channel = size(energies, 1)
	n_events = size(energies, 2)
	cal_energies = Array{T, 2}(undef, n_channel, n_events)
	c_transpose::Array{T, 2} = c'
	for i in axes(energies, 2)
		# evt_energies = zeros(T, n_channel)
		# evt_energies[1] = energies[1, i] * c[1,1]
		# for ichn_tar in 2:n_channel
		# 	iseg_tar = ichn_tar - 1
		# 	for ichn_src in 2:n_channel
		# 		iseg_src = ichn_src - 1
		# 		evt_energies[ichn_tar]  += c[ichn_src, ichn_tar] * energies[ichn_src, i]
		# 	end
		# end
		cal_energies[:, i] = c_transpose * energies[:, i]
	end
	return cal_energies
end
function calibrate_energies(energies::Array{<:Real, 2}, c::Array{<:Real, 2}, T::Type=Float32)::Array{Real, 2}
	return calibrate_energies(T.(energies), c)
end

function get_quick_calibrated_daq_energies(m::Measurement, T::Type=Float32; photon_lines = [609.312, 911.204, 1120.287, 1460.830, 1764.494, 2614.533], nbins = 10000)::Array{T, 2}
	c = determine_calibration_matrix_with_daq_energies(m, photon_lines = photon_lines, nbins = nbins)
	daq_energies = get_daq_energies(m)
	return calibrate_energies(daq_energies, c, T)
end
