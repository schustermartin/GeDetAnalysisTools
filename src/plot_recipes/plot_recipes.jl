@userplot Plot_Energy_Histogram
@recipe function f(gdd::Plot_Energy_Histogram)
	if !( typeof(gdd.args[1]) <: Union{Array{<:Real, 1}, Histogram})
		error("First input argument should be subtype of `Union{Array{<:Real, 1}, Histogram}`.")
	end
	e = gdd.args[1]

	size --> (900, 500)
	legend --> false
	units = "keV"
	xlabel --> "E / $units"
	ylabel --> "Counts"
	# title --> "Energy Histogram"
	seriestype --> :step

	if typeof(e) <: Array{<:Real, 1}
		energyrange = 0:1:3000
		h = fit(Histogram, e, energyrange, closed=:left)
		mean_hight = mean(h.weights)
		if mean_hight < mean(h.weights[1:3]) yscale --> :log10 end
		h
	elseif typeof(e) <: Histogram
		h = e
		mean_hight = mean(h.weights)
		if mean_hight < mean(h.weights[1:3]) yscale --> :log10 end
		h
	end
end


@userplot EnergyHistogramDisplay
@recipe function f(gdd::EnergyHistogramDisplay; edges=0:1:3000)
	if !( typeof(gdd.args[1]) <: Union{Array{<:Real, 2}, Array{<:Histogram, 1}, Measurement})
		error("First input argument should be subtype of `Union{Array{<:Real, 2}, Array{<:Histogram, 1}, Measurement`.")
	end

	if typeof(gdd.args[1]) <: Measurement
		m = gdd.args[1]
		e = get_energies(m)
		n_channel = size(e, 1)

		size --> (1920, 1080)
		# legend --> false
		xlabel --> "E / keV"
		ylabel --> "Counts"
		seriestype := :step

		layout --> if (length(gdd.args) > 1 && typeof(gdd.args[2]) <: Detector)
			@layout [ 	chn1{0.16666w} chn20{0.16666w} 
						chn2 chn3 chn4 chn5 chn6 chn7
						chn8 chn9 chn10 chn11 chn12 chn13
						chn14 chn15 chn16 chn17 chn18 chn19]
		else
			(n_channel)
		end
		channel_order = if (length(gdd.args) > 1 && typeof(gdd.args[2]) <: Detector)
			gdd.args[2].channel_display_order
		else
			Int[chn for chn in 1:n_channel]
		end

		y_max = 1

		for ichn in 1:n_channel
			@series begin
		        subplot := channel_order[ichn]
				h = fit(Histogram, e[ichn, :], edges, closed=:left)
		        if in(10., edges)
					mean_hight = mean(h.weights[StatsBase.binindex(h, 10.):end])
					y_max = maximum(h.weights[StatsBase.binindex(h, 10.):end]) 
					ylims --> (0.5, y_max)
					if mean_hight < mean(h.weights[1:3]) yscale --> :log10 end
				end
				label --> "Chn $(ichn)"
		        h
		    end
		end
	else

		e = gdd.args[1]
		n_channel = minimum(size(e))

		size --> (1920, 1080)
		# legend --> false
		xlabel --> "E / keV"
		ylabel --> "Counts"
		seriestype := :step

		layout --> if (length(gdd.args) > 1 && typeof(gdd.args[2]) <: Detector)
			@layout [ 	chn1{0.16666w} chn20{0.16666w} 
						chn2 chn3 chn4 chn5 chn6 chn7
						chn8 chn9 chn10 chn11 chn12 chn13
						chn14 chn15 chn16 chn17 chn18 chn19]
		else
			(n_channel)
		end
		channel_order = if (length(gdd.args) > 1 && typeof(gdd.args[2]) <: Detector)
			gdd.args[2].channel_display_order
		else
			Int[chn for chn in 1:n_channel]
		end

		y_max = 1

		for ichn in 1:n_channel
			if typeof(e) <: Array{<:Real, 2}
				@series begin
			        subplot := channel_order[ichn]
					h = fit(Histogram, e[ichn, :], edges, closed=:left)
			        if in(10., edges)
						mean_hight = mean(h.weights[StatsBase.binindex(h, 10.):end])
						y_max = maximum(h.weights[StatsBase.binindex(h, 10.):end]) 
						ylims --> (0.5, y_max)
						if mean_hight < mean(h.weights[1:3]) yscale --> :log10 end
					end
					label --> "Chn $(ichn)"
			        h
			    end
			elseif typeof(e) <:  Array{<:Histogram, 1}
				@series begin
					h = e[ichn]
					subplot := channel_order[ichn]
					if in(10., h.edges[1])
						y_max = maximum(h.weights[StatsBase.binindex(h, 10.):end]) 
						mean_hight = mean(h.weights[StatsBase.binindex(h, 10.):end])
						ylims --> (0.5, y_max)
						if mean_hight < mean(h.weights[1:3]) yscale --> :log10 end
					end
					label --> "Chn $(ichn)"
					h
				end
			end
		end
	end
end

@userplot Histogramdisplay
@recipe function f(gdd::Histogramdisplay; nbins=500)
	if !( typeof(gdd.args[1]) <: Union{Array{<:Real, 2}, Array{<:Histogram, 1}})
		error("First input argument should be subtype of `Union{Array{<:Real, 2}, Array{<:Histogram, 1}`.")
	end
	e = gdd.args[1]
	n_channel = minimum(size(e))

	size --> (1920, 1080)
	xlabel --> "E / keV"
	ylabel --> "Counts"
	seriestype := :step

	layout --> if (length(gdd.args) > 1 && typeof(gdd.args[2]) <: Detector)
		@layout [ 	chn1{0.16666w} chn20{0.16666w} 
					chn2 chn3 chn4 chn5 chn6 chn7
					chn8 chn9 chn10 chn11 chn12 chn13
					chn14 chn15 chn16 chn17 chn18 chn19]
	else
		(n_channel)
	end
	channel_order = if (length(gdd.args) > 1 && typeof(gdd.args[2]) <: Detector)
		gdd.args[2].channel_display_order
	else
		Int[chn for chn in 1:n_channel]
	end

	for ichn in 1:n_channel
		if typeof(e) <: Array{<:Real, 2}
			@series begin
		        subplot := channel_order[ichn]
				h = fit(Histogram, e[ichn, :], nbins=nbins, closed=:left)
				label --> "Chn $(ichn)"
		        h
		    end
		elseif typeof(e) <: Array{<:Histogram, 1}
			@series begin
				h = e[ichn]
				subplot := channel_order[ichn]
				label --> "Chn $(ichn)"
				h
			end
		end
	end
end


@userplot Pulsedisplay
@recipe function f(pd::Pulsedisplay)
    if isa(pd.args[1], Array{<:Real, 2})
        pulses = pd.args[1]
    end

    n_samples, n_channel = size(pulses)

    legend --> false
    size --> (1920, 1080)
    println()
    layout --> if (length(pd.args) > 1 && typeof(pd.args[2]) <: Detector)
		@layout [ 	chn1{0.16666w} chn20{0.16666w} 
					chn2 chn3 chn4 chn5 chn6 chn7
					chn8 chn9 chn10 chn11 chn12 chn13
					chn14 chn15 chn16 chn17 chn18 chn19]
	else
		(n_channel)
	end
	channel_order = if (length(pd.args) > 1 && typeof(pd.args[2]) <: Detector)
		pd.args[2].channel_display_order
	else
		Int[chn for chn in 1:n_channel]
	end


    for ichn in 1:n_channel
        @series begin
            subplot := channel_order[ichn]
            title --> "Channel $(ichn)"
            pulses[:, ichn]
        end
    end   
end
