function get_daq_energies(fn::AbstractString)::Array{<:Real,2}
	return h5read(fn, "DAQ_Data/daq_energies")
end
function get_daq_energies(m::Measurement)::Array{<:Real,2}
	inputfiles = gather_absolute_paths_to_hdf5_input_files(m)
	n_events::Int = get_number_of_events(m)
	n_channel::Int = get_number_of_channel(m)
	e1 = get_daq_energies(inputfiles[1])
	T = eltype(e1)
	e = Array{T, 2}(undef, n_channel, n_events)
	last_idx::Int = size(e1, 2)
	e[:, 1:last_idx] = e1
	@inbounds for i in 2:length(inputfiles)
		en = get_daq_energies(inputfiles[i])
		n_new_events::Int = size(en, 2)
		e[:, last_idx+1:last_idx+n_new_events] = en
		last_idx += n_new_events
	end
	return e
end

function get_daq_times(fn::AbstractString)
	return h5read(fn, "DAQ_Data/daq_time")
end
function get_daq_times(m::Measurement)::Array{<:Real,2}
	inputfiles = gather_absolute_paths_to_hdf5_input_files(m)
	n_events::Int = get_number_of_events(m)
	n_channel::Int = get_number_of_channel(m)
	e1 = get_daq_times(inputfiles[1])
	T = eltype(e1)
	e = Array{T, 2}(undef, size(e1, 1), n_events)
	last_idx::Int = size(e1, 2)
	e[:, 1:last_idx] = e1
	@inbounds for i in 2:length(inputfiles)
		en = get_daq_times(inputfiles[i])
		n_new_events::Int = size(en, 2)
		e[:, last_idx+1:last_idx+n_new_events] = en
		last_idx += n_new_events
	end
	return e
end

function get_daq_energy_spectrum(fn::AbstractString, channel::Int = 1; nbins::Int = 10000)::Histogram
	daq_energies = get_daq_energies(fn)
	h = fit(Histogram, daq_energies[channel, :], nbins=nbins, closed=:left)
	return h
end
function get_daq_energy_spectrum(m::Measurement, channel::Int = 1; nbins::Int = 10000)::Histogram
	daq_energies = get_daq_energies(m)
	h = fit(Histogram, daq_energies[channel, :], nbins=nbins, closed=:left)
	return h
end

function get_daq_energy_spectra(fn::AbstractString; nbins::Int = 10000)::Array{Histogram, 1}
	daq_energies = get_daq_energies(fn)
	hists = [fit(Histogram, daq_energies[channel, :], nbins=nbins, closed=:left) for channel in 1:size(daq_energies, 1)]
	return hists
end
function get_daq_energy_spectra(m::Measurement; nbins::Int = 10000)::Array{Histogram, 1}
	daq_energies = get_daq_energies(m)
	hists = [fit(Histogram, daq_energies[channel, :], nbins=nbins, closed=:left) for channel in 1:size(daq_energies, 1)]
	return hists
end


function get_daq_pulse(fn::AbstractString, event_index::Int)::AbstractArray{<:Real, 2}
	h = h5open(fn, "r+")
	try
		g = g_open(h,"DAQ_Data")
		d = d_open(g,"daq_pulses")
		global daq_pulse = d[:, :, event_index][:, :] 
	catch err
		close(h)
		error(err)
	end
	close(h)
	return daq_pulse
end
function get_daq_pulse(m::Measurement, event_index::Int, file_number::Int=1)::AbstractArray{<:Real, 2}
	inputfiles = gather_absolute_paths_to_hdf5_input_files(m)
    return get_daq_pulse(inputfiles[file_number], event_index)
end

function get_daq_livetime(fn::AbstractString; cut = 0.1)::AbstractFloat
	times = get_daq_times(fn)
	??T = if length(size(times)) == 2
		diff(times[1, :])
	elseif length(size(times)) == 1
		diff(times[:])
	else
		error("Strange dimensions of time array.")
	end
	# ??T = diff(get_daq_times(fn)[:])
	T = eltype(??T)
	daqtime::T = 0
	cut::T = T(cut)
	for ?? in ??T
		if ?? < cut
			daqtime += ??
		end
	end
	return daqtime
end
function get_daq_livetime(m::Measurement; cut = 0.1)::AbstractFloat
	inputfiles = gather_absolute_paths_to_hdf5_input_files(m)
	t = 0
	for fn in inputfiles
		t += get_daq_livetime(fn, cut = cut)
	end
	return t
end
