function get_number_of_events(fn::AbstractString)::Int
	h = h5open(fn, "r+")
	g = g_open(h,"DAQ_Data")
	d = d_open(g,"daq_energies")
	s = size(d, 2)
	close(h)
	return s
end
function get_number_of_events(m::Measurement)::Int
	inputfiles = gather_absolute_paths_to_hdf5_input_files(m)
	n_events::Int = 0
	for f in inputfiles
		n_events += get_number_of_events(f)
	end
	return n_events
end

function get_number_of_channel(fn::AbstractString)::Int
	h = h5open(fn, "r+")
	g = g_open(h,"DAQ_Data")
	d = d_open(g,"daq_energies")
	s = size(d, 1)
	close(h)
	return s
end
function get_number_of_channel(m::Measurement)::Int
	inputfiles = gather_absolute_paths_to_hdf5_input_files(m)
	n_channel::Int = get_number_of_channel(inputfiles[1])
	return n_channel
end
