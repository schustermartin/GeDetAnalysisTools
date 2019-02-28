function view_file_structure(fn::AbstractString)::Nothing
	run(Cmd(`h5ls -r "$fn"`))
	nothing
end
function view_file_structure(m::Measurement; file_number::Int = 1)::Nothing
	run(Cmd(`h5ls -r "$(gather_absolute_paths_to_hdf5_input_files(m)[file_number])"`))
	nothing
end

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

function get_eltype_of_dataset(fn::AbstractString, gn::AbstractString, dn::AbstractString)::Type
	h = h5open(fn, "r")
	g = g_open(h, gn)
	d = d_open(g, dn)
	T = eltype(d)
	close(h)
	return T
end
function get_eltype_of_dataset(m::Measurement, gn::AbstractString, dn::AbstractString)::Type
	return get_eltype_of_dataset( gather_absolute_paths_to_hdf5_input_files(m)[1], gn, dn)
end

import HDF5.exists
function exists(fn::AbstractString, path::AbstractString)::Bool
	h = h5open(fn, "r")
	r = exists(h, path)
	close(h)
	return r
end
function exists(m::Measurement, path::AbstractString)::Bool
	exists(gather_absolute_paths_to_hdf5_input_files(m)[1], path)
end

function is_new_pulse_format(fn::AbstractString)::Bool
	new_pulse_format = false
	h = h5open(fn, "r")
	g = g_open(h, "DAQ_Data")
	d = d_open(g, "daq_pulses")
	if size(d, 1) > size(d, 2) new_pulse_format = true end
	close(h)
	return new_pulse_format
end

function is_new_pulse_format(m::Measurement)::Bool
	return is_new_pulse_format(gather_absolute_paths_to_hdf5_input_files(m)[1])
end