function write_daq_calibration_matrix(fn::AbstractString, c::AbstractArray{<:Real, 2})::Nothing
    h5f = h5open(fn, "r+")
    g_results = exists(h5f, "Results") ? g_open(h5f, "Results") : g_create(h5f, "Results")
    if exists(g_results, "daq_calibration_matrix") o_delete(g_results, "daq_calibration_matrix") end
    d_write(g_results, "daq_calibration_matrix", c)
    close(h5f)
    return nothing
end
function write_daq_calibration_matrix(m::Measurement, c::AbstractArray{<:Real, 2})::Nothing
    write_daq_calibration_matrix(gather_absolute_paths_to_hdf5_input_files(m)[1], c)
    return nothing
end

function read_daq_calibration_matrix(fn::AbstractString)::AbstractArray{<:Real, 2} 
	return h5read(fn, "Results/daq_calibration_matrix")
end
function read_daq_calibration_matrix(m::Measurement)::AbstractArray{<:Real, 2}
	return h5read(gather_absolute_paths_to_hdf5_input_files(m)[1], "Results/daq_calibration_matrix")
end

function write_calibration_matrix(fn::AbstractString, c::AbstractArray{<:Real, 2})::Nothing
	h5f = h5open(fn, "r+")
	g_results = exists(h5f, "Results") ? g_open(h5f, "Results") : g_create(h5f, "Results")
	if exists(g_results, "calibration_matrix") o_delete(g_results, "calibration_matrix") end
	d_write(g_results, "calibration_matrix", c)
	close(h5f)
return nothing
end
function write_calibration_matrix(m::Measurement, c::AbstractArray{<:Real, 2})::Nothing
    write_calibration_matrix(gather_absolute_paths_to_hdf5_input_files(m)[1], c)
    return nothing
end

function read_calibration_matrix(fn::AbstractString)::AbstractArray{<:Real, 2} 
	return h5read(fn, "Results/calibration_matrix")
end
function read_calibration_matrix(m::Measurement)::AbstractArray{<:Real, 2}
	return read_calibration_matrix(gather_absolute_paths_to_hdf5_input_files(m)[1])
end



### OLD
function read_core_calibration_factor_from_HDF5(fn::AbstractString)::Real
    return h5read(fn, "Results/core_calibration_factor")
end
function read_core_calibration_factor_from_HDF5(m::Measurement)::Real
    return read_core_calibration_factor_from_HDF5(gather_absolute_paths_to_hdf5_input_files(m)[1])
end
function read_crosstalk_and_calibration_matrix(fn::AbstractString)::AbstractArray{<:Real, 2}
    return h5read(fn, "Results/inverse_calibration_matrix")
end
function read_crosstalk_and_calibration_matrix(m::Measurement)::AbstractArray{<:Real, 2}
    return read_crosstalk_and_calibration_matrix(gather_absolute_paths_to_hdf5_input_files(m)[1])
end

