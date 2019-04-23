function write_analysis_result_dataset(fn::AbstractString, dataset_name::AbstractString, value)::Nothing
    h5f = h5open(fn, "r+")
    g = exists(h5f, "Results") ? g_open(h5f, "Results") : g_create(h5f, "Results")
    if exists(g, dataset_name) o_delete(g, dataset_name) end
    d_write(g, dataset_name, value)
    close(h5f)
    return nothing
end
function write_analysis_result_dataset(m::Measurement, dataset_name::AbstractString, value)::Nothing
    return write_analysis_result_dataset(gather_absolute_paths_to_hdf5_input_files(m)[1], dataset_name, value)
end
function read_analysis_result_dataset(fn::AbstractString, dataset_name::AbstractString)
    return h5read(fn, "Results/$(dataset_name)")
end
function read_analysis_result_dataset(m::Measurement, dataset_name::AbstractString)
    return read_analysis_result_dataset(gather_absolute_paths_to_hdf5_input_files(m)[1], dataset_name)
end

function write_analysis_result_attribute(fn::AbstractString, attr_name::AbstractString, value)::Nothing
    h5f = h5open(fn, "r+")
    g = g_open(h5f, "Results")
    if  exists(attrs(g), attr_name) a_delete(g, attr_name) end
    attrs(g)[attr_name] = value
    close(h5f)
    return nothing
end
function write_analysis_result_attribute(m::Measurement, attr_name::AbstractString, value)::Nothing
    return write_analysis_result_attribute(gather_absolute_paths_to_hdf5_input_files(m)[1], attr_name, value)
end
function read_analysis_result_attribute(fn::AbstractString, attr_name::AbstractString)
    h5f = h5open(fn, "r")
    g = g_open(h5f, "Results")
    data = if exists(attrs(g), attr_name)
        attr = a_open(g, attr_name)
        read(attr)
    else
        -1 # attr not existing
    end
    close(h5f)
    return data
end
function read_analysis_result_attribute(m::Measurement, attr_name::AbstractString)
    return read_analysis_result_attribute(gather_absolute_paths_to_hdf5_input_files(m)[1], attr_name)
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

# for m in dataset
#     c0 = GAT.read_core_calibration_factor_from_HDF5(m)
#     cm = GAT.read_crosstalk_and_calibration_matrix(m)
#     c = zeros(Float32, size(cm, 1)+1, size(cm, 1)+1)
#     c[1, 1] = c0
#     c[2:end, 2:end] = cm
#     GAT.write_analysis_result_dataset(m, "calibration_matrix", c)
# end

function read_decay_time_constants_from_HDF5(fn::AbstractString)
    h5_file =h5open(fn, "r+")
    if exists(h5_file, "Results")
        g_results = g_open(h5_file, "Results")
        # if exists(g_results, "decay_time_constants")
        #     decay_time_constants     = read(d_open(g_results, "decay_time_constants"))
        #     decay_time_constants_err = read(d_open(g_results, "decay_time_constants_err"))
        if exists(g_results, "tau_decay_constants")
            decay_time_constants     = read(d_open(g_results, "tau_decay_constants"))
            decay_time_constants_err = read(d_open(g_results, "tau_decay_constants_err"))
        else
            error("No decay time constants in this file.")
        end
    else
        error("No Results group in this file.")
    end
    close(h5_file)
    return decay_time_constants, decay_time_constants_err
end
function read_decay_time_constants_from_HDF5(m::Measurement)
    return read_decay_time_constants_from_HDF5(gather_absolute_paths_to_hdf5_input_files(m)[1])
end
# for m in dataset
#     tdcs = GAT.read_decay_time_constants_from_HDF5(m)[1]
#     GAT.write_analysis_result_dataset(m, "tau_decay_constants", tdcs)
# end
